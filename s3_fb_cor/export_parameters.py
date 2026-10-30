"""Export saved decisions without fitting, installing or rebuilding data."""
from pathlib import Path
from datetime import datetime
import json
import hashlib
import numpy as np
import pandas as pd

from collections import Counter
import copy
from s3_energy_attenuation.workflow_io import digest

def export_parameters(base, old, review, confirmed, out):
    flagged = {(r["run"], int(r["index"])): r for r in review}
    rows=[];runs=[];positions=[];pixels=[];package={};review_rows=[]
    for run,o in old['runs'].items():
        e=confirmed.get(run);ref=o['reference_run'];active=o['active_rings']
        m=e['models'][e['final']] if e else {'a':o['A_keV'],'b':o['B'],'radial':o['q_position']}
        a=np.asarray(m['a']);b=np.asarray(m['b']);q=np.asarray(m['radial']);assert np.isfinite(a).all() and np.all(b>0) and np.all(q>0)
        actual_source=e['final'] if e else '0229_inherited';coeff=[]
        sm={(r['model'],r['index']):r for r in e['spectral']} if e else {}
        pm={(r['model'],r['index']):r for r in e['channels']} if e else {}
        for i in range(56):
            name=f'ring{i+1:02d}' if i<24 else f'sector{i-23:02d}'
            source=next(c for c in base['runs'][run]['channels'] if c['index']==i)
            refr=next(c for c in base['runs'][ref]['channels'] if c['index']==i)
            original=flagged.get((run,i));production=i<active or i>=24
            row={'run':run,'period':o['period'],'reference_run':ref,'index':i,'channel':name,'production':production,'source':actual_source,
                 'decision':e['decision'] if e else 'outside_targeted_review_preserved_0229','A_keV':a[i],'B':b[i],
                 'p0_keV':a[i],'p1_reference':refr['B'],'k_strip':b[i]/refr['B'],
                 'position_q_min':q[i] if i<24 else min(q[:active]),'position_q_max':q[i] if i<24 else max(q[:active]),
                 'formal_A_keV':source['old_A_keV'],'formal_B':source['old_B'],'A_0228_keV':source['A_keV'],'B_0228':source['B'],
                 'A_0229_keV':o['A_keV'][i],'B_0229':o['B'][i],
                 'original_review_reason':original['channel_confirmation_flags'] if original else '',
                 'confirmation_status':sm.get((actual_source,i),{}).get('status','not_rechecked'),
                 'confirmation_flags':'|'.join(e['confirmation_flags']) if e else ''}
            for model,label in [('formal','formal'),('0228','0228'),('0229','0229'),(actual_source,'final')]:
                s=sm.get((model,i),{});p=pm.get((model,i),{})
                row[f'D_{label}']=s.get('D');row[f'noise_{label}']=s.get('noise_floor');row[f'width_{label}_keV']=p.get('width68_keV');row[f'median_{label}_keV']=p.get('median_keV')
                if label=='final':row['confirmation_prompt_pairs']=p.get('n');row['confirmation_pixels']=len(s.get('pixels',[]))
            rows.append(row);coeff.append(row)
            if original:
                rr=copy.deepcopy(row);rr['previous_reported_D_0228']=float(original['confirmation_D_previous']) if original['confirmation_D_previous'] else None
                rr['previous_reported_D_0229']=float(original['confirmation_D_new']) if original['confirmation_D_new'] else None
                rr['coverage_restored']=original['channel_confirmation_flags']=='insufficient_confirmation_pixel_coverage' and row['confirmation_status']=='comparable'
                review_rows.append(rr)
        for r in range(24):
            positions.append({'run':run,'period':o['period'],'ring':r+1,'q':q[r],'source':actual_source,'production':r<active})
            for s in range(24,56):pixels.append({'run':run,'period':o['period'],'ring':r+1,'sector':s-23,'q_position':q[r],
                                               'A_ring_keV':q[r]*a[r],'B_ring':q[r]*b[r],'A_sector_keV':q[r]*a[s],'B_sector':q[r]*b[s],'production':r<active})
        package[run]={'period':o['period'],'reference_run':ref,'active_rings':active,'source':actual_source,'A_keV':a.tolist(),'B':b.tolist(),'q_position':q.tolist(),
                      'channels':coeff,'requires_ring_context':bool(np.any(q!=1)),'decision':e['decision'] if e else 'outside_targeted_review_preserved_0229'}
        if e:
            d={'run':run,'period':e['period'],'preferred':e['preferred'],'final':e['final'],'decision':e['decision'],
               'confirmation_flags':'|'.join(e['confirmation_flags']),'original_review_rows':sum(r['run']==run for r in review)}
            for model,label in [('formal','formal'),('0228','0228'),('0229','0229'),(e['final'],'final')]:
                d[f'D_{label}']=e['spectral_scores'][model];g=e['global_metrics'][model]['normalized']
                for key in ['median_keV','width68_keV','score_keV','tail2MeV']:d[f'{key}_{label}']=g[key]
            runs.append(d)
    for name,records in [('channel_coefficients',rows),('ring_position_factors',positions),('pixel_coefficients',pixels)]:
        pd.DataFrame(records).to_csv(out/f'{name}.csv',index=False,float_format='%.17g')
    pd.DataFrame(runs).to_csv(out/'run_decisions.csv',index=False,float_format='%.17g')
    pd.DataFrame(review_rows).to_csv(out/'review_decisions.csv',index=False,float_format='%.17g')
    document={'schema':'sm150_s3_targeted_review_candidate_v1','created':datetime.now().astimezone().isoformat(),'runs':package,
              'formula':'Er=q[run,matched_ring]*(A[run,ring]+B[run,ring]*E0r); Es=q[run,matched_ring]*(A[run,sector]+B[run,sector]*E0s)',
              'scope':f'{len(package)}-run candidate; {len(confirmed)} runs rechecked, {len(package)-len(confirmed)} inherited0229. No production adoption. Requires matched-ring context when q differs from1.',
              'reference':'Relative reference energy scale from frozen0228; no independent absolute energy calibration or radiation-damage claim.'}
    (out/'coefficients.json').write_text(json.dumps(document,indent=2,allow_nan=False)+'\n')
    summary={'status':'complete','targeted_runs':len(runs),'review_rows':len(review_rows),'candidate_runs':len(package),
             'channel_rows':len(rows),'position_rows':len(positions),'pixel_rows':len(pixels),
             'run_decisions':dict(Counter(r['decision'] for r in runs)),'models':dict(Counter(r['final'] for r in runs)),
             'old_low_coverage_rows':sum(r['original_review_reason']=='insufficient_confirmation_pixel_coverage' for r in review_rows),
             'coverage_restored':sum(bool(r['coverage_restored']) for r in review_rows),
             'review_rows_final_support':dict(Counter(r['confirmation_status'] for r in review_rows)),
             'new_raw_event_reads_during_export':0,
             'no_production_adoption':True}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({p.name:digest(p) for p in out.iterdir() if p.is_file()},indent=2)+'\n')
    return summary
