"""Export saved decisions without fitting, installing or rebuilding data."""
from pathlib import Path
from datetime import datetime
import json
import hashlib
import numpy as np
import pandas as pd

def export_parameters(sources, previous, evaluation_dir, dest):
    strips=[];radial=[];pixels=[];run_rows=[];pixel_estimates=[];package={}
    for row in sources:
     run=row['run'];e=json.loads((evaluation_dir/f'{run}.json').read_text());f=json.loads(Path(row['conditional']).read_text())
     base=previous['runs'][run];ref=previous['runs'][f['reference_run']];channels=[]
     assert e['conditional_sha256']==hashlib.sha256(Path(row['conditional']).read_bytes()).hexdigest()
     for c,old,refc in zip(f['channels'],base['channels'],ref['channels']):
      i=c['index'];a=e['selected_a'][i];b=e['selected_b'][i];q=e['selected_radial'][i] if i<24 else None
      fixed=f['run']==f['reference_run'];diagnostic=i<24 and i>=e['active_rings']
      status='fixed_reference' if fixed else ('diagnostic_outer_ring' if diagnostic else ('retained_0228' if e['selected']=='0228' else 'run_candidate_passed'))
      spectral={x['model']:x for x in e['spectral_metrics'] if x['fold']==3 and x['channel_index']==i}
      fb={x['model']:x for x in e['channel_pair_metrics'] if x['fold']==3 and x['index']==i}
      d0=spectral.get('0228',{}).get('D');d1=spectral.get(e['selected'],{}).get('D')
      flags=[]
      if d0 is None or d1 is None:flags.append('insufficient_confirmation_pixel_coverage')
      elif d1>d0+.02:flags.append('channel_spectrum_worse_than_0228')
      if not fixed and not diagnostic and e['selected']!='0228' and flags:status='run_update_channel_review'
      item={'run':run,'period':f['period'],'reference_run':f['reference_run'],'channel':c['channel'],'index':i,'model':e['selected'],'status':status,
       'A_keV':a,'B':b,'p0_keV':a,'p1_reference':refc['B'],'k_strip':b/refc['B'],'position_q_for_ring':q,
       'old_A_keV':old['old_A_keV'],'old_B':old['old_B'],'previous_A_keV':old['A_keV'],'previous_B':old['B'],
       'measured_residual_q':c['q'] if c['eligible_pixels'] else None,'training_supported':c['training_supported'],'usable_pixels':len(c['eligible_pixels']),
       'training_net':c['total_training_net'],'training_flags':'|'.join(c['training_flags']),'selection_model':e['preferred'],'confirmation_flags':'|'.join(e['confirmation_flags']),
       'channel_confirmation_flags':'|'.join(flags),'confirmation_D_previous':d0,'confirmation_D_new':d1,
       'confirmation_width_previous_keV':fb.get('0228',{}).get('width68_keV'),'confirmation_width_new_keV':fb.get(e['selected'],{}).get('width68_keV')}
      strips.append(item);channels.append(item)
     for r in range(24):
      q=e['selected_radial'][r];radial.append({'run':run,'period':f['period'],'ring':r+1,'q':q,'model':e['selected'],'production_used':r<e['active_rings']})
      for s in range(24,56):
       pixels.append({'run':run,'period':f['period'],'ring':r+1,'sector':s-23,'q_position':q,'A_ring_keV':q*e['selected_a'][r],'B_ring':q*e['selected_b'][r],
        'A_sector_keV':q*e['selected_a'][s],'B_sector':q*e['selected_b'][s],'production_used':r<e['active_rings']})
     for p in f['pixels']:
      item={'run':run,'period':f['period'],**p};pixel_estimates.append(item)
     package[run]={'period':f['period'],'reference_run':f['reference_run'],'selected_model':e['selected'],'active_rings':e['active_rings'],
      'A_keV':e['selected_a'],'B':e['selected_b'],'q_position':e['selected_radial'],'channels':channels,
      'requires_ring_context':bool(np.any(np.asarray(e['selected_radial'])!=1))}
     metric=e['pair_summary'][e['selected']]['3']['normalized'];prior=e['pair_summary']['0228']['3']['normalized'];formal=e['pair_summary']['formal']['3']['normalized']
     run_rows.append({'run':run,'period':f['period'],'selected':e['selected'],'preferred':e['preferred'],'confirmation_flags':'|'.join(e['confirmation_flags']),
      'sector_coverage':sum(len(c['eligible_pixels'])>=3 and c['total_training_net']>=300 for c in f['channels'][24:]),
      'training_supported_sectors':sum(c['training_supported'] for c in f['channels'][24:]),'training_supported_rings':sum(c['training_supported'] for c in f['channels'][:e['active_rings']]),
      'spectral_D_formal':e['spectral_scores']['formal']['3'],'spectral_D_previous':e['spectral_scores']['0228']['3'],'spectral_D_new':e['spectral_scores'][e['selected']]['3'],
      'fb_width_formal_keV':formal['width68_keV'],'fb_width_previous_keV':prior['width68_keV'],'fb_width_new_keV':metric['width68_keV'],
      'fb_median_formal_keV':formal['median_keV'],'fb_median_previous_keV':prior['median_keV'],'fb_median_new_keV':metric['median_keV'],
      'sample_entries':row['sample_entries'],'source_entries':row['source_entries'],'sample_fraction':row['sample_fraction']})
    pd.DataFrame(strips).to_csv(dest/'channel_coefficients.csv',index=False,float_format='%.17g')
    pd.DataFrame(radial).to_csv(dest/'ring_position_factors.csv',index=False,float_format='%.17g')
    pd.DataFrame(pixels).to_csv(dest/'pixel_coefficients.csv',index=False,float_format='%.17g')
    pd.DataFrame(run_rows).to_csv(dest/'run_summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(pixel_estimates).to_csv(dest/'pixel_estimates.csv',index=False,float_format='%.17g')
    data={'schema':'sm150_s3_conditioned_candidate_v1','status':'candidate_not_production','created':datetime.now().astimezone().isoformat(),
     'formula':'E0 = hardware_gain*(c0+c1*ADC+c2*ADC^2); Er=q_position[ring]*(A[ring]+B[ring]*Er0); Es=q_position[ring]*(A[sector]+B[sector]*Es0)',
     'factorization':'A=p0; B=p1_reference*k_strip. q_position acts on both sides of the same event. A single sector has different total coefficients at different rings when q_position differs.',
     'application':'Full candidate needs paired ring context. Legacy strip tables cannot express the position term; this intermediate schema must pass targeted review before native build export; do not copy only k_strip/p0/p1 into pars.',
     'interpretation':'Relative to frozen reference spectra, conditional on stable physical distributions. Does not uniquely identify radiation damage or restore an independently established absolute energy.',
     'runs':package}
    (dest/'coefficients.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    result={'status':'passed','runs':len(package),'strip_rows':len(strips),'position_rows':len(radial),'pixel_rows':len(pixels),'models':pd.DataFrame(run_rows).groupby(['period','selected']).size().to_dict()}
    assert len(strips)==len(sources)*56 and len(radial)==len(sources)*24 and len(pixels)==len(sources)*24*32
    result['models']={str(k):int(v) for k,v in result['models'].items()}
    for r in package.values():
     assert np.isfinite(r['A_keV']).all() and np.all(np.array(r['B'])>0) and np.all(np.array(r['q_position'])>0)
     if r['reference_run'] in package:
      for c in r['channels']:assert abs(c['p1_reference']*c['k_strip']-c['B'])<1e-12
    (dest/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    (dest/'manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file()},indent=2)+'\n')
    return result
