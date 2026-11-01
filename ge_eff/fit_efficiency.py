"""Fit explicit source spectra and export candidate relative-efficiency curves.

Use Gaussian, left-tail and two-sided-tail hypotheses on identical counts.
Keep every model, including failed or poorly identified fits, in the receipt.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import time
import numpy as np
import pandas as pd
import uproot
from recalibrate import (KEYS, GROUPS, OLD, HERE, save, read_pars, energy,
                         fit_peak, peak_components)


def select_fit(counts, edges, line):
    attempts=[]
    for mode in [False,True,'both']:
        try:
            r=fit_peak(counts,edges,line,tail=mode,search=True)
            n,k=len(r['counts']),r['parameters']
            r['aicc']=r['deviance']+2*k+2*k*(k+1)/(n-k-1)
            r['identifiable_area']=bool(r['success'] and .00001<r['area_se']/r['area']<.20 and 1 not in r['bounds_indices'])
            attempts.append(r)
        except (ValueError,RuntimeError,FloatingPointError) as exc:
            attempts.append(dict(tail_model=mode,error=str(exc),identifiable_area=False))
    valid=[(i,r) for i,r in enumerate(attempts) if r.get('identifiable_area')]
    if not valid:
        usable=[(i,r) for i,r in enumerate(attempts) if 'aicc' in r]
        if not usable:return dict(line=line['id'],energy=line['energy'],status='no_fit',attempts=attempts)
        index,r=min(usable,key=lambda ir:ir[1]['aicc'])
        status='area_not_identified'
    else:
        index,r=min(valid,key=lambda ir:ir[1]['aicc'])
        status='review_peak_shape' if r['deviance_per_dof']>3 else 'conditional_fit'
    same=[a for a in attempts if 'counts' in a]
    assert all(np.array_equal(a['counts'],same[0]['counts']) and np.array_equal(a['edges'],same[0]['edges']) for a in same)
    selected=dict(r)
    selected.update(status=status,selected_model=index,attempts=attempts,
                    model_area_min=min(a['area'] for a in same),model_area_max=max(a['area'] for a in same),
                    model_selection='AICc on identical bins among converged fits with relative conditional area error <20%; all attempts retained')
    return selected


def fit_curve(points, degree=2):
    # The primary curve is supported by gamma lines from 245 to 1408 keV.
    # Low-energy source thresholds/complexes are shown separately, not extrapolated as measured.
    rows=[r for r in points if r['energy']>=240 and r['complete_channel_set'] and r['area']>0 and r['area_se']>0]
    if len(rows)<degree+4:raise RuntimeError('Insufficient common-set efficiency points')
    x=np.log(np.array([r['energy'] for r in rows])/334.)
    y=np.log([r['area']/r['intensity'] for r in rows])
    v=np.array([(r['area_se']/r['area'])**2+(r['intensity_se']/r['intensity'])**2 for r in rows])
    design=np.column_stack([[float(r['isotope']=='Ba133') for r in rows],
                           [float(r['isotope']=='Eu152') for r in rows],*[x**k for k in range(1,degree+1)]])
    weighted=design/np.sqrt(v[:,None]);normal=weighted.T@weighted
    if np.linalg.matrix_rank(normal)<degree+2:raise RuntimeError('Efficiency normalization/shape is not identifiable')
    covariance=np.linalg.inv(normal)
    p=covariance@(design.T@(y/v))
    residual=(y-design@p)/np.sqrt(v)
    return dict(model='exp(sum(c_k*ln(E_keV/334)^k,k=1..degree))',degree=degree,
                reference_energy_keV=334.,coefficients=p[2:].tolist(),parameters=p.tolist(),
                parameter_names=['ln_N_Ba','ln_N_Eu']+[f'c{k}' for k in range(1,degree+1)],
                covariance=covariance.tolist(),shape_covariance=covariance[2:,2:].tolist(),
                log_normalizations={'Ba133':float(p[0]),'Eu152':float(p[1])},
                valid_energy_keV=[min(r['energy'] for r in rows),max(r['energy'] for r in rows)],
                used_lines=[r['line'] for r in rows],chi2=float(residual@residual),ndof=len(rows)-len(p),
                conditional_pull=residual.tolist(),relative_point_residual=np.expm1(y-design@p).tolist(),
                information_condition=float(np.linalg.cond(normal)),
                covariance_scope='conditional local log-area fit; listed photon probability errors independent; no arbitrary chi-square rescaling',
                formal_replacement=False)


def curve_value(model, E):
    x=np.log(np.asarray(E,dtype=float)/334.)
    return np.exp(sum(c*x**(i+1) for i,c in enumerate(model['coefficients'])))


def old_value(group,E):
    a,b,c,d,e,f,g=OLD[group];E=np.asarray(E,dtype=float)
    x=np.log(E/100);y=np.log(E/1000)
    return np.exp(((a+b*x+c*x*x)**(-g)+(d+e*y+f*y*y)**(-g))**(-1/g))


def aggregate(fits,lines,clover=None):
    points=[]
    for group,keys in {'all':KEYS,**GROUPS}.items():
        for line in lines:
            items=[]
            for key in keys:
                if clover is not None and 212<=key<=215:continue
                items.append(fits[(str(key),line['id'])])
            if clover is not None and group in ['all','ring3']:items.append(fits[(clover,line['id'])])
            complete=all('area' in r and r.get('identifiable_area',False) for r in items)
            usable=[r for r in items if 'area' in r]
            points.append(dict(group=group,line=line['id'],energy=line['energy'],isotope=line['isotope'],
               intensity=line['intensity'],intensity_se=line['uncertainty'],
               area=sum(r['area'] for r in usable),area_se=float(np.sqrt(sum(r['area_se']**2 for r in usable))),
               complete_channel_set=complete,components=len(items),peak_review_components=sum(r.get('status')!='conditional_fit' for r in items),
               max_peak_deviance_per_dof=max((r.get('deviance_per_dof',0) for r in items),default=0),
               mode='singles' if clover is None else clover))
    return points


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adc',type=Path,required=True)
    parser.add_argument('--energy',type=Path,required=True)
    parser.add_argument('--seed',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    lines=json.loads((HERE/'source_lines.json').read_text())['lines'];seed=read_pars(args.seed)
    start=time.monotonic();fits={};selected=[]
    save(args.output/'sources.json',dict(created_at=datetime.datetime.now().astimezone().isoformat(),
      files=[dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in
             [args.adc,args.energy,args.seed,HERE/'source_lines.json',HERE/'recalibrate.py',Path(__file__)]]))
    with uproot.open(args.adc) as f:
        for key in KEYS:
            y,adc=f[f'ch{key:05d}/adc_clean'].to_numpy();edges=energy(adc-.5,seed[key])
            for line in lines:
                r=select_fit(y,edges,line);r['detector']=str(key);fits[(str(key),line['id'])]=r;selected.append(r)
            save(args.output/'peak_fits.json',selected)
            print(f'Finished single-crystal Ge{key}',flush=True)
    with uproot.open(args.energy) as f:
        for name in ['clover2_ab50','clover2_ab100','clover2_ab200']:
            y,edges=f[name].to_numpy()
            for line in lines:
                r=select_fit(y,edges,line);r['detector']=name;fits[(name,line['id'])]=r;selected.append(r)
            save(args.output/'peak_fits.json',selected)
            print(f'Finished {name}',flush=True)
    curves={};allpoints=[]
    for clover in [None,'clover2_ab50','clover2_ab100','clover2_ab200']:
        points=aggregate(fits,lines,clover);allpoints.extend(points);mode='singles' if clover is None else clover
        curves[mode]={}
        for group in ['all','ring3','ring4','ring5']:
            rows=[r for r in points if r['group']==group]
            curves[mode][group]={f'degree{d}':fit_curve(rows,d) for d in [2,3]}
    pd.DataFrame(allpoints).to_csv(args.output/'efficiency_points.csv',index=False)
    save(args.output/'efficiency_candidates.json',dict(schema='sm150_ge_source_relative_efficiency_v1',
         status='candidate_requires_peak_shape_and_condition_review',primary_mode='clover2_ab100',primary_degree=2,
         units='relative efficiency; epsilon(334 keV)=1',curves=curves,formal_replacement=False,
         limitations=['Source geometry/activity/live-time records unavailable; no absolute efficiency',
                      'Clover source coarse-time grouping is not exact production timing/grouping',
                      'Peak-shape, model-selection and true-coincidence-summing uncertainties are not included in conditional covariance',
                      'Source channel energy calibration has residual issues in Ge506 and Ge509; retained for review',
                      'No evaluated photon-probability covariance supplied; independent listed line uncertainties assumed']))
    checks=[]
    for group in ['all','ring3','ring4','ring5']:
        m=curves['clover2_ab100'][group]['degree2']
        for E in [334.,439.,586.,737.,1046.,1412.]:
            old=float(old_value(group,E)/old_value(group,334.));new=float(curve_value(m,E))
            checks.append(dict(group=group,energy_keV=E,old_eff_rel334=old,new_eff_rel334=new,change_percent=100*(new/old-1),
                               supported=m['valid_energy_keV'][0]<=E<=m['valid_energy_keV'][1]))
    pd.DataFrame(checks).to_csv(args.output/'old_new_comparison.csv',index=False)
    summary=[]
    for r in selected:
        summary.append({k:r.get(k) for k in ['detector','line','energy','status','area','area_se','sigma','selected_model','deviance_per_dof','identifiable_area','model_area_min','model_area_max']})
    pd.DataFrame(summary).to_csv(args.output/'peak_summary.csv',index=False)
    save(args.output/'completion.json',dict(status='complete',elapsed_seconds=time.monotonic()-start,peak_regions=len(selected),curve_variants=32,
                  formal_replacement=False,candidate_quality='review_required'))
    print('Candidate efficiency fits completed; formal parameters unchanged.',flush=True)


if __name__=='__main__':
    main()
