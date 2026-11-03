#!/usr/bin/env python3
"""Recalibrate saved Si/S3 peak centers without rerunning the source response fit."""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def linear_fit(x,y):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)<2 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.diff(x)>0):
        return None
    centered=x-x.mean();den=np.dot(centered,centered)
    if den<=0:return None
    b=np.dot(centered,y-y.mean())/den;a=y.mean()-b*x.mean()
    if not np.isfinite(a+b) or b<=0:return None
    return np.array([a,b])


def curve(channel,ids):
    peaks=channel.get('peaks',[]);indices=[i for i,p in enumerate(peaks) if p['id'] in ids]
    selected=[peaks[i] for i in indices]
    if len(selected)!=len(ids):return None
    x=np.array([p.get('center',np.nan) for p in selected]);y=np.array([p['energy_keV'] for p in selected])
    coeff=linear_fit(x,y)
    if coeff is None:return None
    old=np.array(channel['seed_coefficients']);old_energy=old[0]+old[1]*x+old[2]*x*x
    residual=coeff[0]+coeff[1]*x-y
    cov=None;joint=channel.get('joint_fit',{})
    if 'covariance' in joint:
        stride=joint.get('parameter_stride',5);ix=np.array(indices)*stride+1
        centers_cov=np.asarray(joint['covariance'])[np.ix_(ix,ix)]
        design=np.column_stack([np.ones(len(x)),x]);normal_inv=np.linalg.inv(design.T@design)
        jac=np.column_stack([normal_inv@np.array([-coeff[1],-r-coeff[1]*xi]) for xi,r in zip(x,residual)])
        candidate=jac@centers_cov@jac.T
        if np.all(np.isfinite(candidate)) and np.min(np.diag(candidate))>=0:cov=candidate.tolist()
    return dict(coefficients=[float(coeff[0]),float(coeff[1]),0.],center_ids=[p['id'] for p in selected],
                centers_adc=x.tolist(),reference_keV=y.tolist(),residual_keV=residual.tolist(),
                old_residual_keV=(old_energy-y).tolist(),anchors=len(x),
                rms_keV=float(np.sqrt(np.mean(residual**2))),max_residual_keV=float(np.max(np.abs(residual))),
                old_rms_keV=float(np.sqrt(np.mean((old_energy-y)**2))),covariance=cov,
                linearity_check='three_point_internal_check' if len(x)>2 else 'not_testable_with_two_points',
                error_scope='Conditional joint-center covariance; excludes source and response-model systematics')


def compare_curves(base,other):
    energies=np.array([min(base['reference_keV']),max(base['reference_keV'])]);a,b,_=base['coefficients']
    adc=(energies-a)/b;aa,bb,_=other['coefficients']
    return float(np.max(np.abs(aa+bb*adc-energies)))


def low_key(key):
    slot=key//100%100
    return key//10000*10000+100*{2:9,9:2,3:10,10:3}.get(slot,slot)+key%100


def run(args):
    out=Path(args.output)
    if out.exists():raise FileExistsError(f'Output already exists: {out}')
    primary_path=Path(args.primary).resolve();window_paths={120:Path(args.window120).resolve(),160:Path(args.window160).resolve()}
    primary=json.loads(primary_path.read_text());windows={w:json.loads(p.read_text()) for w,p in window_paths.items()}
    assert all(set(d['channels'])==set(primary['channels']) for d in windows.values())
    previous=pd.read_csv(args.previous_status,dtype={'run':str})
    previous=previous[previous['run']==args.run].set_index('key')
    channels={};rows=[];peak_rows=[];window_rows=[];start_rows=[]
    for key,c in primary['channels'].items():
        key_int=int(key);peaks=c.get('peaks',[]);joint=c.get('joint_fit',{})
        ids=[p['id'] for p in peaks if p.get('center_usable',False)]
        all_ids=[p['id'] for p in peaks if np.isfinite(p.get('center',np.nan))]
        candidate=curve(c,ids);diagnostic=curve(c,all_ids)
        numeric_ok=bool(candidate and joint.get('success') and joint.get('covariance_status',0)>=2)
        base=candidate or diagnostic
        comparison=[];window_ok=numeric_ok;max_curve=0.;max_center=0.
        for width,wdata in windows.items():
            wc=wdata['channels'][key];chosen=base['center_ids'] if base else []
            alt=curve(wc,chosen);wp={p['id']:p for p in wc.get('peaks',[])}
            reliable=bool(alt and wc.get('joint_fit',{}).get('success') and
                          wc.get('joint_fit',{}).get('covariance_status',0)>=2 and
                          all(wp[p].get('center_usable',False) for p in chosen))
            difference=compare_curves(base,alt) if base and alt else np.nan
            center_difference=float(np.max(np.abs(np.array(base['centers_adc'])-np.array(alt['centers_adc'])))*base['coefficients'][1]) if base and alt else np.nan
            max_curve=max(max_curve,difference) if np.isfinite(difference) else np.inf
            max_center=max(max_center,center_difference) if np.isfinite(center_difference) else np.inf
            window_ok &= reliable
            value=dict(run=args.run,key=key_int,half_width_keV=width,same_primary_anchors=chosen,
                       center_checks_pass=reliable,curve_difference_keV=difference,center_difference_keV=center_difference,
                       shape_bound_peaks=[p for p in chosen if wp.get(p,{}).get('shape_bound')],
                       coefficients=alt['coefficients'] if alt else None,
                       diagnostic_only=not(numeric_ok and reliable))
            comparison.append(value);window_rows.append(value)
        starts=[];spread=0.;center_spread=0.;start_ok=bool(base)
        valid=[(i,a) for i,a in enumerate(joint.get('attempts',[])) if a.get('success') and a.get('covariance_status',0)>=2]
        best=min([joint.get('deviance',np.inf)]+[a['deviance'] for _,a in valid])
        chosen_ix=[i for i,p in enumerate(peaks) if base and p['id'] in base['center_ids']]
        stride=joint.get('parameter_stride',5)
        for index,a in valid:
            if a['deviance']>best+args.near_deviance:continue
            xx=np.array([a['parameters'][stride*i+1] for i in chosen_ix]);yy=np.array(base['reference_keV']) if base else np.array([])
            cp=linear_fit(xx,yy)
            other={'coefficients':[*cp,0.]} if cp is not None else None
            diff=compare_curves(base,other) if base and other else np.nan
            center_diff=float(np.max(np.abs(xx-np.array(base['centers_adc'])))*base['coefficients'][1]) if base and len(xx) else np.nan
            if np.isfinite(diff):spread=max(spread,diff);center_spread=max(center_spread,center_diff)
            else:start_ok=False
            value=dict(run=args.run,key=key_int,attempt=index,deviance=float(a['deviance']),
                       delta_deviance=float(a['deviance']-best),curve_difference_keV=diff,
                       center_difference_keV=center_diff,shape_bound_indices=a.get('bounds_indices',[]),
                       coefficients=[*cp,0.] if cp is not None else None)
            starts.append(value);start_rows.append(value)
        reasons=[]
        if not numeric_ok:reasons.append('fewer_than_two_numerically_usable_centers')
        if numeric_ok and candidate['max_residual_keV']>args.maximum_residual_keV:reasons.append('linear_calibration_residual')
        if numeric_ok and not window_ok:reasons.append('window_center_fit_review')
        if numeric_ok and (max_curve>args.stability_keV or max_center>args.stability_keV):reasons.append('window_center_or_curve_shift')
        if numeric_ok and (not start_ok or spread>args.stability_keV or center_spread>args.stability_keV):reasons.append('near_minimum_start_review')
        stable=numeric_ok and not reasons
        status='window_checked_candidate' if stable else 'candidate_needs_review' if numeric_ok else 'unreliable_centers'
        detector='Si array' if key_int<10900 else 'S3 sector' if key_int<11100 else 'S3 ring'
        cp=candidate['coefficients'] if candidate else [np.nan]*3;dp=diagnostic['coefficients'] if diagnostic else [np.nan]*3
        covariance=candidate['covariance'] if candidate else None
        old=c['seed_coefficients'];old_row=previous.loc[key_int]
        row=dict(run=args.run,key=key_int,low_beam_key=low_key(key_int) if args.run=='0457' else None,
                 detector=detector,status=status,reasons=';'.join(reasons),candidate_available=numeric_ok,
                 window_checked=stable,old_c0=old[0],old_c1=old[1],c0=cp[0],c1=cp[1],c2=cp[2],
                 c0_conditional_se=np.sqrt(covariance[0][0]) if covariance else np.nan,
                 c1_conditional_se=np.sqrt(covariance[1][1]) if covariance else np.nan,
                 all_center_diagnostic_c0=dp[0],all_center_diagnostic_c1=dp[1],
                 fitted_peaks=len(peaks),usable_peaks=len(ids),used_anchors=candidate['anchors'] if candidate else 0,
                 old_rms_keV=candidate['old_rms_keV'] if candidate else np.nan,
                 new_rms_keV=candidate['rms_keV'] if candidate else np.nan,
                 max_residual_keV=candidate['max_residual_keV'] if candidate else np.nan,
                 all_center_rms_keV=diagnostic['rms_keV'] if diagnostic else np.nan,
                 slope_change_percent=100*(cp[1]/old[1]-1) if candidate else np.nan,
                 max_window_curve_difference_keV=max_curve,max_window_center_difference_keV=max_center,
                 window_center_checks_pass=window_ok,near_start_curve_difference_keV=spread,
                 near_start_center_difference_keV=center_spread,
                 local_shape_review_peaks=sum(p.get('deviance_per_dof',np.inf)>8 for p in peaks),
                 previous_status=old_row['status'],previous_stable=bool(old_row['recommended_candidate']))
        rows.append(row)
        for p in peaks:
            xx=p.get('center',np.nan);yy=p['energy_keV'];usable=p.get('center_usable',False)
            peak_rows.append(dict(run=args.run,key=key_int,detector=detector,line=p['id'],reference_keV=yy,
                center_adc=xx,center_se_adc=p.get('center_se'),center_usable=usable,shape_bound=p.get('shape_bound'),
                used_for_new_calibration=numeric_ok and p['id'] in ids,originally_used=p.get('used_for_calibration'),
                local_D_ndof=p.get('deviance_per_dof'),old_residual_keV=old[0]+old[1]*xx+old[2]*xx*xx-yy,
                new_residual_keV=cp[0]+cp[1]*xx-yy if candidate else np.nan,
                all_center_diagnostic_residual_keV=dp[0]+dp[1]*xx-yy if diagnostic else np.nan))
        channels[key]=dict(**row,primary_candidate=candidate,all_center_diagnostic=diagnostic,
                           selected_centers=ids,excluded_centers=[p['id'] for p in peaks if p['id'] not in ids],
                           window_comparisons=comparison,near_minimum_starts=starts)
    table=pd.DataFrame(rows);out.mkdir(parents=True)
    table.to_csv(out/'all_channels.csv',index=False,float_format='%.17g')
    table[table.candidate_available].to_csv(out/'candidate_coefficients.csv',index=False,float_format='%.17g')
    table[table.window_checked].to_csv(out/'window_checked_coefficients.csv',index=False,float_format='%.17g')
    pd.DataFrame(peak_rows).to_csv(out/'peak_positions.csv',index=False,float_format='%.17g')
    pd.DataFrame(window_rows).to_csv(out/'window_comparisons.csv',index=False,float_format='%.17g')
    pd.DataFrame(start_rows).to_csv(out/'start_comparisons.csv',index=False,float_format='%.17g')
    if args.run=='0457':
        remapped=table.copy();remapped['source_key']=remapped.key;remapped['key']=remapped.key.map(low_key)
        remapped.to_csv(out/'low_beam_remapped_channels.csv',index=False,float_format='%.17g')
    sources={str(primary_path):sha(primary_path),**{str(p):sha(p) for p in window_paths.values()}}
    sources[str(Path(args.previous_status).resolve())]=sha(args.previous_status)
    result=dict(schema='sm150_saved_peak_center_linear_calibration_v1',run=args.run,
                generated_at=datetime.now().astimezone().isoformat(),channels=channels,
                formula='E_base_keV = c0 + c1 * original_ADC',
                source_primary=str(primary_path),sources_sha256=sources,program_sha256=sha(__file__),
                settings=dict(maximum_residual_keV=args.maximum_residual_keV,stability_keV=args.stability_keV,
                              near_deviance=args.near_deviance,weights='equal numerically usable centers',
                              local_deviance_veto=False,two_point_linearity_not_testable=True),
                summary=dict(channel_instances=len(table),candidates=int(table.candidate_available.sum()),
                             window_checked=int(table.window_checked.sum()),review=int(((table.candidate_available)&(~table.window_checked)).sum()),
                             unreliable_centers=int((~table.candidate_available).sum()),
                             all_center_diagnostic_curves=int(table.all_center_diagnostic_c0.notna().sum())),
                formal_replacement=False,peak_shapes_refitted=False,
                limitations=['Source/response systematics not included in conditional covariance',
                             'The source interval does not validate in-beam high-energy extrapolation',
                             'Failed-center diagnostic curves are not accepted calibration coefficients'])
    (out/'candidate.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps(dict(formal_replacement=False,partial_candidate_tables=True,
        files={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}),indent=2)+'\n')
    print(json.dumps(result['summary'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['primary','window120','window160','previous-status','output']:parser.add_argument('--'+name,required=True)
    parser.add_argument('--run',required=True)
    parser.add_argument('--stability-keV',dest='stability_keV',type=float,default=5.)
    parser.add_argument('--maximum-residual-keV',dest='maximum_residual_keV',type=float,default=30.)
    parser.add_argument('--near-deviance',type=float,default=2.)
    run(parser.parse_args())
