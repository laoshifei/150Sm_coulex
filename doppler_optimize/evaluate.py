#!/usr/bin/env python3
"""Compare correction candidates on identical weighted spectra, without changing yields."""
from pathlib import Path
import argparse
import copy
import json
import re

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import ndtr
import uproot


def core_fit(edges,counts,variance,energy=334.0,half_width=16):
    """Bin-integrated Gaussian core plus a linear continuum; diagnostic only."""
    centers=(edges[:-1]+edges[1:])/2
    take=(centers>=energy-half_width)&(centers<energy+half_width)
    lo,hi=edges[:-1][take],edges[1:][take]
    y=counts[take];v=variance[take];x=centers[take];dx=hi-lo
    if len(y)<12 or not np.isfinite(y).all() or np.sum(v)<100:return {'valid':False,'reason':'insufficient_variance_or_bins'}
    continuum=np.median(np.r_[y[:8]/dx[:8],y[-8:]/dx[-8:]])
    near=abs(x-energy)<6
    mu0=float(x[near][np.argmax(y[near])]);area0=max(1.,float(np.sum(y-continuum*dx)))
    def prediction(p):
        area,mu,sigma,c,b=p
        return area*(ndtr((hi-mu)/sigma)-ndtr((lo-mu)/sigma))+c*dx+b*((hi-energy)**2-(lo-energy)**2)/2
    error=np.sqrt(np.maximum(v,1.0))
    result=least_squares(lambda p:(prediction(p)-y)/error,[area0,mu0,2.,continuum,0.],
        bounds=([0,energy-6,.3,-np.inf,-np.inf],[np.inf,energy+6,10,np.inf,np.inf]),
        x_scale='jac',max_nfev=400,ftol=1e-9,xtol=1e-9,gtol=1e-9)
    p=result.x;sv=np.linalg.svd(result.jac,compute_uv=False)
    cov=np.linalg.pinv(result.jac.T@result.jac);unc=np.sqrt(np.maximum(np.diag(cov),0))
    chi2=float(np.sum(result.fun**2));dof=len(y)-len(p)
    bound=bool(p[2]<.31 or p[2]>9.99 or abs(p[1]-energy)>5.99)
    sig=float(p[0]/unc[0]) if unc[0]>0 else 0
    valid=bool(result.success and not bound and p[0]>100 and sig>5 and chi2/dof<10 and sv[-1]>sv[0]*1e-12)
    return dict(valid=valid,reason='ok' if valid else 'weak_peak_boundary_or_shape_failure',
        area=float(p[0]),mu=float(p[1]),sigma=float(p[2]),fwhm=float(2.354820045*p[2]),
        mu_error=float(unc[1]),sigma_error=float(unc[2]),area_error=float(unc[0]),
        continuum=float(p[3]),slope=float(p[4]),chi2=chi2,dof=dof,chi2_per_dof=chi2/dof,
        score=float(p[2]**2+(p[1]-energy)**2),bound=bound,success=bool(result.success),nfev=int(result.nfev),
        condition=float(sv[0]/sv[-1]),covariance=cov.tolist())


def read_fit(path,split,peaks=(334.,)):
    rows=[]
    with uproot.open(path) as f:
        m=json.loads(str(f['sm150_doppler_manifest']));configs=m['parameters']
        for key in f.keys(cycle=False):
            hit=re.fullmatch(r'v(\d+)_(high1|low|high2)_(Si|S3)_all',key)
            if not hit:continue
            i,period,family=hit.groups();h=f[key]
            for peak in peaks:
                fit=core_fit(h.axis().edges(),h.values(),h.variances(),peak)
                rows.append(dict(split=split,variant=int(i),label=configs[int(i)]['label'],period=period,family=family,peak=peak,**fit))
    return rows,configs


def summarize(rows):
    frame=pd.DataFrame(rows);base=frame[(frame.variant==0)&(frame.peak==334)]
    comparison=frame[frame.peak==334].merge(base[['split','period','family','score','fwhm','mu','valid']],on=['split','period','family'],suffixes=('','_baseline'))
    comparison['ratio']=comparison.score/comparison.score_baseline
    out=[]
    for (split,family,variant),g in comparison.groupby(['split','family','variant']):
        valid=bool(len(g)==3 and g.valid.all() and g.valid_baseline.all())
        out.append(dict(split=split,family=family,variant=int(variant),label=g.label.iloc[0],
                        valid=valid,mean_score_ratio=float(g.ratio.mean()),worst_score_ratio=float(g.ratio.max()),
                        mean_fwhm=float(g.fwhm.mean()),max_abs_centroid_bias=float(abs(g.mu-334).max())))
    return pd.DataFrame(out)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--task',required=True);args=parser.parse_args();task=Path(args.task)
    rows=[]
    for split in ['train','validation']:
        r,configs=read_fit(task/f'{split}_scan.root',split);rows+=r
    ranked=summarize(rows)
    # Select beta scales from training runs only; held-out runs do not select them.
    candidate=copy.deepcopy(configs[0]);candidate['label']='training_selected_beta';candidate['status']='pilot_candidate_not_formally_adopted'
    choices=[]
    for family in ['Si','S3']:
        options=ranked[(ranked.split=='train')&(ranked.family==family)&ranked.valid&((ranked.variant==0)|ranked.label.str.startswith(f'{family}_energy_scale_'))].sort_values('mean_score_ratio')
        if options.empty:choices.append({'family':family,'reason':'no_valid_training_fit'});continue
        choice=options.iloc[0];candidate[family]['energy_scale']=configs[int(choice.variant)][family]['energy_scale']
        choices.append({'family':family,'variant':int(choice.variant),'training_score_ratio':float(choice.mean_score_ratio),'energy_scale':candidate[family]['energy_scale']})
    combined=[configs[0],candidate]
    # A bounded second stage tests geometry about the selected beta baseline.
    for family in ['Si','S3']:
        for field,values in [('theta_offset_deg',[-2,-1,1,2]),('phi_offset_deg',[-5,5])]:
            for value in values:
                c=copy.deepcopy(candidate);c[family][field]=value;c['label']=f'beta_plus_{family}_{field}_{value}';combined.append(c)
    (task/'combined_parameters.json').write_text(json.dumps(combined,indent=2)+'\n')
    (task/'core_fits.json').write_text(json.dumps(rows,indent=2)+'\n')
    pd.DataFrame([{k:v for k,v in r.items() if k!='covariance'} for r in rows]).to_csv(task/'core_fits.csv',index=False)
    ranked.to_csv(task/'scan_ranking.csv',index=False)
    (task/'training_choices.json').write_text(json.dumps(choices,indent=2)+'\n')
    print(json.dumps(choices,indent=2))
    for family in ['Si','S3']:
        print(family,ranked[(ranked.split=='train')&(ranked.family==family)].sort_values('mean_score_ratio').head(6).to_string(index=False))


if __name__=='__main__':main()
