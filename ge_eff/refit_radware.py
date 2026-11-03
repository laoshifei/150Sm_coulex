"""Refit saved source areas with the original seven-parameter Radware function.

C=0 and G=15 retain the original macro's fixed parameters. No production
parameter is installed. All starts, point selections and limitations are saved.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from recalibrate import OLD, save


def log_radware(energy, pars):
    a,b,c,d,e,f,g = np.asarray(pars, dtype=float)
    x = np.log(np.asarray(energy, dtype=float)/100.)
    y = np.log(np.asarray(energy, dtype=float)/1000.)
    left = a+b*x+c*x*x
    right = d+e*y+f*y*y
    # Clipping is used only inside rejected-domain optimizer trials.
    value = np.exp(-np.logaddexp(-g*np.log(np.maximum(left, 1e-6)),
                                -g*np.log(np.maximum(right, 1e-6)))/g)
    return value


def full_parameters(p):
    return np.array([p[0],p[1],0.,p[2],p[3],p[4],15.])


def branch_minima(pars):
    a,b,c,d,e,f,g = pars
    x = np.log(np.array([50.,2048.])/100.)
    y = np.log(np.array([50.,2048.])/1000.)
    locations = [y[0],y[1]]
    if f > 0 and y[0] < -e/(2*f) < y[1]:
        locations.append(-e/(2*f))
    yy = np.asarray(locations)
    return np.array([np.min(a+b*x), np.min(d+e*yy+f*yy*yy)])


def fit_radware(points, group, normalization='free', minimum_energy=80., starts=8,
                max_nfev=2500):
    rows = points.loc[(points.energy >= minimum_energy) & points.complete_channel_set
                     & (points.area > 0) & (points.area_se > 0)].copy()
    energy = rows.energy.to_numpy()
    obs = (rows.area/rows.intensity).to_numpy()
    se = obs*np.sqrt((rows.area_se/rows.area)**2+
                     (rows.intensity_se/rows.intensity)**2).to_numpy()
    eu = (rows.isotope=='Eu152').to_numpy().astype(float)
    ratio = float((rows.loc[rows.line=='ba356','area'].iloc[0]/rows.loc[rows.line=='ba356','intensity'].iloc[0]) /
                  (rows.loc[rows.line=='eu344','area'].iloc[0]/rows.loc[rows.line=='eu344','intensity'].iloc[0]))
    free = normalization=='free'
    fixed_log_ratio = -np.log(ratio)
    def residual(p):
        pars=full_parameters(p)
        mu=np.exp(np.clip(log_radware(energy,pars)+eu*(p[5] if free else fixed_log_ratio),-50,50))
        domain=np.minimum(branch_minima(pars)-1e-3,0)*1e5
        return np.r_[(mu-obs)/se,domain]
    def jacobian(p):
        pars=full_parameters(p)
        x=np.log(energy/100.);y=np.log(energy/1000.)
        left=pars[0]+pars[1]*x
        right=pars[3]+pars[4]*y+pars[5]*y*y
        lc=np.maximum(left,1e-6);rc=np.maximum(right,1e-6)
        total=np.logaddexp(-15*np.log(lc),-15*np.log(rc))
        m=np.exp(-total/15)
        dl=m*np.exp(-15*np.log(lc)-total)/lc*(left>1e-6)
        dr=m*np.exp(-15*np.log(rc)-total)/rc*(right>1e-6)
        mu=np.exp(np.clip(m+eu*(p[5] if free else fixed_log_ratio),-50,50))
        columns=[dl,dl*x,dr,dr*y,dr*y*y]+([eu] if free else [])
        jac=np.column_stack(columns)*(mu/se)[:,None]
        domain=np.zeros((2,len(p)))
        if branch_minima(pars)[0]<1e-3:
            xx=np.log((50. if pars[1]>=0 else 2048.)/100.)
            domain[0,:2]=[1e5,xx*1e5]
        if branch_minima(pars)[1]<1e-3:
            choices=[np.log(50./1000),np.log(2048./1000)]
            if pars[5]>0 and choices[0]<-pars[4]/(2*pars[5])<choices[1]:
                choices.append(-pars[4]/(2*pars[5]))
            yy=min(choices,key=lambda yy:pars[3]+pars[4]*yy+pars[5]*yy*yy)
            domain[1,2:5]=np.array([1.,yy,yy*yy])*1e5
        return np.vstack([jac,domain])
    lower=np.array([1.,-10.,1.,-15.,-8.]+([-3.] if free else []))
    upper=np.array([45.,20.,45.,15.,8.]+([3.] if free else []))
    seeds=[OLD[group], [12.143,-.614585,0.,20.6668,5.89342,.957826,15.],
           *[OLD[g] for g in OLD], [10.,2.,0.,10.,-.5,.1,15.]]
    rng=np.random.default_rng(20261102)
    attempts=[]
    for i in range(starts):
        source=np.asarray(seeds[i%len(seeds)],dtype=float)
        p=source[[0,1,3,4,5]].copy()
        scale=np.log(float(np.median(obs[eu==0]/np.exp(log_radware(energy[eu==0],source)))))
        p[[0,2]]+=scale
        if i>=len(seeds): p+=rng.normal(0,[.5,.4,.5,.3,.1])
        if free: p=np.r_[p,fixed_log_ratio]
        p=np.clip(p,lower+1e-6,upper-1e-6)
        result=least_squares(residual,p,jac=jacobian,bounds=(lower,upper),x_scale='jac',max_nfev=max_nfev,
                             ftol=1e-11,xtol=1e-11,gtol=1e-8)
        jac=result.jac[:len(rows)]
        _,singular,vt=np.linalg.svd(jac,full_matrices=False)
        threshold=singular[0]*1e-10
        rank=int(np.sum(singular>threshold))
        inverse=np.zeros_like(singular)
        inverse[singular>threshold]=1/singular[singular>threshold]**2
        covariance=(vt.T*inverse)@vt
        touched=np.flatnonzero((result.x-lower<1e-4)|(upper-result.x<1e-4)).tolist()
        attempts.append(dict(start=i,initial=p.tolist(),fit=result.x.tolist(),parameters=full_parameters(result.x).tolist(),
            log_eu_over_ba=float(result.x[5] if free else fixed_log_ratio),chi2=float(np.sum(result.fun[:len(rows)]**2)),
            domain_valid=bool(np.all(branch_minima(full_parameters(result.x))>0)),
            success=bool(result.success),status=int(result.status),message=result.message,
            nfev=result.nfev,optimality=float(result.optimality),active_bounds=touched,
            branch_minima=branch_minima(full_parameters(result.x)).tolist(),
            near_domain_floor=bool(np.any(branch_minima(full_parameters(result.x))<.00101)),
            jacobian_rank=rank,singular_values=singular.tolist(),covariance=covariance.tolist(),
            conditional_pull=(-result.fun[:len(rows)]).tolist()))
    eligible=[a for a in attempts if a['success'] and a['domain_valid']]
    if not eligible: raise RuntimeError('No converged, valid Radware start')
    best=copy.deepcopy(min(eligible,key=lambda a:a['chi2']))
    # A pseudoinverse cannot be interpreted as an ordinary uncertainty when rank deficient.
    best.update(group=group,normalization=normalization,used_lines=rows.line.tolist(),
                npoints=len(rows),ndof=len(rows)-(6 if free else 5),fitted_parameter_names=['A','B','D','E','F']+(['ln_Eu_over_Ba'] if free else []),
                fixed_parameters={'C':0.,'G':15.},lower_bounds=lower.tolist(),upper_bounds=upper.tolist(),
                legacy_ba_eu_factor=ratio,min_energy=minimum_energy,
                fitted_range_keV=[float(energy.min()),float(energy.max())],
                local_covariance_full_rank=best['jacobian_rank']==(6 if free else 5),
                point_error_scope='Peak conditional area plus listed photon probability errors; probabilities treated independently',
                normalization_error_scope='Joint fitted relative source scale' if free else '356/344 factor held fixed; its shared uncertainty is not propagated in this diagnostic',
                attempts=attempts,formal_replacement=False)
    return best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mode',default='clover2_ab100')
    parser.add_argument('--starts',type=int,default=8)
    parser.add_argument('--minimum-energies',type=float,nargs='+',default=[80.,120.,240.])
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    frame=pd.read_csv(args.points);selected=frame.loc[frame['mode']==args.mode]
    started=time.monotonic();models={}
    for minimum in args.minimum_energies:
        for normalization in ['free','legacy_match']:
            for group in ['all','ring3','ring4','ring5']:
                name=f'{group}_{normalization}_emin{minimum:g}'
                model=fit_radware(selected.loc[selected.group==group],group,normalization,minimum,args.starts)
                models[name]=model
                save(args.output/'radware_fits.json',models)
                print(f'{name}: chi2/ndof={model["chi2"]:.3f}/{model["ndof"]}, rank={model["jacobian_rank"]}, bounds={model["active_bounds"]}',flush=True)
    save(args.output/'completion.json',dict(input=str(args.points.resolve()),
         input_sha256=hashlib.sha256(args.points.read_bytes()).hexdigest(),
         code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         mode=args.mode,starts=args.starts,elapsed_seconds=time.monotonic()-started,
         curves=len(models),formal_replacement=False))


if __name__=='__main__':
    main()
