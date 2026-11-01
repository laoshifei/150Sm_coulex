"""Bin-integrated source peak and relative-efficiency fits, without import side effects.

The C++ raw2hist executable supplies complete decoded-source ADC counts. Fits
retain the original bin counts and map ADC bin edges instead of rebinning data.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import os
import time
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import least_squares
from scipy.special import ndtr, log_ndtr, xlogy
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).resolve().parent
GROUPS = {"ring3": [204,205,206,212,213,214,215],
          "ring4": [211,304,305,306,308,507,508,510],
          "ring5": [209,310,311,506,509]}
KEYS = sorted(sum(GROUPS.values(), []))
OLD = {"all":[14.8,4.618,0,10.8,-0.7352,-0.1101,15],
       "ring3":[11.23,2.49,0,9.432,-0.6368,-0.005044,15],
       "ring4":[11.07,1.767,0,9.663,-0.5738,0.0251,15],
       "ring5":[11.06,.2759,0,9.963,-.7139,.1796,15]}


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    os.replace(tmp, path)


def read_pars(path):
    a = np.loadtxt(path, skiprows=1)
    return {int(r[0])*10000+int(r[1])*100+int(r[2]):r[3:6] for r in a}


def energy(x, p):
    return p[0]+p[1]*x+p[2]*x*x


def left_cdf(x, mu, sigma, tau):
    """CDF of a Gaussian minus an independent positive exponential."""
    z = (x-mu)/sigma
    return np.clip(ndtr(z)+np.exp((x-mu)/tau+.5*(sigma/tau)**2+log_ndtr(-z-sigma/tau)), 0., 1.)


def peak_components(edges, p, reference, neighbors, tail):
    # Parameters: ln(area), center, ln(sigma), ln(bg_left), ln(bg_right),
    # ln(step density), optional tail fraction/ln(tau), then ln(neighbor areas).
    A, mu, sig = np.exp(p[0]), p[1], np.exp(p[2])
    f, tau, offset = (p[6], np.exp(p[7]), 8) if tail else (0., 1., 6)
    both = tail == 'both'
    if both: offset = 10
    def shape(center):
        c = (1-f)*ndtr((edges-center)/sig)
        if f:
            fraction_left = p[8] if both else 1.
            c += f*fraction_left*left_cdf(edges, center, sig, tau)
            if both:
                c += f*(1-fraction_left)*(1-left_cdf(2*center-edges,center,sig,np.exp(p[9])))
        return np.maximum(np.diff(c), 0)
    main = A*shape(mu)
    other = np.zeros_like(main)
    for i, e in enumerate(neighbors):
        other += np.exp(p[offset+i])*shape(mu+e-reference)
    width = np.diff(edges)
    mid = (edges[1:]+edges[:-1])/2
    bg = width*(np.exp(p[3])+(np.exp(p[4])-np.exp(p[3]))*(mid-edges[0])/(edges[-1]-edges[0]))
    z = (edges-mu)/sig
    antideriv = (edges-mu)*ndtr(-z)-sig*np.exp(-.5*z*z)/np.sqrt(2*np.pi)
    bg += np.exp(p[5])*np.maximum(np.diff(antideriv), 0)
    return main, other, bg


def poisson_residual(y, model):
    m = np.maximum(model, 1e-200)
    dev = np.maximum(2*(m-y+xlogy(y, y/m)), 0.)
    return np.sign(m-y)*np.sqrt(dev)


def fit_peak(counts, edges, line, tail=True, search=True, window_scale=1., max_nfev=1200):
    reference = line['energy']
    centers = (edges[1:]+edges[:-1])/2
    radius = 2.5 + .006*reference if search else 1.4
    near = np.flatnonzero(np.abs(centers-reference)<radius)
    if len(near)<3:
        raise ValueError('Too few bins near line')
    smooth = gaussian_filter1d(counts.astype(float), .6)
    peakbin = near[np.argmax(smooth[near])]
    guess = float(centers[peakbin])
    shift = guess-reference
    lo = reference + (line['range'][0]-reference)*window_scale + shift
    hi = reference + (line['range'][1]-reference)*window_scale + shift
    select = np.flatnonzero((centers>=lo)&(centers<hi))
    first,last = select[0],select[-1]+1
    x = edges[first:last+1]
    y = counts[first:last].astype(float)
    neighbors = [e for e in line.get('neighbors',[]) if lo<e+shift<hi]
    width = np.diff(x)
    b = max(np.median(np.r_[y[:max(2,len(y)//8)],y[-max(2,len(y)//8):]])/np.median(width), .01)
    a = max(float(y.sum()-b*(x[-1]-x[0])), 10.)
    p = [math.log(a),guess,math.log(.85),math.log(b),math.log(b),math.log(max(b*.1,.01))]
    lower = [math.log(.01),guess-2.,math.log(.20),-18.,-18.,-18.]
    upper = [math.log(max(100.,y.sum()*3)),guess+2.,math.log(4.),math.log(max(100.,y.max()/min(width)*3)),math.log(max(100.,y.max()/min(width)*3)),math.log(max(100.,y.max()/min(width)*3))]
    if tail:
        p += [.04, math.log(1.5)]
        lower += [0.,math.log(.18)]
        upper += [.65,math.log(14.)]
        if tail == 'both':
            p += [.5, math.log(1.5)]
            lower += [0., math.log(.18)]
            upper += [1., math.log(14.)]
    for e in neighbors:
        p.append(math.log(max(2.,a*.005)))
        lower.append(-10.)
        upper.append(math.log(max(100.,y.sum()*3)))
    if len(y)<=len(p)+3:
        raise ValueError('Insufficient bins for peak model')
    def resid(v): return poisson_residual(y,sum(peak_components(x,v,reference,neighbors,tail)))
    result = least_squares(resid,p,bounds=(lower,upper),x_scale='jac',max_nfev=max_nfev,
                           ftol=1e-8,xtol=1e-8,gtol=1e-7)
    main,other,bg = peak_components(x,result.x,reference,neighbors,tail)
    fisher = result.jac.T@result.jac
    covariance = np.linalg.pinv(fisher, rcond=1e-11)
    eig = np.linalg.eigvalsh(fisher)
    area = float(np.exp(result.x[0]));area_se = area*float(np.sqrt(max(covariance[0,0],0)))
    errors = np.sqrt(np.maximum(np.diag(covariance),0))
    boundary = np.flatnonzero((result.x-np.array(lower)<1e-4)|(np.array(upper)-result.x<1e-4)).tolist()
    ndof = len(y)-len(result.x)
    return dict(line=line['id'],energy=reference,isotope=line['isotope'],intensity=line['intensity'],intensity_se=line['uncertainty'],
                area=area,area_se=area_se,center=float(result.x[1]),center_se=float(errors[1]),sigma=float(np.exp(result.x[2])),
                tail_fraction=float(result.x[6]) if tail else 0.,tail_tau=float(np.exp(result.x[7])) if tail else 0.,
                deviance=float(result.fun@result.fun),ndof=ndof,deviance_per_dof=float(result.fun@result.fun)/ndof,
                success=bool(result.success),message=result.message,nfev=result.nfev,bounds_indices=boundary,
                jacobian_rank=int(np.linalg.matrix_rank(result.jac)),parameters=len(result.x),
                information_condition=float(max(eig)/max(min(eig),1e-30)),covariance_method='local deviance-Jacobian information pseudoinverse; unscaled',
                p=result.x.tolist(),lower=lower,upper=upper,covariance=covariance.tolist(),edges=x.tolist(),counts=y.tolist(),
                model=(main+other+bg).tolist(),main=main.tolist(),neighbors=other.tolist(),background=bg.tolist(),
                neighbor_energies=neighbors,tail_model=tail,window_scale=window_scale,residual=(-result.fun).tolist())


def fit_energy(raw, seed, output, lines, energy_model='left'):
    tail_mode = {'gaussian': False, 'left': True, 'both': 'both'}[energy_model]
    fits, rows, parameters = [], [], {}
    with uproot.open(raw) as f:
        for key in KEYS:
            counts, adc = f[f'ch{key:05d}/adc_clean'].to_numpy()
            # A two-channel ADC bin [n,n+2) contains integer ADC n,n+1. Its
            # physical integration boundaries are n-0.5 and n+1.5.
            edges = energy(adc-.5,seed[key])
            if not np.all(np.diff(edges)>0): raise ValueError(f'Nonmonotone seed for {key}')
            local = []
            for line in lines:
                try:
                    r=fit_peak(counts,edges,line,tail=tail_mode,search=True)
                    r['key']=key
                    r['adc_center']=float(np.interp(r['center'],edges,adc-.5))
                    r['used_for_energy'] = bool(r['success'] and r['energy']>=120 and r['area']/max(r['area_se'],1e-20)>20 and r['sigma']<2.8 and 1 not in r['bounds_indices'])
                    # Unresolved 963/964 transitions and weak neighboring components
                    # are retained for area diagnostics but not energy anchoring.
                    if line['id'] in ['eu964','eu411','eu444','ba276']: r['used_for_energy']=False
                    fits.append(r);local.append(r)
                except (ValueError,FloatingPointError) as exc:
                    fits.append(dict(key=key,line=line['id'],energy=line['energy'],error=str(exc),success=False,used_for_energy=False))
            selected=[r for r in local if r['used_for_energy']]
            if len(selected)<5: raise RuntimeError(f'Only {len(selected)} reliable energy anchors for Ge{key}')
            x=np.array([r['adc_center'] for r in selected]); y=np.array([r['energy'] for r in selected])
            # Equal weights prevent high-count but shape-dependent centroids from
            # masquerading as high-precision independent metrological anchors.
            poly=np.polynomial.Polynomial.fit(x,y,2).convert()
            coeff=poly.coef
            residual=poly(x)-y
            maxres=float(max(abs(residual)))
            parameters[str(key)]=dict(coefficients=coeff.tolist(),formula='E_keV=c0+c1*ADC+c2*ADC^2',
                                      anchor_count=len(selected),rms_keV=float(np.sqrt(np.mean(residual**2))),
                                      max_residual_keV=maxres,valid_energy_keV=[float(min(y)),float(max(y))],
                                      status='candidate' if maxres<1 else 'review_required',
                                      weights='equal; diagnostic polynomial, no absolute-energy confidence claim')
            for r in local:
                r['new_center_keV']=float(poly(r['adc_center']))
                r['new_center_residual_keV']=r['new_center_keV']-r['energy']
                rows.append({k:r[k] for k in ['key','line','energy','center','center_se','adc_center','new_center_keV','new_center_residual_keV','used_for_energy','area','area_se','sigma','deviance_per_dof','success']})
            save(output/'energy_peaks.json',fits)
            save(output/'energy_calibration.json',dict(schema='sm150_ge_source_energy_v1',run='0722',channels=parameters,formal_replacement=False))
            pd.DataFrame(rows).to_csv(output/'energy_peaks.csv',index=False)
            print(f'Ge{key}: {len(selected)} anchors, RMS={np.sqrt(np.mean(residual**2)):.3f} keV, max={maxres:.3f} keV',flush=True)
    return parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['energy'])
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--seed',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--energy-model', choices=['gaussian','left','both'], default='left',
                        help='Peak shape for source energy anchors; not the efficiency-area model selection')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    lines=json.loads((HERE/'source_lines.json').read_text())['lines']
    save(args.output/'run.json',dict(stage=args.stage,started_at=datetime.datetime.now().astimezone().isoformat(),input=str(args.input.resolve()),
         input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),seed=str(args.seed.resolve()),
         seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),lines_sha256=hashlib.sha256((HERE/'source_lines.json').read_bytes()).hexdigest(),
         code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),energy_model=args.energy_model))
    fit_energy(args.input,read_pars(args.seed),args.output,lines,args.energy_model)
    save(args.output/'completion.json',dict(status='complete',stage=args.stage,elapsed_seconds=time.monotonic()-start,channels=len(KEYS)))


if __name__=='__main__':
    main()
