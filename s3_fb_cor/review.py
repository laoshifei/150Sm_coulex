"""Matched-event diagnostics with separate distance and counting-support flags."""
from pathlib import Path
import json
import numpy as np
from s3_energy_attenuation.analysis import spectrum_file, interpolate_rows, maps_from_channels, metrics, load_pairs, differences

GRID=np.arange(20000.,90000.1,1000.)

def package_map(row):
    return {'a':np.asarray(row['A_keV'],float),'b':np.asarray(row['B'],float),'radial':np.asarray(row['q_position'],float)}

def candidates(base,previous,fixed):
    result={'0228':maps_from_channels(base['channels']), '0229':package_map(previous)}
    result['formal']={'a':np.array([c['old_A_keV'] for c in base['channels']]),
                      'b':np.array([c['old_B'] for c in base['channels']]),'radial':np.ones(24)}
    for name,m in fixed['models'].items():
        if m['ok']: result[name]=maps_from_channels(m['coefficients'])
    return result

def folded(values,face,region,fold):
    v=values[(face,region)]
    return v.sum(axis=0) if fold<0 else v[fold]

def cumulative_covariance(variance_cdf,x,y):
    """Covariance of two linearly interpolated cumulative histogram counts."""
    x,y=np.broadcast_arrays(np.clip(x/250.,0,800),np.clip(y/250.,0,800))
    low,high=np.minimum(x,y),np.maximum(x,y)
    i=np.minimum(low.astype(int),799);j=np.minimum(high.astype(int),799)
    f=low-i;g=np.where(i==j,high-j,1.)
    rows=np.arange(len(variance_cdf))[:,None]
    return variance_cdf[rows,i]+f*g*(variance_cdf[rows,i+1]-variance_cdf[rows,i])

def cumulative_net(values,face,fold,pixels,a,b):
    pixels=np.asarray(pixels,int);a=np.broadcast_to(np.asarray(a),len(pixels));b=np.broadcast_to(np.asarray(b),len(pixels))
    raw=(GRID[None,:]-a[:,None])/b[:,None]
    pcdf=folded(values,face,'prompt',fold)[pixels];rcdf=folded(values,face,'random',fold)[pixels]
    p=interpolate_rows(pcdf,raw);r=interpolate_rows(rcdf,raw);net=p-.5*r
    total=net[:,-1]-net[:,0];positive=(total>0)&(b>0)&(raw[:,0]>=0)&(raw[:,-1]<=200000)
    f=(net-net[:,0,None])/np.where(total>0,total,1)[:,None]
    v=pcdf+.25*rcdf;lo=raw[:,0,None];hi=raw[:,-1,None]
    vl=cumulative_covariance(v,lo,lo);vh=cumulative_covariance(v,hi,hi);clh=cumulative_covariance(v,lo,hi)
    vx=cumulative_covariance(v,raw,raw);cxl=cumulative_covariance(v,raw,lo);cxh=cumulative_covariance(v,raw,hi)
    vn=vx+vl-2*cxl;vt=vh+vl-2*clh;cnt=cxh-cxl-clh+vl
    vf=np.maximum(0.,(vn+f*f*vt-2*f*cnt)/np.where(total>0,total,1)[:,None]**2)
    return {'F':f,'variance_F':vf,'net':total,'prompt':p[:,-1]-p[:,0],'random':r[:,-1]-r[:,0],
            'positive':positive,'reliable':positive&(total>=15)}

def spectra(values,reference,refmap,maps,conditional,active,fold):
    rows=[]
    for c in conditional['channels']:
        i=c['index'];face='ring' if i<24 else 'sector';eligible=np.asarray(c['eligible_pixels'],int)
        production=i<active or i>=24
        common={'index':i,'channel':c['channel'],'production':production,'training_pixels':len(eligible),'fold':fold}
        if len(eligible):
            ref0=cumulative_net(reference,face,fold,eligible,refmap['a'][i],refmap['b'][i])
            base0=cumulative_net(values,face,fold,eligible,maps['0228']['a'][i],maps['0228']['b'][i])
            use=ref0['reliable']&base0['reliable'];pixels=eligible[use]
        else:pixels=np.array([],int)
        if len(pixels)<3:
            for name in maps:rows.append({**common,'model':name,'D':None,'noise_floor':None,'pixels':pixels.tolist(),
                                        'low_count_pixels':[],'nonpositive_or_outside_pixels':[],'status':'insufficient_common_pixels'})
            continue
        ref=cumulative_net(reference,face,fold,pixels,refmap['a'][i],refmap['b'][i])
        for name,m in maps.items():
            q=m['radial'][pixels//32]
            cur=cumulative_net(values,face,fold,pixels,q*m['a'][i],q*m['b'][i])
            valid=cur['positive'];low=valid&~cur['reliable']
            delta=cur['F'][:,1:-1]-ref['F'][:,1:-1]
            pd=np.sqrt(np.mean(delta**2,axis=1));pn=np.sqrt(np.mean(cur['variance_F'][:,1:-1]+ref['variance_F'][:,1:-1],axis=1))
            d=float(np.sqrt(np.mean(pd**2))) if np.all(valid) else None
            noise=float(np.sqrt(np.mean(pn**2))) if np.all(valid) else None
            rows.append({**common,'model':name,'D':d,'noise_floor':noise,'excess_distance':float(np.sqrt(max(0,d*d-noise*noise))) if d is not None else None,
                         'pixels':pixels.tolist(),'pixel_D':[float(v) if ok else None for v,ok in zip(pd,valid)],'pixel_noise':pn.tolist(),
                         'net':cur['net'].tolist(),'prompt':cur['prompt'].tolist(),'random':cur['random'].tolist(),
                         'low_count_pixels':pixels[low].tolist(),'nonpositive_or_outside_pixels':pixels[~valid].tolist(),
                         'status':'invalid_model_support' if not np.all(valid) else ('low_count_model_support' if np.any(low) else 'comparable')})
    return rows

def masks(v,active,fold):
    oldr=v['r_old_keV'];olds=v['s_old_keV'];dt=v['dt_ns']
    quality=(v['flags']==0)&(v['nr']==1)&(v['ns']==1)&(oldr>=15000)&(oldr<100000)&(olds>=15000)&(olds<100000)
    if fold>=0:quality&=(v['block']%4==fold)
    prompt=quality&(dt>=-80)&(dt<80)
    random=quality&(((dt>=-880)&(dt< -720))|((dt>=720)&(dt<880)))
    return prompt,random,v['r']<active

def trend(x,y):
    points=[]
    for lo in np.arange(15000.,100000.,5000.):
        use=(x>=lo)&(x<lo+5000)
        if use.sum()>=30:points.append({'energy_keV':lo+2500,'n':int(use.sum()),'median_keV':float(np.median(y[use]))})
    return {'points':points,'span_keV':float(np.ptp([p['median_keV'] for p in points])) if len(points)>=4 else None}

def pair_metrics(v,maps,active,fold):
    p,r,production=masks(v,active,fold);global_rows={};rows=[]
    for name,m in maps.items():
        raw,d,_,_=differences(v,m,maps['0228'])
        global_rows[name]={'normalized':metrics(d[p&production]),'raw':metrics(raw[p&production]),'random_pairs':int((r&production).sum())}
        for i in range(56):
            use=p&((v['r']==i) if i<24 else ((v['s']==i)&production))
            rows.append({'model':name,'index':i,'production':i<active or i>=24,'fold':fold,**metrics(d[use]),'energy_trend':trend(v['r_old_keV'][use],d[use])})
    return global_rows,rows

def distance_score(spectral,name):
    values=[r['D'] for r in spectral if r['model']==name and r['production'] and r['D'] is not None]
    return float(np.sqrt(np.mean(np.square(values)))) if values else None

def guard(name,spectral,pairs,channels,confirmation=False):
    flags=[];base=pairs['0228']['normalized'];new=pairs[name]['normalized']
    factor=1.05 if confirmation else 1.03
    if new['n']<100:return ['insufficient_global_pairs']
    for key in ['width68_keV','score_keV']:
        if new[key]>factor*base[key]+1:flags.append('global_'+key+'_worse')
    if new['tail2MeV']>base['tail2MeV']+.01:flags.append('global_tail_increase')
    old={r['index']:r for r in channels if r['model']=='0228'}
    for r in channels:
        if r['model']==name and r['production'] and r['n']>=100 and r['score_keV']>1.10*old[r['index']]['score_keV']+20:
            flags.append(f'front_back_worse:{r["index"]}')
    old={r['index']:r for r in spectral if r['model']=='0228'}
    for r in spectral:
        if r['model']!=name or not r['production']:continue
        if r['status']!='comparable':flags.append(f'spectral_support:{r["index"]}:{r["status"]}')
        elif old[r['index']]['D'] is not None and r['D']>old[r['index']]['D']+.02:flags.append(f'spectral_worse:{r["index"]}')
    d0=distance_score(spectral,'0228');d1=distance_score(spectral,name)
    if d0 is None or d1 is None:flags.append('insufficient_global_spectra')
    elif confirmation:
        if d1>d0+.002:flags.append('global_spectral_worse')
    elif d1>.95*d0 or d0-d1<.001:flags.append('global_spectral_improvement_below_threshold')
    return flags

def model_json(m):return {k:np.asarray(v).tolist() for k,v in m.items()}
