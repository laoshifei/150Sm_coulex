"""Numerical utilities for fixed-pixel spectra and frozen-candidate checks."""
from pathlib import Path
import json
import numpy as np
import uproot
ENERGY_GRID=np.arange(20000.,90000.1,1000.)
def spectrum_file(path):
 with uproot.open(path,handler=uproot.source.file.MemmapSource) as f:
  meta=json.loads(str(f['metadata']));values={}
  for face in ['ring','sector']:
   for region in ['prompt','random']:
    # Array axes are energy,pixel,fold. Keep each actual fold separate.
    h=f[f'{face}_{region}'].values().transpose(2,1,0)
    values[(face,region)]=np.concatenate([np.zeros((*h.shape[:2],1)),np.cumsum(h,axis=2)],axis=2)
 return meta,values
def interpolate_rows(cdf,raw):
 raw=np.broadcast_to(raw,(len(cdf),raw.shape[-1]));scaled=np.clip(raw/250.,0,800)
 i=np.minimum(scaled.astype(int),799);fraction=scaled-i
 rows=np.arange(len(cdf))[:,None]
 return cdf[rows,i]*(1-fraction)+cdf[rows,i+1]*fraction
def normalized_net(values,face,fold,pixels,a,b,grid=ENERGY_GRID):
 a=np.broadcast_to(np.asarray(a),len(pixels));b=np.broadcast_to(np.asarray(b),len(pixels))
 raw=(grid[None,:]-a[:,None])/b[:,None]
 p=interpolate_rows(values[(face,'prompt')][fold,pixels],raw)
 r=interpolate_rows(values[(face,'random')][fold,pixels],raw)
 net=p-.5*r;total=net[:,-1]-net[:,0]
 valid=(total>=15)&(raw[:,0]>=0)&(raw[:,-1]<=200000)&(b>0)
 normalized=(net-net[:,0,None])/np.where(total>0,total,1)[:,None]
 return normalized,total,valid,p[:,-1]-p[:,0],r[:,-1]-r[:,0]
def maps_from_channels(rows):
 a=np.zeros(56);b=np.ones(56)
 for row in rows:a[row['index']]=row['A_keV'];b[row['index']]=row['B']
 return {'a':a,'b':b,'radial':np.ones(24)}
def model_maps(base,conditional,coupled):
 result={'0228':maps_from_channels(base['channels'])}
 result['formal']={'a':np.array([c['old_A_keV'] for c in base['channels']]),'b':np.array([c['old_B'] for c in base['channels']]),'radial':np.ones(24)}
 result['direct']={k:v.copy() for k,v in result['0228'].items()}
 for c in conditional['channels']:
  if c['training_supported']:
   result['direct']['a'][c['index']]=c['direct_A_keV'];result['direct']['b'][c['index']]=c['direct_B']
 for name,model in coupled['models'].items():
  if model['ok']:result[name]=maps_from_channels(model['coefficients'])
 result['radial']={k:v.copy() for k,v in result['0228'].items()}
 for c in conditional['channels'][:24]:
  if c['training_supported']:result['radial']['radial'][c['index']]=c['q']
 return result
def spectral_metrics(values,reference,refmaps,maps,conditional,active,fold):
 rows=[]
 for c in conditional['channels']:
  index=c['index'];face='ring' if index<24 else 'sector';pixels=np.array(c['eligible_pixels'],int)
  if len(pixels)<3:continue
  ref,_,rv,_,_=normalized_net(reference,face,fold,pixels,refmaps['a'][index],refmaps['b'][index])
  base,_,bv,_,_=normalized_net(values,face,fold,pixels,maps['0228']['a'][index],maps['0228']['b'][index])
  use=rv&bv
  if use.sum()<3:continue
  pixels=pixels[use];ref=ref[use]
  for name,model in maps.items():
   scale=model['radial'][pixels//32]
   cur,total,valid,p,r=normalized_net(values,face,fold,pixels,scale*model['a'][index],scale*model['b'][index])
   d=np.sqrt(np.mean((cur[:,1:-1]-ref[:,1:-1])**2,axis=1));d[~valid]=1.
   rows.append({'channel_index':index,'model':name,'fold':fold,'production':index<active or index>=24,'pixels':pixels.tolist(),
    'D':float(np.sqrt(np.mean(d*d))),'pixel_D':d.tolist(),'net':total.tolist(),'prompt':p.tolist(),'random':r.tolist(),'invalid_pixels':int((~valid).sum())})
 return rows
def metrics(x):
 if len(x)<20:return {'n':len(x),'median_keV':None,'width68_keV':None,'score_keV':None,'tail2MeV':None}
 q=np.quantile(x,[.16,.5,.84]);w=(q[2]-q[0])/2
 return {'n':len(x),'median_keV':float(q[1]),'width68_keV':float(w),'score_keV':float(np.hypot(q[1],w)),'tail2MeV':float(np.mean(np.abs(x)>2000))}
def load_pairs(path):
 branches=['block','r','s','nr','ns','flags','dt_ns','r_e0_keV','s_e0_keV','r_old_keV','s_old_keV']
 with uproot.open(path,handler=uproot.source.file.MemmapSource) as f:return f['pairs'].arrays(branches,library='np'),json.loads(str(f['metadata']))
def pair_masks(v,active,fold):
 r=v['r'];x=v['r_old_keV'];y=v['s_old_keV'];dt=v['dt_ns']
 base=(v['block']%4==fold)&(v['flags']==0)&(v['nr']==1)&(v['ns']==1)&(x>=15000)&(x<100000)&(y>=15000)&(y<100000)
 prompt=(dt>=-80)&(dt<80);random=((dt>=-880)&(dt< -720))|((dt>=720)&(dt<880))
 return base,prompt,random,r<active
def differences(v,model,baseline):
 r=v['r'];s=v['s'];q=model['radial'][r]
 er=q*(model['a'][r]+model['b'][r]*v['r_e0_keV']);es=q*(model['a'][s]+model['b'][s]*v['s_e0_keV'])
 # Keep a common energy unit while testing front/back matching. A common gain
 # cannot count as an intrinsic improvement in the difference width.
 scale=q*.5*(model['b'][r]/baseline['b'][r]+model['b'][s]/baseline['b'][s])
 return es-er,(es-er)/scale,er,es
