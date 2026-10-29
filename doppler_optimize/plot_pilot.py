"""Reproduce the Doppler pilot figures from preserved pair files and fit receipts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot
from uproot.source.file import MemmapSource
from evaluate import core_fit

PERIODS={'high1':('High intensity I',530),'low':('Low intensity',647),'high2':('High intensity II',700)}
COLORS={'raw':'0.6','baseline':'#2366ad','candidate':'#e57920'}

class PilotPlots:
    def __init__(self,task):
        self.task=Path(task);self.data={};self.narrow={}
        self.decisions=json.loads((self.task/'decisions.json').read_text())
        self.metrics=pd.read_csv(self.task/'multi_peak_fits.csv')
        self.ranking=pd.read_csv(self.task/'scan_ranking.csv')
        fields=['ge','reference','random_side','weight','e_raw','e_dc_r','particle_MeV','expected_particle_MeV','ring','sector','theta_ge_deg','phi_ge_deg','beta','cos_psi','exclusive_max_energy']
        for period,(_,run) in PERIODS.items():
            states={}
            for state,folder in [('baseline','pilot'),('candidate','candidate')]:
                path=self.task/f'{folder}/data{run:04d}_doppler_wide.root'
                with uproot.open(path,handler=MemmapSource,decompression_executor=uproot.TrivialExecutor(),interpretation_executor=uproot.TrivialExecutor()) as f:
                    states[state]=f['pairs'].arrays(fields,library='np')
            self.data[period]=states
            with uproot.open(self.task/f'candidate/data{run:04d}_doppler_narrow.root') as f:
                self.narrow[period]={family:(f[f'{family}_corrected_net'].axis().edges(),f[f'{family}_corrected_net'].values(),f[f'{family}_corrected_net'].variances()) for family in ['Si','S3']}
        self.ge=sorted(set(np.concatenate([s['baseline']['ge'] for s in self.data.values()])))

    def hist(self,period,reference,state,edges,ge=None,exclusive=False):
        a=self.data[period]['baseline' if state=='raw' else state]
        mask=a['reference']==reference
        if ge is not None:mask &= a['ge']==ge
        if exclusive:mask &= a['exclusive_max_energy']
        energy=a['e_raw'] if state=='raw' else a['e_dc_r']
        return np.histogram(energy[mask],edges,weights=a['weight'][mask])[0],np.histogram(energy[mask],edges,weights=a['weight'][mask]**2)[0]

    def spectrum_grid(self,ge=None,limits=(300,800)):
        fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True)
        edges=np.arange(limits[0],limits[1]+1,1.)
        for i,(period,(name,run)) in enumerate(PERIODS.items()):
            for ref in [0,1]:
                ax=axes[i,ref]
                minimum,maximum=0.,0.
                for state in ['raw','baseline','candidate']:
                    y,v=self.hist(period,ref,state,edges,ge)
                    minimum=min(minimum,float(y.min()));maximum=max(maximum,float(y.max()))
                    ax.stairs(y,edges,baseline=None,color=COLORS[state],lw=1.1,ls='--' if state=='candidate' else '-',label=state)
                for e in [334,406.5,439.4,737]:
                    if limits[0]<e<limits[1]:ax.axvline(e,color='0.2',alpha=.3,lw=.7)
                ax.set_yscale('symlog',linthresh=2);ax.set_title(f'({chr(97+2*i+ref)}) {name}, run {run:04d} | '+('S3' if ref else 'Si array'))
                ax.set_ylim(min(-2,minimum*2),max(5,maximum*2))
                ax.set_ylabel('Net pairs / 1 keV');ax.grid(alpha=.15);ax.set_xlim(*limits)
        axes[0,0].legend(fontsize=8);axes[-1,0].set_xlabel('Gamma energy (keV)');axes[-1,1].set_xlabel('Gamma energy (keV)')
        fig.suptitle(('All Ge' if ge is None else f'Ge{ge}')+' | raw, baseline, tested candidate | held-out runs')
        fig.tight_layout(rect=(0,0,1,.96));return fig

    def scan_curves(self):
        rows=json.loads((self.task/'core_fits.json').read_text())+json.loads((self.task/'combined_core_fits.json').read_text())
        fig,axes=plt.subplots(1,3,figsize=(13,3.8))
        for ax,(period,(name,run)) in zip(axes,PERIODS.items()):
            for split,marker in [('train','o'),('validation','s')]:
                r=[x for x in rows if x['split']==split and x['period']==period and x['family']=='S3' and (x['label'].startswith('S3_energy_scale_') or x['label'].startswith('bracket_S3_energy_scale_') or x['variant']==0)]
                values={}
                for x in r:
                    scale=.5 if x['variant']==0 else float(x['label'].split('_')[-1]);values[scale]=x['score']
                xs=sorted(values);ys=np.array([values[x] for x in xs])/values[.5]
                ax.plot(xs,ys,marker=marker,ms=4,label='Parameter-selection run' if split=='train' else 'Held-out run')
            ax.axhline(1,color='0.5',ls=':');ax.axvline(.5,color='0.5',ls=':');ax.set_title(name);ax.set_xlabel('Recoil exit-energy scale');ax.grid(alpha=.2)
        axes[0].set_ylabel('334-keV score / baseline score');axes[0].legend(fontsize=8)
        fig.suptitle('S3 beta scan: width and centroid jointly scored');fig.tight_layout(rect=(0,0,1,.94));return fig

    def energy_maps(self,period,by_ring=False):
        name,run=PERIODS[period];fig,axes=plt.subplots(2,3,figsize=(14,7),sharey=True)
        fig.subplots_adjust(top=.9,bottom=.13,hspace=.44,wspace=.28,left=.07,right=.86)
        for ref in [0,1]:
            data=[];a=self.data[period]['baseline'];s=(a['reference']==ref)&(a['random_side']==0)
            x=a['ring'] if by_ring else np.searchsorted(self.ge,a['ge'])
            xb=np.arange(.5,(24 if ref else 8)+1.5) if by_ring else np.arange(-.5,len(self.ge)+.5)
            yb=np.arange(310.,356.,1.)
            for state in ['raw','baseline','candidate']:
                e=a['e_raw'] if state=='raw' else self.data[period][state]['e_dc_r']
                h,_,_=np.histogram2d(x[s],e[s],bins=(xb,yb));data.append(h)
            vmax=max(2,max(h.max() for h in data))
            for j,(state,h) in enumerate(zip(['raw','baseline','candidate'],data)):
                ax=axes[ref,j];im=ax.pcolormesh(xb,yb,np.ma.masked_less_equal(h.T,0),norm=LogNorm(1,vmax),cmap='viridis',rasterized=True)
                ax.axhline(334,color='red',lw=.8);ax.set_title(f'({chr(97+ref*3+j)}) '+('S3' if ref else 'Si array')+' | '+state)
                if not by_ring:ax.set_xticks(range(len(self.ge)),self.ge,rotation=90,fontsize=7)
                ax.set_xlabel('Particle ring' if by_ring else 'Ge channel');ax.set_ylabel('Gamma energy (keV)')
            position=axes[ref,2].get_position()
            color_axis=fig.add_axes([.89,position.y0,.012,position.height])
            fig.colorbar(im,cax=color_axis,label='Prompt pairs / bin')
        fig.suptitle(f'{name}, run {run:04d} | same prompt pairs, 1-keV energy bins')
        return fig

    def particle_energy(self):
        fig,axes=plt.subplots(3,2,figsize=(12,10))
        for i,(period,(name,run)) in enumerate(PERIODS.items()):
            a=self.data[period]['baseline']
            for ref in [0,1]:
                ax=axes[i,ref];s=(a['reference']==ref)&(a['random_side']==0)&(abs(a['e_dc_r']-334)<6)
                xb=np.arange(.5,(24 if ref else 8)+1.5);yb=np.arange(0,101,1.)
                h,_,_=np.histogram2d(a['ring'][s],a['particle_MeV'][s],bins=(xb,yb))
                ax.pcolormesh(xb,yb,np.ma.masked_less_equal(h.T,0),norm=LogNorm(1,max(2,h.max())),cmap='viridis')
                rings=sorted(set(a['ring'][s]));med=[np.median(a['particle_MeV'][s&(a['ring']==r)]) for r in rings]
                expected=[np.median(a['expected_particle_MeV'][s&(a['ring']==r)]) for r in rings]
                ax.plot(rings,expected,'r.-',label='Old exit-energy table');ax.plot(rings,med,'w.--',label='Measured median')
                ax.set_title(f'{name}, {run:04d} | '+('S3' if ref else 'Si array'));ax.set_xlabel('Particle ring');ax.set_ylabel('Measured energy (MeV)')
        axes[0,0].legend(fontsize=8);fig.suptitle('Particle-energy check | prompt, baseline gamma gate 334 +/- 6 keV')
        fig.tight_layout(rect=(0,0,1,.96));return fig

    def prompt_random(self):
        fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True)
        edges=np.arange(40,2001,1.)
        for i,(period,(name,run)) in enumerate(PERIODS.items()):
            a=self.data[period]['baseline']
            for ref in [0,1]:
                ax=axes[i,ref];s=a['reference']==ref
                p=np.histogram(a['e_dc_r'][s&(a['random_side']==0)],edges)[0]
                r=np.histogram(a['e_dc_r'][s&(a['random_side']!=0)],edges,weights=-a['weight'][s&(a['random_side']!=0)])[0]
                ax.stairs(p,edges,color='0.65',lw=.6,label='Prompt');ax.stairs(r,edges,color='#884090',lw=.6,label='Scaled random');ax.stairs(p-r,edges,color='#238850',lw=.7,label='Net')
                ax.set_yscale('symlog',linthresh=2);ax.set_title(f'{name}, {run:04d} | '+('S3' if ref else 'Si array'));ax.set_ylabel('Pairs / 1 keV');ax.grid(alpha=.15)
        axes[0,0].legend(fontsize=8);axes[-1,0].set_xlabel('Recoil-corrected gamma energy (keV)');axes[-1,1].set_xlabel('Recoil-corrected gamma energy (keV)')
        fig.suptitle('Baseline correction | per-pair random normalization');fig.tight_layout(rect=(0,0,1,.96));return fig

    def window_comparison(self):
        fig,axes=plt.subplots(3,2,figsize=(12,8),sharex=True)
        for i,(period,(name,run)) in enumerate(PERIODS.items()):
            for ref,family in enumerate(['Si','S3']):
                ax=axes[i,ref];edges,y,v=self.narrow[period][family]
                yw,vw=self.hist(period,ref,'candidate',edges)
                ax.stairs(y,edges,label='Narrow',color='#2366ad');ax.stairs(yw,edges,label='Wide',color='#e57920')
                ax.set_xlim(315,353);ax.set_yscale('symlog',linthresh=2);ax.set_title(f'{name}, {run:04d} | {family}');ax.set_ylabel('Net pairs / 1 keV');ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=8);axes[-1,0].set_xlabel('Candidate-corrected energy (keV)');axes[-1,1].set_xlabel('Candidate-corrected energy (keV)')
        fig.suptitle('Window cross-check | same candidate parameters, no area renormalization');fig.tight_layout(rect=(0,0,1,.95));return fig

    def ge_metrics(self):
        result=[];edges=np.arange(318.,350.5,.5)
        for period in PERIODS:
            for ge in self.ge:
                for ref,family in enumerate(['Si','S3']):
                    for state in ['baseline','candidate']:
                        y,v=self.hist(period,ref,state,edges,ge)
                        fit=core_fit(edges,y,v)
                        result.append(dict(period=period,ge=int(ge),family=family,state=state,**fit))
        return pd.DataFrame(result)

    def parameter_tables(self):
        out=[];ge=[]
        for period in PERIODS:
            b=self.data[period]['baseline'];a=self.data[period]['candidate']
            for ref in [0,1]:
                for r in sorted(set(a['ring'][a['reference']==ref])):
                    s=(a['reference']==ref)&(a['ring']==r)
                    out.append(dict(period=period,family='S3' if ref else 'Si',ring=int(r),beta_baseline=float(np.median(b['beta'][s])),beta_candidate=float(np.median(a['beta'][s]))))
            for key in self.ge:
                s=a['ge']==key
                if s.any():ge.append(dict(period=period,Ge=int(key),theta_deg=float(a['theta_ge_deg'][s][0]),phi_deg=float(a['phi_ge_deg'][s][0]),theta_change_deg=0.,phi_change_deg=0.))
        return pd.DataFrame(out),pd.DataFrame(ge)

    def fit_diagnostics(self):
        from scipy.special import ndtr
        fig,axes=plt.subplots(2,3,figsize=(14,7),gridspec_kw={'height_ratios':[3,1]})
        cases=[('low',1,334.),('low',1,439.4),('high2',0,334.)]
        for col,(period,reference,peak) in enumerate(cases):
            width=16 if peak==334 else 12
            full_edges=np.arange(280.,800.5,.5)
            centers=(full_edges[1:]+full_edges[:-1])/2
            indices=np.flatnonzero((centers>=peak-width)&(centers<peak+width))
            edges=full_edges[indices[0]:indices[-1]+2];x=(edges[1:]+edges[:-1])/2
            for state in ['baseline','candidate']:
                y,v=self.hist(period,reference,state,edges)
                f=core_fit(edges,y,v,peak,width)
                if 'area' not in f:continue
                model=f['area']*(ndtr((edges[1:]-f['mu'])/f['sigma'])-ndtr((edges[:-1]-f['mu'])/f['sigma']))+f['continuum']*np.diff(edges)+f['slope']*((edges[1:]-peak)**2-(edges[:-1]-peak)**2)/2
                axes[0,col].stairs(y,edges,baseline=None,color=COLORS[state],lw=.7,alpha=.7)
                axes[0,col].plot(x,model,color=COLORS[state],ls='--',label=f"{state}: chi2/ndf={f['chi2_per_dof']:.1f}")
                axes[1,col].plot(x,(y-model)/np.sqrt(np.maximum(v,1)),color=COLORS[state],lw=.8)
            axes[0,col].set_title(f'{PERIODS[period][0]}, {PERIODS[period][1]:04d}\n'+('S3' if reference else 'Si array')+f', {peak:g} keV')
            axes[0,col].legend(fontsize=8);axes[0,col].set_ylabel('Net pairs / 0.5 keV')
            axes[1,col].axhline(0,color='0.3',lw=.7);axes[1,col].set_ylabel('(data-model)/error');axes[1,col].set_xlabel('Corrected gamma energy (keV)')
        fig.suptitle('Core-fit diagnostics | narrowing, centroid shifts, and a failed simple peak model')
        fig.tight_layout(rect=(0,0,1,.94));return fig
