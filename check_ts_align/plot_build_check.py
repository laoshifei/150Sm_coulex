"""Plot the C++ timing-check histograms; never read raw events or refit timing."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

PERIODS = ['high1', 'low', 'high2']
PERIOD_LABELS = {'high1': 'High beam I', 'low': 'Low beam', 'high2': 'High beam II'}
CONTROL_NAMES = ['s3_ring_minus_sector', 'si_minus_s3', 'si_low_channel_minus_high_channel']
CONTROL_LABELS = ['S3 ring - sector', 'Si array - S3 ring', 'Si: lower channel - higher channel']

class TimingCheck:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.hist = {}; self.meta = {}
        for period in PERIODS:
            for candidate in ['narrow', 'wide']:
                key = (period, candidate)
                source = self.directory / f'{period}_{candidate}.root'
                receipt = json.loads((self.directory / f'{period}_{candidate}_receipt.json').read_text())
                assert hashlib.sha256(source.read_bytes()).hexdigest() == receipt['output_sha256'], source
                with uproot.open(source) as f:
                    self.meta[key] = json.loads(str(f['sm150_timing_check_manifest']))
                    assert self.meta[key]['common']['group'] == period
                    assert self.meta[key]['common']['candidate'] == ('A' if candidate == 'narrow' else 'B')
                    assert self.meta[key]['common']['time_alignment_profile'] == 'joint_zero_20261028'
                    self.hist[key] = {name: f[name].to_numpy() for name in f.keys(cycle=False) if name != 'sm150_timing_check_manifest'}
        self.ges = sorted({int(g['ge']) for m in self.meta.values() for g in m['common']['windows']})
        self.edges = self.hist['low', 'wide'][f'ge{self.ges[0]}_si_prompt'][1]
        self.centers = (self.edges[1:] + self.edges[:-1])/2

    def h1(self, period, candidate, name):
        return self.hist[period, candidate].get(name, (np.zeros(len(self.centers)), self.edges))[0]

    def gate(self, period, candidate, ge, reference):
        return next(g for g in self.meta[period, candidate]['common']['windows'] if g['ge'] == ge and g['reference'] == reference)

    def mode(self, y):
        if not np.sum(y): return np.nan
        return self.centers[np.argmax(np.convolve(y, [1,2,3,2,1], mode='same'))]

    def tables(self):
        coverage = []; ge_rows = []; run_rows = []; channel_rows = []
        for period in PERIODS:
            for candidate in ['narrow', 'wide']:
                m = self.meta[period, candidate]
                coverage.append(dict(period=period, candidate=candidate, runs=len(m['runs']),
                    prompt_rows=sum(r['tr_event']['events'] for r in m['runs']),
                    random_rows=sum(r['tr_bg']['events'] for r in m['runs']),
                    si_pairs=sum(r['tr_event']['si_pairs'] for r in m['runs']),
                    s3_pairs=sum(r['tr_event']['s3_pairs'] for r in m['runs'])))
                for ge in self.ges:
                    for reference, family in [(0, 'si'), (1, 's3')]:
                        base = f'ge{ge}_{family}'
                        gate = self.gate(period, candidate, ge, reference)
                        for state, suffix in [('all',''),('normal','_normal'),('failed','_failed')]:
                            y = self.h1(period,candidate,base+'_prompt'+suffix)
                            ge_rows.append(dict(period=period,candidate=candidate,ge=ge,reference=family,cfd=state,
                                prompt_pairs=int(y.sum()),mode_ns=self.mode(y),lower_ns=gate['lower_ns'],upper_ns=gate['upper_ns'],
                                mode_edge_distance_ns=min(self.mode(y)-gate['lower_ns'],gate['upper_ns']-self.mode(y)),
                                random_left_pairs=int(self.h1(period,candidate,base+'_random_left'+suffix).sum()),
                                random_right_pairs=int(self.h1(period,candidate,base+'_random_right'+suffix).sum()),
                                alpha=gate['alpha']))
                            for r in m['runs']:
                                point = r['per_ge_prompt'].get(base+suffix,{})
                                run_rows.append(dict(period=period,candidate=candidate,run=r['run'],ge=ge,reference=family,cfd=state,
                                    prompt_pairs=point.get('pairs',0),mode_ns=point.get('local_mode_ns')))
            mapping = self.meta[period,'wide']['runs'][0]['channel_mapping']
            for key, item in sorted(mapping.items(),key=lambda pair:int(pair[0])):
                if item['kind'] not in ['si','spider','s3_ring','s3_sector']: continue
                kind = item['kind']
                name = f'channel{key}_ge_prompt' if kind != 's3_sector' else f'channel{key}_ring_minus_sector'
                row = dict(period=period,channel=int(key),kind=kind,ring=item['ring'],sector=item['sector'],
                           difference='ring - sector' if kind=='s3_sector' else 'particle - Ge')
                for candidate in ['narrow','wide']:
                    y = self.h1(period,candidate,name)
                    row[candidate+'_pairs']=int(y.sum()); row[candidate+'_mode_ns']=self.mode(y)
                channel_rows.append(row)
        return {name:pd.DataFrame(rows) for name,rows in [('coverage',coverage),('ge_summary',ge_rows),('per_run',run_rows),('channels',channel_rows)]}

    def line(self, ax, y, label, color, **kwargs):
        ax.stairs(np.where(y>0,y,np.nan),self.edges,label=label,color=color,baseline=None,**kwargs)

    @staticmethod
    def decorate(ax, xlabel='Corrected time difference (ns)', xlim=(-550,550)):
        ax.set(xlabel=xlabel,ylabel='Retained pairs / 2 ns',xlim=xlim,yscale='log',ylim=(0.7,None))
        ax.axvline(0,color='0.4',lw=0.8,ls=':'); ax.grid(alpha=.18,which='major')

    def overview(self):
        fig, axes = plt.subplots(3,2,figsize=(13,10),layout='constrained')
        for row,period in enumerate(PERIODS):
            for col,family in enumerate(['si','s3']):
                ax = axes[row,col]
                for candidate,color in [('wide','#167d63'),('narrow','#426bb8')]:
                    y = sum((self.h1(period,candidate,f'ge{ge}_{family}_prompt') for ge in self.ges),np.zeros(len(self.centers)))
                    self.line(ax,y,f'{candidate} prompt',color,lw=1.2)
                for side,color in [('left','#c36e22'),('right','#a04179')]:
                    y = sum((self.h1(period,'wide',f'ge{ge}_{family}_random_{side}') for ge in self.ges),np.zeros(len(self.centers)))
                    self.line(ax,y,'wide random '+side,color,lw=1)
                self.decorate(ax);ax.set_title(f'({chr(97+row*2+col)}) {PERIOD_LABELS[period]} | '+('Si array - Ge' if family=='si' else 'S3 ring - Ge'))
                ax.legend(fontsize=8,ncol=2)
        return fig

    def ge_figure(self,period,ge):
        fig,axes=plt.subplots(3,2,figsize=(13,12),layout='constrained')
        for col,family in enumerate(['si','s3']):
            base=f'ge{ge}_{family}';ax=axes[0,col]
            gw=self.gate(period,'wide',ge,col);gn=self.gate(period,'narrow',ge,col)
            ax.axvspan(gw['lower_ns'],gw['upper_ns'],color='#167d63',alpha=.07)
            for interval in gw['random']: ax.axvspan(*interval,color='#ba7a21',alpha=.08)
            for x in [gn['lower_ns'],gn['upper_ns']]:ax.axvline(x,color='#426bb8',lw=.9,ls=':')
            self.line(ax,self.h1(period,'wide',base+'_prompt'),'wide prompt (all CFD)','#167d63',lw=1.5)
            self.line(ax,self.h1(period,'narrow',base+'_prompt'),'narrow prompt (all CFD)','#426bb8',lw=1)
            for state,color in [('normal','#777777'),('failed','#d36523')]:
                self.line(ax,self.h1(period,'wide',base+'_prompt_'+state),'wide '+state,color,lw=.8,ls='--')
            for side,color in [('left','#a46b21'),('right','#a04179')]:
                self.line(ax,self.h1(period,'wide',base+'_random_'+side),'random '+side,color,lw=.9)
            self.decorate(ax);ax.set_title(f'({chr(97+col)}) '+('Si array - Ge' if family=='si' else 'S3 ring - Ge'))
            ax.legend(fontsize=7,ncol=2)
            data=[self.hist[period,'wide'][base+'_energy_'+s] for s in ['normal','failed']]
            vmax=max(2,max(v[0].max() for v in data))
            for row,(state,(counts,xe,ye)) in enumerate(zip(['normal','failed'],data),start=1):
                ax=axes[row,col]
                image=ax.pcolormesh(xe,ye,np.ma.masked_less_equal(counts.T,0),norm=LogNorm(vmin=1,vmax=vmax),cmap='viridis',rasterized=True)
                lower=min(-60,gw['lower_ns']-12);upper=max(60,gw['upper_ns']+12)
                ax.axhspan(lower,gw['lower_ns'],color='0.88');ax.axhspan(gw['upper_ns'],upper,color='0.88')
                ax.axhline(0,color='tomato',lw=.9,ls='--')
                for y in [gw['lower_ns'],gw['upper_ns']]:ax.axhline(y,color='0.35',lw=.8,ls=':')
                ax.set(xlim=(40,2000),ylim=(lower,upper),xlabel='Ge energy (keV)',ylabel='Particle - Ge (ns)',
                       title=f'({chr(97+row*2+col)}) Wide prompt: {state} CFD | N={int(counts.sum()):,}')
                fig.colorbar(image,ax=ax,label='Pairs / (20 keV x 4 ns)',pad=.02)
        fig.suptitle(f'{PERIOD_LABELS[period]} | Ge{ge} | stored corrected times',fontsize=14)
        return fig

    def controls(self,period):
        fig,axes=plt.subplots(1,3,figsize=(15,4.3),layout='constrained')
        for i,(name,title) in enumerate(zip(CONTROL_NAMES,CONTROL_LABELS)):
            ax=axes[i]
            for candidate,color in [('wide','#167d63'),('narrow','#426bb8')]:
                y=self.h1(period,candidate,name)
                self.line(ax,y,f'{candidate}; N={int(y.sum()):,}',color)
            if i==0:
                ax.axvspan(-80,80,color='#167d63',alpha=.06)
                for x in [-80,80]:ax.axvline(x,color='0.4',ls=':',lw=1)
            self.decorate(ax,xlim=(-250,250));ax.set_title(f'({chr(97+i)}) {title}');ax.legend(fontsize=8)
        fig.suptitle(PERIOD_LABELS[period]+' | unique particle pairs retained in prompt records',fontsize=12)
        return fig

    def channel_figure(self,period):
        mapping=self.meta[period,'wide']['runs'][0]['channel_mapping']
        groups=[('si','_ge_prompt','Si array - Ge'),('s3_ring','_ge_prompt','S3 ring - Ge'),('s3_sector','_ring_minus_sector','S3 ring - sector (by sector)')]
        fig,axes=plt.subplots(1,3,figsize=(14,13),layout='constrained')
        for i,(kind,suffix,label) in enumerate(groups):
            keys=sorted(int(k) for k,v in mapping.items() if v['kind']==kind or (kind=='si' and v['kind']=='spider'))
            y=np.array([self.h1(period,'wide',f'channel{k}'+suffix) for k in keys])
            ax=axes[i]
            if not keys:ax.text(.5,.5,'No mapped channels',transform=ax.transAxes);continue
            im=ax.pcolormesh(self.edges,np.arange(len(keys)+1),np.ma.masked_less_equal(y,0),norm=LogNorm(vmin=1,vmax=max(2,y.max())),cmap='viridis',rasterized=True)
            ax.set(yticks=np.arange(len(keys))+.5,yticklabels=[f'{k}  ({int(n):,})' for k,n in zip(keys,y.sum(axis=1))],
                   xlim=(-150,150),xlabel='Time difference (ns)',title=f'({chr(97+i)}) {label}')
            ax.tick_params(axis='y',labelsize=6.5);ax.invert_yaxis();ax.axvline(0,color='tomato',ls='--',lw=.8)
            fig.colorbar(im,ax=ax,label='Pairs / 2 ns',pad=.01,shrink=.65)
        fig.suptitle(PERIOD_LABELS[period]+' | wide prompt | channel (retained pair count)',fontsize=13)
        return fig

    def run_modes(self,table):
        fig,axes=plt.subplots(3,2,figsize=(15,10),layout='constrained')
        for row,period in enumerate(PERIODS):
            for col,family in enumerate(['si','s3']):
                d=table[(table.period==period)&(table.candidate=='wide')&(table.reference==family)&(table.cfd=='all')]
                mode=d.pivot(index='ge',columns='run',values='mode_ns').reindex(self.ges)
                count=d.pivot(index='ge',columns='run',values='prompt_pairs').reindex(self.ges)
                values=np.ma.masked_where((count.values<100)|~np.isfinite(mode.values.astype(float)),mode.values.astype(float))
                ax=axes[row,col];im=ax.imshow(values,aspect='auto',interpolation='nearest',cmap='coolwarm',vmin=-20,vmax=20)
                positions=np.unique(np.linspace(0,len(mode.columns)-1,min(8,len(mode.columns))).round().astype(int))
                ax.set(xticks=positions,xticklabels=[f'{mode.columns[j]:04d}' for j in positions],yticks=range(len(self.ges)),yticklabels=self.ges,
                       xlabel='Run',ylabel='Ge',title=f'({chr(97+row*2+col)}) {PERIOD_LABELS[period]} | {family} - Ge')
                ax.tick_params(axis='y',labelsize=7);fig.colorbar(im,ax=ax,label='Smoothed prompt mode (ns)',pad=.02,extend='both')
        fig.suptitle('Wide prompt, all CFD | 2 ns bins | blank if fewer than 100 retained pairs',fontsize=13)
        return fig
