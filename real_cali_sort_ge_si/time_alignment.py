"""Read the project's active run/Ge/CFD timing profile without changing data.

New processing should instantiate Calibration() and use the complete offsets
with raw decoded timestamps. Historical additive profiles remain explicit.
"""
from pathlib import Path
import json
import numpy as np

class Calibration:
    def __init__(self,profile=None,pars=None):
        pars=Path(pars) if pars is not None else Path(__file__).resolve().parents[1]/'pars'
        if profile is None:
            active=pars/'time_alignment/active.txt'
            profile=active.read_text().strip() if active.exists() else 'legacy'
        self.profile=profile;self.data=None;self.held=set();self.total=False
        if profile=='legacy':return
        directory=pars/'time_alignment/versions'/profile
        self.metadata=json.loads((directory/'metadata.json').read_text())
        self.data=json.loads((directory/'calibration.json').read_text())
        # Complete profiles need no historical additive-branch exclusion list.
        self.held={(x['run'],x['ge'],x['state']) for x in self.metadata.get('held_back',[])}
        self.total=self.data.get('correction_mode')=='total'

    def _require_additive(self):
        if self.total:
            raise ValueError('Complete profile: use total_ge_offset_ns / total_particle_offset_ns on raw decoded times; do not add legacy offsets')

    def _total_run(self,run):
        if not self.total:raise ValueError('A complete timing profile is required')
        run=f'{int(run):04d}'
        if run not in self.data['groups']:raise ValueError('Run outside complete profile; use explicit legacy processing')
        return run

    def total_ge_offset_ns(self,run,ge,energy_keV,failed=False):
        """Complete offset: add once to native integer CFD/coarse Ge time."""
        run=self._total_run(run);E=np.asarray(energy_keV,dtype=float)
        base=self.data['channel_offsets'][self.data['groups'][run]].get(str(int(ge)))
        if base is None:raise ValueError('Unmapped Ge timing channel')
        result=np.full_like(E,base)
        curve=self.data['ge_curves'].get(f'{run}:{int(ge)}:{int(bool(failed))}')
        if curve is None:return result
        valid=np.isfinite(E)&(E>=40)&(E<2000)
        values=np.interp(np.log(np.clip(np.where(np.isfinite(E),E,40.),40.,2000.)),np.log(curve['knots_keV']),curve['offset_ns'])
        if failed:
            models=self.data.get('failed_center_models',{})
            model=models.get(f'{run}:{int(ge)}',models.get(f"{self.data['groups'][run]}:{int(ge)}"))
            if model is not None:
                values=values+self._failed_center_value(model,np.where(valid,E,40.))
        return np.where(valid,values,result)

    @staticmethod
    def _failed_center_value(model,energy_keV):
        """Exact fitted family, clamped to measured support; units ns/keV."""
        e=np.clip(np.asarray(energy_keV,dtype=float),*model['support_keV'])
        c=np.asarray(model['coefficients_ns']);u=np.log(e/334.)
        family=model['family']
        if family=='constant':return np.full_like(e,c[0])
        if family=='log_linear':return c[0]+c[1]*u
        if family=='log_quadratic':return c[0]+c[1]*u+c[2]*u*u
        if family=='inverse_sqrt':return c[0]+c[1]*(np.sqrt(334./e)-1.)
        if family=='inverse':return c[0]+c[1]*(334./e-1.)
        if family=='piecewise_log_linear':return np.interp(np.log(e),np.log(model['knots_keV']),c)
        raise ValueError('Unsupported failed-CFD center model: '+family)

    def total_particle_offset_ns(self,run,key,energy_keV):
        """Complete particle offset, physical energy in keV; no reference guess."""
        run=self._total_run(run);E=np.asarray(energy_keV,dtype=float)
        if np.any(~np.isfinite(E)) or np.any(E>=1e9):raise ValueError('Invalid physical particle timing energy')
        E=np.maximum(E,0.);pieces=self.data['particle_tables'][run].get(str(int(key)))
        if pieces is None:raise ValueError('Unmapped particle timing channel')
        out=np.zeros_like(E)
        for lo,hi,slope,intercept in pieces:
            out=np.where((E>=lo)&(E<hi),slope*E+intercept,out)
        return out

    def total_particle_offset_from_adc_ns(self,run,key,adc):
        """Apply the physical ADC calibration in every beam period.

        Older high-beam profiles without ADC coefficients are compatible only
        when the entire timing function is constant, as in the C++ reader.
        """
        run=self._total_run(run);adc=np.asarray(adc,dtype=float)
        if np.any(~np.isfinite(adc)):raise ValueError('Non-finite particle ADC')
        key=str(int(key))
        c=self.data.get('particle_calibration',{}).get(run,{}).get(key)
        if c is None:
            pieces=self.data['particle_tables'][run].get(key)
            if self.data['groups'][run]=='low' or not pieces:
                raise ValueError('Unmapped physical particle calibration')
            if any(slope!=0. or intercept!=pieces[0][3]
                   for lo,hi,slope,intercept in pieces):
                raise ValueError('Energy-dependent particle timing requires ADC calibration')
            return self.total_particle_offset_ns(run,key,np.zeros_like(adc))
        return self.total_particle_offset_ns(run,key,c[0]+c[1]*adc+c[2]*adc*adc)

    def correct_raw_dt(self,run,ge,energy_keV,dt_raw_decoded_ns,particle_key,particle_energy_keV,failed=False):
        """Real-valued dt from decoded clocks BEFORE ANY channel/time-walk offset.

        For production integer results correct the two absolute timestamps
        separately with correct_raw_timestamps(); dt alone loses quantization.
        """
        return (np.asarray(dt_raw_decoded_ns,dtype=float)
                +self.total_particle_offset_ns(run,particle_key,particle_energy_keV)
                -self.total_ge_offset_ns(run,ge,energy_keV,failed))

    def correct_raw_timestamps(self,run,ge,energy_keV,ge_decoded_ns,particle_decoded_ns,
                               particle_key,particle_adc,failed=False):
        """One calibration addition per raw integer timestamp, prior to time dither.

        Ge input includes the native CFD decoding/int64 conversion, but includes
        no channel offset or time-walk correction. Returns (Ge, particle).
        """
        gt=np.asarray(ge_decoded_ns,dtype=np.float64)+self.total_ge_offset_ns(run,ge,energy_keV,failed)
        pt=np.asarray(particle_decoded_ns,dtype=np.float64)+self.total_particle_offset_from_adc_ns(run,particle_key,particle_adc)
        return gt.astype(np.int64),pt.astype(np.int64)

    def delta_ge_ns(self,run,ge,energy_keV,failed=False):
        self._require_additive()
        E=np.asarray(energy_keV,dtype=float);answer=np.zeros_like(E)
        if np.any(~np.isfinite(E)):raise ValueError('Energy must be finite')
        run=f'{int(run):04d}';key=(run,int(ge),int(bool(failed)))
        if self.data is None or f'{run}:{int(ge)}:{int(bool(failed))}' not in self.data['models'] or key in self.held:return answer
        m=self.data['models'][f'{run}:{int(ge)}:{int(bool(failed))}']
        valid=(E>=40)&(E<2000)
        def interp(c):return np.interp(np.log(np.clip(E,40,2000)),np.log(c['knots_keV']),c['offset_ns'])
        if m.get('group')=='low':
            return np.where(valid,interp(self.data['low_curves'][m['low_curve_key']]),answer)
        p=self.data['period_curves'][m['period_key']]
        return np.where(valid,m['normal_anchor_ns']+interp(p)+interp(m['deviation']),answer)

    def reference_shift_ns(self,run,reference='Si array'):
        self._require_additive()
        if reference not in ['Si array','S3 ring','S3 sector']:raise ValueError('Unknown particle reference')
        run=f'{int(run):04d}'
        if self.data is None or run not in self.data['runs']:return 0.
        if run in self.data.get('low_particle_tables',{}):
            raise ValueError('Low-beam references need particle_key and particle_energy_keV; use delta_particle_ns or correct_dt with both arguments')
        if reference=='Si array':return 0.
        return float(self.data['models'][f'{run}:{self.data["ge"][0]}:0']['reference_shift_ns'][1])

    def delta_particle_ns(self,run,key,energy_keV,reference='Si array'):
        """Additional particle timestamp term; physical particle energy is keV."""
        self._require_additive()
        E=np.asarray(energy_keV,dtype=float)
        if np.any(~np.isfinite(E)):raise ValueError('Particle energy must be finite')
        run=f'{int(run):04d}'
        if self.data is None or run not in self.data['runs']:return np.zeros_like(E)
        if run not in self.data.get('low_particle_tables',{}):
            return np.full_like(E,-self.reference_shift_ns(run,reference))
        if np.any(E>=1e9):raise ValueError('Particle energy exceeds physical ADC domain')
        E=np.maximum(E,0.)
        pieces=self.data['low_particle_tables'][run].get(str(int(key)))
        if pieces is None:raise ValueError('Unmapped low-beam particle channel')
        result=np.zeros_like(E)
        for lo,hi,slope,intercept in pieces:
            select=(E>=lo)&(E<hi)
            result=np.where(select,slope*E+intercept,result)
        return result

    def delta_particle_from_adc_ns(self,run,key,adc,reference='Si array'):
        self._require_additive()
        run=f'{int(run):04d}';adc=np.asarray(adc,dtype=float)
        if self.data is None or run not in self.data.get('low_particle_calibration',{}):
            return self.delta_particle_ns(run,key,adc,reference)
        c=self.data['low_particle_calibration'][run].get(str(int(key)))
        if c is None:raise ValueError('Unmapped low-beam particle calibration')
        return self.delta_particle_ns(run,key,c[0]+c[1]*adc+c[2]*adc*adc,reference)

    def correct_dt(self,run,ge,energy_keV,dt_native_ns,failed=False,reference='Si array',
                   particle_key=None,particle_energy_keV=None):
        self._require_additive()
        normalized=f'{int(run):04d}'
        if self.data is not None and normalized in self.data.get('low_particle_tables',{}):
            if particle_key is None or particle_energy_keV is None:
                raise ValueError('Low-beam correct_dt requires particle_key and particle_energy_keV')
            return np.asarray(dt_native_ns,dtype=float)+self.delta_particle_ns(run,particle_key,particle_energy_keV,reference)-self.delta_ge_ns(run,ge,energy_keV,failed)
        return np.asarray(dt_native_ns,dtype=float)-self.reference_shift_ns(run,reference)-self.delta_ge_ns(run,ge,energy_keV,failed)
