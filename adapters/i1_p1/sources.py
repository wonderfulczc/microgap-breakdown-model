"""Ideal contact-separation geometry feeding the existing frozen-charge source."""
from dataclasses import dataclass
import math
import numpy as np
from adapters.i1_teng_source.source import TENGFrozenPositionSource
from .device import value, need

EPS0=8.8541878128e-12  # SI numerical reference; no fitted dielectric constant.


@dataclass(frozen=True)
class TENGContactSeparationSource:
    contract: dict

    def validate(self):
        c=self.contract
        need(c.get('model_id')=='IDEAL_PARALLEL_PLATE_CONTACT_SEPARATION','SOURCE_MODEL_ID')
        need(c.get('mechanical_mode')=='FROZEN_POSITION_AND_SIGMA' and 'x_of_t' not in c,
             'DYNAMIC_MECHANICAL_SOURCE_NOT_IMPLEMENTED')
        need(c.get('leakage_mode')=='EXPLICIT_DISABLED_REFERENCE','LEAKAGE_MODEL_UNAVAILABLE')
        need(c.get('contact_model')=='IDEAL_SEPARATION_NO_CONTACT_CHARGING_HISTORY','CONTACT_HISTORY_UNAVAILABLE')
        area=value(c['A_overlap'],'m2',True)
        upper=value(c['upper_area'],'m2',True);lower=value(c['lower_area'],'m2',True)
        need(area<=min(upper,lower),'OVERLAP_EXCEEDS_ACTIVE_SURFACE')
        need(c.get('exposed_copper_counted') is False,'EXPOSED_COPPER_NOT_TRIBO_AREA')
        value(c['sigma_eff'],'C/m2')
        need(set(x['material'] for x in c['dielectrics'])>={'PET','FEP'},'DIELECTRIC_PATH_INCOMPLETE')
        for layer in c['dielectrics']:
            value(layer['thickness'],'m',True);value(layer['epsilon_r'],'1',True)
            need(type(layer.get('in_metal_to_metal_path')) is bool,'ADHESIVE_PATH_UNKNOWN')
            need(bool(layer.get('path_basis')) and bool(layer.get('frequency_domain')),'DIELECTRIC_DOMAIN_MISSING')
        need(bool(c.get('model_limitations')),'MODEL_ERROR_IDENTITY_REQUIRED')
        return self

    @property
    def d_eff_m(self):
        self.validate()
        return sum(value(l['thickness'],'m')/value(l['epsilon_r'],'1')
                   for l in self.contract['dielectrics'] if l['in_metal_to_metal_path'])

    def at(self, x_m):
        need(type(x_m) in (float,int) and math.isfinite(x_m) and x_m>=0,'INVALID_TENG_SEPARATION')
        d=self.d_eff_m
        need(d>0,'ZERO_DIELECTRIC_THICKNESS')
        area=value(self.contract['A_overlap'],'m2')
        sigma=value(self.contract['sigma_eff'],'C/m2')
        ct=EPS0*area/(x_m+d);voc=sigma*x_m/EPS0
        return dict(x_T_m=x_m,d_eff_m=d,C_T_F=ct,Voc_V=voc,Qsc_C=ct*voc)

    def freeze(self, x_quantity, charge_quantity):
        x=value(x_quantity,'m')
        qt=value(charge_quantity,'C')
        a=self.at(x)
        return TENGFrozenPositionSource(a['C_T_F'],a['Voc_V'],qt,x,
             'IDEAL_CONTACT_SEPARATION_FROZEN_GEOMETRY','DERIVED_CONDITIONAL_UNCALIBRATED')


class HVAmplifierSource:
    """Finite-band waveform and explicit series R+L impedance family only."""
    def __init__(self, c):
        need(c.get('source_kind')=='HV_AMPLIFIER','AMPLIFIER_IDENTITY')
        self.identity=c.get('qualification')
        need(self.identity in ('IDEAL_AMPLIFIER_LIMIT','SOURCE_BOUND_RL_OUTPUT'),'AMPLIFIER_IMPEDANCE_UNKNOWN')
        self.times=np.asarray(c.get('times_s'),float);self.volts=np.asarray(c.get('voltage_V'),float)
        need(self.times.ndim==1 and self.times.size>=2 and self.volts.shape==self.times.shape
             and np.all(np.isfinite(self.times)) and np.all(np.isfinite(self.volts))
             and np.all(np.diff(self.times)>0),'AMPLIFIER_SIGNAL_DOMAIN')
        need(c.get('signal_unit')=='V' and c.get('time_unit')=='s' and c.get('signal_source'),
             'AMPLIFIER_SIGNAL_SOURCE')
        from .archive import load
        evidence=load(c['signal_source'])
        need(evidence['times_s']==c['times_s'] and evidence['voltage_V']==c['voltage_V'],'AMPLIFIER_SIGNAL_HASH')
        if self.identity=='IDEAL_AMPLIFIER_LIMIT':
            need(c.get('Rout') is None and c.get('Lout') is None,'IDEAL_LIMIT_MUST_NOT_HIDE_IMPEDANCE')
            self.Rout=0.;self.Lout=0.
        else:
            self.Rout=value(c['Rout'],'ohm');self.Lout=value(c['Lout'],'H')
            need(self.Rout>=0 and self.Lout>=0 and c.get('impedance_domain'),'OUTPUT_IMPEDANCE_DOMAIN')
        need(c.get('interpolation')=='LINEAR_NO_EXTRAPOLATION','AMPLIFIER_INTERPOLATION')

    def voltage(self,t):
        need(self.times[0]<=t<=self.times[-1],'AMPLIFIER_SIGNAL_OUT_OF_DOMAIN')
        return float(np.interp(t,self.times,self.volts))

    def interval_voltage(self,t0,t1):
        need(not np.any((self.times>t0)&(self.times<t1)),'AMPLIFIER_BREAKPOINT_REQUIRES_STEP_SPLIT')
        return .5*(self.voltage(t0)+self.voltage(t1))
