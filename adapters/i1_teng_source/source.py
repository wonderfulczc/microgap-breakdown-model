"""Frozen-position V-Q-x interface; parameters are never inferred from a gap."""
from dataclasses import dataclass
import math

from adapters.i1_fp3_port_circuit.ports import PortError


def quantity(spec, unit, name):
    if not isinstance(spec,dict) or spec.get('unit')!=unit or not spec.get('identity'):
        raise PortError('MISSING_OR_INCOMPATIBLE_SOURCE_QUANTITY:'+name)
    value=spec.get('value')
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise PortError('NONFINITE_OR_MISSING_SOURCE_VALUE:'+name)
    return float(value)


@dataclass(frozen=True)
class TENGFrozenPositionSource:
    C_T_F: float
    V_oc_V: float
    q_T0_C: float
    x0_m: float
    source_identity: str
    parameter_identity: str
    source_kind: str = 'TENG_VQX_FROZEN_POSITION'

    def __post_init__(self):
        if (not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
                    for v in (self.C_T_F,self.V_oc_V,self.q_T0_C,self.x0_m))
                or self.C_T_F<=0):
            raise PortError('INVALID_SOURCE_PARAMETERS')
        if (not all(isinstance(v,str) and v.strip() for v in (self.source_identity,self.parameter_identity))
                or self.source_kind!='TENG_VQX_FROZEN_POSITION'):
            raise PortError('MISSING_OR_INVALID_SOURCE_IDENTITY')

    @classmethod
    def from_mapping(cls, spec):
        if not isinstance(spec,dict):
            raise PortError('SOURCE_PARAMETERS_REQUIRED')
        if spec.get('mechanical_mode')!='FROZEN_POSITION' or 'x_of_t' in spec:
            raise PortError('DYNAMIC_MECHANICAL_SOURCE_NOT_IMPLEMENTED')
        if spec.get('topology')!='CONNECTED_SINGLE_SERIES_LOOP':
            raise PortError('ALTERNATIVE_DISCONNECTED_SOURCE_NOT_SELECTED')
        if ('q_T0' not in spec or spec['q_T0'] is None or
                isinstance(spec['q_T0'],dict) and spec['q_T0'].get('value') is None):
            raise PortError('INITIAL_SOURCE_CHARGE_UNKNOWN')
        for key in ('C_T','V_oc','x0','source_identity','parameter_identity'):
            if key not in spec:
                raise PortError('SOURCE_PARAMETER_MISSING:'+key)
        return cls(quantity(spec['C_T'],'F','C_T'),quantity(spec['V_oc'],'V','V_oc'),
                   quantity(spec['q_T0'],'C','q_T0'),quantity(spec['x0'],'m','x0'),
                   spec['source_identity'],spec['parameter_identity'])

    @property
    def inverse_capacitance(self):
        return 1/self.C_T_F

    def voltage(self,q_T_C):
        if not math.isfinite(q_T_C):
            raise PortError('NONFINITE_SOURCE_CHARGE')
        return self.V_oc_V-q_T_C/self.C_T_F

    def coordinate_energy(self,q_T_C):
        return .5*q_T_C*q_T_C/self.C_T_F

    def require_mechanical_mode(self, mode):
        if mode!='FROZEN_POSITION':
            raise PortError('DYNAMIC_MECHANICAL_SOURCE_NOT_IMPLEMENTED')


@dataclass(frozen=True)
class IdealVoltageLimit:
    """Explicit mathematical C_T -> infinity limit, NOT a calibrated TENG."""
    V_oc_V: float
    q_T0_C: float = 0.
    source_identity: str = 'SYNTHETIC_IDEAL_VOLTAGE_LIMIT'
    source_kind: str = 'SYNTHETIC_IDEAL_VOLTAGE_LIMIT'
    inverse_capacitance: float = 0.

    def __post_init__(self):
        if not math.isfinite(self.V_oc_V) or not math.isfinite(self.q_T0_C):
            raise PortError('NONFINITE_SOURCE_PARAMETERS')
        if self.inverse_capacitance!=0:
            raise PortError('INVALID_IDEAL_SOURCE_LIMIT')

    def voltage(self,q_T_C):
        return self.V_oc_V

    def coordinate_energy(self,q_T_C):
        return 0.
