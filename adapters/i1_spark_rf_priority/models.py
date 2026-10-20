"""Thin timed-resistor wrapper over unchanged, accepted source/port protocols."""
from copy import copy
from dataclasses import replace
import math

from adapters.i1_fp3_port_circuit.ports import PortError
from adapters.i1_teng_source.circuit import ReferenceParallelRCPort


def resistance(curve, time_s):
    if not math.isfinite(time_s) or time_s < 0:
        raise PortError('OUT_OF_MODEL_DOMAIN')
    if curve in ('R20', 'R60', 'R200'):
        return float(curve[1:])
    if curve == 'RAMP20_200':
        return 20.+180.*min(time_s/25e-9, 1.)
    raise PortError('UNREGISTERED_RESISTANCE_CURVE')


class TimedResistorPort(ReferenceParallelRCPort):
    def __init__(self, curve, **kwargs):
        self.curve = curve
        super().__init__(.2e-12, resistance(curve, 0.), **kwargs)
        self.descriptor = replace(self.descriptor, provider_id='SYNTHETIC_PRESCRIBED_TIME_R_PARALLEL_CP',
                                  operator_id='MIDPOINT_RESISTOR_WORK_V1',
                                  model_identity=(('curve', curve), ('Cp_F', self.Cp)),
                                  validity_domain='SYNTHETIC_EVENT_NOT_CALIBRATED_SPARK')

    def evaluate_trial(self, request):
        # A shallow proxy only changes a scalar, preserving the immutable accepted token.
        proxy = copy(self)
        proxy.R_test = resistance(self.curve, request.t_n_s+.5*request.dt_s)
        result = ReferenceParallelRCPort.evaluate_trial(proxy, request)
        return replace(result, provenance='SYNTHETIC_PRESCRIBED_TIME_R_NOT_LTE_OR_TOEPLER_RW')


def toepler_guard(K_T, ell, q0, signed_charge, parameter_identity):
    if not parameter_identity or not all(math.isfinite(v) for v in (K_T, ell, q0, signed_charge)):
        raise PortError('MISSING_OR_INVALID_MODEL_PARAMETERS')
    if K_T <= 0 or ell <= 0 or q0 <= 0 or q0+signed_charge <= 0:
        raise PortError('OUT_OF_MODEL_DOMAIN')
    return K_T*ell/(q0+signed_charge)


def rw_guard(K_RW, ell, S0, squared_current_integral, parameter_identity):
    if not parameter_identity or not all(math.isfinite(v) for v in (K_RW, ell, S0, squared_current_integral)):
        raise PortError('MISSING_OR_INVALID_MODEL_PARAMETERS')
    if K_RW <= 0 or ell <= 0 or S0 <= 0 or squared_current_integral < 0:
        raise PortError('OUT_OF_MODEL_DOMAIN')
    return K_RW*ell/math.sqrt(S0+squared_current_integral)
