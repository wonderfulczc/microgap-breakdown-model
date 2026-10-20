"""Pure implicit-midpoint formation trials on the accepted FP3 capacitor protocol.

Q/S are explicitly dimensioned history variables, never stored in x_V or counted
as thermal energy. Existing TENGSeriesCircuit owns the joint source/circuit commit.
"""
from dataclasses import asdict, dataclass, replace
import math

from adapters.i1_fp3_port_circuit.ports import AcceptedState, IdealCapacitorPort, PortError
from adapters.i1_teng_source.source import quantity


def real(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise PortError('MISSING_OR_NONFINITE:'+name)
    if positive and value <= 0:
        raise PortError('OUT_OF_MODEL_DOMAIN:'+name)
    return float(value)


@dataclass(frozen=True)
class FormationParameters:
    law: str
    K: float
    ell_m: float
    initial_value: float
    valid_end_s: float
    parameter_identity: str
    initialization_identity: str
    parameter_status: str = 'SYNTHETIC_NOT_CALIBRATED'

    def __post_init__(self):
        if self.law not in ('TOEPLER', 'RW'):
            raise PortError('UNSUPPORTED_FORMATION_LAW')
        for k in ('K','ell_m','initial_value','valid_end_s'):
            real(getattr(self,k),k,positive=True)
        if not all(isinstance(x,str) and x.strip() for x in (self.parameter_identity,self.initialization_identity)):
            raise PortError('SOURCE_IDENTITY_REQUIRED')
        if self.parameter_status != 'SYNTHETIC_NOT_CALIBRATED':
            raise PortError('DEVICE_PHYSICAL_RUN_NOT_AUTHORIZED_PARAMETER_PACKAGE_REQUIRED')
        real(self.K*self.ell_m,'K_times_ell',positive=True)
        self.resistance(self.initial_value)

    @property
    def state_unit(self):
        return 'C' if self.law == 'TOEPLER' else 'A^2*s'

    @classmethod
    def from_mapping(cls, spec):
        if not isinstance(spec,dict):
            raise PortError('FORMATION_PARAMETER_PACKAGE_REQUIRED')
        if set(spec) & {'recovery_tau','R_floor','abs_current','recovery_law'}:
            raise PortError('RECOVERY_OR_MODIFIED_LAW_NOT_IMPLEMENTED')
        required = ('law','K','ell','initial_state','valid_end','parameter_identity',
                    'initialization_identity','parameter_status')
        for key in required:
            if key not in spec:
                raise PortError('FORMATION_PARAMETER_MISSING:'+key)
        if spec['law'] not in ('TOEPLER','RW'):
            raise PortError('UNSUPPORTED_FORMATION_LAW')
        toepler = spec['law']=='TOEPLER'
        return cls(spec['law'],quantity(spec['K'],'V*s/m' if toepler else 'V*sqrt(s)/m','K'),
                   quantity(spec['ell'],'m','ell'),
                   quantity(spec['initial_state'],'C' if toepler else 'A^2*s','initial_state'),
                   quantity(spec['valid_end'],'s','valid_end'),spec['parameter_identity'],
                   spec['initialization_identity'],spec['parameter_status'])

    def resistance(self, value):
        value = real(value,'formation_state',positive=True)
        result = self.K*self.ell_m/(value if self.law=='TOEPLER' else math.sqrt(value))
        return real(result,'resistance',positive=True)

    def exact_current_interval(self, value, current_A, t0_s, t1_s):
        """Formula-only signed-current test, not a physical bipolar port advance."""
        value = real(value,'formation_state',positive=True)
        current_A = real(current_A,'I_col_A')
        t0_s,t1_s = real(t0_s,'t0_s'),real(t1_s,'t1_s')
        if t0_s < 0 or t1_s <= t0_s or t1_s > self.valid_end_s:
            raise PortError('OUT_OF_MODEL_TIME_DOMAIN')
        delta = (t1_s-t0_s)*(current_A if self.law=='TOEPLER' else current_A**2)
        result = real(value+delta,'formation_state',positive=True)
        self.resistance(result)
        return result


@dataclass(frozen=True)
class FormationAccepted(AcceptedState):
    formation_value: float = 0.
    formation_unit: str = ''
    formation_law: str = ''


class FormationPort(IdealCapacitorPort):
    def __init__(self, parameters, Cp_F, **kwargs):
        if not isinstance(parameters,FormationParameters):
            raise PortError('FORMATION_PARAMETER_PACKAGE_REQUIRED')
        self.parameters = parameters
        super().__init__(Cp_F,**kwargs)
        if self.accepted.voltage_V < 0:
            raise PortError('FORMATIVE_UNIPOLAR_ONLY')
        self._accepted = FormationAccepted(**asdict(self.accepted),
            formation_value=parameters.initial_value,formation_unit=parameters.state_unit,
            formation_law=parameters.law)
        self.descriptor = replace(self.descriptor,provider_id=parameters.law+'_FORMATION_PORT',
            operator_id='IMPLICIT_MIDPOINT_FORMATION_V1',
            validity_domain='FORMATIVE_UNIPOLAR_ONLY_SYNTHETIC_NOT_CALIBRATED',
            model_identity=(('law',parameters.law),('K',parameters.K),('ell_m',parameters.ell_m),
                            ('state_unit',parameters.state_unit),('valid_end_s',parameters.valid_end_s),
                            ('parameter_identity',parameters.parameter_identity)),
            capabilities=self.descriptor.capabilities+('formation_state_transactions',))

    def evaluate_trial(self, request):
        base = super().evaluate_trial(request)
        s,p = self.accepted,self.parameters
        dt = request.dt_s
        if request.t_n_s < 0 or request.t_n_s+dt > p.valid_end_s:
            raise PortError('FORMATION_TIME_LIMIT_REACHED')
        if request.voltage_n_V < 0 or request.voltage_n1_V < 0:
            raise PortError('FORMATIVE_UNIPOLAR_ONLY')
        vm = .5*(request.voltage_n_V+request.voltage_n1_V)
        a = p.K*p.ell_m
        value = real(s.formation_value,'accepted_formation_state',positive=True)
        if p.law == 'TOEPLER':
            denominator = a-.5*dt*vm
            if denominator <= 0:
                raise PortError('TRIAL_INTERVAL_OUT_OF_MODEL_DOMAIN')
            j = vm*value/denominator
            dj_dv1 = .5*value*a/denominator**2
            value1 = value+dt*j
        else:
            denominator = a*a-.5*dt*vm*vm
            if denominator <= 0:
                raise PortError('TRIAL_INTERVAL_OUT_OF_MODEL_DOMAIN')
            j = vm*math.sqrt(value)/math.sqrt(denominator)
            dj_dv1 = .5*math.sqrt(value)*a*a/denominator**1.5
            value1 = value+dt*j*j
        value1 = real(value1,'candidate_formation_state',positive=True)
        p.resistance(value1)
        dq,work = dt*j,vm*dt*j
        real(j,'I_col_mid_A')
        real(work,'column_work_J')
        real(dj_dv1,'dIcol_dVp1')
        candidate_data = asdict(base.candidate)
        candidate_data.update(cumulative_charge_C=base.candidate.cumulative_charge_C+dq,
            cumulative_port_work_J=base.candidate.cumulative_port_work_J+work,
            cumulative_component_work_J=base.candidate.cumulative_component_work_J+work,
            cumulative_dissipation_J=base.candidate.cumulative_dissipation_J+work)
        candidate = FormationAccepted(**candidate_data,formation_value=value1,
                                      formation_unit=p.state_unit,formation_law=p.law)
        return replace(base,candidate=candidate,delta_q_transfer_C=dq,
            delta_q_port_C=base.delta_q_port_C+dq,
            port_work_J=base.port_work_J+work,component_work_J=base.component_work_J+work,
            resistor_dissipation_J=work,transfer_current_mid_A=j,
            d_deltaq_port_d_voltage_n1_F=base.d_deltaq_port_d_voltage_n1_F+dt*dj_dv1,
            component_energy_residual_J=base.component_energy_residual_J,
            total_energy_residual_J=base.total_energy_residual_J,
            provenance='FORMATION_HISTORY_NOT_HEAT_NOT_RECOVERY_NOT_DEVICE_CALIBRATION')


class ToeplerFormationPort(FormationPort):
    def __init__(self,parameters,Cp_F,**kwargs):
        if parameters.law != 'TOEPLER':
            raise PortError('FORMATION_LAW_MISMATCH')
        super().__init__(parameters,Cp_F,**kwargs)


class RWFormationPort(FormationPort):
    def __init__(self,parameters,Cp_F,**kwargs):
        if parameters.law != 'RW':
            raise PortError('FORMATION_LAW_MISMATCH')
        super().__init__(parameters,Cp_F,**kwargs)
