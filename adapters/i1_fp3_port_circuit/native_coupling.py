"""Bounded full-trial charge residual, only LOCAL_NUMERICAL_TEST scope."""
from dataclasses import dataclass

from .circuit import CircuitParameters, SyntheticSource
from .ports import PortError


@dataclass(frozen=True)
class LocalCircuitState:
    version: int
    provider_token: str
    time_s: float
    i_A: float
    Vc_V: float
    cumulative_source_work_J: float = 0.
    cumulative_Rs_loss_J: float = 0.


class LocalNumericalFullCircuit:
    def __init__(self, port, *, operation_scope, i0_A=.01, Vc0_V=0.,
                 parameters=CircuitParameters(), source=SyntheticSource(500.03)):
        if operation_scope != "LOCAL_NUMERICAL_TEST":
            raise PortError("UNAVAILABLE_CAPABILITY_GLOBAL_PHYSICAL_PORT")
        port.require("provide_local_numerical_port_charge")
        self.port, self.p, self.source = port, parameters, source
        self.state = LocalCircuitState(0, port.accepted.token, port.accepted.time_s, i0_A, Vc0_V)
        self.trials, self.records = [], []

    def advance(self, dt_s=20e-15):
        if self.state.version >= 3 or dt_s != 20e-15:
            raise PortError("LOCAL_MACRO_BUDGET_OR_SCOPE")
        s, p = self.state, self.p
        if s.provider_token != self.port.accepted.token:
            raise PortError("PROVIDER_CIRCUIT_DESYNCHRONIZED")
        v0 = self.port.accepted.voltage_V
        denominator = 2*p.L_H/dt_s**2+p.Rs_ohm/dt_s+1/(2*p.Cext_F)
        def circuit_q(v1):
            return (self.source.voltage_V+2*p.L_H*s.i_A/dt_s-s.Vc_V-.5*(v0+v1))/denominator
        evaluations = []
        def evaluate(v):
            if len(evaluations) >= 8 or not 499 <= v <= 501:
                raise PortError("BOUNDED_FULL_ROOT_LIMIT")
            trial = self.port.evaluate_trial(self.port.request(dt_s, v))
            # Explicit circuit-owned Cp in SPLIT, component total in TOTAL.
            from .ports import CapacitanceOwner
            q = trial.delta_Q_selected_C+trial.delta_q_transfer_C
            if self.port.owner is CapacitanceOwner.SPLIT:
                q += self.port.Cp*(v-v0)
            residual = q-circuit_q(v)
            row = dict(macro=s.version+1, evaluation=len(evaluations)+1,
                       owner=self.port.owner.value, Vp1_trial_V=v,
                       q_FULL_C=q, q_circuit_C=circuit_q(v), residual_C=residual,
                       residual_over_Qref=abs(residual)/2e-16,
                       native_trial_id=trial.native_trial_id)
            self.trials.append(row);evaluations.append((v,residual,trial))
            return residual, trial
        fa, ta = evaluate(499.)
        fb, tb = evaluate(501.)
        if fa*fb > 0:
            raise PortError("NO_CREDIBLE_AUTHORIZED_BRACKET")
        lo, hi = 499., 501.
        winner = None
        for _ in range(6):
            v = (lo*fb-hi*fa)/(fb-fa)
            if not lo <= v <= hi:
                raise PortError("SECANT_LEFT_AUTHORIZED_BRACKET")
            f, trial = evaluate(v)
            if abs(f)/2e-16 <= 1e-5:
                winner = (v,f,trial);break
            if fa*f <= 0:
                hi, fb = v, f
            else:
                lo, fa = v, f
        if winner is None:
            raise PortError("FULL_CHARGE_ROOT_NOT_CONVERGED")
        v1, residual, trial = winner
        m = dict(trial.native_metrics)
        q = trial.delta_q_port_C
        i1 = 2*q/dt_s-s.i_A
        vc1 = s.Vc_V+q/p.Cext_F
        im = .5*(s.i_A+i1)
        source_work = dt_s*self.source.voltage_V*im
        loss = dt_s*p.Rs_ohm*im*im
        dWcirc = .5*p.L_H*(i1+s.i_A)*(i1-s.i_A)+.5*p.Cext_F*(vc1+s.Vc_V)*(vc1-s.Vc_V)
        row = dict(owner=self.port.owner.value, macro=s.version+1,
            time_n_s=s.time_s, time_n1_s=s.time_s+dt_s, dt_s=dt_s,
            Vp0_V=v0, Vp1_V=v1, i0_A=s.i_A, i1_A=i1, Vc0_V=s.Vc_V, Vc1_V=vc1,
            q_port_C=q, charge_residual_C=residual, residual_over_Qref=abs(residual)/2e-16,
            residual_over_actual_charge=abs(residual)/abs(q) if q != 0 else None,
            residual_over_unsigned_transport=abs(residual)/m["unsigned_native_transport_C"] if m["unsigned_native_transport_C"] != 0 else None,
            native_trials=len(evaluations), native_substeps=m["substeps"],
            source_work_J=source_work, Rs_loss_J=loss, circuit_storage_change_J=dWcirc,
            Wport_circuit_J=.5*(v0+v1)*q, Wport_native_J=m["Wport_native_J"],
            circuit_energy_residual_J=source_work-loss-dWcirc-.5*(v0+v1)*q,
            DeltaU_E_J=m["DeltaU_E_J"], Wcond_face_J=m["Wcond_face_J"],
            Wcorrection_J=m["Wcorrection_J"], Wouter_J=m["Wouter_J"],
            native_electrical_energy_residual_J=m["energy_residual_J"],
            component_U_n1_J=trial.U_n1_J, circuit_Cp_U_n1_J=trial.external_Cp_U_n1_J,
            native_physical_hash=trial.candidate.physical_hash, native_history_hash=trial.candidate.history_hash)
        # Prepare both candidates before the unique commit; no PDE reevaluation here.
        candidate = LocalCircuitState(s.version+1, "PENDING_PROVIDER_COMMIT", s.time_s+dt_s,
            i1, vc1, s.cumulative_source_work_J+source_work, s.cumulative_Rs_loss_J+loss)
        accepted = self.port.commit(trial, converged=True)
        from dataclasses import replace
        self.state = replace(candidate, provider_token=accepted.token)
        self.records.append(row)
        return row


class IncrementalLocalFullCircuit(LocalNumericalFullCircuit):
    """Opt-in weak-signal solver. Commit KVL increments, not charge-derived i."""
    def advance(self, dt_s=20e-15):
        import math
        from dataclasses import replace
        if not self.port.stable_evaluation or self.state.version>=3 or dt_s!=20e-15:
            raise PortError("INCREMENTAL_LOCAL_SCOPE")
        s,p=self.state,self.p
        if s.provider_token!=self.port.accepted.token:
            raise PortError("PROVIDER_CIRCUIT_DESYNCHRONIZED")
        v0=self.port.accepted.voltage_V;vs=self.source.voltage_V
        alpha=p.L_H/dt_s+p.Rs_ohm/2+dt_s/(4*p.Cext_F)
        def di_for(dv):
            return math.fsum([vs,-p.Rs_ohm*s.i_A,-s.Vc_V,-v0,-dv/2,-dt_s*s.i_A/(2*p.Cext_F)])/alpha
        evaluations=[]
        def evaluate(v1):
            if len(evaluations)>=12 or not 499<=v1<=501:
                raise PortError("BOUNDED_INCREMENT_ROOT_LIMIT")
            dv=v1-v0;di=di_for(dv)
            t=self.port.evaluate_trial(self.port.request(dt_s,v1));m=dict(t.native_metrics)
            rq=math.fsum([self.port.Cp*dv,-dt_s*s.i_A,m['stable_q_pl_C'],-dt_s*di/2])
            self.trials.append(dict(owner=self.port.owner.value,macro=s.version+1,
                evaluation=len(evaluations)+1,Vp1_trial_V=v1,delta_V_represented_V=dv,
                delta_i_A=di,residual_q_C=rq,native_trial_id=t.native_trial_id))
            evaluations.append((v1,rq,t,di,m));return rq
        fa=evaluate(499.);fb=evaluate(501.)
        if fa*fb>0:raise PortError('NO_CREDIBLE_AUTHORIZED_BRACKET')
        lo,hi=499.,501.;winner=None
        for _ in range(10):
            v=lo-fa*(hi-lo)/(fb-fa)
            f=evaluate(v);last=evaluations[-1];m=last[4]
            if abs(f)<=1e-7*m['unsigned_native_transport_C']:
                winner=last;break
            if f*fa<=0:hi,fb=v,f
            else:lo,fa=v,f
        if winner is None:raise PortError('INCREMENTAL_CHARGE_ROOT_PRECISION_LIMIT')
        v1,rq,t,di,m=winner;dv=v1-v0;im=s.i_A+di/2;qc=dt_s*im
        i1=s.i_A+di;vc1=s.Vc_V+qc/p.Cext_F;vm=v0+dv/2
        rv=math.fsum([p.L_H*di/dt_s,p.Rs_ohm*im,s.Vc_V+qc/(2*p.Cext_F),vm,-vs])
        source=dt_s*vs*im;loss=dt_s*p.Rs_ohm*im*im
        # Small increments are kept before adding to the accepted large state.
        dw=p.L_H*(s.i_A+di/2)*di+(s.Vc_V+qc/(2*p.Cext_F))*qc
        wp=vm*t.delta_q_port_C
        energy=math.fsum([source,-loss,-dw,-wp])
        identity=-vm*rq-dt_s*im*rv
        row=dict(owner=self.port.owner.value,macro=s.version+1,time_n_s=s.time_s,
            Vp0_V=v0,Vp1_V=v1,i0_A=s.i_A,i1_A=i1,Vc0_V=s.Vc_V,Vc1_V=vc1,
            delta_i_A=di,delta_V_V=dv,q_port_C=t.delta_q_port_C,q_circuit_C=qc,
            R_q_C=rq,R_V_V=rv,circuit_energy_residual_J=energy,
            energy_residual_from_constraints_J=identity,energy_identity_difference_J=energy-identity,
            R_q_over_Qref=abs(rq)/2e-16,R_q_over_Qtr=abs(rq)/m['unsigned_native_transport_C'],
            circuit_energy_over_signed_work=abs(energy)/abs(m['stable_Wcond_J']),
            native_charge_over_Qtr=abs(m['stable_charge_residual_C'])/m['unsigned_native_transport_C'],
            native_energy_over_signed_work=abs(m['stable_energy_residual_J'])/abs(m['stable_Wcond_J']),
            native_energy_over_abs_work=abs(m['stable_energy_residual_J'])/m['stable_Wcond_abs_J'],
            Wport_J=wp,source_work_J=source,Rs_loss_J=loss,DeltaW_circuit_J=dw,
            Wcond_signed_J=m['stable_Wcond_J'],Wcond_abs_J=m['stable_Wcond_abs_J'],
            native_energy_residual_J=m['stable_energy_residual_J'],
            native_charge_residual_C=m['stable_charge_residual_C'],
            native_trials=len(evaluations),native_physical_hash=t.candidate.physical_hash,
            native_history_hash=t.candidate.history_hash,
            component_U_n1_J=t.U_n1_J,circuit_Cp_U_n1_J=t.external_Cp_U_n1_J)
        if abs(rv)>1e-6 or row['R_q_over_Qtr']>1e-3 or row['circuit_energy_over_signed_work']>1e-3:
            raise PortError('INCREMENTAL_OUTER_ACCEPTANCE_FAILED:'+str(row))
        if max(row['native_charge_over_Qtr'],row['native_energy_over_signed_work'],row['native_energy_over_abs_work'])>1e-3:
            raise PortError('STABLE_WEAK_PORT_ACCEPTANCE_FAILED:'+str(row))
        candidate=LocalCircuitState(s.version+1,'PENDING_PROVIDER_COMMIT',s.time_s+dt_s,
            i1,vc1,s.cumulative_source_work_J+source,s.cumulative_Rs_loss_J+loss)
        accepted=self.port.commit(t,converged=True)
        self.state=replace(candidate,provider_token=accepted.token)
        self.records.append(row);return row
