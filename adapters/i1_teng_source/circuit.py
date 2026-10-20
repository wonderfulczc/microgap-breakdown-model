"""Incremental frozen-source extension of the accepted FP3 series protocol.

No thermal/full provider is imported. The parent bounded Newton scheduler and
immutable trial dataclasses are reused; the residual and source energy differ.
"""
from dataclasses import replace
import numpy as np

from adapters.i1_fp3_port_circuit.circuit import SeriesPortCircuit, CircuitAccepted
from adapters.i1_fp3_port_circuit.ports import IdealCapacitorPort, PortError, CapacitanceOwner


class ReferenceParallelRCPort(IdealCapacitorPort):
    """Mathematical constant resistor in parallel with a synthetic Cp."""
    def __init__(self, Cp_F, R_test, **kwargs):
        if not np.isfinite(R_test) or R_test<=0:
            raise PortError('INVALID_REFERENCE_RESISTOR')
        super().__init__(Cp_F,**kwargs)
        self.R_test=R_test
        self.descriptor=replace(self.descriptor,provider_id='REFERENCE_FROZEN_R_PARALLEL_CP',
                                model_identity=(('R_test_ohm',R_test),('Cp_F',Cp_F)))

    def evaluate_trial(self,request):
        result=super().evaluate_trial(request)
        vm=.5*(request.voltage_n_V+request.voltage_n1_V)
        dq=request.dt_s*vm/self.R_test
        work=vm*dq
        candidate=replace(result.candidate,
            cumulative_charge_C=result.candidate.cumulative_charge_C+dq,
            cumulative_port_work_J=result.candidate.cumulative_port_work_J+work,
            cumulative_component_work_J=result.candidate.cumulative_component_work_J+work,
            cumulative_dissipation_J=result.candidate.cumulative_dissipation_J+work)
        port_work=result.port_work_J+work
        component_work=result.component_work_J+work
        return replace(result,candidate=candidate,delta_q_transfer_C=dq,
            delta_q_port_C=result.delta_q_port_C+dq,
            port_work_J=port_work,component_work_J=component_work,
            resistor_dissipation_J=work,transfer_current_mid_A=vm/self.R_test,
            d_deltaq_port_d_voltage_n1_F=result.d_deltaq_port_d_voltage_n1_F+request.dt_s/(2*self.R_test),
            component_energy_residual_J=component_work-(result.U_n1_J-result.U_n_J)-work,
            total_energy_residual_J=port_work-(result.U_n1_J+result.external_Cp_U_n1_J-
                                            result.U_n_J-result.external_Cp_U_n_J)-work,
            provenance='REFERENCE_FROZEN_R_NOT_MEASURED_GAP')


class TENGSeriesCircuit(SeriesPortCircuit):
    def __init__(self,params,port,source,*,i0_A,Vc0_V,**kwargs):
        self._q_T_C=source.q_T0_C
        super().__init__(params,port,source,i0_A=i0_A,Vc0_V=Vc0_V,**kwargs)

    @property
    def q_T_C(self):
        return self._q_T_C

    def total_energy(self,i=None,vc=None,vp=None,x=None,q_T_C=None):
        q=self._q_T_C if q_T_C is None else q_T_C
        return super().total_energy(i,vc,vp,x)+self.source.coordinate_energy(q)

    def _residual(self,dt,i1,vp1):
        s,c=self.accepted,self.params
        if s.port_token!=self.port.accepted.token:
            raise PortError('COMPONENT_CIRCUIT_DESYNCHRONIZED')
        self.trial_evaluations+=1
        trial=self.port.evaluate_trial(self.port.request(dt,vp1))
        charge=trial.delta_Q_selected_C+trial.delta_q_transfer_C
        if self.port.descriptor.capacitance_owner is CapacitanceOwner.SPLIT:
            charge+=self.port.Cp*(vp1-self.port.accepted.voltage_V)
        im=.5*(s.i_A+i1)
        vc1=s.Vc_V+dt*im/c.Cext_F
        source_bar=self.source.voltage(self._q_T_C+.5*dt*im)
        r=np.array([(c.L_H*(i1-s.i_A)/dt+c.Rs_ohm*im+.5*(s.Vc_V+vc1)
                     +.5*(self.port.accepted.voltage_V+vp1)-source_bar)/self.scales.V,
                    (im-charge/dt)/self.scales.I])
        jac=np.array([[(c.L_H/dt+.5*c.Rs_ohm+dt/(4*c.Cext_F)+
                        dt*self.source.inverse_capacitance/4)/self.scales.V,.5/self.scales.V],
                      [.5/self.scales.I,-trial.d_deltaq_port_d_voltage_n1_F/dt/self.scales.I]])
        return r,jac,trial,vc1

    def commit(self,trial):
        if trial.base_token!=self.accepted.token:
            raise PortError('STALE_OR_FOREIGN_CIRCUIT_TRIAL')
        residual,_,verified,vc1=self._residual(trial.dt_s,trial.i1_A,trial.port_trial.candidate.voltage_V)
        if (np.max(abs(residual))>1e-10 or verified!=trial.port_trial or vc1!=trial.Vc1_V):
            raise PortError('UNCONVERGED_OR_ALTERED_SOURCE_TRIAL')
        self.port.validate_commit(verified,converged=True)
        s,p,c=self.accepted,self.port.accepted,self.params
        im=.5*(s.i_A+trial.i1_A)
        dq=trial.dt_s*im
        q0=self._q_T_C
        q1=q0+dq
        vm=.5*(p.voltage_V+verified.candidate.voltage_V)
        Woc=self.source.V_oc_V*dq
        Wterminal=self.source.voltage(.5*(q0+q1))*dq
        Rs_loss=trial.dt_s*c.Rs_ohm*im*im
        energy0=self.total_energy()
        energy1=self.total_energy(trial.i1_A,vc1,verified.candidate.voltage_V,q_T_C=q1)
        source_energy0=self.source.coordinate_energy(q0)
        source_energy1=self.source.coordinate_energy(q1)
        candidate=CircuitAccepted(self._token(s.version+1),s.version+1,verified.candidate.token,
            verified.candidate.time_s,trial.i1_A,vc1,s.source_work_J+Woc,s.series_loss_J+Rs_loss)
        row=dict(time_s=candidate.time_s,dt_s=trial.dt_s,i_A=candidate.i_A,Vc_V=vc1,
            Vp_V=verified.candidate.voltage_V,q_T_C=q1,V_T_V=self.source.voltage(q1),
            Q_Cext_C=c.Cext_F*vc1,Q_Cp_C=self.port.Cp*verified.candidate.voltage_V,
            Icol_mid_A=verified.transfer_current_mid_A,I_port_mid_A=im,
            I_Cp_mid_A=self.port.Cp*(verified.candidate.voltage_V-p.voltage_V)/trial.dt_s,
            delta_q_port_C=verified.delta_q_port_C,delta_q_column_C=verified.delta_q_transfer_C,
            delta_Q_selected_C=verified.delta_Q_selected_C,delta_q_external_Cp_C=verified.delta_q_external_Cp_C,
            E_total_J=energy1,delta_E_total_J=energy1-energy0,
            E_source_coordinate_J=source_energy1,
            E_Cp_J=.5*self.port.Cp*verified.candidate.voltage_V**2,
            E_component_J=verified.U_n1_J,E_external_Cp_J=verified.external_Cp_U_n1_J,
            W_oc_J=Woc,W_terminal_source_J=Wterminal,W_Rs_J=Rs_loss,
            W_column_J=verified.resistor_dissipation_J,
            W_port_J=verified.port_work_J,W_component_J=verified.component_work_J,
            source_coordinate_residual_J=source_energy1-source_energy0-Woc+Wterminal,
            source_charge_residual_C=q1-q0-dq,
            Cext_charge_residual_C=c.Cext_F*(vc1-s.Vc_V)-dq,
            port_charge_residual_C=dq-verified.delta_q_port_C,
            kvl_residual_V=float(residual[0]*self.scales.V),
            energy_residual_J=energy1-energy0+Rs_loss+verified.resistor_dissipation_J-Woc,
            W_oc_cumulative_J=candidate.source_work_J,W_Rs_cumulative_J=candidate.series_loss_J,
            W_column_cumulative_J=verified.candidate.cumulative_dissipation_J,
            accepted_version=candidate.version,quadrature_id=verified.quadrature_id)
        if not all(np.isfinite(v) for v in row.values() if isinstance(v,(int,float))):
            raise PortError('NONFINITE_SOURCE_CIRCUIT_RESULT')
        # All validations precede the joint port/source/circuit commit.
        self.port.commit(verified,converged=True)
        self._q_T_C=q1
        self._accepted=candidate
        self.records.append(row)
        return row
