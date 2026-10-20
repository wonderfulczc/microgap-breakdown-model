"""Cext shunt leakage on the original FP3 Newton/port transaction scheduler."""
import numpy as np
from adapters.i1_fp3_port_circuit.circuit import SeriesPortCircuit, CircuitAccepted
from adapters.i1_fp3_port_circuit.ports import CapacitanceOwner, PortError


class AmplifierSeriesCircuit(SeriesPortCircuit):
    def __init__(self,params,port,source,*,Rbleed_ohm,i0_A,Vc0_V,**kwargs):
        if not np.isfinite(Rbleed_ohm) or Rbleed_ohm<=0:
            raise PortError('RBLEED_REQUIRED_POSITIVE_NOT_RS')
        self.Rbleed=Rbleed_ohm
        self.losses=dict(Rs=0.,Rout=0.,bleed=0.)
        super().__init__(params,port,source,i0_A=i0_A,Vc0_V=Vc0_V,**kwargs)

    def total_energy(self,i=None,vc=None,vp=None,x=None):
        current=self.accepted.i_A if i is None else i
        return super().total_energy(i,vc,vp,x)+.5*self.source.Lout*current**2

    def _residual(self,dt,i1,vp1):
        s,c=self.accepted,self.params
        if s.port_token!=self.port.accepted.token: raise PortError('COMPONENT_CIRCUIT_DESYNCHRONIZED')
        self.trial_evaluations+=1
        p=self.port.evaluate_trial(self.port.request(dt,vp1))
        q=p.delta_Q_selected_C+p.delta_q_transfer_C
        if self.port.descriptor.capacitance_owner is CapacitanceOwner.SPLIT:
            q+=self.port.Cp*(vp1-self.port.accepted.voltage_V)
        im=.5*(s.i_A+i1)
        a=dt/(2*c.Cext_F*self.Rbleed)
        vc1=(s.Vc_V*(1-a)+dt*im/c.Cext_F)/(1+a)
        v=self.source.interval_voltage(s.time_s,s.time_s+dt)
        r=np.array([((c.L_H+self.source.Lout)*(i1-s.i_A)/dt+(c.Rs_ohm+self.source.Rout)*im+
                     .5*(s.Vc_V+vc1)+.5*(self.port.accepted.voltage_V+vp1)-v)/self.scales.V,
                    (im-q/dt)/self.scales.I])
        jac=np.array([[((c.L_H+self.source.Lout)/dt+.5*(c.Rs_ohm+self.source.Rout)+
                       dt/(4*c.Cext_F*(1+a)))/self.scales.V,.5/self.scales.V],
                      [.5/self.scales.I,-p.d_deltaq_port_d_voltage_n1_F/dt/self.scales.I]])
        return r,jac,p,vc1

    def commit(self,trial):
        if trial.base_token!=self.accepted.token: raise PortError('STALE_OR_FOREIGN_CIRCUIT_TRIAL')
        r,_,p,vc1=self._residual(trial.dt_s,trial.i1_A,trial.port_trial.candidate.voltage_V)
        if max(abs(r))>1e-10 or p!=trial.port_trial or vc1!=trial.Vc1_V:
            raise PortError('UNCONVERGED_OR_ALTERED_AMPLIFIER_TRIAL')
        self.port.validate_commit(p,converged=True)
        s,c=self.accepted,self.params
        dt=trial.dt_s;im=.5*(s.i_A+trial.i1_A);vcm=.5*(s.Vc_V+vc1)
        v=self.source.interval_voltage(s.time_s,s.time_s+dt)
        ws=dt*v*im
        qr=dt*c.Rs_ohm*im**2;qo=dt*self.source.Rout*im**2;qb=dt*vcm**2/self.Rbleed
        e0=self.total_energy();e1=self.total_energy(trial.i1_A,vc1,p.candidate.voltage_V)
        losses={k:self.losses[k]+q for k,q in [('Rs',qr),('Rout',qo),('bleed',qb)]}
        candidate=CircuitAccepted(self._token(s.version+1),s.version+1,p.candidate.token,p.candidate.time_s,
            trial.i1_A,vc1,s.source_work_J+ws,s.series_loss_J+qr+qo+qb)
        row=dict(time_s=candidate.time_s,dt_s=dt,i_A=candidate.i_A,Vc_V=vc1,Vp_V=p.candidate.voltage_V,
            V_T_V=v-self.source.Rout*im-
                  self.source.Lout*(candidate.i_A-s.i_A)/dt,
            Vamp_endpoint_V=self.source.voltage(candidate.time_s),
            q_T_C=None,Q_Cext_C=c.Cext_F*vc1,Q_Cp_C=self.port.Cp*p.candidate.voltage_V,
            Icol_mid_A=p.transfer_current_mid_A,I_port_mid_A=im,I_bleed_mid_A=vcm/self.Rbleed,
            I_Cext_mid_A=c.Cext_F*(vc1-s.Vc_V)/dt,
            I_Cp_mid_A=self.port.Cp*(p.candidate.voltage_V-self.port.accepted.voltage_V)/dt,
            delta_q_port_C=p.delta_q_port_C,delta_q_column_C=p.delta_q_transfer_C,
            delta_Q_selected_C=p.delta_Q_selected_C,delta_q_external_Cp_C=p.delta_q_external_Cp_C,
            E_total_J=e1,delta_E_total_J=e1-e0,E_source_coordinate_J=None,
            E_Cp_J=.5*self.port.Cp*p.candidate.voltage_V**2,
            E_Lout_J=.5*self.source.Lout*candidate.i_A**2,
            E_component_J=p.U_n1_J,E_external_Cp_J=p.external_Cp_U_n1_J,
            W_oc_J=ws,W_Rs_J=qr,W_Rout_J=qo,W_bleed_J=qb,
            W_terminal_source_J=dt*im*(v-self.source.Rout*im-self.source.Lout*(candidate.i_A-s.i_A)/dt),
            W_column_J=p.resistor_dissipation_J,W_port_J=p.port_work_J,
            Cext_charge_residual_C=c.Cext_F*(vc1-s.Vc_V)+dt*vcm/self.Rbleed-dt*im,
            port_charge_residual_C=dt*im-p.delta_q_port_C,kvl_residual_V=float(r[0]*self.scales.V),
            energy_residual_J=e1-e0+qr+qo+qb+p.resistor_dissipation_J-ws,
            W_oc_cumulative_J=candidate.source_work_J,W_Rs_cumulative_J=losses['Rs'],
            W_Rout_cumulative_J=losses['Rout'],W_bleed_cumulative_J=losses['bleed'],
            W_column_cumulative_J=p.candidate.cumulative_dissipation_J,accepted_version=candidate.version,
            quadrature_id=p.quadrature_id,source_coordinate_status='NOT_APPLICABLE_AMPLIFIER',
            voltage_time_identity='V_T is interval midpoint including output R/L; Vamp_endpoint_V is endpoint')
        if not all(np.isfinite(v) for v in row.values() if isinstance(v,(float,int))):
            raise PortError('NONFINITE_AMPLIFIER_RESULT')
        self.port.commit(p,converged=True)
        self._accepted=candidate;self.losses=losses;self.records.append(row)
        return row
