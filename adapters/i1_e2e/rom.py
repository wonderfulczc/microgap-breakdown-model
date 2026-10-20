"""Conservative axial species ROM, same-network Poisson charge and FP3 trials.

This is not a reduction certified against 2-D FULL. Native electron coefficients
are reused through the prior parity-checked Python interface; ion data are input.
"""
from dataclasses import replace
import time
import numpy as np
from adapters.i1_conductive_channel.kinetics import Gas,transport,sources,E_CHARGE
from adapters.i1_fp3_port_circuit.ports import IdealCapacitorPort,CapacitanceOwner,PortError,TrialResult
from adapters.i1_teng_source.circuit import TENGSeriesCircuit

EPS=8.8541878128e-12


class AxialParticlePort(IdealCapacitorPort):
    def __init__(self,c,geometry,v0,deadline):
        self.c=c;self.deadline=deadline;self.n=c['cells']
        if type(self.n) is not int or not 3<=self.n<=256:raise ValueError('ROM_CELL_COUNT')
        self.length=geometry.length;self.dz=self.length/self.n
        self.z=(np.arange(self.n)+.5)*self.dz
        factor=np.asarray(c.get('area_factors',[1.]*self.n),float)
        if factor.shape!=(self.n,) or np.any(factor<=0) or not np.all(np.isfinite(factor)):
            raise ValueError('EFFECTIVE_AREA_REQUIRED')
        self.area=geometry.area*factor;self.volume=self.area*self.dz
        self.af=np.r_[self.area[0],2*self.area[:-1]*self.area[1:]/(self.area[:-1]+self.area[1:]),self.area[-1]]
        self.distance=np.r_[self.dz/2,np.full(self.n-1,self.dz),self.dz/2]
        self.links=EPS*self.af/self.distance
        self.K=np.diag(self.links[:-1]+self.links[1:])+np.diag(-self.links[1:-1],1)+np.diag(-self.links[1:-1],-1)
        self.Kinv=np.linalg.inv(self.K)
        b=np.zeros(self.n);b[0]=self.links[0]
        self.psi=self.Kinv@b;cp=1/np.sum(1/self.links)
        self.gas=Gas(c['pressure_Pa'],c['temperature_K'])
        for key in ('ion_mobility_m2_V_s','ion_diffusion_m2_s','gamma_i'):
            if key not in c or not np.isfinite(c[key]) or c[key]<0:raise ValueError('EXPLICIT_ION_SURFACE_PARAMETER_REQUIRED')
        if c['parameter_identity']!='SYNTHETIC_ION_SURFACE_NOT_DEVICE':raise ValueError('PARAMETER_SCOPE')
        self.mu_p=c['ion_mobility_m2_V_s'];self.Dp=c['ion_diffusion_m2_s'];self.gamma=c['gamma_i']
        self.state=np.array(c['initial_species_m3'],float)
        if self.state.shape==(3,):self.state=np.tile(self.state[:,None],(1,self.n))
        if self.state.shape!=(3,self.n) or not np.all(np.isfinite(self.state)) or np.any(self.state<0):
            raise ValueError('INITIAL_SPECIES')
        if not np.array_equal(self.state[1],self.state[0]+self.state[2]):raise ValueError('EXPLICIT_NEUTRAL_INITIAL_STATE')
        super().__init__(cp,owner=CapacitanceOwner.TOTAL,voltage0_V=v0,
            geometry_id='AXIAL_EFFECTIVE_GEOMETRY:'+str(self.length),operator_id='ROM_AXIAL_FV_POISSON_V1')
        self.descriptor=replace(self.descriptor,provider_id='PHYSICS_INFORMED_CHANNEL_ROM',
            capabilities=('local_rom_charge','conservative_species_flux','field_energy','transactions'),
            validity_domain='POISSON_1D_APPROXIMATION_NO_MATERIAL_ENERGY_NO_SP3_NO_NEGATIVE_ION_TRANSPORT',
            source_identity='NATIVE_PARITY_ELECTRON_COEFFICIENTS_SYNTHETIC_ION_PARAMETERS')
        self.cache={};self.active=None;self.last=None;self.cell_audits=[]
        self.trial_substeps=0;self.accepted_substeps=0
        self.energy=self.field(self.state,v0)[2]
        self.sum_Wcond=0.;self.sum_Wport_native=0.;self.sum_time_work_difference=0.

    def field(self,y,V):
        q=E_CHARGE*(y[1]-y[0]-y[2])*self.volume
        b=q.copy();b[0]+=self.links[0]*V
        phi=np.r_[V,self.Kinv@b,0.]
        ef=(phi[:-1]-phi[1:])/self.distance
        energy=.5*np.dot(self.links,(phi[:-1]-phi[1:])**2)
        return phi,ef,float(energy)

    def flux(self,n,v,D):
        f=np.zeros(self.n+1)
        f[1:-1]=self.af[1:-1]*(v[1:-1]*np.where(v[1:-1]>=0,n[:-1],n[1:])-
                                      D[1:-1]*(n[1:]-n[:-1])/self.dz)
        f[0]=-self.af[0]*(max(-v[0],0)+D[0]/self.distance[0])*n[0]
        f[-1]=self.af[-1]*(max(v[-1],0)+D[-1]/self.distance[-1])*n[-1]
        outleft=self.af[:-1]*(np.maximum(-v[:-1],0)+D[:-1]/self.distance[:-1])
        outright=self.af[1:]*(np.maximum(v[1:],0)+D[1:]/self.distance[1:])
        return f,(outleft+outright)/self.volume

    def rhs(self,y,V):
        phi,E,U=self.field(y,V);ec=.5*(E[:-1]+E[1:])
        coeff=[transport(float(e),self.gas) for e in ec]
        mu=np.array([a['mobility'] for a in coeff]);D=np.array([a['diffusion'] for a in coeff])
        muf=np.r_[mu[0],.5*(mu[:-1]+mu[1:]),mu[-1]]
        df=np.r_[D[0],.5*(D[:-1]+D[1:]),D[-1]]
        fe,le=self.flux(y[0],-muf*E,df)
        fp,lp=self.flux(y[1],self.mu_p*E,np.full(self.n+1,self.Dp))
        emitted=np.zeros(2)
        if V>0:emitted[1]=self.gamma*max(fp[-1],0);fe[-1]-=emitted[1]
        elif V<0:emitted[0]=self.gamma*max(-fp[0],0);fe[0]+=emitted[0]
        S=np.empty_like(y);mechanisms=np.empty((6,self.n))
        for j,a in enumerate(coeff):S[:,j],mechanisms[:,j]=sources(y[:,j],a,self.gas,0.)
        rhs=S.copy();rhs[0]-=np.diff(fe)/self.volume;rhs[1]-=np.diff(fp)/self.volume
        rates=np.zeros_like(y)
        # Actual removal-rate bounds, without a density floor.
        for k in range(3):
            nz=y[k]>0
            rates[k,nz]=np.maximum(-S[k,nz],0)/y[k,nz]
        rates[0]+=le;rates[1]+=lp
        rate=max(float(rates.max()),max(a['ionization_frequency'] for a in coeff))
        return rhs,fe,fp,emitted,mechanisms,phi,E,U,rate,mu

    def integrate(self,request):
        y=self.state.copy();elapsed=0.;dt=request.dt_s;v0=request.voltage_n_V;v1=request.voltage_n1_V
        qHg=0.;wc=0.;wp=0.;reaction_work=0.;max_rel=0.;audits=[];sub=0
        qfaces=np.zeros(self.n+1);species_boundary=np.zeros((2,2));created=np.zeros(6);emission=np.zeros(2)
        while elapsed<dt:
            if time.monotonic()>=self.deadline:raise TimeoutError('ROM_RESOURCE_LIMIT')
            if sub>=self.c['max_substeps_per_trial']:raise ValueError('ROM_STIFFNESS_LIMIT_NO_MODEL_ADJUSTMENT')
            v=v0+(v1-v0)*elapsed/dt
            rhs,fe,fp,emit,mechanisms,phi,E,U,rate,mu=self.rhs(y,v)
            step=min(dt-elapsed,.15/rate if rate>0 else dt)
            if step<dt*1e-12:raise ValueError('ROM_UNRESOLVED_TIMESTEP')
            vn=v0+(v1-v0)*min(elapsed+step,dt)/dt
            yn=y+step*rhs
            if np.any(yn<0) or not np.all(np.isfinite(yn)):raise ValueError('ROM_NEGATIVITY_NO_FLOOR')
            phin,_,Un=self.field(yn,vn)
            dq=E_CHARGE*((yn[1]-y[1])-(yn[0]-y[0])-(yn[2]-y[2]))*self.volume
            faceq=E_CHARGE*(fp-fe)*step
            # Physical reaction stoichiometry independently records roundoff defect.
            netS=np.array([sources(y[:,j],transport(float(.5*(E[j]+E[j+1])),self.gas),self.gas)[0] for j in range(self.n)]).T
            qr=E_CHARGE*(netS[1]-netS[0]-netS[2])*self.volume*step
            residual=dq+np.diff(faceq)-qr
            scale=E_CHARGE*np.sum(abs(yn-y),axis=0)*self.volume+abs(faceq[:-1])+abs(faceq[1:])+abs(qr)
            nz=scale>0
            rel=float(np.max(abs(residual[nz])/scale[nz])) if np.any(nz) else 0.
            if np.any(residual[~nz]!=0) or rel>1e-8:raise ValueError('ROM_CELL_CHARGE_CONTINUITY')
            max_rel=max(max_rel,rel)
            dqh=self.Cp*(vn-v)-np.dot(self.psi,dq)
            portq=dqh+faceq[0]
            bar=.5*(phi+phin)
            wc+=float(np.dot(faceq,bar[:-1]-bar[1:]))
            wp+=.5*(v+vn)*portq
            reaction_work+=float(np.dot(bar[1:-1],qr))
            qHg+=faceq[0];qfaces+=faceq
            species_boundary+=step*np.array([[-fe[0],fe[-1]],[-fp[0],fp[-1]]])
            emission+=step*emit
            created+=step*(mechanisms@self.volume)
            audits.append(dict(substep=sub,dt_s=step,max_cell_relative=rel,
                global_residual_C=float(np.sum(residual)),L1_residual_C=float(np.sum(abs(residual))),
                unsigned_cell_exchange_C=float(np.sum(scale))))
            y=yn;elapsed+=step;sub+=1;self.trial_substeps+=1
        dq=E_CHARGE*((y[1]-self.state[1])-(y[0]-self.state[0])-(y[2]-self.state[2]))*self.volume
        dqH=self.Cp*(v1-v0)-float(np.dot(self.psi,dq))
        q=dqH+qHg;u1=self.field(y,v1)[2]
        return dict(y=y,dqH=dqH,qHg=qHg,q=q,U=u1,Wcond=wc,Wport_native=wp,
            reaction_work=reaction_work,native_energy_residual=wp+reaction_work-(u1-self.energy)-wc,
            max_cell_relative=max_rel,substeps=sub,audits=audits,faceq=qfaces,
            species_boundary=species_boundary,emission=emission,mechanisms=created)

    def evaluate_trial(self,request):
        base=super().evaluate_trial(request)
        key=(request.dt_s,request.voltage_n1_V)
        if key not in self.cache:
            info=self.integrate(request)
            h=1e-5*max(1.,abs(request.voltage_n1_V))
            shifted=self.integrate(replace(request,voltage_n1_V=request.voltage_n1_V+h))
            derivative=(shifted['q']-info['q'])/h
            vm=.5*(request.voltage_n_V+request.voltage_n1_V)
            work=vm*info['q']
            candidate=replace(base.candidate,cumulative_charge_C=self.accepted.cumulative_charge_C+info['q'],
                cumulative_port_work_J=self.accepted.cumulative_port_work_J+work,
                cumulative_component_work_J=self.accepted.cumulative_component_work_J+work,
                cumulative_dissipation_J=self.accepted.cumulative_dissipation_J+info['Wcond'])
            result=replace(base,candidate=candidate,delta_Q_selected_C=info['dqH'],delta_q_external_Cp_C=0.,
                delta_q_transfer_C=info['qHg'],delta_q_port_C=info['q'],
                U_n_J=self.energy,U_n1_J=info['U'],external_Cp_U_n_J=0.,external_Cp_U_n1_J=0.,
                port_work_J=work,component_work_J=work,resistor_dissipation_J=info['Wcond'],
                transfer_current_mid_A=info['qHg']/request.dt_s,d_deltaq_port_d_voltage_n1_F=derivative,
                component_energy_residual_J=work-(info['U']-self.energy)-info['Wcond'],
                total_energy_residual_J=work-(info['U']-self.energy)-info['Wcond'],
                quadrature_id='ROM_EXPLICIT_FV_SUBCYCLES_AND_OUTER_MIDPOINT',
                provenance='SIGNED_FACE_ELECTRICAL_WORK_NOT_RESISTOR_HEAT')
            self.cache[key]=(result,info)
        result,info=self.cache[key];self.active=(result,info)
        return result

    def commit(self,result,*,converged):
        self.validate_commit(result,converged=converged)
        info=self.active[1]
        self.state=info['y'].copy();self.energy=info['U'];self._accepted=result.candidate
        self.last=info;self.accepted_substeps+=info['substeps']
        self.sum_Wcond+=info['Wcond'];self.sum_Wport_native+=info['Wport_native']
        self.sum_time_work_difference+=result.port_work_J-info['Wport_native']
        self.cell_audits.extend(dict(accepted_step=result.candidate.version,**r) for r in info['audits'])
        self.cache={};return self._accepted

    def summary(self):
        phi,E,U=self.field(self.state,self.accepted.voltage_V)
        mu=np.array([transport(float(e),self.gas)['mobility'] for e in .5*(E[:-1]+E[1:])])
        sigma=E_CHARGE*(mu*self.state[0]+self.mu_p*self.state[1])
        resistance=float(np.sum(self.dz/(sigma*self.area))) if np.all(sigma>0) else None
        last=self.last
        return dict(ne_mean_m3=float(np.dot(self.state[0],self.volume)/self.volume.sum()),
            np_mean_m3=float(np.dot(self.state[1],self.volume)/self.volume.sum()),
            nn_mean_m3=float(np.dot(self.state[2],self.volume)/self.volume.sum()),
            total_electrons=float(np.dot(self.state[0],self.volume)),
            total_positive_ions=float(np.dot(self.state[1],self.volume)),
            net_gas_charge_C=float(E_CHARGE*np.dot(self.state[1]-self.state[0]-self.state[2],self.volume)),
            G_channel_S=1/resistance if resistance is not None else 0.,R_channel_ohm=resistance,
            G_scope='FROZEN_DRIFT_COLUMN_NOT_PORT_R',E_min_V_m=float(E.min()),E_max_V_m=float(E.max()),
            field_energy_J=U,Wcond_native_cumulative_J=self.sum_Wcond,
            Wport_native_cumulative_J=self.sum_Wport_native,
            macro_substep_work_difference_J=self.sum_time_work_difference,
            native_energy_residual_J=last['native_energy_residual'],
            max_cell_charge_relative=last['max_cell_relative'],substeps=last['substeps'],
            impact_events=last['mechanisms'][0],attachment2_events=last['mechanisms'][1],
            attachment3_events=last['mechanisms'][2],electron_recombination_events=last['mechanisms'][3],
            ion_recombination_events=last['mechanisms'][4],photo_events=last['mechanisms'][5],
            H_electron_net_collection=last['species_boundary'][0,0],
            Return_electron_net_collection=last['species_boundary'][0,1],
            H_ion_collection=last['species_boundary'][1,0],Return_ion_collection=last['species_boundary'][1,1],
            secondary_emitted_count=float(sum(last['emission'])),
            accepted_substeps=self.accepted_substeps,trial_substeps=self.trial_substeps,
            model='POISSON_1D_APPROXIMATION_NOT_FULL',material_energy='UNAVAILABLE')


class ROMSeriesCircuit(TENGSeriesCircuit):
    """Reuse original residual/Newton/commit; replace only total field inventory."""
    def total_energy(self,i=None,vc=None,vp=None,x=None,q_T_C=None):
        base=super().total_energy(i,vc,vp,x,q_T_C)
        if vp is None:
            return base+self.port.energy-.5*self.port.Cp*self.port.accepted.voltage_V**2
        result,info=self.port.active
        if result.candidate.voltage_V!=vp:raise PortError('TRIAL_ENERGY_IDENTITY')
        return base+info['U']-.5*self.port.Cp*vp**2

    def commit(self,trial):
        row=super().commit(trial)
        last=self.port.last
        row.update(Icond_H_mid_A=last['qHg']/trial.dt_s,
            Idisplacement_H_mid_A=last['dqH']/trial.dt_s,
            Idisplacement_space_charge_mid_A=(last['dqH']-self.port.Cp*(
                trial.port_trial.input.voltage_n1_V-trial.port_trial.input.voltage_n_V))/trial.dt_s,
            field_energy_J=last['U'],Wcond_face_J=last['Wcond'],Wport_native_J=last['Wport_native'],
            Wport_circuit_minus_native_J=trial.port_trial.port_work_J-last['Wport_native'],
            reaction_charge_defect_work_J=last['reaction_work'],
            electrical_energy_native_residual_J=last['native_energy_residual'],
            W_column_identity='Compatibility column is signed face work, not gas heat or lumped resistor loss')
        return row
