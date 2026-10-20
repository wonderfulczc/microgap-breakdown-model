"""SI transcription of local native MorrowLowke/ReactionModel; no spatial closure."""
from dataclasses import dataclass,asdict
import math
import numpy as np

E_CHARGE=1.602176634e-19
KB=1.380649e-23

@dataclass(frozen=True)
class Gas:
    pressure_Pa: float=101325.
    temperature_K: float=300.
    @property
    def neutral_density_m3(self):
        if not all(math.isfinite(v) and v>0 for v in [self.pressure_Pa,self.temperature_K]):
            raise ValueError('INVALID_GAS')
        return self.pressure_Pa/(KB*self.temperature_K)


def transport(E_V_m,gas=Gas()):
    if not math.isfinite(E_V_m):raise ValueError('NONFINITE_FIELD')
    N=gas.neutral_density_m3; E=abs(E_V_m)
    c=dict(mobility=0.,diffusion=0.,ionization_frequency=0.,attachment_two_body_frequency=0.,
           attachment_three_body_frequency=0.,ionization_townsend=0.,attachment_two_body_townsend=0.,
           attachment_three_body_townsend=0.,drift_speed=0.)
    if E==0:
        c['mobility']=6.87e24/N
        return c
    x=E/N*1e4; Nc=N/1e6
    alpha=2e-16*math.exp(-7.248e-15/x) if x>1.5e-15 else 6.619e-17*math.exp(-5.593e-15/x)
    eta2=max(0.,8.889e-5*x+2.567e-19 if x>1.05e-15 else 6.089e-4*x-2.893e-19)
    eta3=4.7778e-59*x**(-1.2749)
    if x>2e-15:w=7.4e21*x+7.1e6
    elif x>=1e-16:w=1.03e22*x+1.3e6
    elif x>=2.6e-17:
        x0,x1=2.6e-17,1e-16;t=(x-x0)/(x1-x0)
        y0=6.87e22*x0+3.38e4;y1=1.03e22*x1+1.3e6
        m0=6.87e22*(x1-x0);m1=1.03e22*(x1-x0)
        w=(2*t**3-3*t*t+1)*y0+(t**3-2*t*t+t)*m0+(-2*t**3+3*t*t)*y1+(t**3-t*t)*m1
    else:w=6.87e22*x+3.38e4
    c.update(drift_speed=w/100,mobility=w/100/E,diffusion=.3341e9*x**.54069*w/E*1e-2,
             ionization_townsend=alpha*Nc*100,attachment_two_body_townsend=eta2*Nc*100,
             attachment_three_body_townsend=eta3*Nc*Nc*100)
    for kind in ['ionization','attachment_two_body','attachment_three_body']:
        c[kind+'_frequency']=c[kind+'_townsend']*c['drift_speed']
    return c


def coefficients(c,gas=Gas()):
    Te=E_CHARGE*c['diffusion']/(c['mobility']*KB)
    return dict(electron_temperature_K=Te,beta_ep_m3_s=1.138e-11*max(Te,1.)**(-.7),
                beta_np_m3_s=2e-13*math.sqrt(300/gas.temperature_K))


def sources(y,c,gas=Gas(),photo_source_m3_s=0.):
    ne,np_,nn=map(float,y)
    if not all(math.isfinite(x) and x>=0 for x in [ne,np_,nn,photo_source_m3_s]):
        raise ValueError('NEGATIVE_OR_NONFINITE_STATE_OR_SOURCE')
    b=coefficients(c,gas)
    v=np.array([c['ionization_frequency']*ne,c['attachment_two_body_frequency']*ne,
        c['attachment_three_body_frequency']*ne,b['beta_ep_m3_s']*ne*np_,b['beta_np_m3_s']*nn*np_,photo_source_m3_s])
    # Columns: ionization, attach2, attach3, e-p recombination, ion recombination, photo.
    return STOICH@v,v


STOICH=np.array([[1,-1,-1,-1,0,1],[1,0,0,-1,-1,1],[0,1,1,0,-1,0]],dtype=float)


def density_jacobian(y,c,gas=Gas()):
    ne,np_,nn=y;b=coefficients(c,gas);ep=b['beta_ep_m3_s'];ip=b['beta_np_m3_s']
    grad=np.array([[c['ionization_frequency'],0,0],[c['attachment_two_body_frequency'],0,0],
        [c['attachment_three_body_frequency'],0,0],[ep*np_,ep*ne,0],[0,ip*nn,ip*np_],[0,0,0]])
    return STOICH@grad,grad


def column_demand(mu,R_ohm,radius_m,ell_m=70e-6):
    if not all(math.isfinite(x) and x>0 for x in [mu,R_ohm,radius_m,ell_m]):
        raise ValueError('INVALID_COLUMN_ASSUMPTION')
    A=math.pi*radius_m**2
    number=ell_m**2/(E_CHARGE*mu*R_ohm)
    return dict(area_m2=A,volume_m3=A*ell_m,N_e_required=number,n_e_required_m3=number/(A*ell_m))


def column_conductance(ne,E_V_m,radius_m,ell_m=70e-6):
    if not math.isfinite(ne) or ne<0:raise ValueError('INVALID_DENSITY')
    return E_CHARGE*transport(E_V_m)['mobility']*ne*math.pi*radius_m**2/ell_m


def conventional_drift_current_density(ne,E_V_m,gas=Gas()):
    return E_CHARGE*transport(E_V_m,gas)['mobility']*ne*E_V_m
