import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
import numpy as np
from scipy.integrate import quad
from streamer_rf.lifecycle import current_moment
from streamer_rf.fourier import *
P=dict(I0=.44,alpha=1.7e9,beta=7.5e7,T0=6e-9)
def test_zero_frequency_matches_integral():
 z=current_moment_transform_omega(0,**P); q=quad(lambda t:current_moment(t,**P),-1e-6,1e-6,points=[P['T0']],epsabs=1e-13)[0]; assert np.isclose(z.real,q,rtol=1e-8)
def test_conjugate_symmetry():
 w=2*np.pi*1e9; assert np.allclose(current_moment_transform_omega(-w,**P),np.conj(current_moment_transform_omega(w,**P)),rtol=1e-13)
def test_quad_validation_representative():
 w=2*np.pi*3e8;q=P['alpha']+P['beta'];a=P['alpha']/q;k=w/q;L=-max(80,38/a);R=max(80,38/(1-a));base=lambda y:np.exp(a*y-np.logaddexp(0,y));re=quad(base,L,R,weight='cos',wvar=k,epsabs=2e-12,limit=1000)[0];im=-quad(base,L,R,weight='sin',wvar=k,epsabs=2e-12,limit=1000)[0];z=P['I0']/q*(re+1j*im)*np.exp(-1j*w*P['T0']);assert abs(current_moment_transform_omega(w,**P)-z)/abs(z)<1e-7
def test_equation5():
 w=2*np.pi*1e9; assert np.isclose(derivative_spectrum_omega(w,**P),w*abs(current_moment_transform_omega(w,**P)),rtol=1e-13)
def test_high_frequency_stable():
 f=np.logspace(6,np.log10(3e10),500); y=derivative_spectrum_hz(f,**P); assert np.all(np.isfinite(y)) and np.all(y>0)
 # At 100 GHz the true log-amplitude is below binary64 range, so an explicit zero is legitimate.
 assert derivative_spectrum_hz(1e11,**P)==0.0
def test_T0_phase_only():
 w=2*np.pi*3e9; p2=P|{'T0':8e-9}; assert np.isclose(abs(current_moment_transform_omega(w,**P)),abs(current_moment_transform_omega(w,**p2)),rtol=1e-13)
def test_I0_linear_scaling():
 w=2*np.pi*1e9; p2=P|{'I0':.88}; assert np.isclose(abs(current_moment_transform_omega(w,**p2))/abs(current_moment_transform_omega(w,**P)),2,rtol=1e-13)
def test_fft_trusted_band():
 t=build_fft_time_grid(**P,max_frequency_hz=30e9); fg,z=discrete_fft_derivative(t,**P,method='analytic'); f=np.array([30e6,100e6,300e6,1e9,3e9,10e9,30e9]); ana=derivative_spectrum_hz(f,**P);rel=np.abs(np.interp(f,fg,np.abs(z))-ana)/ana;trusted=(f>=5*fg[1])&(f<=fg[-1]/10)&(ana>=np.max(ana)*1e-6);assert np.max(rel[trusted])<.02;assert not trusted[-1]
