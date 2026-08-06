import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
import numpy as np
from scipy.integrate import simpson
from streamer_rf.esd import *
P=dict(I0=.44,alpha=1.7e9,beta=7.5e7,T0=6e-9)
def test_esd_nonnegative(): assert np.all(esd_per_hz(np.logspace(6,11,100),**P)>=0)
def test_I0_quadratic():
 f=1e9;p2=P|{'I0':.88};assert np.isclose(esd_per_hz(f,**p2)/esd_per_hz(f,**P),4,rtol=1e-13)
def test_omega_hz_integrals_match():
 f=np.logspace(6,11,100000);w=2*np.pi*f; a=simpson(esd_per_hz(f,**P),x=f);b=2*simpson(esd_per_angular_frequency(w,**P),x=w);assert abs(a-b)/a<1e-10
def test_band_energy_monotonic():
 f=np.logspace(6,10,10000);e=esd_per_hz(f,**P);assert integrate_band_energy(f,e,3e7,3e9)>=integrate_band_energy(f,e,3e7,3e8)
def test_no_extra_two_pi():
 f=1e9;w=2*np.pi*f;assert np.isclose(esd_per_hz(f,**P),4*np.pi*esd_per_angular_frequency(w,**P),rtol=1e-13)
