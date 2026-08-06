"""Energy spectral density in angular-frequency and hertz conventions."""
import numpy as np
from scipy.constants import epsilon_0,c
from scipy.integrate import simpson
from .fourier import derivative_spectrum_omega,derivative_spectrum_hz

def esd_per_angular_frequency(omega,I0,alpha,beta,T0):
    s=derivative_spectrum_omega(omega,I0,alpha,beta,T0); return s*s/(6*np.pi*epsilon_0*c**3)
def esd_per_hz(f,I0,alpha,beta,T0):
    s=derivative_spectrum_hz(f,I0,alpha,beta,T0); return 2*s*s/(3*epsilon_0*c**3)
def integrate_band_energy(f,esd,f_low,f_high):
    f=np.asarray(f,float); esd=np.asarray(esd,float)
    if f.ndim!=1 or f.shape!=esd.shape or np.any(np.diff(f)<=0): raise ValueError('f and esd must be matching increasing 1D arrays')
    tol=32*np.finfo(float).eps*max(abs(f[0]),abs(f[-1]),abs(f_low),abs(f_high)); f_low=max(f_low,f[0]) if f_low>=f[0]-tol else f_low; f_high=min(f_high,f[-1]) if f_high<=f[-1]+tol else f_high
    if f_low< f[0] or f_high>f[-1] or f_low>=f_high: raise ValueError('band outside grid or invalid')
    inner=(f>f_low)&(f<f_high); x=np.r_[f_low,f[inner],f_high]; y=np.interp(x,f,esd)
    return float(simpson(y,x=x))

