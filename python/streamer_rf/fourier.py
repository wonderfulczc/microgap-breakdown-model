"""Continuous Fourier transform and discrete FFT validation for A2-Lifecycle."""
import numpy as np
from .lifecycle import current_moment,current_moment_derivative,peak_time,peak_value

def _p(I0,alpha,beta,T0):
    if not all(np.isfinite([I0,alpha,beta,T0])) or min(I0,alpha,beta)<=0: raise ValueError('finite positive I0, alpha, beta required')

def _log_sinh_abs(x):
    x=np.asarray(x,float); y=np.empty_like(x); small=x<20
    with np.errstate(divide='ignore',invalid='ignore'): y[small]=np.log(np.sinh(x[small]))
    y[~small]=x[~small]-np.log(2)+np.log1p(-np.exp(-2*x[~small])); return y

def current_moment_transform_magnitude_omega(omega,I0,alpha,beta,T0):
    _p(I0,alpha,beta,T0); w=np.asarray(omega,float)
    if not np.all(np.isfinite(w)): raise ValueError('omega must be finite')
    q=alpha+beta; aa=np.pi*alpha/q; b=np.pi*np.abs(w)/q
    lsin2=2*np.log(abs(np.sin(aa))); lsinh2=2*_log_sinh_abs(b)
    logden=.5*np.logaddexp(lsin2,lsinh2)
    y=np.exp(np.log(np.pi*I0/q)-logden)
    return float(y) if w.ndim==0 else y

def current_moment_transform_omega(omega,I0,alpha,beta,T0):
    _p(I0,alpha,beta,T0); w=np.asarray(omega,float); q=alpha+beta; aa=np.pi*alpha/q; bb=np.pi*w/q; ab=np.abs(bb)
    # scaled sin(aa-i bb), avoiding cosh/sinh overflow
    e2=np.exp(-2*ab); real=np.sin(aa)*.5*(1+e2); imag=-np.sign(bb)*np.cos(aa)*.5*(1-e2)
    inv_scaled=(real-1j*imag)/(real*real+imag*imag)
    z=(np.pi*I0/q)*np.exp(-ab)*inv_scaled*np.exp(-1j*w*T0)
    if not np.all(np.isfinite(z)): raise FloatingPointError('non-finite transform')
    return complex(z) if w.ndim==0 else z

def derivative_spectrum_omega(omega,I0,alpha,beta,T0): return np.abs(np.asarray(omega,float))*current_moment_transform_magnitude_omega(omega,I0,alpha,beta,T0)
def derivative_spectrum_hz(f,I0,alpha,beta,T0): return derivative_spectrum_omega(2*np.pi*np.asarray(f,float),I0,alpha,beta,T0)

def build_fft_time_grid(I0,alpha,beta,T0,max_frequency_hz=30e9,edge_ratio=1e-14,max_points=2**23):
    _p(I0,alpha,beta,T0)
    if max_frequency_hz<=0 or not 0<edge_ratio<1: raise ValueError('invalid grid controls')
    tp=peak_time(I0,alpha,beta,T0); ip=peak_value(I0,alpha,beta,T0); left=tp; right=tp; dl=1/alpha; dr=1/beta
    while current_moment(left,I0,alpha,beta,T0)/ip>edge_ratio: left-=dl; dl*=1.4
    while current_moment(right,I0,alpha,beta,T0)/ip>edge_ratio: right+=dr; dr*=1.4
    dtmax=1/(20*max_frequency_hz); nmin=int(np.ceil((right-left)/dtmax))+1; n=1<<(nmin-1).bit_length()
    if n>max_points: raise ValueError(f'FFT grid requires {n} points, exceeding max_points={max_points}')
    dt=(right-left)/(n-1); t=left+np.arange(n)*dt
    if np.max(np.abs(np.diff(t)-dt))>max(1e-25,abs(dt)*1e-9): raise RuntimeError('time grid is not uniform')
    return t

def discrete_fft_current_moment(t,I0,alpha,beta,T0,n_fft=None):
    t=np.asarray(t,float); dt=np.diff(t)
    if t.ndim!=1 or len(t)<8 or not np.allclose(dt,dt[0],rtol=1e-10,atol=0): raise ValueError('uniform 1D grid required')
    n=len(t) if n_fft is None else int(n_fft)
    if n<len(t): raise ValueError('n_fft cannot truncate data')
    y=current_moment(t,I0,alpha,beta,T0); f=np.fft.rfftfreq(n,dt[0]); z=dt[0]*np.fft.rfft(y,n=n)*np.exp(-2j*np.pi*f*t[0])
    return f,z

def discrete_fft_derivative(t,I0,alpha,beta,T0,method='analytic',n_fft=None):
    t=np.asarray(t,float); dt=t[1]-t[0]
    if method=='analytic': y=current_moment_derivative(t,I0,alpha,beta,T0)
    elif method=='five_point':
        x=current_moment(t,I0,alpha,beta,T0); y=np.gradient(x,dt,edge_order=2); y[2:-2]=(x[:-4]-8*x[1:-3]+8*x[3:-1]-x[4:])/(12*dt)
    elif method=='spectral':
        f,z=discrete_fft_current_moment(t,I0,alpha,beta,T0,n_fft); return f,1j*2*np.pi*f*z
    else: raise ValueError('method must be analytic, five_point, or spectral')
    n=len(t) if n_fft is None else int(n_fft); f=np.fft.rfftfreq(n,dt); z=dt*np.fft.rfft(y,n=n)*np.exp(-2j*np.pi*f*t[0]); return f,z

def compare_fft_to_analytic(t,frequencies,I0,alpha,beta,T0):
    fg,z=discrete_fft_current_moment(t,I0,alpha,beta,T0); frequencies=np.asarray(frequencies,float); mag=np.interp(frequencies,fg,np.abs(z)); ana=current_moment_transform_magnitude_omega(2*np.pi*frequencies,I0,alpha,beta,T0); rel=np.abs(mag-ana)/ana
    resolution=fg[1]-fg[0]; nyquist=fg[-1]; trusted=(frequencies>=5*resolution)&(frequencies<=nyquist/10)&(ana>=np.max(ana)*1e-6)
    return mag,ana,rel,trusted,resolution,nyquist

