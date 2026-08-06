"""Numerically stable Shi et al. (2019) Eq. (4) lifecycle model."""
import numpy as np
from scipy.special import expit

def _validate(t,I0,alpha,beta,T0):
    for name,val in [('I0',I0),('alpha',alpha),('beta',beta),('T0',T0)]:
        if not np.isfinite(val): raise ValueError(f'{name} must be finite')
    if I0<=0 or alpha<=0 or beta<=0: raise ValueError('I0, alpha and beta must be positive')
    a=np.asarray(t,dtype=float)
    if not np.all(np.isfinite(a)): raise ValueError('t must contain only finite values')
    return a

def _scalar(a,y): return float(y) if np.ndim(a)==0 else y

def current_moment(t,I0,alpha,beta,T0):
    a=_validate(t,I0,alpha,beta,T0); x=a-T0; q=alpha+beta
    logI=np.log(I0)+alpha*x-np.logaddexp(0.0,q*x)
    y=np.exp(logI)
    if not np.all(np.isfinite(y)): raise FloatingPointError('non-finite current moment')
    return _scalar(a,y)

def current_moment_derivative(t,I0,alpha,beta,T0):
    a=_validate(t,I0,alpha,beta,T0); q=alpha+beta; z=q*(a-T0)
    sigmoid=expit(z)
    y=np.asarray(current_moment(a,I0,alpha,beta,T0))*(alpha-q*sigmoid)
    if not np.all(np.isfinite(y)): raise FloatingPointError('non-finite derivative')
    return _scalar(a,y)

def peak_time(I0,alpha,beta,T0):
    _validate(T0,I0,alpha,beta,T0); return T0+np.log(alpha/beta)/(alpha+beta)

def peak_value(I0,alpha,beta,T0):
    _validate(T0,I0,alpha,beta,T0); q=alpha+beta
    return float(np.exp(np.log(I0)+np.log(beta/q)+(alpha/q)*np.log(alpha/beta)))

