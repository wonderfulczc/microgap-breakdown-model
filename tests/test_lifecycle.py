import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
import numpy as np
from scipy.optimize import minimize_scalar
from streamer_rf.lifecycle import *
P=dict(I0=.44,alpha=1.7e9,beta=7.5e7,T0=6e-9)
def test_tails_tend_to_zero():
 assert current_moment(-1e-6,**P)<1e-100; assert current_moment(1e-5,**P)<1e-100
def test_derivative_matches_complex_step_like_central():
 t=np.linspace(2e-9,50e-9,20); h=1e-15; num=(current_moment(t+h,**P)-current_moment(t-h,**P))/(2*h); ana=current_moment_derivative(t,**P); assert np.allclose(num,ana,rtol=2e-7,atol=1e-3)
def test_peak_derivative_zero(): assert abs(current_moment_derivative(peak_time(**P),**P))<1e-5*P['alpha']*P['I0']
def test_peak_value_matches_optimization():
 x=minimize_scalar(lambda t:-current_moment(t,**P),bounds=(-1e-8,2e-7),method='bounded',options={'xatol':1e-18}); assert np.isclose(-x.fun,peak_value(**P),rtol=1e-10)
def test_scalar_array_consistency():
 t=1e-8; assert current_moment(t,**P)==current_moment(np.array(t),**P); assert current_moment_derivative(t,**P)==current_moment_derivative(np.array(t),**P)
def test_extreme_input_stable():
 y=current_moment(np.array([-1e6,0,1e6]),**P); assert np.all(np.isfinite(y))
def test_validation_rejects_bad_input():
 import pytest
 with pytest.raises(ValueError): current_moment(0,-1,P['alpha'],P['beta'],P['T0'])
