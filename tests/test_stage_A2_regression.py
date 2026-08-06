import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
import numpy as np
from scipy.optimize import minimize_scalar
from streamer_rf.lifecycle import *
from streamer_rf.fourier import derivative_spectrum_hz
from streamer_rf.esd import esd_per_hz,integrate_band_energy
P=dict(I0=.44,alpha=1.7e9,beta=7.5e7,T0=6e-9)
def test_primary_peak_regression():
 assert np.isclose(peak_time(**P),7.758250938877745e-9,rtol=1e-12);assert np.isclose(peak_value(**P),0.3693459484952068,rtol=1e-12)
def test_derivative_peak_regression():
 tp=peak_time(**P);span=25/P['beta'];pos=minimize_scalar(lambda t:-current_moment_derivative(t,**P),bounds=(tp-span,tp),method='bounded',options={'xatol':1e-18});neg=minimize_scalar(lambda t:current_moment_derivative(t,**P),bounds=(tp,tp+span),method='bounded',options={'xatol':1e-18});assert np.isclose(-pos.fun,1.800363132966759e8,rtol=2e-8);assert np.isclose(neg.fun,-2.4152203608836766e7,rtol=2e-8)
def test_spectrum_regression():
 f=np.array([30e6,300e6,1e9,3e9]);expected=np.array([4.02495473e-1,1.04561251e-1,1.44865062e-4,9.52322643e-14]);assert np.allclose(derivative_spectrum_hz(f,**P),expected,rtol=5e-9)
def test_band_energy_regression():
 f=np.logspace(np.log10(3e7),np.log10(3e10),50000);e=esd_per_hz(f,**P);vals=[integrate_band_energy(f,e,3e7,3e8),integrate_band_energy(f,e,3e8,3e9),integrate_band_energy(f,e,3e9,3e10)];expected=[6.072019861646134e-8,1.8445438555603278e-9,1.1741415386837378e-33];assert np.allclose(vals,expected,rtol=2e-10)
