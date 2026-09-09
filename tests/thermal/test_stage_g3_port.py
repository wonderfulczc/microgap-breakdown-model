from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"python"))
from streamer_rf.circuit.model import CircuitParameters, SeriesRLCGapCircuit
from streamer_rf.circuit.signals import ConductanceProfile, PiecewiseLinearSignal
from streamer_rf.thermal.port import dynamic_impedance, port_frame, energy_metrics, PortTransientContract


def waveform(g=.01,c=1e-6,n=1001):
    t=np.linspace(0,1e-3,n)
    v=2+100*t
    ce=2e-6
    series=g*v+(c+ce)*100
    p=CircuitParameters(1e-6,ce,c,1.,ConductanceProfile.constant(g),PiecewiseLinearSignal.constant(0.))
    return port_frame(SeriesRLCGapCircuit(p).evaluate_solution(t,series,v))


@pytest.mark.parametrize("g,c",[(.01,0.),(0.,1e-6),(.01,1e-6)])
def test_R_C_RC_and_energy(g,c):
    f=waveform(g,c)
    np.testing.assert_allclose(f.I_conduction_A,g*f.V_port_V,atol=1e-16)
    np.testing.assert_allclose(f.I_gap_displacement_A,c*100,atol=1e-16)
    np.testing.assert_allclose(f.I_port_A,g*f.V_port_V+c*100,atol=1e-16)
    assert energy_metrics(f)["energy_residual_relative"]<1e-12
    if c==0:
        np.testing.assert_allclose(f.Z_dynamic_ohm,1/g)
    if g==0:
        assert np.all(f.I_conduction_A==0)


def test_signed_power_zero_crossing_and_reference_reversal():
    v=np.array([2.,2.,0.,-2.])
    i=np.array([1.,-1.,1.,0.])
    z,valid,_=dynamic_impedance(v,i)
    np.testing.assert_array_equal(v*i,[2.,-2.,0.,0.])
    np.testing.assert_allclose(z[:3],[2.,-2.,0.])
    assert not valid[-1] and np.isnan(z[-1])
    reversed_z,_,_=dynamic_impedance(v,-i)
    np.testing.assert_allclose(reversed_z[:3],-z[:3])
    # Reverse both terminal voltage and entering-current convention: power invariant.
    np.testing.assert_allclose((-v)*(-i),v*i)


def test_relative_gate_and_no_inf():
    z,valid,gate=dynamic_impedance(np.ones(4),[1.,1e-9,0.,-1.])
    assert gate==1e-8
    np.testing.assert_array_equal(valid,[True,False,False,True])
    assert not np.any(np.isinf(z))


@pytest.mark.parametrize("v,i",[([np.nan],[1.]),([1.],[np.inf]),([] ,[]),([1.],[1.,2.])])
def test_nonfinite_or_invalid_waveform(v,i):
    with pytest.raises(ValueError):
        dynamic_impedance(v,i)


def test_contract_and_frozen_status():
    f=waveform()
    c=PortTransientContract("synthetic",f,1e-6).validate()
    assert c["production_calibration_status"]=="PRODUCTION_PORT_TRANSIENT_NOT_CALIBRATED"
    assert c["G0_status"]=="HANDOFF_CALIBRATION_PENDING"
    assert "H_GAP_CAPACITANCE_PARTITION_REQUIRED" in c["interface_constraints"]
    assert "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED" in c["interface_constraints"]


@pytest.mark.parametrize("change",["nan","time","kcl","direction","impedance"])
def test_invalid_contract(change):
    f=waveform()
    if change=="nan": f.loc[0,"V_port_V"]=np.nan
    if change=="time": f.loc[1,"time_s"]=0.
    if change=="kcl": f.loc[0,"I_conduction_A"]+=1
    if change=="direction": f.loc[0,"I_series_A"]+=1
    if change=="impedance": f.loc[0,"Z_dynamic_ohm"]+=100
    with pytest.raises(ValueError):
        PortTransientContract("synthetic",f,1e-6).validate()


def test_resampling_preserves_energy_and_peaks():
    dense,sparse=waveform(n=2001),waveform(n=1001)
    a,b=energy_metrics(dense),energy_metrics(sparse)
    for key in ("port_energy_J","conduction_energy_J"):
        assert a[key]==pytest.approx(b[key],rel=1e-8)
    for key in ("V_port_V","I_port_A"):
        assert dense[key].max()==sparse[key].max()


def test_invalid_energy_time_and_power():
    f=waveform()
    f.loc[2,"time_s"]=-1
    with pytest.raises(ValueError): energy_metrics(f)
    f=waveform()
    f.loc[2,"P_port_W"]=np.inf
    with pytest.raises(ValueError): energy_metrics(f)
