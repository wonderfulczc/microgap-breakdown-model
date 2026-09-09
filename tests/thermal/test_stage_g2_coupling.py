from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"python"))

from streamer_rf.circuit.model import CircuitParameters, SeriesRLCGapCircuit
from streamer_rf.circuit.signals import ConductanceProfile, PiecewiseLinearSignal
from streamer_rf.thermal.coupling import (ThermalCircuitCoupling, circuit_interval,
    require_g0_initialization, resistance_conductance)
from streamer_rf.thermal.radial import RadialGrid, RadialThermalSolver, joule_heating
from streamer_rf.thermal.properties import ATM


def params(g=.001, source=0.):
    return CircuitParameters(1e-6,.8e-12,.2e-12,2.,ConductanceProfile.constant(g),
                            PiecewiseLinearSignal.constant(source))


def thermal():
    return RadialThermalSolver(RadialGrid(8,5e-4),1e-3,10000.,ATM,boundary="closed",conduction=False)


@pytest.mark.parametrize("varying",[False,True])
def test_canonical_fixed_and_prescribed_reproduction(varying):
    p = params()
    if varying:
        p = replace(p, conductance=ConductanceProfile.from_samples(
            np.array([0.,1e-9,2e-9]),np.array([.001,.002,.003])))
    adapted = circuit_interval(p,[0.,100.],0.,2e-9)
    old = SeriesRLCGapCircuit(p).solve(t_span=(0.,2e-9),initial_state=(0.,100.),
        output_dt_s=1e-12,max_step_s=1e-12,rtol=2e-10,atol=1e-13)
    np.testing.assert_allclose(adapted["state"],[old.I_L_A[-1],old.Vgap_V[-1]],rtol=2e-8)
    assert adapted["Qsp"] == pytest.approx(old.gap_conductive_loss_J[-1],rel=2e-6)
    assert adapted["kcl_max"] < 1e-14
    assert adapted["kvl_max"] < 1e-12


def test_constant_current_energy_and_nonzero_source_accounting():
    p = params(.001,102.)
    # Exact steady circuit state: I_L=G*V=.1 A, Vs=V+Rs*I.
    p = replace(p,source_voltage=PiecewiseLinearSignal.constant(100.2))
    c = circuit_interval(p,[.1,100.],0.,1e-9)
    assert c["Qsp"] == pytest.approx(.1**2*1000*1e-9,rel=1e-13)
    assert c["Qfixed"] == pytest.approx(.1**2*2*1e-9,rel=1e-13)
    assert c["Qsp"]+c["Qfixed"] == pytest.approx(c["Wsource"],rel=1e-13)


def test_zero_current_thermal_energy():
    t = thermal()
    s = ThermalCircuitCoupling(t,params(),[0.,0.],1e-11)
    initial=t.U.copy()
    s.run(3e-11)
    assert s.Qsp == s.thermal.QJ == 0.
    np.testing.assert_allclose(s.thermal.U,initial,rtol=1e-13,atol=1e-15)


def test_independent_power_and_accepted_energy():
    s = ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11)
    s.run(3e-11)
    assert s.Qsp == pytest.approx(s.thermal.QJ,rel=2e-14)
    assert sum(r["P_sp_circuit_W"]*r["interval_dt_s"] for r in s.rows) == pytest.approx(s.Qsp,rel=1e-14)
    for row in s.rows[1:]:
        assert row["P_sp_circuit_W"] == pytest.approx(row["P_sp_thermal_W"],rel=2e-14)
    assert s.rows[-1]["circuit_energy_residual_J"] == pytest.approx(0.,abs=1e-19)
    # Initial caller-owned G1 state is not mutated by the adapter.
    assert s.thermal.time == pytest.approx(3e-11)


def test_no_feedback_and_monotonic_conductivity():
    s = ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11,feedback=False)
    s.run(2e-11)
    assert all(r["Rsp_circuit_ohm"] == s.fixed_R for r in s.rows)
    area=RadialGrid(8,5e-4).volumes
    a=joule_heating(np.ones(8)*100,area,1e-3,1.)
    b=joule_heating(np.ones(8)*200,area,1e-3,1.)
    assert b["Rsp_ohm"] == pytest.approx(a["Rsp_ohm"]/2)


def test_resistance_changes_actual_parallel_branch_current():
    a=circuit_interval(params(.001),[0.,100.],0.,1e-11)
    b=circuit_interval(params(.002),[0.,100.],0.,1e-11)
    assert b["peak_current"] > a["peak_current"]
    assert b["state"][1] < a["state"][1]


@pytest.mark.parametrize("r",[0.,-1.,np.nan,np.inf,1e-320])
def test_invalid_resistance(r):
    with pytest.raises(ValueError):
        resistance_conductance(r)


def test_unresolved_g0():
    with pytest.raises(ValueError,match="G0_THERMAL_INITIALIZATION_UNRESOLVED"):
        require_g0_initialization({"G1_handoff_candidate_contract":{
            "state_valid":False,"thermal_profile_status":"NOT_AVAILABLE"}})


def test_rejected_interval_transaction(monkeypatch):
    original=RadialThermalSolver.step
    calls=[]
    def shortened(self,dt):
        calls.append(dt)
        return original(self,dt/2 if len(calls)==1 else dt)
    t=thermal()
    original_state=t.U.copy()
    s=ThermalCircuitCoupling(t,params(),[0.,100.],1e-11)
    monkeypatch.setattr(RadialThermalSolver,"step",shortened)
    dt=s.step()
    assert dt == pytest.approx(5e-12)
    assert s.discarded_intervals == 1
    assert s.thermal.steps == 1
    assert s.Qsp == pytest.approx(s.thermal.QJ,rel=2e-14)
    assert t.time == 0.
    np.testing.assert_array_equal(t.U,original_state)
    expected=circuit_interval(replace(params(),conductance=resistance_conductance(s.fixed_R)),[0.,100.],0.,dt)
    assert s.Qsp == pytest.approx(expected["Qsp"],rel=1e-14)


def test_failed_trial_leaves_state_and_budgets_unchanged(monkeypatch):
    s=ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11)
    def failed(self,dt):
        raise ValueError("SYNTHETIC_FAILURE")
    monkeypatch.setattr(RadialThermalSolver,"step",failed)
    with pytest.raises(ValueError):
        s.step()
    assert s.thermal.time == s.Qsp == s.thermal.QJ == 0.
    np.testing.assert_array_equal(s.state,[0.,100.])


def test_coupling_refinement():
    a=ThermalCircuitCoupling(thermal(),params(),[0.,100.],2e-11)
    b=ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11)
    a.run(1e-10)
    b.run(1e-10)
    assert a.Qsp == pytest.approx(b.Qsp,rel=.01)
    assert a.rows[-1]["Rsp_ohm"] == pytest.approx(b.rows[-1]["Rsp_ohm"],rel=.01)


@pytest.mark.parametrize("dt",[0.,-1.,np.nan,np.inf])
def test_invalid_interval(dt):
    with pytest.raises(ValueError):
        ThermalCircuitCoupling(thermal(),params(),[0.,100.],dt)
    s=ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11)
    with pytest.raises(ValueError):
        s.step(dt)


def test_feedback_uses_updated_thermal_resistance():
    s=ThermalCircuitCoupling(thermal(),params(),[0.,100.],1e-11)
    s.step()
    previous=s.rows[-1]["Rsp_ohm"]
    s.step()
    assert s.rows[-1]["interval_Rsp_ohm"] == previous
    assert s.rows[-1]["Rsp_ohm"] != s.fixed_R
    assert s.rows[-1]["I_A"] == pytest.approx(s.rows[-1]["V_channel_V"]/s.rows[-1]["Rsp_ohm"])
