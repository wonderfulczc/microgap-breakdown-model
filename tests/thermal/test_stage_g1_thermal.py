from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import h, k as kb, m_e, electron_volt
from scipy.special import j0, jn_zeros

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.thermal.properties import ATM, COEFFICIENT_FILE, EquilibriumAirProperties, gaussian, sigmoid, xi, phi
from streamer_rf.thermal.radial import (RadialGrid, RadialThermalSolver, conduction_power,
    joule_heating, prescribed_current_from_csv, required_thermal_energy,
    thermal_initialization_from_g0, validate_calibrated_radial_state)


@pytest.fixture(scope="module")
def air():
    return EquilibriumAirProperties()


def test_provenance_and_continuation_coefficients(air):
    record = json.loads((ROOT / "thermal/g1/dangola_property_provenance.json").read_text())
    assert record["doi"] == "10.1140/epjd/e2007-00305-4"
    assert record["coefficient_sha256"] == air.coefficient_sha256
    assert len(air.tables["22"]["terms"][6][0]) == 10
    assert air.tables["22"]["terms"][6][0][6:] == [6.598926e-3,-2.119755e-4,-1.369506e-4,-8.311253e-6]
    assert len(air.tables["23"]["terms"][2][0]) == 7
    assert len(air.tables["26"]["terms"]) == 10
    assert air.tables["26"]["terms"][9][2] == [-1.374699,2.577156e-2,-1.763376e-3]
    assert air.coefficient_sha256 == hashlib.sha256(COEFFICIENT_FILE.read_bytes()).hexdigest()


def test_base_functions():
    assert gaussian(5.,5.,2.) == 1.
    assert sigmoid(5.,5.,2.) == .5
    assert sigmoid(1e9,0.,1.) == 1.
    assert xi(2.,3.,4.,2.,1.) == pytest.approx(3-4/np.e)
    assert phi(3.,2.,2.) == 18.


@pytest.mark.parametrize("T,p", [(49,ATM),(60001,ATM),(300,.009*ATM),(300,101*ATM),(np.nan,ATM),(300,np.inf)])
def test_range_rejection(air,T,p):
    result = air.evaluate(T,p)
    assert not result["valid"]
    assert result["status"] == "PROPERTY_OUT_OF_RANGE"
    assert np.isnan(result["rho_kg_m3"])


def test_exact_property_bounds_and_positive_domain(air):
    result = air.evaluate(np.geomspace(50,60000,71)[:,None],np.geomspace(.01,100,7)*ATM)
    assert np.all(result["valid"])


def test_selected_atmospheric_values_and_continuity(air):
    T = np.array([300,1000,3000,5000,10000,15000,30000])
    a, b = air.evaluate(T,ATM), air.evaluate(T+1e-3,ATM)
    assert np.all(a["valid"])
    # Rounded transcription regression anchors, NOT independent experimental data.
    np.testing.assert_allclose(a["sigma_S_m"][2:], [.0235425743,24.7192793,2991.45598,7332.13134,12060.4154],rtol=2e-8)
    assert a["k_W_mK"][0] == pytest.approx(.02611087, rel=3e-7)
    for name in ("rho_kg_m3","h_J_kg","cp_J_kgK","k_W_mK","sigma_S_m"):
        assert np.all(np.isfinite(a[name]))
        np.testing.assert_allclose(a[name],b[name],rtol=.001)
    np.testing.assert_allclose(a["e_J_kg"],a["h_J_kg"]-ATM/a["rho_kg_m3"],rtol=1e-14)


@pytest.mark.parametrize("T,p", [(300,1),(1000,.1),(5000,10),(15000,1),(30000,50)])
def test_thermodynamic_roundtrip(air,T,p):
    state = air.evaluate(T,p*ATM)
    tr,pr = air.invert(state["rho_kg_m3"],state["e_J_kg"],T*.8,p*ATM*1.2)
    assert tr == pytest.approx(T,rel=3e-8)
    assert pr == pytest.approx(p*ATM,rel=3e-8)


def test_inversion_invalid_target(air):
    with pytest.raises(ValueError):
        air.invert(-1,1e5)
    with pytest.raises(ValueError):
        air.invert(1,np.nan)
    with pytest.raises(ValueError):
        air.invert(1,1e6,np.nan,ATM)


def test_inversion_far_guess_and_unreachable_state(air):
    state=air.evaluate(15000,ATM)
    T,p=air.invert(state["rho_kg_m3"],state["e_J_kg"],300,ATM)
    assert T==pytest.approx(15000,rel=2e-8)
    assert p==pytest.approx(ATM,rel=2e-8)
    with pytest.raises(ValueError):
        air.invert(1,-1e9)


def test_cp_fit_consistent_with_enthalpy_derivative_and_sound(air):
    T=np.array([300.,1000.,3000.,5000.,10000.,15000.,30000.])
    result=air.evaluate(T,ATM)
    dh=(air.evaluate(T+.001,ATM)["h_J_kg"]-air.evaluate(T-.001,ATM)["h_J_kg"])/.002
    np.testing.assert_allclose(result["cp_J_kgK"],dh,rtol=.01)
    sound,cv=air.acoustic_properties(T,ATM)
    assert np.all(sound>0) and np.all(cv>0)


def test_radial_geometry():
    g = RadialGrid(100,.001)
    assert sum(g.volumes) == pytest.approx(np.pi*1e-6)
    assert g.areas[0] == 0
    assert g.areas[-1] == pytest.approx(2*np.pi*.001)
    np.testing.assert_allclose(g.volumes, 2*np.pi*g.centers*g.dr)
    with pytest.raises(ValueError):
        RadialGrid(1,1.)


def test_uniform_stationary_and_axis_symmetry(air):
    s = RadialThermalSolver(RadialGrid(12,.001),.001,300.,ATM,properties=air)
    initial = s.U.copy()
    for _ in range(4):
        s.step(1e-10)
    np.testing.assert_allclose(s.U,initial,rtol=1e-13,atol=1e-15)
    assert s.QJ == 0.
    assert abs(s.diagnostics()["energy_residual_J"]) < 1e-18


def test_conducting_cylinder_joule_identity_and_zero_current():
    g=RadialGrid(40,.002)
    d=joule_heating(np.full(40,17.),g.volumes,.003,2.)
    assert d["Rsp_ohm"] == pytest.approx(.003/(17*np.pi*.002**2))
    assert d["PJ_volume_W"] == pytest.approx(d["PJ_circuit_W"],rel=1e-14)
    assert d["joule_closure_relative"] < 1e-14
    zero=joule_heating(np.full(40,17.),g.volumes,.003,0.)
    assert np.all(zero["qJ_W_m3"] == 0)
    with pytest.raises(ValueError):
        joule_heating(np.zeros(40),g.volumes,.003,1.)


def test_constant_power_energy_and_closed_conservation(air):
    g=RadialGrid(8,.001)
    s=RadialThermalSolver(g,.001,10000.,ATM,properties=air,boundary="closed")
    power=.5
    # Test-only feedback selects prescribed power, without circuit dynamics.
    s.current=lambda t: np.sqrt(power*float(np.dot(air.evaluate(s.T,s.p)["sigma_S_m"],g.volumes))/.001)
    for _ in range(5):
        s.step(1e-10)
    assert s.QJ == pytest.approx(power*s.time,rel=1e-14)
    assert s.total_energy()-s.initial_energy == pytest.approx(s.QJ,rel=1e-8)
    assert abs(s.diagnostics()["mass_residual_kg"]) < 1e-22


def test_conduction_bessel_decay_and_spatial_convergence():
    # Synthetic constant diffusivity, closed cylindrical Neumann eigenmode J0.
    root=jn_zeros(1,1)[0]
    errors=[]
    for n in (40,80):
        g=RadialGrid(n,1.)
        T=300+j0(root*g.centers)
        initial=T.copy()
        t=0.
        while t < .005:
            dt=min(.1*g.dr**2,.005-t)
            heat,faces=conduction_power(g,T,np.ones(n))
            assert faces[0] == 0.
            assert faces[-1] == 0.
            assert abs(np.dot(heat,g.volumes)) < 1e-10
            T+=dt*heat
            t+=dt
        exact=300+(initial-300)*np.exp(-root**2*t)
        errors.append(np.linalg.norm(T-exact))
    assert errors[1] < .45*errors[0]


def test_invalid_states_and_waveform(tmp_path,air):
    with pytest.raises(ValueError):
        RadialThermalSolver(RadialGrid(5,.001),.001,-300.,ATM,properties=air)
    with pytest.raises(ValueError):
        joule_heating(np.array([-1.]),np.array([1.]),1.,1.)
    path=tmp_path/"current.csv"
    path.write_text("time_s,current_A\n0,0\n1,2\n2,0\n")
    signal=prescribed_current_from_csv(path)
    assert signal(.5)==1.
    with pytest.raises(ValueError):
        signal(3.)
    with pytest.raises(ValueError):
        signal(np.nan)
    path.write_text("time_s,current_A\n0,1\n")
    signal=prescribed_current_from_csv(path)
    with pytest.raises(ValueError):
        signal(1.)


def test_g0_unresolved_and_future_profile_contract(air):
    summary=json.loads((ROOT/"circuit/g0_handoff/g0_handoff_summary.json").read_text())
    assert thermal_initialization_from_g0(summary)["status"] == "G0_THERMAL_INITIALIZATION_UNRESOLVED"
    g=RadialGrid(4,.001)
    with pytest.raises(ValueError):
        validate_calibrated_radial_state(g,np.full(4,300.),np.full(4,ATM),"","",air)
    data=validate_calibrated_radial_state(g,np.full(4,300.),np.full(4,ATM),"test","synthetic",air)
    assert data["source_id"] == "test"


def test_required_energy_has_explicit_target_only():
    assert required_thermal_energy([1.,2.],[100.,100.],[200.,300.]) == 500.
    with pytest.raises(ValueError):
        required_thermal_energy(1.,100.,50.)


def test_single_ionization_saha_sanity_only():
    # Ground-state hydrogen: 2 U_ion/U_neutral = 1. NIST H I energy.
    # University of Alabama AY521 notes 2018-09-17, equilibrium equation.
    T=np.array([3000.,10000.,30000.])
    n=1e22
    rhs=(2*np.pi*m_e*kb*T/h**2)**1.5*np.exp(-13.598434599702*electron_volt/(kb*T))
    y=rhs/n
    fraction=2/(1+np.sqrt(1+4/y))
    np.testing.assert_allclose(n*fraction**2/(1-fraction),rhs,rtol=1e-10)
    assert np.all(np.diff(fraction)>0)
    assert fraction[0] < 1e-8 and fraction[-1] > .99


def test_hot_channel_has_outward_motion(air):
    g=RadialGrid(20,5e-4)
    f=np.exp(-(g.centers/5e-5)**2)
    s=RadialThermalSolver(g,.001,300+9700*f,ATM*(1+f),properties=air)
    initial=s.p.copy()
    for _ in range(3):
        s.step(1e-10)
    d=s.diagnostics()
    assert d["outward_velocity_max_m_s"] > 0
    assert np.max(np.abs(s.p-initial)) > 0
    assert d["T_min_K"] > 0 and d["rho_min_kg_m3"] > 0 and d["p_min_Pa"] > 0
    assert abs(d["energy_residual_J"])/s.initial_energy < 1e-12


def test_rejected_trial_does_not_accumulate_energy(monkeypatch):
    air=EquilibriumAirProperties()
    solver=RadialThermalSolver(RadialGrid(5,.001),.001,10000.,ATM,properties=air,
                              boundary="closed",current=lambda t:.1)
    invert=air.invert
    calls=0
    def fail_once(*args):
        nonlocal calls
        calls+=1
        if calls==1:
            raise ValueError("test trial rejection")
        return invert(*args)
    monkeypatch.setattr(air,"invert",fail_once)
    power=solver.diagnostics()["PJ_volume_W"]
    dt=solver.step(1e-11)
    assert solver.rejected_steps==1 and solver.steps==1
    assert dt==5e-12
    assert solver.QJ==pytest.approx(power*dt,rel=1e-14)
