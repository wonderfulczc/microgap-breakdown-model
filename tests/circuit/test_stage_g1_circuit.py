from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.circuit import CircuitParameters, ConductanceProfile, PiecewiseLinearSignal, SeriesRLCGapCircuit, coupling_error  # noqa: E402
from streamer_rf.circuit.validation import estimate_damping_alpha, estimate_oscillation_frequency, rlc_theory  # noqa: E402


def make_circuit(*, G: float = 0.0, Cgap: float = 0.0, Rs: float = 10.0) -> SeriesRLCGapCircuit:
    params = CircuitParameters(
        L_H=1.0e-6,
        Cext_F=1.0e-12,
        Cgap_F=Cgap,
        Rs_ohm=Rs,
        conductance=ConductanceProfile.constant(G),
        source_voltage=PiecewiseLinearSignal.constant(0.0),
    )
    return SeriesRLCGapCircuit(params)


def solve(circuit: SeriesRLCGapCircuit, *, t1: float = 80e-9, dt: float = 1e-11):
    return circuit.solve(t_span=(0.0, t1), initial_state=(0.0, 100.0), output_dt_s=dt, max_step_s=dt, rtol=2e-10, atol=1e-13)


def test_constant_rlc_frequency_and_damping() -> None:
    circuit = make_circuit(G=0.0, Cgap=0.0, Rs=10.0)
    sol = solve(circuit, t1=100e-9)
    theory = rlc_theory(1.0e-6, 1.0e-12, 10.0)
    omega = estimate_oscillation_frequency(sol.time_s, sol.Vgap_V)
    alpha = estimate_damping_alpha(sol.time_s, sol.Vgap_V)
    assert abs(omega - theory["omega_d_rad_s"]) / theory["omega_d_rad_s"] < 5e-3
    assert abs(alpha - theory["alpha_1_s"]) / theory["alpha_1_s"] < 1e-2


def test_open_gap_displacement_current_only_and_residuals() -> None:
    circuit = make_circuit(G=0.0, Cgap=2e-13, Rs=0.0)
    sol = solve(circuit, t1=20e-9)
    assert np.max(np.abs(sol.currents.I_gap_cond_A)) == 0.0
    np.testing.assert_allclose(sol.currents.I_gap_disp_A, circuit.params.Cgap_F * sol.dVgap_dt_Vps)
    assert np.max(np.abs(sol.KCL_residual_A)) < 1e-12
    assert np.max(np.abs(sol.KVL_residual_V)) < 1e-10
    assert np.max(np.abs(sol.gap_conductive_loss_J)) == 0.0


def test_finite_G_conductive_current_and_power_integration() -> None:
    G = 2e-4
    circuit = make_circuit(G=G, Cgap=1e-13, Rs=2.0)
    sol = solve(circuit, t1=30e-9)
    np.testing.assert_allclose(sol.currents.I_gap_cond_A, G * sol.Vgap_V)
    direct_loss = np.trapezoid(G * sol.Vgap_V**2, sol.time_s)
    assert abs(direct_loss - sol.gap_conductive_loss_J[-1]) / direct_loss < 1e-12
    assert sol.gap_conductive_loss_J[-1] > 0.0


def test_kcl_kvl_and_energy_balance_are_controlled() -> None:
    circuit = make_circuit(G=1e-4, Cgap=1e-13, Rs=3.0)
    sol = solve(circuit, t1=30e-9)
    assert np.max(np.abs(sol.KCL_residual_A)) < 1e-12
    assert np.max(np.abs(sol.KVL_residual_V)) < 1e-10
    e0 = sol.stored_energy_J[0]
    assert np.max(np.abs(sol.energy_balance_residual_J)) / e0 < 2e-3


def test_kcl_residual_explicit_identity() -> None:
    circuit = make_circuit(G=5e-5, Cgap=2e-13, Rs=1.0)
    sol = solve(circuit, t1=10e-9)
    recomputed = sol.I_L_A - sol.currents.I_Cext_A - sol.currents.I_gap_cond_A - sol.currents.I_gap_disp_A
    np.testing.assert_allclose(sol.KCL_residual_A, recomputed, atol=1e-16)


def test_kvl_residual_explicit_identity() -> None:
    circuit = make_circuit(G=5e-5, Cgap=2e-13, Rs=1.0)
    sol = solve(circuit, t1=10e-9)
    recomputed = sol.source_voltage_V - circuit.params.Rs_ohm * sol.I_L_A - circuit.params.L_H * sol.dI_L_dt_Aps - sol.Vgap_V
    np.testing.assert_allclose(sol.KVL_residual_V, recomputed, atol=1e-12)


def test_energy_balance_components_are_nonnegative_losses() -> None:
    circuit = make_circuit(G=1e-4, Cgap=1e-13, Rs=3.0)
    sol = solve(circuit, t1=20e-9)
    assert np.all(np.diff(sol.gap_conductive_loss_J) >= -1e-24)
    assert np.all(np.diff(sol.resistor_loss_J) >= -1e-24)
    assert sol.gap_conductive_loss_J[-1] > 0.0
    assert sol.resistor_loss_J[-1] > 0.0


def test_time_varying_G_interpolation_and_stability() -> None:
    conductance = ConductanceProfile.from_samples(
        np.array([0.0, 5e-9, 15e-9, 30e-9]),
        np.array([1e-9, 1e-9, 1e-3, 1e-3]),
        extrapolation="invalid",
    )
    circuit = SeriesRLCGapCircuit(
        CircuitParameters(
            L_H=8e-7,
            Cext_F=1e-12,
            Cgap_F=1e-13,
            Rs_ohm=5.0,
            conductance=conductance,
            source_voltage=PiecewiseLinearSignal(np.array([0.0, 30e-9]), np.array([0.0, 0.0]), extrapolation="invalid"),
        )
    )
    assert np.isclose(conductance(10e-9), 5.000005e-4)
    sol = circuit.solve(t_span=(0.0, 30e-9), initial_state=(0.0, 100.0), output_dt_s=1e-11)
    assert np.all(np.isfinite(sol.to_frame().select_dtypes(include=[float, int]).to_numpy()))
    assert np.max(np.abs(sol.currents.I_gap_total_A)) > 0.0


def test_voltage_waveform_interpolation() -> None:
    sig = PiecewiseLinearSignal(np.array([0.0, 1.0, 2.0]), np.array([0.0, 10.0, 0.0]), extrapolation="invalid")
    assert sig(0.5) == 5.0
    assert sig(np.array([0.5, 1.5])).tolist() == [5.0, 5.0]


def test_hold_extrapolation_is_explicit() -> None:
    sig = PiecewiseLinearSignal(np.array([0.0, 1.0]), np.array([2.0, 4.0]), extrapolation="hold")
    assert sig(-1.0) == 2.0
    assert sig(2.0) == 4.0


def test_negative_conductance_rejected() -> None:
    with pytest.raises(ValueError):
        ConductanceProfile.constant(-1.0)
    with pytest.raises(ValueError):
        ConductanceProfile.from_samples(np.array([0.0, 1.0]), np.array([0.0, -1.0]))


def test_csv_Gb_and_Rb_inputs(tmp_path: Path) -> None:
    g_csv = tmp_path / "g.csv"
    pd.DataFrame({"time_s": [0.0, 1.0], "Gb_S": [1.0, 2.0]}).to_csv(g_csv, index=False)
    assert ConductanceProfile.from_csv(g_csv)(0.5) == 1.5
    r_csv = tmp_path / "r.csv"
    pd.DataFrame({"time_s": [0.0, 1.0], "Rb_ohm": [2.0, np.inf]}).to_csv(r_csv, index=False)
    assert ConductanceProfile.from_csv(r_csv)(0.0) == 0.5
    assert ConductanceProfile.from_csv(r_csv)(1.0) == 0.0


def test_invalid_time_and_extrapolation_handling() -> None:
    with pytest.raises(ValueError):
        PiecewiseLinearSignal(np.array([0.0, 0.0]), np.array([1.0, 2.0]))
    sig = PiecewiseLinearSignal(np.array([0.0, 1.0]), np.array([0.0, 1.0]), extrapolation="invalid")
    with pytest.raises(ValueError):
        sig(1.1)
    circuit = SeriesRLCGapCircuit(
        CircuitParameters(
            L_H=1e-6,
            Cext_F=1e-12,
            Cgap_F=0.0,
            Rs_ohm=1.0,
            conductance=ConductanceProfile.from_samples(np.array([0.0, 1e-9]), np.array([0.0, 0.0]), extrapolation="invalid"),
            source_voltage=sig,
        )
    )
    with pytest.raises(ValueError):
        circuit.solve(t_span=(0.0, 2e-9), initial_state=(0.0, 1.0), output_dt_s=1e-11)


def test_invalid_output_dt_rejected() -> None:
    circuit = make_circuit()
    with pytest.raises(ValueError):
        circuit.solve(t_span=(0.0, 1e-9), initial_state=(0.0, 1.0), output_dt_s=0.0)


def test_coupling_error_metric() -> None:
    t = np.linspace(0.0, 1.0, 10)
    ref = np.ones_like(t) * 100.0
    got = ref + 2.0
    assert coupling_error(t, got, ref) == 0.02
    assert coupling_error(t, np.ones_like(t), np.zeros_like(t)) == 1.0
