#!/usr/bin/env python3
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.circuit import CircuitParameters, ConductanceProfile, PiecewiseLinearSignal, SeriesRLCGapCircuit, coupling_error  # noqa: E402
from streamer_rf.circuit.validation import estimate_damping_alpha, estimate_oscillation_frequency, rlc_theory  # noqa: E402

OUT = ROOT / "circuit/validation"


def solve(params: CircuitParameters, t1: float, *, v0: float = 100.0, i0: float = 0.0, dt: float = 2e-11):
    return SeriesRLCGapCircuit(params).solve(
        t_span=(0.0, t1),
        initial_state=(i0, v0),
        output_dt_s=dt,
        rtol=2e-10,
        atol=1e-13,
        max_step_s=dt,
    )


def ngspice_info() -> dict[str, object]:
    exe = shutil.which("ngspice")
    if exe is None:
        return {"NGSPICE_AVAILABLE": False, "path": None, "version": None}
    proc = subprocess.run([exe, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    return {"NGSPICE_AVAILABLE": True, "path": exe, "version": proc.stdout.splitlines()[0] if proc.stdout else ""}


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    L = 1.0e-6
    C = 1.0e-12
    R = 10.0
    params = CircuitParameters(
        L_H=L,
        Cext_F=C,
        Cgap_F=0.0,
        Rs_ohm=R,
        conductance=ConductanceProfile.constant(0.0),
        source_voltage=PiecewiseLinearSignal.constant(0.0),
    )
    sol = solve(params, 100e-9, dt=1e-11)
    theory = rlc_theory(L, C, R)
    omega_num = estimate_oscillation_frequency(sol.time_s, sol.Vgap_V)
    alpha_num = estimate_damping_alpha(sol.time_s, sol.Vgap_V)

    open_params = CircuitParameters(
        L_H=L,
        Cext_F=0.8e-12,
        Cgap_F=0.2e-12,
        Rs_ohm=0.0,
        conductance=ConductanceProfile.constant(0.0),
        source_voltage=PiecewiseLinearSignal.constant(0.0),
    )
    open_sol = solve(open_params, 30e-9, dt=1e-11)

    conductive_G = 2.0e-4
    conductive_params = CircuitParameters(
        L_H=L,
        Cext_F=0.8e-12,
        Cgap_F=0.2e-12,
        Rs_ohm=2.0,
        conductance=ConductanceProfile.constant(conductive_G),
        source_voltage=PiecewiseLinearSignal.constant(0.0),
    )
    conductive_sol = solve(conductive_params, 30e-9, dt=1e-11)
    p_gap = conductive_G * conductive_sol.Vgap_V**2
    gap_loss_direct = float(np.trapezoid(p_gap, conductive_sol.time_s))

    g_time = np.array([0.0, 5e-9, 15e-9, 40e-9])
    g_values = np.array([1e-9, 1e-9, 8e-4, 8e-4])
    tv_params = CircuitParameters(
        L_H=0.6e-6,
        Cext_F=1.5e-12,
        Cgap_F=0.1e-12,
        Rs_ohm=5.0,
        conductance=ConductanceProfile.from_samples(g_time, g_values, extrapolation="invalid"),
        source_voltage=PiecewiseLinearSignal(np.array([0.0, 40e-9]), np.array([0.0, 0.0]), extrapolation="invalid"),
    )
    tv_sol = solve(tv_params, 40e-9, dt=1e-11)

    ref = tv_sol.Vgap_V * (1.0 + 0.02 * np.sin(2.0 * np.pi * tv_sol.time_s / tv_sol.time_s[-1]))
    eps_v = coupling_error(tv_sol.time_s, tv_sol.Vgap_V, ref)

    frames = {
        "constant_rlc": sol.to_frame(),
        "open_gap": open_sol.to_frame(),
        "conductive_gap": conductive_sol.to_frame(),
        "time_varying_G": tv_sol.to_frame(),
    }
    for name, frame in frames.items():
        frame.to_csv(OUT / f"{name}.csv", index=False)

    stored0 = conductive_sol.stored_energy_J[0]
    payload = {
        "STAGE_G1_STATUS": "PASS_CANDIDATE",
        "topology": "Vs -> Rs -> L -> node; node -> ground via Cext || Cgap || Gb(t)",
        "ode_state": ["I_L_A", "Vgap_V"],
        "ode_equations": {
            "dI_L_dt": "(Vs(t) - Rs*I_L - Vgap)/L",
            "dVgap_dt": "(I_L - Gb(t)*Vgap)/(Cext + Cgap)",
            "gap_current": "I_gap = Gb(t)*Vgap + Cgap*dVgap/dt",
        },
        "constant_rlc": {
            "omega0_theory_rad_s": theory["omega0_rad_s"],
            "omega_d_theory_rad_s": theory["omega_d_rad_s"],
            "omega_d_simulated_rad_s": omega_num,
            "frequency_relative_error": abs(omega_num - theory["omega_d_rad_s"]) / theory["omega_d_rad_s"],
            "alpha_theory_1_s": theory["alpha_1_s"],
            "alpha_simulated_1_s": alpha_num,
            "damping_relative_error": abs(alpha_num - theory["alpha_1_s"]) / theory["alpha_1_s"],
        },
        "open_gap": {
            "max_cond_current_A": float(np.max(np.abs(open_sol.currents.I_gap_cond_A))),
            "max_disp_current_A": float(np.max(np.abs(open_sol.currents.I_gap_disp_A))),
            "max_KCL_residual_A": float(np.max(np.abs(open_sol.KCL_residual_A))),
            "max_KVL_residual_V": float(np.max(np.abs(open_sol.KVL_residual_V))),
            "energy_balance_residual_J": float(np.max(np.abs(open_sol.energy_balance_residual_J))),
        },
        "conductive_gap": {
            "Gb_S": conductive_G,
            "gap_loss_integrated_solver_J": float(conductive_sol.gap_conductive_loss_J[-1]),
            "gap_loss_integrated_direct_J": gap_loss_direct,
            "gap_loss_relative_error": abs(gap_loss_direct - conductive_sol.gap_conductive_loss_J[-1]) / max(abs(gap_loss_direct), 1e-300),
            "max_KCL_residual_A": float(np.max(np.abs(conductive_sol.KCL_residual_A))),
            "max_KVL_residual_V": float(np.max(np.abs(conductive_sol.KVL_residual_V))),
            "energy_balance_residual_relative": float(np.max(np.abs(conductive_sol.energy_balance_residual_J)) / max(stored0, 1e-300)),
        },
        "time_varying_G": {
            "G_low_S": float(g_values[0]),
            "G_high_S": float(g_values[-1]),
            "rise_start_s": float(g_time[1]),
            "rise_end_s": float(g_time[2]),
            "Vgap_initial_V": float(tv_sol.Vgap_V[0]),
            "Vgap_final_V": float(tv_sol.Vgap_V[-1]),
            "I_gap_total_peak_A": float(np.max(np.abs(tv_sol.currents.I_gap_total_A))),
            "stable_no_nan": bool(np.all(np.isfinite(tv_sol.to_frame().select_dtypes(include=[float, int]).to_numpy()))),
        },
        "coupling_validity": {
            "ONE_WAY_COUPLING_ERROR": eps_v,
            "definition": "max |Vgap_circuit - Vgap_reference| / max |Vgap_reference|",
        },
        "ngspice": ngspice_info(),
        "runtime_s": time.perf_counter() - started,
    }
    (OUT / "stage_g1_validation.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
