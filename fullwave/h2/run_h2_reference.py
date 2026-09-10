"""Canonical two-dipole H2 receiver fixture; never uses the G3 waveform."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.ports import CurvePort

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.foundation import mesh_report  # noqa: E402
from streamer_rf.fullwave.receiver import near_far_classification  # noqa: E402


def _integrate(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def run(output, case, spacing_mm, domain_scale):
    started = time.perf_counter()
    if output.exists():
        raise ValueError("EXISTING_RAW_CASE_REQUIRES_EXPLICIT_REUSE")
    output.mkdir(parents=True)

    f_start, f_stop, f_center = 0.7e9, 1.3e9, 1.0e9
    arm_mm, separation_mm = 75.0, 300.0
    fdtd = openEMS(NrTS=30000, EndCriteria=1e-5)
    fdtd.SetGaussExcite(f_center, 0.5 * (f_stop - f_start))
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)

    tx_x, rx_x = -separation_mm / 2, separation_mm / 2
    inner_half = np.array([300.0, 150.0, 225.0]) * domain_scale
    grid.AddLine("x", [-inner_half[0], tx_x, 0, rx_x, inner_half[0]])
    grid.AddLine("y", [-inner_half[1], 0, inner_half[1]])
    grid.AddLine("z", [-inner_half[2], -arm_mm, 0, arm_mm, inner_half[2]])
    grid.SmoothMeshLines("all", spacing_mm, 1.4)
    for axis in "xyz":
        lines = np.asarray(grid.GetLines(axis), dtype=float)
        extra_left = lines[0] - np.arange(8, 0, -1) * (lines[1] - lines[0])
        extra_right = lines[-1] + np.arange(1, 9) * (lines[-1] - lines[-2])
        grid.AddLine(axis, np.r_[extra_left, extra_right])
    axes_mm = [np.asarray(grid.GetLines(axis), dtype=float) for axis in "xyz"]
    mesh = mesh_report([axis * 1e-3 for axis in axes_mm], f_stop)

    tx = CurvePort(
        csx, 1, R=50.0,
        start=[tx_x, 0, -arm_mm], stop=[tx_x, 0, arm_mm], excite=1,
    )
    rx = CurvePort(
        csx, 2, R=50.0,
        start=[rx_x, 0, -arm_mm], stop=[rx_x, 0, arm_mm], excite=0,
    )
    nf2ff = fdtd.CreateNF2FFBox(name="h2_nf2ff", frequency=[f_center])
    fdtd.Run(str(output), cleanup=False, numThreads=2)

    frequency = np.linspace(f_start, f_stop, 301)
    tx.CalcPort(str(output), frequency, ref_impedance=50.0)
    rx.CalcPort(str(output), frequency, ref_impedance=50.0)
    s11 = tx.uf_ref / tx.uf_inc
    s21 = rx.uf_ref / tx.uf_inc
    zin = tx.uf_tot / tx.if_tot
    phase = np.unwrap(np.angle(s21))
    group_delay = -np.gradient(phase, 2 * np.pi * frequency)
    loaded_voltage_transfer = rx.uf_tot / tx.uf_inc
    # Identical geometry and exchange symmetry give S22=S11 for this fixture.
    receiver_zin = 50.0 * (1 + s11) / (1 - s11)
    open_circuit_voltage_transfer = loaded_voltage_transfer * (receiver_zin + 50.0) / 50.0

    columns = np.column_stack([
        frequency, s11.real, s11.imag, s21.real, s21.imag,
        20 * np.log10(np.maximum(abs(s11), np.finfo(float).tiny)),
        20 * np.log10(np.maximum(abs(s21), np.finfo(float).tiny)),
        phase, group_delay, zin.real, zin.imag,
        loaded_voltage_transfer.real, loaded_voltage_transfer.imag,
        open_circuit_voltage_transfer.real, open_circuit_voltage_transfer.imag,
    ])
    np.savetxt(
        output / "sparameters.csv", columns, delimiter=",",
        header=(
            "frequency_Hz,S11_real,S11_imag,S21_real,S21_imag,S11_dB,S21_dB,"
            "S21_phase_unwrapped_rad,group_delay_s,Zin_real_ohm,Zin_imag_ohm,"
            "Vrx_50ohm_per_tx_inc_real,Vrx_50ohm_per_tx_inc_imag,"
            "Vrx_open_circuit_equivalent_per_tx_inc_real,"
            "Vrx_open_circuit_equivalent_per_tx_inc_imag"
        ), comments="",
    )

    theta = np.arange(0.0, 181.0, 15.0)
    phi = [0.0, 90.0]
    far = nf2ff.CalcNF2FF(str(output), f_center, theta, phi, radius=1.0)
    e_norm = np.asarray(far.E_norm[0], dtype=float)
    np.savetxt(
        output / "nf2ff_cut.csv",
        np.column_stack([theta, e_norm[:, 0], e_norm[:, 1]]),
        delimiter=",", header="theta_deg,E_norm_phi0,E_norm_phi90", comments="",
    )

    idx_center = int(np.argmin(abs(frequency - f_center)))
    idx_resonance = int(np.argmin(abs(s11)))
    minus10 = 20 * np.log10(np.maximum(abs(s11), np.finfo(float).tiny)) <= -10
    if np.any(minus10):
        bandwidth = [float(frequency[minus10][0]), float(frequency[minus10][-1])]
    else:
        bandwidth = None
    result = {
        "case": case,
        "raw_directory": str(output),
        "spacing_mm": spacing_mm,
        "domain_scale": domain_scale,
        "PML": "PML_8",
        **mesh,
        "frequency_start_Hz": f_start,
        "frequency_stop_Hz": f_stop,
        "frequency_samples": len(frequency),
        "tx_rx_distance_m": separation_mm * 1e-3,
        "dipole_total_length_m": 2 * arm_mm * 1e-3,
        "Z0_ohm": 50.0,
        "S11_min_dB": float(np.min(columns[:, 5])),
        "S11_min_frequency_Hz": float(frequency[idx_resonance]),
        "S11_minus10dB_band_Hz": bandwidth,
        "Zin_at_1GHz_ohm": [float(zin[idx_center].real), float(zin[idx_center].imag)],
        "S21_at_1GHz": [float(s21[idx_center].real), float(s21[idx_center].imag)],
        "S21_at_1GHz_dB": float(columns[idx_center, 6]),
        "S21_at_1GHz_phase_rad": float(phase[idx_center]),
        "group_delay_at_1GHz_s": float(group_delay[idx_center]),
        "response_finite": bool(np.all(np.isfinite(columns))),
        "near_far_at_1GHz": near_far_classification(0.15, 0.3, 1e9),
        "nf2ff": {
            "frequency_Hz": f_center,
            "Prad": float(np.asarray(far.Prad).ravel()[0]),
            "Dmax": float(np.asarray(far.Dmax).ravel()[0]),
            "E_norm_max": float(np.max(e_norm)),
            "E_norm_axis": float(np.max(e_norm[[0, -1], :])),
            "E_norm_broadside": float(np.max(e_norm[theta == 90, :])),
            "finite": bool(np.all(np.isfinite(e_norm))),
        },
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "raw_bytes": sum(path.stat().st_size for path in output.rglob("*") if path.is_file()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--spacing-mm", type=float, required=True)
    parser.add_argument("--domain-scale", type=float, default=1.0)
    args = parser.parse_args()
    run(args.output, args.case, args.spacing_mm, args.domain_scale)
