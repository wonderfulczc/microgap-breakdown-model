#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.axisymmetric import rotate_axisymmetric_series  # noqa: E402
from streamer_rf.rf.jefimenko.constants import C0, K_B, K_E  # noqa: E402
from streamer_rf.rf.jefimenko.diagnostics import current_moment_radiation_approx, near_far_audit  # noqa: E402
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_observer, evaluate_waveform  # noqa: E402
from streamer_rf.rf.jefimenko.synthetic import (  # noqa: E402
    dipole_far_field_E,
    dipole_radiation_series,
    static_gaussian_charge_series,
    steady_current_element_series,
)
from streamer_rf.rf.source.adapters import load_afivo_stage_e_source_csv, source_record_from_petsc_axisymmetric  # noqa: E402
from streamer_rf.rf.source.schema import SourceMetadata, SourceSeries  # noqa: E402


AFIVO_SHA = "a50b5508775086e90dfe423455fb58d812578410"


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def retarded_time(pos: np.ndarray, tr: float) -> float:
    return tr + float(np.linalg.norm(pos)) / C0


def static_charge_metrics() -> dict[str, float]:
    q = 1e-12
    series = static_gaussian_charge_series(total_charge_C=q, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x_far", 0.2, 0.0, 0.0)
    sample = evaluate_observer(series, obs, retarded_time(obs.position, 2e-9))
    expected = K_E * q / 0.2**2
    radii = np.array([0.08, 0.14, 0.25])
    amps = []
    for R in radii:
        o = Observer(f"R{R}", R, 0.0, 0.0)
        amps.append(np.linalg.norm(evaluate_observer(series, o, retarded_time(o.position, 2e-9)).E_rho_near))
    slope = float(np.polyfit(np.log(radii), np.log(amps), 1)[0])
    return {
        "E_expected_Vpm": expected,
        "E_simulated_Vpm": float(sample.E_total[0]),
        "relative_error": rel(float(sample.E_total[0]), expected),
        "B_norm_T": float(np.linalg.norm(sample.B_total)),
        "near_field_slope": slope,
    }


def magnetic_metrics() -> dict[str, float]:
    mz = 2e-6
    series = steady_current_element_series(current_moment_Am=mz, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x_far", 0.1, 0.0, 0.0)
    sample = evaluate_observer(series, obs, retarded_time(obs.position, 2e-9))
    expected = K_B * mz / 0.1**2
    return {
        "By_expected_T": expected,
        "By_simulated_T": float(sample.B_total[1]),
        "relative_error": rel(float(sample.B_total[1]), expected),
        "B_dJ_norm_T": float(np.linalg.norm(sample.B_dJ_radiation)),
    }


def radiation_metrics() -> dict[str, object]:
    omega = 2.0 * math.pi * 1.0e9
    p0 = 1e-18
    tr = 0.25e-9
    series = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 121))
    radii = np.array([0.2, 0.3, 0.5])
    e_perp = []
    eb_errors = []
    for R in radii:
        obs = Observer(f"R{R}", R, 0.0, 0.0)
        sample = evaluate_observer(series, obs, retarded_time(obs.position, tr))
        rhat = obs.position / np.linalg.norm(obs.position)
        evec = sample.E_total - np.dot(sample.E_total, rhat) * rhat
        e_perp.append(np.linalg.norm(evec))
        eb_errors.append(rel(np.linalg.norm(evec) / np.linalg.norm(sample.B_total), C0))
    slope = float(np.polyfit(np.log(radii), np.log(e_perp), 1)[0])

    angles = np.deg2rad([30.0, 60.0, 90.0])
    amps = []
    for theta in angles:
        pos = np.array([0.3 * math.sin(theta), 0.0, 0.3 * math.cos(theta)])
        sample = evaluate_observer(series, Observer(f"theta{theta}", *pos), retarded_time(pos, tr))
        rhat = pos / np.linalg.norm(pos)
        amps.append(np.linalg.norm(sample.E_total - np.dot(sample.E_total, rhat) * rhat))
    ratios = np.asarray(amps) / amps[-1]
    angular_error = float(np.max(np.abs(ratios - np.sin(angles))))
    pddot = -p0 * omega**2 * math.sin(omega * tr)
    expected_90 = dipole_far_field_E(math.pi / 2.0, 0.3, pddot)

    R1, R2 = 0.2, 0.35
    delay_error = abs((retarded_time(np.array([R2, 0.0, 0.0]), tr) - retarded_time(np.array([R1, 0.0, 0.0]), tr)) - (R2 - R1) / C0)

    coarse = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 61))
    obs = Observer("interp", 0.3, 0.0, 0.0)
    ec = np.linalg.norm(evaluate_observer(coarse, obs, retarded_time(obs.position, tr)).E_total)
    ef = np.linalg.norm(evaluate_observer(series, obs, retarded_time(obs.position, tr)).E_total)

    dMdt = np.array([0.0, 0.0, pddot])
    cm_errors = []
    for R in (0.12, 0.4):
        o = Observer(f"cm{R}", R, 0.0, 0.0)
        full = evaluate_observer(series, o, retarded_time(o.position, tr)).E_total
        approx = current_moment_radiation_approx(dMdt, o.position)
        cm_errors.append(float(np.linalg.norm(full - approx) / np.linalg.norm(full)))

    return {
        "radiation_1_over_R_slope": slope,
        "E_over_B_max_relative_error": float(max(eb_errors)),
        "angular_sin_theta_max_abs_error": angular_error,
        "E90_expected_Vpm": expected_90,
        "E90_simulated_Vpm": float(amps[-1]),
        "E90_relative_error": rel(float(amps[-1]), expected_90),
        "propagation_delay_error_s": float(delay_error),
        "interpolation_coarse_relative_error": rel(float(ec), expected_90),
        "interpolation_fine_relative_error": rel(float(ef), expected_90),
        "current_moment_error_R1": cm_errors[0],
        "current_moment_error_R2": cm_errors[1],
        "near_far_audit_R2": near_far_audit(series, np.array([0.4, 0.0, 0.0]), 1.0 / omega),
    }


def axisymmetric_metrics() -> dict[str, float]:
    r = np.linspace(0.2e-3, 1.0e-3, 5)
    z = np.linspace(-0.4e-3, 0.4e-3, 5)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    rho = 1e-7 * np.exp(-((rr - 0.6e-3) ** 2 + zz**2) / (2 * (0.25e-3) ** 2))
    zeros = np.zeros_like(rho)
    meta = SourceMetadata("axisym", "PETSc-2D", "test", 0.0, "axisymmetric_rz", 1.0, 1.0, "ring", "none", "off", "static axisymmetric")
    base = source_record_from_petsc_axisymmetric(r, z, r[1] - r[0], z[1] - z[0], rho, zeros, zeros, meta)
    series = SourceSeries(
        tuple(base.with_source_columns(base.columns["rho"], np.zeros((base.n_cells, 3)), time_s=t) for t in (0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    )
    obs = Observer("off_axis", 0.01, 0.003, 0.002)
    vals = {}
    for nphi in (8, 16, 32, 64, 128):
        vals[nphi] = evaluate_observer(rotate_axisymmetric_series(series, nphi), obs, retarded_time(obs.position, 2e-9)).E_total
    ref = vals[128]
    return {
        "Nphi_8_relative_error": float(np.linalg.norm(vals[8] - ref) / np.linalg.norm(ref)),
        "Nphi_16_relative_error": float(np.linalg.norm(vals[16] - ref) / np.linalg.norm(ref)),
        "Nphi_32_relative_error": float(np.linalg.norm(vals[32] - ref) / np.linalg.norm(ref)),
        "Nphi_reference": 128,
    }


def real_afivo_pipeline_smoke() -> dict[str, object]:
    manifest = json.loads((ROOT / "rf/source/audit/rf_source_manifest.json").read_text())
    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    files = [
        ROOT / f"solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_00000{i}.csv"
        for i in (2, 3, 4, 5)
    ]
    if not all(p.exists() for p in files):
        return {
            "ran": False,
            "SCIENTIFIC_RF_VALID": False,
            "reason": "No three-snapshot real Stage E source series is available.",
        }
    selected = []
    with files[-1].open(newline="") as handle:
        reader = csv.DictReader(handle)
        for i, row in enumerate(reader):
            activity = (
                abs(float(row["rho_Cpm3"]))
                + abs(float(row["Jx_Apm2"]))
                + abs(float(row["Jy_Apm2"]))
                + abs(float(row["Jz_Apm2"]))
                + 1e-30 * abs(float(row["ne_m3"]))
            )
            if activity > 0.0:
                selected.append(i)
            if len(selected) >= 2048:
                break
    if len(selected) < 128:
        return {
            "ran": False,
            "SCIENTIFIC_RF_VALID": False,
            "reason": "Could not find enough active rows in real Stage E source snapshots.",
        }
    records = [
        load_afivo_stage_e_source_csv(
            path,
            cfg,
            case_id="E3_aligned_pipeline_smoke",
            solver_version=AFIVO_SHA,
            geometry_id="E3_aligned_needle_pair",
            voltage_state="500 V constant",
            photoionization="off",
            row_indices=np.asarray(selected, dtype=int),
        )
        for path in files
    ]
    series = SourceSeries(tuple(records))
    obs = Observer("pipeline_observer", 2e-4, 0.0, 0.0)
    xyz = np.column_stack(
        (
            records[0].columns["x_center"],
            records[0].columns["y_center"],
            records[0].columns["z_center"],
        )
    )
    median_delay = float(np.median(np.linalg.norm(obs.position[None, :] - xyz, axis=1)) / C0)
    times = np.array([3.5e-12 + median_delay])
    samples = evaluate_waveform(series, [obs], times, source_manifest_id="PIPELINE_ONLY", chunk_size=4096)
    rows = [sample.to_row() for sample in samples]
    out_csv = ROOT / "rf/jefimenko/validation/real_afivo_pipeline_smoke.csv"
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return {
        "ran": True,
        "csv": str(out_csv),
        "sample": samples[0].to_dict(),
        "source_rows_used_per_snapshot": len(selected),
        "SCIENTIFIC_RF_VALID": False,
        "reason": "Stage E source J is drift-only and this is an older sparse 5 ps development source window.",
        "F1_scientific_spectrum_valid": manifest["usable_for_vhf_uhf_scientific_spectrum"],
    }


def main() -> None:
    out = ROOT / "rf/jefimenko/validation/stage_f2_validation.json"
    data = {
        "stage": "F2",
        "implemented_equations": "Full Jefimenko E/B with rho, drho/dt, dJ/dt, J near-field terms and native SourceRecord volume integration.",
        "static_charge": static_charge_metrics(),
        "steady_current_magnetic": magnetic_metrics(),
        "radiation": radiation_metrics(),
        "axisymmetric_phi": axisymmetric_metrics(),
        "real_afivo_pipeline_smoke": real_afivo_pipeline_smoke(),
    }
    out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
