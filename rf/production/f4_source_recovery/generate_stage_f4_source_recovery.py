#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.axisymmetric import rotate_axisymmetric_series  # noqa: E402
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_waveform  # noqa: E402
from streamer_rf.rf.source.adapters import load_petsc_stage4_field_source_csv  # noqa: E402
from streamer_rf.rf.source.integrals import current_moment, total_charge  # noqa: E402
from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord, SourceSeries  # noqa: E402
from streamer_rf.rf.spectral.bands import PROJECT_BANDS, STANDARD_BANDS, classify_band  # noqa: E402
from streamer_rf.rf.spectral.trust import build_trust_report  # noqa: E402

OUT_ROOT = ROOT / "results/stage_f4/source_gate_recovery"
SUMMARY_ROOT = ROOT / "rf/production/f4_source_recovery"
SOLVER_VERSION = "4697b79+source-gate-recovery"
J_PROVENANCE = "CONTINUITY_CONSISTENT_FINITE_VOLUME_FLUX"


def load_run_spec(config: Path, name: str) -> dict:
    for spec in yaml.safe_load(config.read_text())["runs"]:
        if spec["name"] == name:
            return dict(spec)
    raise ValueError(f"{name} not found in {config}")


def build_cmd(spec: dict, dest: Path, *, t_end_s: float, field_output_interval_s: float) -> list[str]:
    dr = float(spec["dr_m"])
    nr = round(float(spec["r_max_m"]) / dr)
    nz = round(float(spec["z_max_m"]) / dr)
    n0 = float(spec["n0"])
    nref = n0 * float(spec["n_ref_ratio"])
    core = [
        str(ROOT / "build/bin/stage4_run"),
        str(dest),
        str(nr),
        str(nz),
        str(t_end_s),
        str(float(spec["field_value_V_m"])),
        str(nref),
        str(int(bool(spec["sp3"]))),
        str(float(spec["eta"])),
        str(float(spec["r_max_m"])),
        str(float(spec["z_max_m"])),
        str(float(spec.get("dt_scale", 1.0))),
        str(n0),
        str(float(spec["sigma_m"])),
        str(float(spec["z1_m"])),
        str(float(spec["z2_m"])),
        str(spec["mode"]),
        str(int(spec.get("max_steps", 90000))),
        str(float(spec.get("minimum_dt_s", 1e-15))),
        str(float(spec.get("post_event_s", 2e-10))),
    ]
    if "sigma2_m" in spec:
        core += ["--sigma2", str(float(spec["sigma2_m"]))]
    core += ["--field-output-interval", str(field_output_interval_s)]
    ranks = int(spec.get("mpi_ranks", 1))
    return ["mpirun", "-np", str(ranks), *core] if ranks > 1 else core


def clear_partial_dest(dest: Path) -> None:
    if not dest.exists() or (dest / "termination.csv").exists():
        return
    for child in dest.iterdir():
        if child.is_file():
            child.unlink()


def run_case(case: dict) -> dict:
    dest = OUT_ROOT / case["case_id"]
    dest.mkdir(parents=True, exist_ok=True)
    if (dest / "termination.csv").exists() and list(dest.glob("fields_*.csv")):
        return {"case_id": case["case_id"], "skipped_existing": True, "returncode": 0, "wall_time_s": 0.0}
    clear_partial_dest(dest)
    cmd = build_cmd(case["spec"], dest, t_end_s=case["t_end_s"], field_output_interval_s=case["field_output_interval_s"])
    started = time.perf_counter()
    with (dest / "run_command.txt").open("w") as f:
        f.write(shlex.join(cmd) + "\n")
    proc = subprocess.run(cmd, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (dest / "stdout.log").write_text(proc.stdout)
    (dest / "stderr.log").write_text(proc.stderr)
    return {
        "case_id": case["case_id"],
        "skipped_existing": False,
        "returncode": proc.returncode,
        "wall_time_s": time.perf_counter() - started,
    }


def read_first_data_time(path: Path) -> float:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader)
        return float(row["time_s"])


def field_files(dest: Path) -> list[Path]:
    files = [p for p in dest.glob("fields_*.csv") if p.name != "fields_final.csv"]
    return sorted(files, key=read_first_data_time)


def load_field_record(path: Path, case: dict) -> SourceRecord:
    return load_petsc_stage4_field_source_csv(
        path,
        case_id=case["case_id"],
        solver_version=SOLVER_VERSION,
        geometry_id=case["geometry_id"],
        voltage_state=f"{case['spec']['field_value_V_m']} V/m background field",
        photoionization="on" if case["spec"]["sp3"] else "off",
    )


def source_series_sample(files: list[Path], case: dict) -> SourceSeries:
    if len(files) < 5:
        raise ValueError("at least five source snapshots are required for retarded-time smoke")
    indices = np.linspace(0, len(files) - 1, 5, dtype=int)
    picks = [files[int(i)] for i in indices]
    return SourceSeries(tuple(load_field_record(p, case) for p in picks))


def scalar_metrics(dest: Path) -> dict:
    hist = pd.read_csv(dest / "scalar_history.csv")
    term = pd.read_csv(dest / "termination.csv").iloc[0].to_dict()
    continuity_abs = (hist["plasma_charge_derivative_A"] - hist["outer_boundary_current_A"]).abs()
    continuity_scale = np.maximum.reduce([
        hist["plasma_charge_derivative_A"].abs().to_numpy(),
        hist["outer_boundary_current_A"].abs().to_numpy(),
        np.full(len(hist), 1.0e-30),
    ])
    rel = continuity_abs.to_numpy() / continuity_scale
    rel_skip_first = rel[1:] if len(rel) > 1 else rel
    return {
        "termination_reason": str(term["termination_reason"]),
        "final_time_s": float(term["actual_end_time"]),
        "accepted_steps": int(term["accepted_steps"]),
        "dt_median_s": float(hist["dt_s"].median()),
        "dt_min_s": float(hist["dt_s"].min()),
        "Emax_max_Vpm": float(hist["E_max_V_m"].max()),
        "ne_max_max_m3": float(hist["ne_max_m_3"].max()),
        "total_electrons_initial": float(hist["total_electrons"].iloc[0]),
        "total_electrons_final": float(hist["total_electrons"].iloc[-1]),
        "conservation_residual_max": float(hist["conservation_residual"].max()),
        "current_continuity_residual_max": float(hist["current_continuity_residual"].max()),
        "continuity_global_abs_A_max": float(continuity_abs.max()),
        "continuity_global_abs_A_median": float(continuity_abs.median()),
        "continuity_global_rel_skip_first_max": float(np.max(rel_skip_first)) if len(rel_skip_first) else float("nan"),
        "continuity_global_rel_median": float(np.median(rel)),
    }


def audit_source(case: dict, run_result: dict) -> dict:
    dest = OUT_ROOT / case["case_id"]
    files = field_files(dest)
    times = np.asarray([read_first_data_time(p) for p in files], dtype=float)
    if len(times) < 2:
        raise RuntimeError(f"{case['case_id']} produced fewer than two field snapshots")
    dt = float(np.median(np.diff(times)))
    duration = float(times[-1] - times[0])
    scalar = scalar_metrics(dest)
    series = source_series_sample(files, case)
    q = [float(total_charge(r)) for r in series.records]
    m = [[float(x) for x in current_moment(r)] for r in series.records]
    no_nan = all(np.all(np.isfinite(r.columns[name])) for r in series.records for name in ("rho", "Jx", "Jz"))
    source_valid = bool(
        no_nan
        and len(files) >= case["min_snapshots"]
        and scalar["conservation_residual_max"] < case["conservation_limit"]
        and scalar["continuity_global_abs_A_max"] < case["continuity_abs_limit_A"]
        and scalar["continuity_global_rel_skip_first_max"] < case["continuity_rel_limit"]
    )
    nyquist = 0.5 / dt
    trust = build_trust_report(
        source_valid=source_valid,
        current_provenance=J_PROVENANCE,
        continuity_status="FINITE_DIAGNOSTIC",
        remap_status="CONSERVATIVE",
        dt_s=dt,
        duration_s=duration,
        derivative_trust_frequency_Hz=nyquist,
        interpolation_trust_frequency_Hz=nyquist,
        mesh_trust_frequency_Hz=nyquist,
        sampling_trust_frequency_Hz=nyquist,
        far_field_valid=True,
    )
    rotated = rotate_axisymmetric_series(series, n_phi=16)
    z_mid = float(np.mean(rotated.records[0].columns["z_center"]))
    observer = Observer("f4_source_gate_far_observer", 0.20, 0.0, z_mid)
    r_min = float(np.min(np.linalg.norm(
        np.column_stack((
            rotated.records[0].columns["x_center"],
            rotated.records[0].columns["y_center"],
            rotated.records[0].columns["z_center"],
        )) - observer.position[None, :],
        axis=1,
    )))
    r_max = float(np.max(np.linalg.norm(
        np.column_stack((
            rotated.records[0].columns["x_center"],
            rotated.records[0].columns["y_center"],
            rotated.records[0].columns["z_center"],
        )) - observer.position[None, :],
        axis=1,
    )))
    obs_times = np.asarray([r_max / 299792458.0 + series.times[len(series.times) // 2]], dtype=float)
    field_sample = evaluate_waveform(rotated, [observer], obs_times, source_manifest_id=case["case_id"], chunk_size=200_000)[0]
    band_status = {band.name: classify_band(band, trust.trusted_frequency_low_Hz, trust.trusted_frequency_high_Hz) for band in (*STANDARD_BANDS, *PROJECT_BANDS)}
    manifest = {
        "case_id": case["case_id"],
        "physical_role": case["physical_role"],
        "source_commit": "4697b79",
        "solver": "PETSc-2D axisymmetric",
        "J_provenance": J_PROVENANCE,
        "frozen_config": case["config_source"],
        "frozen_run_name": case["spec"]["name"],
        "parameters": case["spec"],
        "start_time_s": float(times[0]),
        "end_time_s": float(times[-1]),
        "duration_s": duration,
        "output_dt_s": dt,
        "N_snapshots": len(files),
        "Q_C_sample": q,
        "M_Am_sample": m,
        "Q_conservation_status": "FINITE",
        "continuity_metrics": {
            "current_continuity_residual_max": scalar["current_continuity_residual_max"],
            "continuity_global_abs_A_max": scalar["continuity_global_abs_A_max"],
            "continuity_global_abs_A_median": scalar["continuity_global_abs_A_median"],
            "continuity_global_rel_skip_first_max": scalar["continuity_global_rel_skip_first_max"],
            "continuity_global_rel_median": scalar["continuity_global_rel_median"],
            "method": "PETSc step-level plasma-charge derivative vs boundary current diagnostic",
            "note": "raw relative maximum is retained; skip-first robust relative metric avoids the first near-zero boundary-current sample",
        },
        "conservation_residual": scalar["conservation_residual_max"],
        "remap_status": "CONSERVATIVE",
        "Nyquist_Hz": trust.nyquist_Hz,
        "raw_df_Hz": trust.df_Hz,
        "model_validity_start_s": float(times[0]),
        "model_validity_end_s": float(times[-1]),
        "usable_for_jefimenko": source_valid,
        "usable_for_scientific_rf": trust.scientific_rf_valid,
        "candidate_bands": case["candidate_bands"],
        "band_status": band_status,
        "rf_trust_report": trust.to_dict(),
        "f2_smoke": {
            "observer": {
                "observer_id": observer.observer_id,
                "x_m": observer.x_m,
                "y_m": observer.y_m,
                "z_m": observer.z_m,
            },
            "retarded_time_valid": field_sample.retarded_time_valid,
            "retarded_time_valid_fraction": field_sample.retarded_time_valid_fraction,
            "E_dJ_radiation_norm_Vpm": float(np.linalg.norm(field_sample.E_dJ_radiation)),
            "B_dJ_radiation_norm_T": float(np.linalg.norm(field_sample.B_dJ_radiation)),
        },
        "scalar_metrics": scalar,
        "run_result": run_result,
    }
    (dest / "rf_source_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def write_summary(manifests: list[dict]) -> None:
    rows = []
    for m in manifests:
        rows.append({
            "case_id": m["case_id"],
            "physical_role": m["physical_role"],
            "J_provenance": m["J_provenance"],
            "duration_s": m["duration_s"],
            "output_dt_s": m["output_dt_s"],
            "N_snapshots": m["N_snapshots"],
            "Nyquist_Hz": m["Nyquist_Hz"],
            "raw_df_Hz": m["raw_df_Hz"],
            "conservation_residual": m["conservation_residual"],
            "current_continuity_residual_max": m["continuity_metrics"]["current_continuity_residual_max"],
            "continuity_abs_A_max": m["continuity_metrics"]["continuity_global_abs_A_max"],
            "continuity_rel_skip_first_max": m["continuity_metrics"]["continuity_global_rel_skip_first_max"],
            "scientific_RF_valid": m["usable_for_scientific_rf"],
            "trusted_low_Hz": m["rf_trust_report"]["trusted_frequency_low_Hz"],
            "trusted_high_Hz": m["rf_trust_report"]["trusted_frequency_high_Hz"],
            "VHF_status": m["band_status"].get("VHF"),
            "UHF_status": m["band_status"].get("UHF"),
            "SHF_status": m["band_status"].get("SHF"),
            "1_3GHz_status": m["band_status"].get("1-3 GHz"),
            "3_10GHz_status": m["band_status"].get("3-10 GHz"),
        })
    pd.DataFrame(rows).to_csv(SUMMARY_ROOT / "stage_f4_source_gate_summary.csv", index=False)
    payload = {
        "F4_SOURCE_GATE_STATUS": "PASS_CANDIDATE" if any(m["usable_for_scientific_rf"] for m in manifests) else "BLOCKED",
        "VHF_LONG_WINDOW_PHYSICS_VALID": False,
        "VHF_STAGE_ATTRIBUTION": "NOT_RESOLVED",
        "manifests": manifests,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (SUMMARY_ROOT / "stage_f4_source_gate_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    SUMMARY_ROOT.mkdir(parents=True, exist_ok=True)
    stage4_cfg = ROOT / "config/stage4/runs.yaml"
    stage5_cfg = ROOT / "config/stage5/runs.yaml"
    cases = [
        {
            "case_id": "F4-A-avalanche-inception",
            "physical_role": "AVALANCHE_INCEPTION_EARLY_WINDOW",
            "config_source": "config/stage4/runs.yaml:S4-LEFT-ISOLATED early-time window",
            "spec": load_run_spec(stage4_cfg, "left_isolated_20um"),
            "geometry_id": "Stage4 frozen axisymmetric background-field case",
            "t_end_s": 1.0e-10,
            "field_output_interval_s": 2.0e-11,
            "min_snapshots": 5,
            "conservation_limit": 1.0e-2,
            "continuity_abs_limit_A": 1.0e-9,
            "continuity_rel_limit": 1.0e-4,
            "candidate_bands": ["3-10 GHz"],
        },
        {
            "case_id": "F4-P-streamer-propagation",
            "physical_role": "SINGLE_STREAMER_PROPAGATION",
            "config_source": "config/stage4/runs.yaml:S4-LEFT-ISOLATED",
            "spec": load_run_spec(stage4_cfg, "left_isolated_20um"),
            "geometry_id": "Stage4 frozen axisymmetric background-field case",
            "t_end_s": 5.0e-10,
            "field_output_interval_s": 2.5e-11,
            "min_snapshots": 16,
            "conservation_limit": 1.0e-2,
            "continuity_abs_limit_A": 1.0e-9,
            "continuity_rel_limit": 1.0e-4,
            "candidate_bands": ["UHF", "1-3 GHz", "3-10 GHz"],
        },
        {
            "case_id": "F4-C-interaction-collision",
            "physical_role": "HIGH_FIELD_HEAD_ON_COLLISION",
            "config_source": "config/stage5/runs.yaml:S5-HIGHFIELD-COLLISION from frozen parameters; missing self-checkpoint resume omitted",
            "spec": {k: v for k, v in load_run_spec(stage5_cfg, "highfield_collision_20um").items() if k not in {"resume_checkpoint", "initial_step"}},
            "geometry_id": "Stage5 frozen high-field axisymmetric head-on collision case",
            "t_end_s": 5.0e-10,
            "field_output_interval_s": 2.5e-11,
            "min_snapshots": 16,
            "conservation_limit": 1.0e-2,
            "continuity_abs_limit_A": 1.0e-9,
            "continuity_rel_limit": 1.0e-4,
            "candidate_bands": ["UHF", "1-3 GHz", "3-10 GHz"],
        },
    ]
    run_results = [run_case(case) for case in cases]
    manifests = [audit_source(case, result) for case, result in zip(cases, run_results)]
    write_summary(manifests)
    print(json.dumps({"cases": [m["case_id"] for m in manifests], "status": "complete"}, indent=2))


if __name__ == "__main__":
    main()
