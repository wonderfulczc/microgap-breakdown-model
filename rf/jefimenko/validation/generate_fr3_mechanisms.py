from __future__ import annotations

import csv
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.mechanisms import (  # noqa: E402
    classify_mechanism_status,
    compute_mechanism_decomposition,
    relative_change,
)


SOURCE_ROOT = ROOT / "results/stage_f4/source_gate_recovery"
F4_ROOT = ROOT / "rf/production/f4_attribution"
FR2 = ROOT / "rf/jefimenko/validation/fr2_em_applicability/fr2_em_applicability_summary.json"
OUT = ROOT / "rf/jefimenko/validation/fr3_mechanisms"
PRIMARY_DERIVATIVE_RADIUS = 5
ALTERNATIVE_DERIVATIVE_RADIUS = 2
DERIVATIVE_DEGREE = 2


DATASETS = [
    {
        "case_id": "F4-P-stage4-left-isolated",
        "source_dir": SOURCE_ROOT / "F4-P-streamer-propagation",
        "labels": F4_ROOT / "F4-P-stage4-left-isolated_stage_labels.csv",
        "trust": F4_ROOT / "F4-P-stage4-left-isolated_rf_trust_report.json",
        "summary": F4_ROOT / "F4-P-stage4-left-isolated_summary.json",
        "threshold": 0.2,
        "full_maxwell_note": "LOW_EM_FEEDBACK_RISK",
    },
    {
        "case_id": "F4-C-stage5-highfield-collision",
        "source_dir": SOURCE_ROOT / "F4-C-interaction-collision",
        "labels": F4_ROOT / "F4-C-stage5-highfield-collision_stage_labels.csv",
        "trust": F4_ROOT / "F4-C-stage5-highfield-collision_rf_trust_report.json",
        "summary": F4_ROOT / "F4-C-stage5-highfield-collision_summary.json",
        "threshold": 0.2,
        "full_maxwell_note": "FULL_MAXWELL_REFERENCE_PENDING",
    },
]


def first_time(path: Path) -> float:
    with path.open(newline="") as handle:
        return float(next(csv.DictReader(handle))["time_s"])


def field_files(path: Path) -> list[Path]:
    files = [p for p in path.glob("fields_*.csv") if p.name != "fields_final.csv"]
    return sorted(files, key=first_time)


def _axisymmetric_volume(df: pd.DataFrame) -> np.ndarray:
    r = df["r_m"].to_numpy(dtype=float)
    z = df["z_m"].to_numpy(dtype=float)
    ru = np.unique(np.round(r, decimals=18))
    zu = np.unique(np.round(z, decimals=18))
    dr = float(np.min(np.diff(ru))) if ru.size > 1 else 1.0
    dz = float(np.min(np.diff(zu))) if zu.size > 1 else 1.0
    return 2.0 * math.pi * r * dr * dz


def extract_primary_head(path: Path, *, threshold: float = 0.2, requested_polarity: int = 0) -> dict[str, object]:
    df = pd.read_csv(path)
    time_s = float(df["time_s"].iloc[0])
    rho = df["rho_C_m_3"].to_numpy(dtype=float)
    E = df["E_V_m"].to_numpy(dtype=float)
    if not np.all(np.isfinite(rho)) or not np.all(np.isfinite(E)):
        return {"time_s": time_s, "head_valid": False, "head_status": "NONFINITE_INPUT"}
    if not (0.0 < threshold <= 1.0):
        return {"time_s": time_s, "head_valid": False, "head_status": "INVALID_SEGMENTATION_PARAMETER"}
    if requested_polarity not in {-1, 0, 1}:
        return {"time_s": time_s, "head_valid": False, "head_status": "INVALID_SEGMENTATION_PARAMETER"}

    if requested_polarity == 0:
        peak_idx = int(np.argmax(np.abs(rho)))
        peak_abs = float(abs(rho[peak_idx]))
        polarity = 1 if rho[peak_idx] >= 0.0 else -1
    else:
        signed = requested_polarity * rho
        peak_idx = int(np.argmax(signed))
        peak_abs = float(signed[peak_idx])
        polarity = requested_polarity
    if not (peak_abs > 0.0):
        return {"time_s": time_s, "head_valid": False, "head_status": "INSUFFICIENT_CHARGE"}

    selected = polarity * rho >= threshold * peak_abs
    if not bool(selected[peak_idx]):
        return {"time_s": time_s, "head_valid": False, "head_status": "INSUFFICIENT_CHARGE"}
    ij_to_row = {(int(row.i), int(row.j)): idx for idx, row in df[["i", "j"]].reset_index(drop=True).iterrows()}
    peak_key = (int(df["i"].iloc[peak_idx]), int(df["j"].iloc[peak_idx]))
    selected_keys = {
        (int(i), int(j))
        for i, j in zip(df.loc[selected, "i"].to_numpy(), df.loc[selected, "j"].to_numpy(), strict=True)
    }
    stack = [peak_key]
    visited = {peak_key}
    comp_rows: list[int] = []
    while stack:
        key = stack.pop()
        comp_rows.append(ij_to_row[key])
        i, j = key
        for nb in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
            if nb in visited or nb not in selected_keys:
                continue
            visited.add(nb)
            stack.append(nb)
    comp = np.asarray(comp_rows, dtype=int)
    volume = _axisymmetric_volume(df)
    w = np.abs(rho[comp]) * volume[comp]
    wsum = float(np.sum(w, dtype=np.float64))
    if not (wsum > 0.0):
        return {"time_s": time_s, "head_valid": False, "head_status": "INSUFFICIENT_CHARGE"}
    r = df["r_m"].to_numpy(dtype=float)
    z = df["z_m"].to_numpy(dtype=float)
    q = float(np.sum(rho[comp] * volume[comp], dtype=np.float64))
    zc = float(np.sum(z[comp] * w, dtype=np.float64) / wsum)
    rmean = float(np.sum(r[comp] * w, dtype=np.float64) / wsum)
    rrms = float(math.sqrt(np.sum(r[comp] * r[comp] * w, dtype=np.float64) / wsum))
    peak_E = float(np.max(E[comp]))
    peak_rho = float(rho[peak_idx])
    return {
        "time_s": time_s,
        "source_file": str(path.relative_to(ROOT)),
        "head_valid": True,
        "head_status": "VALID",
        "head_polarity": polarity,
        "head_cell_count": int(comp.size),
        "head_charge_C": q,
        "head_z_m": zc,
        "head_r_mean_m": rmean,
        "head_r_rms_m": rrms,
        "head_peak_rho_C_m3": peak_rho,
        "head_peak_E_V_m": peak_E,
        "segmentation_fraction": float(comp.size / max(len(df), 1)),
        "segmentation_parameter": threshold,
        "position_method": "C-R3 same-polarity rho connected-component centroid from F4 source snapshots",
    }


def head_timeseries(source_dir: Path, *, threshold: float) -> pd.DataFrame:
    return pd.DataFrame([extract_primary_head(p, threshold=threshold) for p in field_files(source_dir)])


def label_window(path: Path) -> tuple[pd.DataFrame, float, float]:
    labels = pd.read_csv(path)
    return labels, float(labels["start_time_s"].min()), float(labels["end_time_s"].max())


def current_moment_on_head_grid(source_dir: Path, head: pd.DataFrame, start: float, end: float) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    cm = pd.read_csv(source_dir / "current_moment.csv")
    t_cm = cm["time"].to_numpy(dtype=float)
    m_cm = cm["I_CM_electron"].to_numpy(dtype=float)
    h = head[(head["time_s"] >= start) & (head["time_s"] <= end) & (head["head_valid"])].copy()
    th = h["time_s"].to_numpy(dtype=float)
    in_range = (th >= t_cm[0]) & (th <= t_cm[-1])
    h = h[in_range].copy()
    th = h["time_s"].to_numpy(dtype=float)
    if th.size == 0:
        return th, np.empty((0, 3)), {"head_frame": h, "interpolation_fraction": math.nan, "exact_common_count": 0}
    Mz = np.interp(th, t_cm, m_cm)
    exact = np.zeros(th.size, dtype=bool)
    for i, tt in enumerate(th):
        exact[i] = bool(np.any(np.isclose(t_cm, tt, rtol=0.0, atol=max(1e-24, 1e-12 * abs(tt)))))
    M = np.column_stack((np.zeros(th.size), np.zeros(th.size), Mz))
    return th, M, {
        "head_frame": h,
        "interpolation_fraction": float(1.0 - np.mean(exact)),
        "exact_common_count": int(np.count_nonzero(exact)),
    }


def vector_waveform_metrics(a: np.ndarray, b: np.ndarray, valid: np.ndarray) -> dict[str, float | int]:
    mask = valid & np.all(np.isfinite(a), axis=1) & np.all(np.isfinite(b), axis=1)
    if np.count_nonzero(mask) < 3:
        return {"sample_count": int(np.count_nonzero(mask)), "correlation": math.nan, "normalized_L2_error": math.nan}
    aa = a[mask].reshape(-1)
    bb = b[mask].reshape(-1)
    corr = float(np.corrcoef(aa, bb)[0, 1]) if np.std(aa) > 0.0 and np.std(bb) > 0.0 else math.nan
    err = float(np.linalg.norm(aa - bb) / max(np.linalg.norm(aa), 1e-300))
    return {"sample_count": int(np.count_nonzero(mask)), "correlation": corr, "normalized_L2_error": err}


def rel_metrics(primary: dict[str, object], candidate: dict[str, object]) -> dict[str, float]:
    pm = primary["metrics"]
    cm = candidate["metrics"]
    return {
        "rms_head_charge_evolution": relative_change(pm.rms_head_charge_evolution, cm.rms_head_charge_evolution),
        "rms_head_acceleration": relative_change(pm.rms_head_acceleration, cm.rms_head_acceleration),
        "rms_redistribution": relative_change(pm.rms_redistribution, cm.rms_redistribution),
        "normalized_rms_closure": relative_change(pm.normalized_rms_closure, cm.normalized_rms_closure),
    }


def build_timeseries_output(case_id: str, head: pd.DataFrame, decomp: dict[str, object], path: Path) -> None:
    n = len(head)
    out = pd.DataFrame(
        {
            "time_s": head["time_s"].to_numpy(dtype=float),
            "q_head_C": head["head_charge_C"].to_numpy(dtype=float),
            "qdot_head_C_s": decomp["qdot_C_s"][:, 0],
            "head_position_z_m": head["head_z_m"].to_numpy(dtype=float),
            "head_velocity_z_m_s": decomp["velocity_m_s"][:, 0, 2],
            "head_acceleration_z_m_s2": decomp["acceleration_m_s2"][:, 0, 2],
            "M_total_x_Am": decomp["M_total_Am"][:, 0],
            "M_total_y_Am": decomp["M_total_Am"][:, 1],
            "M_total_z_Am": decomp["M_total_Am"][:, 2],
            "M_head_x_Am": decomp["M_head_Am"][:, 0],
            "M_head_y_Am": decomp["M_head_Am"][:, 1],
            "M_head_z_Am": decomp["M_head_Am"][:, 2],
            "M_redis_x_Am": decomp["M_redis_Am"][:, 0],
            "M_redis_y_Am": decomp["M_redis_Am"][:, 1],
            "M_redis_z_Am": decomp["M_redis_Am"][:, 2],
            "dM_total_x_Am_s": decomp["dM_total_Am_s"][:, 0],
            "dM_total_y_Am_s": decomp["dM_total_Am_s"][:, 1],
            "dM_total_z_Am_s": decomp["dM_total_Am_s"][:, 2],
            "dM_charge_evolution_x_Am_s": decomp["dM_charge_evolution_Am_s"][:, 0],
            "dM_charge_evolution_y_Am_s": decomp["dM_charge_evolution_Am_s"][:, 1],
            "dM_charge_evolution_z_Am_s": decomp["dM_charge_evolution_Am_s"][:, 2],
            "dM_head_acceleration_x_Am_s": decomp["dM_head_acceleration_Am_s"][:, 0],
            "dM_head_acceleration_y_Am_s": decomp["dM_head_acceleration_Am_s"][:, 1],
            "dM_head_acceleration_z_Am_s": decomp["dM_head_acceleration_Am_s"][:, 2],
            "dM_redistribution_x_Am_s": decomp["dM_redistribution_Am_s"][:, 0],
            "dM_redistribution_y_Am_s": decomp["dM_redistribution_Am_s"][:, 1],
            "dM_redistribution_z_Am_s": decomp["dM_redistribution_Am_s"][:, 2],
            "closure_residual_x_Am_s": decomp["closure_residual_Am_s"][:, 0],
            "closure_residual_y_Am_s": decomp["closure_residual_Am_s"][:, 1],
            "closure_residual_z_Am_s": decomp["closure_residual_Am_s"][:, 2],
            "head_valid": np.ones(n, dtype=bool),
            "derivative_valid": decomp["derivative_valid"],
            "trust_valid": np.ones(n, dtype=bool),
            "head_representation": "PRIMARY_HEAD_ONLY",
        }
    )
    out.to_csv(path, index=False)


def analyze_dataset(spec: dict[str, object], fr2: dict[str, object]) -> dict[str, object]:
    labels, start, end = label_window(spec["labels"])
    trust = json.loads(Path(spec["trust"]).read_text())
    f4_summary = json.loads(Path(spec["summary"]).read_text())
    source_manifest = json.loads((Path(spec["source_dir"]) / "rf_source_manifest.json").read_text())

    head = head_timeseries(Path(spec["source_dir"]), threshold=float(spec["threshold"]))
    t, M, common = current_moment_on_head_grid(Path(spec["source_dir"]), head, start, end)
    h = common["head_frame"].reset_index(drop=True)
    positions = np.column_stack((np.zeros(t.size), np.zeros(t.size), h["head_z_m"].to_numpy(dtype=float)))
    q = h["head_charge_C"].to_numpy(dtype=float)
    primary = compute_mechanism_decomposition(
        t,
        M,
        q,
        positions,
        derivative_radius=PRIMARY_DERIVATIVE_RADIUS,
        derivative_degree=DERIVATIVE_DEGREE,
    )
    alt = (
        compute_mechanism_decomposition(
            t,
            M,
            q,
            positions,
            derivative_radius=ALTERNATIVE_DERIVATIVE_RADIUS,
            derivative_degree=DERIVATIVE_DEGREE,
        )
        if t.size >= 2 * ALTERNATIVE_DERIVATIVE_RADIUS + 5
        else primary
    )

    threshold_metrics: dict[str, dict[str, object]] = {}
    for thr in (0.1, 0.2, 0.3):
        th = head_timeseries(Path(spec["source_dir"]), threshold=thr)
        tt, MM, cc = current_moment_on_head_grid(Path(spec["source_dir"]), th, start, end)
        hh = cc["head_frame"].reset_index(drop=True)
        if tt.size >= 9:
            rr = np.column_stack((np.zeros(tt.size), np.zeros(tt.size), hh["head_z_m"].to_numpy(dtype=float)))
            dec = compute_mechanism_decomposition(
                tt,
                MM,
                hh["head_charge_C"].to_numpy(dtype=float),
                rr,
                derivative_radius=PRIMARY_DERIVATIVE_RADIUS,
                derivative_degree=DERIVATIVE_DEGREE,
            )
            threshold_metrics[str(thr)] = dec["metrics"].to_dict()

    primary_csv = OUT / f"{spec['case_id']}_mechanisms.csv"
    build_timeseries_output(str(spec["case_id"]), h, primary, primary_csv)
    recon = primary["dM_charge_evolution_Am_s"] + primary["dM_head_acceleration_Am_s"] + primary["dM_redistribution_Am_s"]
    mechanism_sum_cross = vector_waveform_metrics(primary["dM_total_Am_s"], recon, primary["derivative_valid"])
    deriv_change = rel_metrics(primary, alt)
    seg_change = {}
    if "0.1" in threshold_metrics and "0.3" in threshold_metrics and "0.2" in threshold_metrics:
        for key in ("rms_head_charge_evolution", "rms_head_acceleration", "rms_redistribution", "normalized_rms_closure"):
            ref = float(threshold_metrics["0.2"][key])
            seg_change[key] = max(relative_change(ref, float(threshold_metrics["0.1"][key])), relative_change(ref, float(threshold_metrics["0.3"][key])))

    status = classify_mechanism_status(
        primary["metrics"],
        derivative_relative_changes=deriv_change,
        segmentation_relative_changes=seg_change,
        field_correlation=f4_summary["M_dMdt_vs_jefimenko"]["correlation"],
    )
    return {
        "case_id": spec["case_id"],
        "source_dir": str(Path(spec["source_dir"]).relative_to(ROOT)),
        "J_RF_definition": source_manifest.get("CURRENT_SOURCE_PROVENANCE", source_manifest.get("current_provenance", "CONTINUITY_CONSISTENT_FINITE_VOLUME_FLUX")),
        "trusted_time_window_s": [start, end],
        "trusted_frequency_low_Hz": trust["trusted_frequency_low_Hz"],
        "trusted_frequency_high_Hz": trust["trusted_frequency_high_Hz"],
        "common_sample_count": int(t.size),
        "head_valid_coverage": float(t.size / max(len(head[(head["time_s"] >= start) & (head["time_s"] <= end)]), 1)),
        "interpolation_fraction": common["interpolation_fraction"],
        "exact_common_count": common["exact_common_count"],
        "head_representation": "PRIMARY_HEAD_ONLY",
        "head_segmentation": "C-R3 net-rho same-polarity connected component, relative threshold 0.2",
        "primary_derivative_operator": {
            "type": "local polynomial derivative",
            "degree": DERIVATIVE_DEGREE,
            "radius_samples": PRIMARY_DERIVATIVE_RADIUS,
            "validity": "centered stencil only; invalid across insufficient history or invalid head samples",
        },
        "alternative_derivative_operator": {
            "type": "local polynomial derivative",
            "degree": DERIVATIVE_DEGREE,
            "radius_samples": ALTERNATIVE_DERIVATIVE_RADIUS,
        },
        "primary_metrics": primary["metrics"].to_dict(),
        "alternative_derivative_metrics": alt["metrics"].to_dict(),
        "derivative_relative_changes": deriv_change,
        "segmentation_threshold_metrics": threshold_metrics,
        "segmentation_relative_changes": seg_change,
        "current_moment_full_jefimenko_cross_check": f4_summary["M_dMdt_vs_jefimenko"],
        "mechanism_sum_current_moment_cross_check": mechanism_sum_cross,
        "mechanism_status": status.to_dict(),
        "full_maxwell_note": spec["full_maxwell_note"],
        "output_csv": str(primary_csv.relative_to(ROOT)),
        "stage_labels": labels.to_dict("records"),
    }


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    fr2 = json.loads(FR2.read_text())
    datasets = [analyze_dataset(spec, fr2) for spec in DATASETS]
    summary = {
        "F_R3_STATUS": "PASS_CANDIDATE",
        "scope": "time-domain reduced-order current-moment mechanism decomposition; no frequency attribution",
        "F_R2_checkpoint_commit": "020729f",
        "F_R2_overall_decision": fr2["overall_architecture_decision"],
        "datasets": datasets,
    }
    summary["runtime_s"] = time.perf_counter() - started
    summary["peak_rss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    summary["new_output_bytes"] = sum(p.stat().st_size for p in OUT.glob("*") if p.is_file())
    (OUT / "fr3_mechanism_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
