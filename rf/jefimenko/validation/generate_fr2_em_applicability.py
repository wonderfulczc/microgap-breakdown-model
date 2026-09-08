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

from streamer_rf.rf.applicability import (  # noqa: E402
    current_moment_tau_M,
    current_weighted_source_scale,
    em_applicability_audit,
)
from streamer_rf.rf.source.adapters import load_petsc_stage4_field_source_csv  # noqa: E402


SOURCE_ROOT = ROOT / "results/stage_f4/source_gate_recovery"
F4_ROOT = ROOT / "rf/production/f4_attribution"
OUT = ROOT / "rf/jefimenko/validation/fr2_em_applicability"
SOLVER_VERSION = "eb1dd95+796336d+fr2"


def first_time(path: Path) -> float:
    with path.open(newline="") as handle:
        return float(next(csv.DictReader(handle))["time_s"])


def field_files(path: Path) -> list[Path]:
    files = [p for p in path.glob("fields_*.csv") if p.name != "fields_final.csv"]
    return sorted(files, key=first_time)


def load_record(path: Path, *, case_id: str, geometry_id: str, voltage_state: str) -> object:
    return load_petsc_stage4_field_source_csv(
        path,
        case_id=case_id,
        solver_version=SOLVER_VERSION,
        geometry_id=geometry_id,
        voltage_state=voltage_state,
        photoionization="on",
    )


def spatial_timeseries(source_dir: Path, *, case_id: str, geometry_id: str, voltage_state: str) -> pd.DataFrame:
    rows = []
    for path in field_files(source_dir):
        record = load_record(path, case_id=case_id, geometry_id=geometry_id, voltage_state=voltage_state)
        scale = current_weighted_source_scale(record)
        row = scale.to_dict()
        row["source_file"] = str(path.relative_to(ROOT))
        rows.append(row)
    return pd.DataFrame(rows)


def window_mask(times: np.ndarray, start: float, end: float) -> np.ndarray:
    return (times >= start) & (times <= end)


def summarize_spatial(spatial: pd.DataFrame, labels: pd.DataFrame) -> dict[str, object]:
    ok = spatial["status"] == "OK"
    out: dict[str, object] = {
        "valid_scale_samples": int(ok.sum()),
        "total_scale_samples": int(len(spatial)),
        "trusted_coverage": float(ok.mean()) if len(spatial) else 0.0,
        "all_window": _stats(spatial[ok]),
        "stages": {},
    }
    for row in labels.to_dict("records"):
        mask = ok & (spatial["time_s"] >= row["start_time_s"]) & (spatial["time_s"] <= row["end_time_s"])
        out["stages"][row["stage"]] = _stats(spatial[mask])
    return out


def _stats(frame: pd.DataFrame) -> dict[str, float | int]:
    if frame.empty:
        return {
            "sample_count": 0,
            "L95_median_m": math.nan,
            "L95_p95_m": math.nan,
            "L95_max_m": math.nan,
            "Lbounding_median_m": math.nan,
            "Lbounding_p95_m": math.nan,
        }
    return {
        "sample_count": int(len(frame)),
        "L95_median_m": float(frame["L95_m"].median()),
        "L95_p95_m": float(frame["L95_m"].quantile(0.95)),
        "L95_max_m": float(frame["L95_m"].max()),
        "Lbounding_median_m": float(frame["Lbounding_m"].median()),
        "Lbounding_p95_m": float(frame["Lbounding_m"].quantile(0.95)),
    }


def current_moment_window(cm: pd.DataFrame, start: float, end: float) -> tuple[np.ndarray, np.ndarray]:
    t = cm["time"].to_numpy(dtype=float)
    vals = cm["I_CM_electron"].to_numpy(dtype=float)
    mask = window_mask(t, start, end)
    if np.count_nonzero(mask) < 3:
        return t[mask], vals[mask]
    return t[mask], vals[mask]


def conservative_tau_from_derivative(cm: pd.DataFrame, start: float, end: float) -> float:
    t, m = current_moment_window(cm, start, end)
    if t.size < 3:
        return math.nan
    dm = np.gradient(m, t, edge_order=2)
    amp = float(np.quantile(np.abs(m), 0.95))
    rate = float(np.quantile(np.abs(dm), 0.95))
    if amp <= 0.0 or rate <= 0.0:
        return math.nan
    return amp / rate


def sampling_sensitivity(cm: pd.DataFrame, stats: dict[str, float], start: float, end: float) -> dict[str, object]:
    t, m = current_moment_window(cm, start, end)
    full = current_moment_tau_M(t, m)
    half = current_moment_tau_M(t[::2], m[::2]) if t.size >= 6 else current_moment_tau_M(np.array([0.0]), np.array([0.0]))
    full_em = em_applicability_audit(
        L95_p95_m=float(stats["L95_p95_m"]),
        L95_max_m=float(stats["L95_max_m"]),
        Lbounding_p95_m=float(stats["Lbounding_p95_m"]),
        tau=full,
        conservative_tau_s=conservative_tau_from_derivative(cm, start, end),
    )
    half_em = em_applicability_audit(
        L95_p95_m=float(stats["L95_p95_m"]),
        L95_max_m=float(stats["L95_max_m"]),
        Lbounding_p95_m=float(stats["Lbounding_p95_m"]),
        tau=half,
        conservative_tau_s=conservative_tau_from_derivative(pd.DataFrame({"time": t[::2], "I_CM_electron": m[::2]}), float(t[::2][0]), float(t[::2][-1])) if t.size >= 6 else math.nan,
    )
    return {
        "full_tau_M_s": full.tau_M_s,
        "every_second_tau_M_s": half.tau_M_s,
        "tau_relative_change": abs(half.tau_M_s - full.tau_M_s) / max(abs(full.tau_M_s), abs(half.tau_M_s), 1e-300)
        if math.isfinite(full.tau_M_s) and math.isfinite(half.tau_M_s)
        else math.nan,
        "full_epsilon_EM": full_em.epsilon_EM,
        "every_second_epsilon_EM": half_em.epsilon_EM,
        "epsilon_relative_change": abs(half_em.epsilon_EM - full_em.epsilon_EM)
        / max(abs(full_em.epsilon_EM), abs(half_em.epsilon_EM), 1e-300)
        if math.isfinite(full_em.epsilon_EM) and math.isfinite(half_em.epsilon_EM)
        else math.nan,
        "full_decision": full_em.decision,
        "every_second_decision": half_em.decision,
    }


def analyze_dataset(
    *,
    case_id: str,
    source_dir: Path,
    labels_path: Path,
    trust_path: Path,
    geometry_id: str,
    voltage_state: str,
) -> dict[str, object]:
    spatial = spatial_timeseries(source_dir, case_id=case_id, geometry_id=geometry_id, voltage_state=voltage_state)
    spatial.to_csv(OUT / f"{case_id}_timeseries.csv", index=False)
    labels = pd.read_csv(labels_path)
    cm = pd.read_csv(source_dir / "current_moment.csv")
    trust = json.loads(trust_path.read_text())

    start = float(labels["start_time_s"].min())
    end = float(labels["end_time_s"].max())
    spatial_window = spatial[(spatial["time_s"] >= start) & (spatial["time_s"] <= end)].copy()
    t, m = current_moment_window(cm, start, end)
    tau = current_moment_tau_M(t, m)
    spatial_summary = summarize_spatial(spatial_window, labels)
    all_stats = spatial_summary["all_window"]
    conservative_tau = conservative_tau_from_derivative(cm, start, end)
    em = em_applicability_audit(
        L95_p95_m=float(all_stats["L95_p95_m"]),
        L95_max_m=float(all_stats["L95_max_m"]),
        Lbounding_p95_m=float(all_stats["Lbounding_p95_m"]),
        tau=tau,
        conservative_tau_s=conservative_tau,
    )
    return {
        "case_id": case_id,
        "source_dir": str(source_dir.relative_to(ROOT)),
        "trusted_time_start_s": start,
        "trusted_time_end_s": end,
        "trusted_sample_count": int(t.size),
        "source_snapshot_count": int(len(spatial)),
        "spatial": spatial_summary,
        "tau_M": tau.to_dict(),
        "conservative_tau_s": conservative_tau,
        "em_audit": em.to_dict(),
        "temporal_sampling_sensitivity": sampling_sensitivity(cm, all_stats, start, end),
        "trusted_frequency_context": {
            "trusted_frequency_low_Hz": trust["trusted_frequency_low_Hz"],
            "trusted_frequency_high_Hz": trust["trusted_frequency_high_Hz"],
            "band_status": trust["band_status"],
            "VHF_UHF_1_3GHz_interpretation": "NOT_RESOLVED remains unchanged",
        },
    }


def overall_decision(datasets: list[dict[str, object]]) -> str:
    decisions = [d["em_audit"]["decision"] for d in datasets]
    if "SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED" in decisions:
        return "SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED"
    if "REVIEW_SELECTED_FULL_MAXWELL_REFERENCE" in decisions:
        return "SELECTED_FULL_MAXWELL_REFERENCE_RECOMMENDED"
    return "FULL_MAXWELL_REFERENCE_NOT_REQUIRED_FOR_CURRENT_TRUSTED_CASES"


def main() -> None:
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    stage4 = analyze_dataset(
        case_id="F4-P-stage4-left-isolated",
        source_dir=SOURCE_ROOT / "F4-P-streamer-propagation",
        labels_path=F4_ROOT / "F4-P-stage4-left-isolated_stage_labels.csv",
        trust_path=F4_ROOT / "F4-P-stage4-left-isolated_rf_trust_report.json",
        geometry_id="S4-LEFT-ISOLATED frozen 20 um axisymmetric propagation case",
        voltage_state="4.8e6 V/m frozen background field",
    )
    stage5 = analyze_dataset(
        case_id="F4-C-stage5-highfield-collision",
        source_dir=SOURCE_ROOT / "F4-C-interaction-collision",
        labels_path=F4_ROOT / "F4-C-stage5-highfield-collision_stage_labels.csv",
        trust_path=F4_ROOT / "F4-C-stage5-highfield-collision_rf_trust_report.json",
        geometry_id="S5-HIGHFIELD-COLLISION frozen 20 um axisymmetric collision case",
        voltage_state="6.4e6 V/m frozen background field",
    )
    summary = {
        "purpose": "quasi-static/full-Maxwell source-evolution applicability audit",
        "epsilon_EM_definition": "L_source/(c*tau_M)",
        "L_source_choice": "temporal p95 of current-weighted L95(t)",
        "engineering_bands": {
            "epsilon_EM <= 0.1": "LOW_EM_FEEDBACK_RISK",
            "0.1 < epsilon_EM <= 0.3": "REVIEW_SELECTED_FULL_MAXWELL_REFERENCE",
            "epsilon_EM > 0.3": "SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED",
            "scope": "engineering audit bands, not universal physical constants",
        },
        "datasets": [stage4, stage5],
        "overall_architecture_decision": overall_decision([stage4, stage5]),
        "runtime_s": time.perf_counter() - started,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUT / "fr2_em_applicability_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": summary["overall_architecture_decision"], "runtime_s": summary["runtime_s"]}, indent=2))


if __name__ == "__main__":
    main()
