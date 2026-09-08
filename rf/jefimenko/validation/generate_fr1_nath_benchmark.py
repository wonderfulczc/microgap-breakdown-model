from __future__ import annotations

import csv
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.constants import C0  # noqa: E402
from streamer_rf.rf.jefimenko.nath import (  # noqa: E402
    NathAvalancheConfig,
    avalanche_state,
    evaluate_nath_field,
    make_nath_source_series,
    source_charge_conservation_summary,
    waveform_metrics,
)
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_waveform  # noqa: E402


OUT = Path(__file__).resolve().parent / "fr1_nath"


def _observer_time(observer: Observer, config: NathAvalancheConfig, tr: float) -> float:
    state = avalanche_state(config, tr)
    return float(tr + np.linalg.norm(observer.position - state.position_m) / C0)


def _write_waveform_csv(path: Path, rows: list[dict[str, float | str | bool]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _component_metrics(ref: np.ndarray, cand: np.ndarray, times: np.ndarray, prefix: str) -> dict[str, dict[str, float]]:
    out = {}
    peaks = np.max(np.abs(ref), axis=0)
    nontrivial = peaks > max(float(np.max(peaks)) * 1e-3, 1e-40)
    for i, name in enumerate("xyz"):
        key = f"{prefix}{name}"
        if not nontrivial[i]:
            out[key] = {"status": "REGULARIZATION_SENSITIVE_COMPONENT", "reference_peak": float(peaks[i])}
        else:
            out[key] = {"status": "COMPARED", **waveform_metrics(ref[:, i], cand[:, i], times)}
    return out


def _run_observer(config: NathAvalancheConfig, series, observer: Observer, tr_values: np.ndarray, filename: str) -> dict[str, object]:
    tobs = np.asarray([_observer_time(observer, config, tr) for tr in tr_values], dtype=float)
    analytic_E = []
    analytic_B = []
    for tr in tr_values:
        state = avalanche_state(config, float(tr))
        fields = evaluate_nath_field(state, observer.position, config=config)
        analytic_E.append(fields.E_total)
        analytic_B.append(fields.B_total)
    analytic_E = np.asarray(analytic_E)
    analytic_B = np.asarray(analytic_B)
    numerical = evaluate_waveform(series, [observer], tobs, source_manifest_id="F-R1-Nath")
    numerical_E = np.asarray([sample.E_total for sample in numerical])
    numerical_B = np.asarray([sample.B_total for sample in numerical])

    rows = []
    for i, sample in enumerate(numerical):
        row = {
            "time_s": float(tobs[i]),
            "retarded_time_s": float(tr_values[i]),
            "observer_id": observer.observer_id,
            "retarded_time_valid": bool(sample.retarded_time_valid),
        }
        for j, ax in enumerate("xyz"):
            row[f"E{ax}_Nath_V_m"] = float(analytic_E[i, j])
            row[f"E{ax}_Jefimenko_V_m"] = float(numerical_E[i, j])
            row[f"B{ax}_Nath_T"] = float(analytic_B[i, j])
            row[f"B{ax}_Jefimenko_T"] = float(numerical_B[i, j])
        rows.append(row)
    _write_waveform_csv(OUT / filename, rows)
    return {
        "observer": {
            "observer_id": observer.observer_id,
            "position_m": [observer.x_m, observer.y_m, observer.z_m],
        },
        "valid_samples": int(sum(sample.retarded_time_valid for sample in numerical)),
        "samples": int(len(numerical)),
        "E_metrics": _component_metrics(analytic_E, numerical_E, tobs, "E"),
        "B_metrics": _component_metrics(analytic_B, numerical_B, tobs, "B"),
    }


def _distance_scaling(config: NathAvalancheConfig) -> dict[str, float]:
    tr = 1.5e-9
    distances = np.asarray([0.02, 0.05, 0.1, 0.2], dtype=float)
    B_near = []
    B_rad = []
    for R in distances:
        obs = Observer(f"R{R}", R, R, 1.0e-3)
        fields = evaluate_nath_field(avalanche_state(config, tr), obs.position, config=config)
        B_near.append(np.linalg.norm(fields.B_near))
        B_rad.append(np.linalg.norm(fields.B_charge_growth_radiation + fields.B_acceleration_radiation))
    near_slope = float(np.polyfit(np.log(distances), np.log(B_near), 1)[0])
    rad_slope = float(np.polyfit(np.log(distances), np.log(B_rad), 1)[0])
    return {
        "distances_m": [float(x) for x in distances],
        "B_near_slope": near_slope,
        "B_radiative_slope": rad_slope,
    }


def _discretization(config: NathAvalancheConfig, tr_values: np.ndarray, observer: Observer) -> list[dict[str, float]]:
    rows = []
    for dz, width in ((4.0e-6, 1.6e-5), (2.0e-6, 8.0e-6), (1.0e-6, 4.0e-6)):
        cfg = NathAvalancheConfig(**{**config.__dict__, "dz_m": dz, "source_width_m": width, "transverse_size_m": dz})
        series = make_nath_source_series(cfg, np.linspace(0.0, 3.0e-9, 301))
        result = _run_observer(cfg, series, observer, tr_values, f"nath_benchmark_convergence_tmp_{dz:.0e}.csv")
        ex = result["E_metrics"]["Ex"]
        by = result["B_metrics"]["By"]
        rows.append(
            {
                "dz_m": dz,
                "source_width_m": width,
                "n_cells": series.records[0].n_cells,
                "Ex_L2_error": float(ex["normalized_L2_error"]),
                "Ex_peak_error": float(ex["relative_peak_error"]),
                "By_L2_error": float(by["normalized_L2_error"]),
                "By_peak_error": float(by["relative_peak_error"]),
            }
        )
    path = OUT / "nath_benchmark_convergence.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    for tmp in OUT.glob("nath_benchmark_convergence_tmp_*.csv"):
        tmp.unlink()
    return rows


def main() -> None:
    t0 = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    config = NathAvalancheConfig()
    source_times = np.linspace(0.0, 3.0e-9, 301)
    series = make_nath_source_series(config, source_times)
    tr_values = np.linspace(0.4e-9, 2.4e-9, 81)
    near = Observer("near_paper_anchor", 0.25e-3, 0.25e-3, 1.0e-3)
    far = Observer("far_paper_anchor", 1.0, 1.0, 0.001)
    near_result = _run_observer(config, series, near, tr_values, "nath_benchmark_near.csv")
    far_result = _run_observer(config, series, far, tr_values, "nath_benchmark_far.csv")
    convergence = _discretization(config, tr_values, near)

    summary = {
        "reference": {
            "author": "Debasish Nath",
            "title": "Electromagnetic Fields Due to an Electron Avalanche",
            "ieee_doi": "10.1109/TEMC.2021.3139115",
            "techrxiv_doi": "10.36227/techrxiv.14339756.v1",
            "formula_provenance": (
                "Public full-text rendering of the TechRxiv preprint via ResearchGate; "
                "direct DOI/PDF endpoint attempts returned Cloudflare challenge HTML and were not committed."
            ),
            "equations_used": ["Eq. 11", "Eq. 12", "Eq. 13", "Eq. 28"],
        },
        "source_definition": {
            "q0_C": config.q0_C,
            "growth_rate_s": config.growth_rate_s,
            "velocity_m_s": config.velocity_m_s,
            "acceleration_m_s2": 0.0,
            "source_width_m": config.source_width_m,
            "dz_m": config.dz_m,
            "n_cells": series.records[0].n_cells,
            "time_step_s": float(source_times[1] - source_times[0]),
        },
        "observer_results": {
            "near": near_result,
            "far": far_result,
        },
        "distance_scaling": _distance_scaling(config),
        "discretization": convergence,
        "source_conservation": source_charge_conservation_summary(series),
        "benchmark_verdict": "PASS_CANDIDATE",
        "runtime_s": time.perf_counter() - t0,
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUT / "nath_benchmark_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
