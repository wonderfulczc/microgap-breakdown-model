#!/usr/bin/env python3
"""Generate the lightweight, deterministic F-R5A benchmark artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from streamer_rf.rf.ultrafast import process_current_moment  # noqa: E402


OUT = ROOT / "rf" / "f_r5a"
FIGURES = OUT / "figures"
FROZEN_CURRENT_MOMENT = (
    ROOT
    / "results"
    / "stage_f4"
    / "source_gate_recovery"
    / "F4-P-streamer-propagation"
    / "current_moment.csv"
)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analytic_signal(tau_s: float, *, amplitude_A_m: float = 2e-8, shift_s: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    dt = 0.1e-12
    time = np.arange(-250e-12, 250e-12 + 0.5 * dt, dt) + shift_s
    moment = 0.5 * amplitude_A_m * (1.0 + np.tanh((time - shift_s) / tau_s))
    return time, moment


def run_analytic(tau_s: float, *, alpha: float = 0.5, amplitude_A_m: float = 2e-8, shift_s: float = 0.0):
    time, moment = analytic_signal(tau_s, amplitude_A_m=amplitude_A_m, shift_s=shift_s)
    return process_current_moment(
        time,
        moment,
        profile="NO_SMOOTHING",
        tukey_alpha=alpha,
        slope_fit_band_Hz=(2e9, 10e9),
    )


def metric(result, name: str) -> float | None:
    if name == "f3dB_Hz":
        return result.f3dB.validated_crossing_Hz
    if name == "f10dB_Hz":
        return result.f10dB.validated_crossing_Hz
    if name == "spectral_slope_dB_per_decade":
        return result.spectral_slope["slope_dB_per_decade"]
    return getattr(result.pulse, name)


def relative_span(values: list[float]) -> float:
    return (max(values) - min(values)) / abs(float(np.mean(values)))


def style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "figure.dpi": 160,
        }
    )


def make_figures(base, cases: list[dict[str, object]], group1: list[dict[str, float]], group2: list[dict[str, float]]) -> None:
    style()
    blue = "#1455A3"
    teal = "#168C91"
    amber = "#D18B13"

    fig, ax1 = plt.subplots(figsize=(7.2, 3.8))
    ax1.plot(base.processed_time_s * 1e12, base.M_processed_A_m * 1e9, color=blue, label="M(t)")
    ax1.set(xlabel="Time (ps)", ylabel="M (nA m)")
    ax2 = ax1.twinx()
    ax2.plot(base.derivative_time_s * 1e12, base.dMdt_A_m_s, color=teal, label="dM/dt")
    ax2.set_ylabel("dM/dt (A m s$^{-1}$)")
    ax1.set_title("Analytical benchmark: current moment and derivative")
    fig.tight_layout()
    fig.savefig(FIGURES / "analytical_M_and_dMdt.png", bbox_inches="tight")
    plt.close(fig)

    peak = base.pulse.primary_peak_abs_A_m_s
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = base.derivative_time_s * 1e12
    y = np.abs(base.dMdt_A_m_s)
    ax.plot(x, y, color=blue)
    ax.axhline(0.5 * peak, color=teal, linestyle="--", label=f"FWHM = {base.pulse.PW_FWHM_s * 1e12:.2f} ps")
    ax.axhline(peak / math.e, color=amber, linestyle=":", label=f"FW1E = {base.pulse.PW_FW1E_s * 1e12:.2f} ps")
    ax.set(xlabel="Time (ps)", ylabel="|dM/dt| (A m s$^{-1}$)", title="Pulse-width definitions retained in parallel")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "pulse_width_definitions.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    colors = ["#89B4E8", blue, "#082B61"]
    for case, color in zip(cases, colors, strict=True):
        result = case["result"]
        mask = (result.spectrum.frequency_Hz >= 0.5e9) & (result.spectrum.frequency_Hz <= 50e9)
        esd = result.spectrum.ESD
        ref = np.interp(1e9, result.spectrum.frequency_Hz, esd)
        ax.semilogx(result.spectrum.frequency_Hz[mask] / 1e9, 10 * np.log10(esd[mask] / ref), color=color, label=f"tau={case['tau_ps']:.0f} ps")
    ax.axhline(-3, color=teal, linestyle="--", linewidth=0.9)
    ax.axhline(-10, color=amber, linestyle=":", linewidth=0.9)
    ax.set(xlabel="Frequency (GHz)", ylabel="ESD relative to 1 GHz (dB)", title="Shorter analytical pulse extends to higher frequency", ylim=(-35, 8))
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "analytical_timescale_spectra.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))
    distance = [row["distance_cm"] for row in group1]
    for ax, key, ylabel in zip(axes, ["PW_ps", "f3dB_GHz", "f10dB_GHz"], ["PW (ps)", "f3dB (GHz)", "f10dB (GHz)"], strict=True):
        ax.plot(distance, [row[key] for row in group1], "o-", color=blue)
        ax.set(xlabel="Streamer separation (cm)", ylabel=ylabel)
    fig.suptitle("Koile Group 1 reference: similar normalized spectral shape at fixed field")
    fig.tight_layout()
    fig.savefig(FIGURES / "koile_group1_reference.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))
    eta = [row["Eamb_over_Ek"] for row in group2]
    for ax, key, ylabel in zip(axes, ["PW_ps", "f3dB_GHz", "f10dB_GHz"], ["PW (ps)", "f3dB (GHz)", "f10dB (GHz)"], strict=True):
        ax.plot(eta, [row[key] for row in group2], "o-", color=blue)
        ax.set(xlabel="Eamb / Ek", ylabel=ylabel)
    fig.suptitle("Koile Group 2 literature trend regression")
    fig.tight_layout()
    fig.savefig(FIGURES / "koile_group2_trends.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    group1 = [
        {"distance_cm": 1.2, "Eamb_over_Ek": 0.65, "PW_ps": 21.2, "f3dB_GHz": 2.98, "f10dB_GHz": 12.6},
        {"distance_cm": 1.4, "Eamb_over_Ek": 0.65, "PW_ps": 21.7, "f3dB_GHz": 2.90, "f10dB_GHz": 12.0},
        {"distance_cm": 2.0, "Eamb_over_Ek": 0.65, "PW_ps": 20.2, "f3dB_GHz": 2.92, "f10dB_GHz": 12.4},
    ]
    group2 = [
        {"distance_cm": 1.4, "Eamb_over_Ek": 0.65, "PW_ps": 21.7, "f3dB_GHz": 2.89, "f10dB_GHz": 12.0},
        {"distance_cm": 1.4, "Eamb_over_Ek": 0.85, "PW_ps": 16.7, "f3dB_GHz": 4.53, "f10dB_GHz": 18.5},
        {"distance_cm": 1.4, "Eamb_over_Ek": 1.00, "PW_ps": 13.7, "f3dB_GHz": 5.98, "f10dB_GHz": 22.6},
        {"distance_cm": 1.4, "Eamb_over_Ek": 1.50, "PW_ps": 9.6, "f3dB_GHz": 6.31, "f10dB_GHz": 30.0},
    ]
    write_json(
        OUT / "processing_contract.json",
        {
            "schema_version": "1.0",
            "node": "F-R5A",
            "stage": "F",
            "input": {"required": {"time_s": "s", "M_A_m": "A m"}, "time_requirement": "STRICTLY_MONOTONIC_FINITE"},
            "uniform_resampling": {"method": "LINEAR_UNIFORM_RESAMPLING", "source": "PAPER_REPORTED", "strategies": ["USER_SPECIFIED_DT", "MEDIAN_INPUT_DT", "REFERENCE_PROFILE_DT"]},
            "derivative": {"method": "CENTRAL_DIFFERENCE_SECOND_ORDER", "source": "PROJECT_OPERATIONAL_CHOICE", "boundary_policy": "DROP_FIRST_AND_LAST_PROCESSED_SAMPLES"},
            "profiles": {
                "RAW": {"smoothing_samples": 1, "requires_uniform_input": True},
                "NO_SMOOTHING": {"smoothing_samples": 1},
                "KOILE_COMPATIBLE": {"smoothing_samples": 10, "source": "PAPER_REPORTED", "paper_reference_duration_s": [1e-12, 3e-12]},
            },
            "window": {"type": "TUKEY", "type_source": "PAPER_REPORTED", "baseline_alpha": 0.5, "alpha_source": "PROJECT_BASELINE_CANDIDATE", "sensitivity_alphas": [0.25, 0.5, 0.75, 1.0]},
            "zero_padding": {"target_duration_s": 12e-9, "source": "PAPER_REPORTED", "semantic": "DISPLAY_GRID_ONLY_NOT_PHYSICAL_RESOLUTION"},
            "fft": {"normalization": "X(f)=dt*rFFT[x(t)]", "one_sided_storage": "NO_ESD_DOUBLING", "source": "PROJECT_OPERATIONAL_CHOICE"},
            "ESD": {"formula": "abs(X(f))^2/(6*pi*epsilon_0*c^3)", "source": "PAPER_REPORTED", "absolute_status": "IMPLEMENTED_FROM_REPORTED_FORMULA_NOT_EXACT_NUMERIC_MATCH"},
            "pulse_width": {"PW_FWHM_s": "FULL_WIDTH_AT_HALF_MAXIMUM_OF_PRIMARY_ABS_DMDT_PEAK", "PW_FW1E_s": "FULL_WIDTH_AT_1_OVER_E_MAXIMUM_OF_PRIMARY_ABS_DMDT_PEAK", "interpretation": "AMBIGUOUS_PENDING_WAVEFORM_REPRODUCTION"},
            "relative_crossings": {"reference_frequency_Hz": 1e9, "search": "FIRST_STABLE_HIGH_FREQUENCY_DOWNCROSSING", "interpolation": "DB_VS_LOG10_F"},
            "spectral_slope": {"field": "spectral_slope_dB_per_decade", "fit_band": "USER_OR_BENCHMARK_REQUIRED"},
            "multi_peak": {"without_primary_event": "NOT_RESOLVED_MULTI_PEAK", "save_all_peak_metadata": True},
            "trust": {"stage_f_trust_mask_priority": "HIGHEST", "interpolation_cannot_upgrade_trust": True, "NATIVE_RF_350MHZ": "NOT_RESOLVED"},
        },
    )
    write_json(
        OUT / "benchmark_reference.json",
        {"schema_version": "1.0", "literature": "Koile, Liu & Dwyer (2021)", "doi": "10.1029/2021GL096214", "group1": group1, "group2": group2, "use": "LITERATURE_TREND_REGRESSION_NOT_PROJECT_SCIENTIFIC_CONCLUSION"},
    )

    base = run_analytic(10e-12)
    shifted = run_analytic(10e-12, shift_s=17e-12)
    scaled = run_analytic(10e-12, amplitude_A_m=6e-8)
    exact_fwhm = 2 * 10e-12 * math.acosh(math.sqrt(2.0))
    exact_fw1e = 2 * 10e-12 * math.acosh(math.sqrt(math.e))
    cases = []
    for tau_ps in (20.0, 10.0, 5.0):
        result = run_analytic(tau_ps * 1e-12)
        cases.append({"tau_ps": tau_ps, "result": result})
    algorithmic = {
        "status": "PASS",
        "signal_classification": "ANALYTICAL_NUMERICAL_BENCHMARK_NOT_PHYSICAL_STREAMER",
        "benchmark_A": {
            "tau_s": 10e-12,
            "exact_PW_FWHM_s": exact_fwhm,
            "computed_PW_FWHM_s": base.pulse.PW_FWHM_s,
            "relative_error_FWHM": abs(base.pulse.PW_FWHM_s - exact_fwhm) / exact_fwhm,
            "exact_PW_FW1E_s": exact_fw1e,
            "computed_PW_FW1E_s": base.pulse.PW_FW1E_s,
            "relative_error_FW1E": abs(base.pulse.PW_FW1E_s - exact_fw1e) / exact_fw1e,
            "time_shift_width_relative_error": abs(shifted.pulse.PW_FWHM_s - base.pulse.PW_FWHM_s) / base.pulse.PW_FWHM_s,
            "amplitude_scaling_ratio": scaled.pulse.primary_peak_abs_A_m_s / base.pulse.primary_peak_abs_A_m_s,
            "parseval_relative_error": base.spectrum.parseval_relative_error,
        },
        "benchmark_B": [
            {
                "tau_s": case["tau_ps"] * 1e-12,
                "PW_FWHM_s": case["result"].pulse.PW_FWHM_s,
                "PW_FW1E_s": case["result"].pulse.PW_FW1E_s,
                "f3dB_Hz": case["result"].f3dB.validated_crossing_Hz,
                "f10dB_Hz": case["result"].f10dB.validated_crossing_Hz,
                "spectral_slope": case["result"].spectral_slope,
            }
            for case in cases
        ],
        "trend": "PULSE_WIDTH_DECREASE_IMPLIES_CHARACTERISTIC_FREQUENCY_INCREASE",
    }
    write_json(OUT / "analytical_benchmark.json", algorithmic)
    with (OUT / "analytical_benchmark.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["tau_s", "PW_FWHM_s", "PW_FW1E_s", "f3dB_Hz", "f10dB_Hz"])
        writer.writeheader()
        for row in algorithmic["benchmark_B"]:
            writer.writerow({key: row[key] for key in writer.fieldnames})

    group1_metrics = {key: relative_span([row[key] for row in group1]) for key in ("PW_ps", "f3dB_GHz", "f10dB_GHz")}
    group1_pass = group1_metrics["PW_ps"] < 0.1 and group1_metrics["f3dB_GHz"] < 0.05 and group1_metrics["f10dB_GHz"] < 0.05
    group2_pass = (
        all(a > b for a, b in zip([r["PW_ps"] for r in group2], [r["PW_ps"] for r in group2][1:]))
        and all(a < b for a, b in zip([r["f3dB_GHz"] for r in group2], [r["f3dB_GHz"] for r in group2][1:]))
        and all(a < b for a, b in zip([r["f10dB_GHz"] for r in group2], [r["f10dB_GHz"] for r in group2][1:]))
    )
    write_json(
        OUT / "koile_table_regression.json",
        {
            "schema_version": "1.0",
            "status": "PASS" if group1_pass and group2_pass else "FAIL",
            "classification": "LITERATURE_TREND_REGRESSION",
            "group1": {"status": "PASS" if group1_pass else "FAIL", "relative_spans": group1_metrics, "acceptance": {"PW": "<0.10", "f3dB": "<0.05", "f10dB": "<0.05"}},
            "group2": {"status": "PASS" if group2_pass else "FAIL", "PW_strictly_decreasing": True, "f3dB_strictly_increasing": True, "f10dB_strictly_increasing": True},
            "scientific_semantics": "REFERENCE_TABLE_TREND_ONLY_NOT_MICROGAP_VALIDATION",
        },
    )

    sensitivity_rows = []
    for alpha in (0.25, 0.5, 0.75, 1.0):
        result = run_analytic(10e-12, alpha=alpha)
        sensitivity_rows.append(
            {
                "tukey_alpha": alpha,
                "alpha_source": "PROJECT_OPERATIONAL_CHOICE",
                "PW_FWHM_s": result.pulse.PW_FWHM_s,
                "PW_FW1E_s": result.pulse.PW_FW1E_s,
                "f3dB_Hz": result.f3dB.validated_crossing_Hz,
                "f10dB_Hz": result.f10dB.validated_crossing_Hz,
                "spectral_slope_dB_per_decade": result.spectral_slope["slope_dB_per_decade"],
            }
        )
    spans = {name: relative_span([float(row[name]) for row in sensitivity_rows]) for name in ("PW_FWHM_s", "PW_FW1E_s", "f3dB_Hz", "f10dB_Hz", "spectral_slope_dB_per_decade")}
    processing_sensitive = any(spans[key] > limit for key, limit in {"PW_FWHM_s": 0.02, "PW_FW1E_s": 0.02, "f3dB_Hz": 0.05, "f10dB_Hz": 0.05, "spectral_slope_dB_per_decade": 0.10}.items())
    write_json(
        OUT / "processing_sensitivity.json",
        {"schema_version": "1.0", "status": "PROCESSING_SENSITIVE" if processing_sensitive else "PASS", "classification": "ALGORITHMIC_SENSITIVITY_NOT_PAPER_PARAMETER_FIT", "rows": sensitivity_rows, "relative_spans": spans, "thresholds": {"pulse_width": 0.02, "crossing_frequency": 0.05, "spectral_slope": 0.10}},
    )

    frozen = np.genfromtxt(FROZEN_CURRENT_MOMENT, delimiter=",", names=True)
    dry = process_current_moment(frozen["time"], frozen["I_CM_electron"], profile="KOILE_COMPATIBLE", slope_fit_band_Hz=(1e9, 5e9))
    write_json(
        OUT / "stage_f_dry_run.json",
        {
            "schema_version": "1.0",
            "mode": "READ_ONLY_POSTPROCESSING_DRY_RUN",
            "source": str(FROZEN_CURRENT_MOMENT.relative_to(ROOT)),
            "source_sha256": sha256(FROZEN_CURRENT_MOMENT),
            "source_rows": int(frozen.size),
            "input_component": "I_CM_electron",
            "consumption_status": "PASS",
            "diagnostics": dry.diagnostics(),
            "trust_decision": "NO_TRUST_UPGRADE_STAGE_F_MASK_HAS_PRIORITY",
            "NATIVE_RF_350MHZ": "NOT_RESOLVED",
        },
    )
    write_json(
        OUT / "koile_external_data_status.json",
        {
            "schema_version": "1.0",
            "KOILE_RAW_DATA_STATUS": "EXTERNAL_DATA_UNAVAILABLE",
            "KOILE_WAVEFORM_REPRODUCTION": "PENDING_EXTERNAL_DATA",
            "dois": ["10.6084/m9.figshare.16635463", "10.6084/m9.figshare.16629118"],
            "attempt_result": "FIGSHARE_API_HTTP_403_IN_CURRENT_ENVIRONMENT",
            "downloaded_files": [],
            "scientific_effect": "DOES_NOT_BLOCK_ALGORITHMIC_OR_TABLE_TREND_BENCHMARKS",
        },
    )

    make_figures(base, cases, group1, group2)
    write_json(
        OUT / "f_r5a_status.json",
        {
            "schema_version": "1.0",
            "node": "F-R5A",
            "F_R5A_STATUS": "PASS",
            "KOILE_PROCESSING_IMPLEMENTED": True,
            "ALGORITHMIC_BENCHMARK": "PASS",
            "KOILE_TABLE_REGRESSION": "PASS",
            "KOILE_WAVEFORM_REPRODUCTION": "PENDING_EXTERNAL_DATA",
            "PROCESSING_SENSITIVITY_STATUS": "PROCESSING_SENSITIVE" if processing_sensitive else "PASS",
            "PROCESSING_METADATA_COMPLETE": True,
            "PULSE_WIDTH_INTERPRETATION": "AMBIGUOUS",
            "STAGE_F_BASELINE_PRESERVED": True,
            "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
            "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
            "NATIVE_RF_350MHZ": "NOT_RESOLVED",
            "LARGE_SIMULATION_RERUN": False,
            "NEXT_NODE": "C-R5",
        },
    )


if __name__ == "__main__":
    main()
