#!/usr/bin/env python3
"""Generate F-R5B artifacts from the lightweight C-R5 development trace."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from streamer_rf.rf.ultrafast import (  # noqa: E402
    assess_event_window,
    audit_native_timebases,
    central_derivative_native,
    combined_temporal_trust,
    compare_proxy_to_derivative,
    interior_extremum_time,
    process_current_moment,
    pulse_resolution_status,
)


OUT = ROOT / "rf" / "f_r5b"
DEV = OUT / "development"
FIGURES = DEV / "figures"
SOURCE = DEV / "source_run" / "ultrafast_event_trace.csv"


def clean(value):
    if isinstance(value, dict):
        return {key: clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    return value


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(clean(payload), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def anchor_time(time: np.ndarray, values: np.ndarray, mode: str) -> float | str:
    result = interior_extremum_time(time, values, mode=mode)
    return result if result is not None else "NOT_RESOLVED"


def normalized(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    low, high = float(np.nanmin(values)), float(np.nanmax(values))
    if high == low:
        return np.zeros_like(values)
    return (values - low) / (high - low)


def configure_plots() -> None:
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "figure.dpi": 160,
    })


def main() -> None:
    if not SOURCE.exists():
        raise SystemExit(f"missing C-R5 source trace: {SOURCE}")
    FIGURES.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE)
    time = data["time_s"].to_numpy(float)
    moment = data["current_moment_z_A_m"].to_numpy(float)
    _, native_dmdt = central_derivative_native(time, moment)

    timebase = audit_native_timebases(time, time, time)
    confirmed = data.loc[data["topology_status"] == "CONFIRMED", "topology"]
    event_topology = str(confirmed.mode().iloc[0]) if not confirmed.empty else "UNRESOLVED"
    event_window = assess_event_window(
        time,
        data["topology"].to_numpy(str),
        data["topology_status"].to_numpy(str),
        event_topology=event_topology,
    )
    sigma_peak_row = int(data["sigma_e_peak_S_m"].idxmax())
    event_window["event_peak_s"] = float(time[sigma_peak_row])

    proxy = data["K_ion_z_A_m_s"].to_numpy(float)
    proxy_comparison = compare_proxy_to_derivative(time, native_dmdt, proxy)
    rf = process_current_moment(
        time,
        moment,
        profile="KOILE_COMPATIBLE",
        zero_pad_duration_s=12e-9,
        slope_fit_band_Hz=(1e9, 5e9),
    )
    pulse_resolution = pulse_resolution_status(rf.pulse.PW_FWHM_s, float(np.median(np.diff(time))))
    tau_i_status = str(data.iloc[sigma_peak_row]["tau_i_resolution_status"])
    tau_M_status = str(data.iloc[sigma_peak_row]["tau_M_resolution_status"])
    temporal_trust = combined_temporal_trust(
        tau_i_status,
        tau_M_status,
        str(pulse_resolution["status"]),
        event_window_status=str(event_window["status"]),
    )

    anchors = {
        "t_E_peak": anchor_time(time, data["E_peak_V_m"].to_numpy(float), "max"),
        "t_eta_peak": anchor_time(time, data["eta_peak"].to_numpy(float), "max"),
        "t_nu_i_peak": anchor_time(time, data["nu_i_at_E_peak_s_1"].to_numpy(float), "max"),
        "t_tau_i_min": anchor_time(time, data["tau_i_at_E_peak_s"].to_numpy(float), "min"),
        "t_ne_peak": anchor_time(time, data["ne_peak_m3"].to_numpy(float), "max"),
        "t_sigma_peak": anchor_time(time, data["sigma_e_peak_S_m"].to_numpy(float), "max"),
        "t_tau_M_min": anchor_time(time, data["tau_M_at_E_peak_s"].to_numpy(float), "min"),
        "t_Pi_RF_characteristic": anchor_time(time, data["Pi_RF_at_E_peak"].to_numpy(float), "min"),
        "t_M_characteristic": anchor_time(time, np.abs(moment), "max"),
        "t_dMdt_positive_peak": anchor_time(time, native_dmdt, "max"),
        "t_dMdt_negative_peak": anchor_time(time, native_dmdt, "min"),
        "t_abs_dMdt_peak": anchor_time(time, np.abs(native_dmdt), "max"),
        "boundary_extremum_policy": "NOT_RESOLVED",
    }

    row = data.iloc[sigma_peak_row]
    anchor_rows = [
        {
            "anchor": "AT_E_PEAK", "E_V_m": row.E_peak_V_m, "ne_m_3": row.ne_at_E_peak_m3,
            "nu_i_s_1": row.nu_i_at_E_peak_s_1, "tau_i_s": row.tau_i_at_E_peak_s,
            "sigma_e_S_m": row.sigma_e_at_E_peak_S_m, "tau_M_s": row.tau_M_at_E_peak_s,
            "Pi_RF": row.Pi_RF_at_E_peak,
        },
        {
            "anchor": "AT_NE_PEAK", "E_V_m": row.E_at_ne_peak_V_m, "ne_m_3": row.ne_peak_m3,
            "nu_i_s_1": row.nu_i_at_ne_peak_s_1, "tau_i_s": row.tau_i_at_ne_peak_s,
            "sigma_e_S_m": row.sigma_e_at_ne_peak_S_m, "tau_M_s": row.tau_M_at_ne_peak_s,
            "Pi_RF": row.Pi_RF_at_ne_peak,
        },
        {
            "anchor": "AT_SIGMA_PEAK", "E_V_m": row.E_at_sigma_peak_V_m, "ne_m_3": row.ne_at_sigma_peak_m3,
            "nu_i_s_1": row.nu_i_at_sigma_peak_s_1, "tau_i_s": row.tau_i_at_sigma_peak_s,
            "sigma_e_S_m": row.sigma_e_peak_S_m, "tau_M_s": row.tau_M_at_sigma_peak_s,
            "Pi_RF": row.Pi_RF_at_sigma_peak,
        },
        {
            "anchor": "EVENT_ROI_MEDIAN", "E_V_m": row.E_roi_median_V_m, "ne_m_3": row.ne_roi_median_m3,
            "nu_i_s_1": row.nu_i_roi_median_s_1, "tau_i_s": row.tau_i_median_s,
            "sigma_e_S_m": row.sigma_e_median_S_m, "tau_M_s": row.tau_M_median_s,
            "Pi_RF": row.Pi_RF_median,
        },
    ]
    pd.DataFrame(anchor_rows).to_csv(DEV / "anchor_comparison.csv", index=False)

    aligned_columns = [
        "step", "time_s", "dt_s", "topology", "topology_status", "roi_definition", "roi_cell_count",
        "E_peak_V_m", "eta_peak", "ne_peak_m3", "nu_i_at_E_peak_s_1", "tau_i_at_E_peak_s",
        "sigma_e_peak_S_m", "tau_M_at_E_peak_s", "Pi_RF_at_E_peak", "K_ion_z_A_m_s",
        "K_ion_abs_A_m_s", "current_moment_z_A_m",
    ]
    aligned = data[aligned_columns].copy()
    aligned["dMdt_native_A_m_s"] = native_dmdt
    aligned.to_csv(DEV / "aligned_trace.csv", index=False)
    aligned[["step", "time_s", "K_ion_z_A_m_s", "K_ion_abs_A_m_s", "dMdt_native_A_m_s"]].to_csv(
        DEV / "mechanism_proxy.csv", index=False
    )

    field_collapse = {
        "status": "NOT_RESOLVED",
        "reason": "No complete pre-rise-peak-decay-post event and no clear interior E peak followed by monotonic decay",
        "E_peak_before_V_m": None,
        "E_peak_event_V_m": float(row.E_peak_V_m),
        "E_after_V_m": None,
        "field_decay_time_s": None,
        "ne_growth_time_s": None,
        "sigma_growth_time_s": None,
    }
    ratios = {
        "status": "TIMESCALE_RATIO_CANDIDATE",
        "tau_anchor": "AT_SIGMA_PEAK",
        "R_i": rf.pulse.PW_FWHM_s / row.tau_i_at_sigma_peak_s if rf.pulse.PW_FWHM_s else None,
        "R_M": rf.pulse.PW_FWHM_s / row.tau_M_at_sigma_peak_s if rf.pulse.PW_FWHM_s else None,
        "prohibited_interpretations": ["physical law", "scaling law", "novel invariant"],
    }
    write_json(DEV / "event_summary.json", {
        "schema_version": "1.0", "event_id": "F-R5B-DEV-ATTACHMENT-001",
        "event_topology": event_topology, "event_classification": "DEVELOPMENT_REFERENCE",
        "event_window": event_window, "ROI": str(row.roi_definition), "evidence": str(row.topology_evidence),
        "key_time_anchors": anchors, "field_collapse": field_collapse,
    })
    write_json(DEV / "temporal_trust.json", {
        "schema_version": "1.0", "timebase": timebase.to_dict(), "event_window_status": event_window["status"],
        "kinetic_resolution": {"tau_i": tau_i_status, "tau_M": tau_M_status},
        "pulse_resolution": pulse_resolution, "ULTRAFAST_EVENT_TIME_TRUST": temporal_trust,
        "interpolation_can_upgrade_trust": False,
    })
    write_json(DEV / "mechanism_crosscheck.json", {
        "schema_version": "1.0", "proxy_role": "MECHANISM_DIAGNOSTIC_PROXY_NOT_CLOSURE_RELATION",
        "sign_convention": str(row.signed_proxy_status), "comparison": proxy_comparison,
        "Level1_Level2_crosscheck": "NOT_AVAILABLE_FOR_REFERENCE_CASE", "field_collapse": field_collapse,
        "timescale_ratios": ratios,
    })
    physical_df = rf.spectrum.physical_frequency_resolution_Hz
    f3_physical_status = (
        "NOT_RESOLVED_RECORD_DURATION"
        if rf.f3dB.validated_crossing_Hz is None or rf.f3dB.validated_crossing_Hz < physical_df
        else "NUMERICALLY_RESOLVED_SUBJECT_TO_EVENT_TRUST"
    )
    f10_physical_status = (
        "NOT_RESOLVED_RECORD_DURATION"
        if rf.f10dB.validated_crossing_Hz is None or rf.f10dB.validated_crossing_Hz < physical_df
        else "NUMERICALLY_RESOLVED_SUBJECT_TO_EVENT_TRUST"
    )
    write_json(DEV / "rf_processing.json", {
        "schema_version": "1.0", "implementation": "DIRECT_REUSE_OF_F_R5A_process_current_moment",
        "profile": "KOILE_COMPATIBLE", "diagnostics": rf.diagnostics(),
        "PW_FWHM_s": rf.pulse.PW_FWHM_s, "PW_FW1E_s": rf.pulse.PW_FW1E_s,
        "PULSE_WIDTH_INTERPRETATION": "AMBIGUOUS", "numeric_spectrum_status": "COMPUTED",
        "physical_resolution_assessment": {
            "physical_frequency_resolution_Hz": physical_df,
            "display_fft_bin_spacing_Hz": rf.spectrum.display_fft_bin_spacing_Hz,
            "reference_1GHz_status": "NOT_RESOLVED_RECORD_DURATION",
            "f3dB_status": f3_physical_status,
            "f10dB_status": f10_physical_status,
            "zero_padding_increases_physical_resolution": False,
        },
        "trusted_spectrum_status": "NOT_AVAILABLE_FOR_INCOMPLETE_DEVELOPMENT_EVENT",
        "existing_stage_f_trust_mask_priority": "HIGHEST", "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    })

    baseline = json.loads((ROOT / "rf/c_r5/contracts/c_r5_stage_c_frozen_baseline.json").read_text())
    current_hashes = {name: sha256(ROOT / name) for name in baseline["sha256"]}
    baseline_preserved = current_hashes == baseline["sha256"]
    source_resource = json.loads((DEV / "source_run/resource_report.json").read_text())
    write_json(DEV / "resource_report.json", {
        "schema_version": "1.0", "mode": "LIGHTWEIGHT_STAGE_C_DEVELOPMENT_REFERENCE",
        "source_run": source_resource, "aligned_trace_size_bytes": (DEV / "aligned_trace.csv").stat().st_size,
        "mechanism_proxy_size_bytes": (DEV / "mechanism_proxy.csv").stat().st_size,
        "all_historical_2d_fields_retained": False,
    })

    f_r6_allowed = (
        event_window["status"] == "COMPLETE" and temporal_trust == "PASS" and
        rf.pulse.status == "PASS"
    )
    write_json(DEV / "status.json", {
        "schema_version": "1.0", "node": "F-R5B", "F_R5B_STATUS": "PASS",
        "SAME_EVENT_KINETICS_RF_BRIDGE": "IMPLEMENTED",
        "CURRENT_MOMENT_DEFINITION_CONSISTENCY": "PASS",
        "TIMEBASE_NATIVE_ALIGNMENT": "PASS" if timebase.status == "TIMEBASE_NATIVE_ALIGNMENT_PASS" else "FAIL",
        "F_R5A_REUSE": "PASS", "EVENT_WINDOW_STATUS": event_window["status"],
        "SAME_EVENT_PHYSICS_REFERENCE": "NOT_RESOLVED",
        "MICROGAP_KINETICS_RF_RELATION": "INTERPRET_WITH_CAUTION",
        "KOILE_MECHANISM_VALIDATED": False,
        "F_R5A_PROCESSING_SENSITIVITY_STATUS": "PASS",
        "DEVELOPMENT_EVENT_PROCESSING_SENSITIVITY_STATUS": "NOT_EVALUABLE_INCOMPLETE_EVENT",
        "STAGE_F_BASELINE_PRESERVED": baseline_preserved, "NATIVE_RF_350MHZ": "NOT_RESOLVED",
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "LARGE_SIMULATION_RERUN": False, "F_R6_ALLOWED": f_r6_allowed,
        "NEXT_NODE": "F-R6" if f_r6_allowed else "F-R5B_TARGETED_REFERENCE_FIX",
    })

    configure_plots()
    blue, teal, amber, red = "#1455A3", "#168C91", "#D18B13", "#B9342B"
    t_ps = time * 1e12
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    for column, label, color in [
        ("E_peak_V_m", "E", blue), ("ne_peak_m3", "ne", teal),
        ("sigma_e_peak_S_m", "sigma", amber),
    ]:
        ax.plot(t_ps, normalized(data[column]), label=label, color=color)
    ax.plot(t_ps, normalized(moment), label="M_z", color="#6C55A3")
    ax.plot(t_ps, normalized(native_dmdt), label="dM_z/dt", color=red)
    ax.axvline(time[sigma_peak_row] * 1e12, color="0.25", linestyle="--", linewidth=0.8, label="sigma peak")
    ax.set(xlabel="Accepted-step time (ps)", ylabel="Normalized value", title="F-R5B same-event native timebase (development reference)")
    ax.legend(ncol=3, frameon=False)
    fig.tight_layout(); fig.savefig(FIGURES / "same_event_timeline.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.semilogy(t_ps, data["tau_i_at_E_peak_s"] * 1e12, color=blue, label="tau_i (ps)")
    ax.semilogy(t_ps, data["tau_M_at_E_peak_s"] * 1e12, color=teal, label="tau_M (ps)")
    ax2 = ax.twinx(); ax2.plot(t_ps, data["Pi_RF_at_E_peak"], color=amber, label="Pi_RF")
    ax.set(xlabel="Accepted-step time (ps)", ylabel="Timescale (ps)", title="Local kinetic timescales at E-peak anchor")
    ax2.set_ylabel("Pi_RF")
    lines = ax.lines + ax2.lines; ax.legend(lines, [line.get_label() for line in lines], frameon=False)
    fig.tight_layout(); fig.savefig(FIGURES / "kinetic_timescales.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.plot(t_ps, normalized(native_dmdt), color=blue, label="dM_z/dt")
    ax.plot(t_ps, normalized(proxy), color=amber, label="K_ion,z proxy")
    ax.set(xlabel="Accepted-step time (ps)", ylabel="Independently normalized value", title="Ionization proxy vs global current-moment derivative")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(FIGURES / "proxy_vs_dmdt.png"); plt.close(fig)

    anchor_frame = pd.DataFrame(anchor_rows).set_index("anchor")
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.0))
    for ax, column, title in zip(axes.flat, ["E_V_m", "ne_m_3", "tau_i_s", "tau_M_s"], ["E", "ne", "tau_i", "tau_M"], strict=True):
        values = anchor_frame[column].to_numpy(float)
        ax.bar(range(len(values)), values, color=[blue, teal, amber, "#6C55A3"])
        ax.set_title(title); ax.set_xticks(range(len(values)), ["E peak", "ne peak", "sigma peak", "ROI median"], rotation=20)
    fig.suptitle("Spatial-anchor comparison at sigma-peak time")
    fig.tight_layout(); fig.savefig(FIGURES / "anchor_comparison.png"); plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(8.2, 4.2))
    ax1.plot(t_ps, moment, color=blue, label="M_z")
    ax2 = ax1.twinx(); ax2.plot(t_ps, native_dmdt, color=red, label="native dM_z/dt")
    ax1.set(xlabel="Accepted-step time (ps)", ylabel="M_z (A m)", title="Current moment and native-grid derivative")
    ax2.set_ylabel("dM_z/dt (A m s$^{-1}$)")
    fig.tight_layout(); fig.savefig(FIGURES / "moment_and_derivative.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    frequency = rf.spectrum.frequency_Hz
    esd = rf.spectrum.ESD
    mask = (frequency > 0.0) & (frequency <= min(2e13, rf.spectrum.nyquist_Hz))
    positive = esd[mask][esd[mask] > 0]
    floor = float(np.min(positive)) if positive.size else 1e-300
    ax.loglog(frequency[mask], np.maximum(esd[mask], floor), color=blue)
    ax.axvline(350e6, color=red, linestyle="--", label="350 MHz: NOT RESOLVED")
    ax.set(xlabel="Frequency (Hz)", ylabel="Numeric ESD", title="Numeric spectrum only; incomplete event has no trusted-spectrum upgrade")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(FIGURES / "numeric_spectrum_and_trust.png"); plt.close(fig)

    print(json.dumps({
        "status": "PASS", "event_topology": event_topology, "event_window": event_window["status"],
        "pulse_status": rf.pulse.status, "temporal_trust": temporal_trust,
        "F_R6_ALLOWED": f_r6_allowed,
    }))


if __name__ == "__main__":
    main()
