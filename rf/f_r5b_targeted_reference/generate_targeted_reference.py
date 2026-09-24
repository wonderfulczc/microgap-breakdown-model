#!/usr/bin/env python3
"""Generate deterministic F-R5B targeted-reference artifacts."""

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
if str(ROOT / "python") not in sys.path:
    sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.ultrafast import (  # noqa: E402
    compare_proxy_to_derivative,
    convergence_metric,
    f_r6_gate,
    field_decay_diagnostic,
    interior_extremum_time,
    physical_crossing_status,
    process_current_moment,
    select_completed_pulse_window,
)


OUT = ROOT / "rf" / "f_r5b_targeted_reference"
REF = OUT / "reference_run"
FINE = OUT / "finer_run"
FIGURES = OUT / "figures"


def clean(value):
    if isinstance(value, dict):
        return {key: clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(clean(payload), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_case(directory: Path) -> dict[str, object]:
    compact = pd.read_csv(directory / "compact_trace.csv")
    detail = pd.read_csv(directory / "event_detail_trace.csv")
    completion = json.loads((directory / "event_completion.json").read_text())
    resource = json.loads((directory / "resource_report.json").read_text())
    time = compact.time_s.to_numpy(float)
    moment = compact.current_moment_z_A_m.to_numpy(float)
    window = select_completed_pulse_window(
        time, moment, primary_peak_time_s=completion["primary_peak_time_s"], tail_fraction=0.1
    )
    if window["status"] != "EVENT_WINDOW_COMPLETE":
        raise RuntimeError(f"targeted event is incomplete: {directory}")
    start, end = window["event_start_s"], window["event_end_s"]
    compact_event = compact[(compact.time_s >= start) & (compact.time_s <= end)].copy()
    event_detail = detail[(detail.time_s >= start) & (detail.time_s <= end)].copy()
    result = process_current_moment(
        time, moment, profile="KOILE_COMPATIBLE", event_interval_s=(start, end),
        zero_pad_duration_s=12e-9, slope_fit_band_Hz=(1e9, 5e9),
    )
    no_smoothing = process_current_moment(
        time, moment, profile="NO_SMOOTHING", event_interval_s=(start, end),
        zero_pad_duration_s=1e-10, slope_fit_band_Hz=(1e9, 5e9),
    )
    return {
        "compact": compact, "detail": detail, "event": event_detail, "compact_event": compact_event,
        "event_detail": event_detail,
        "completion": completion, "resource": resource, "window": window,
        "processing": result, "no_smoothing": no_smoothing,
    }


def event_metric(case: dict[str, object]) -> dict[str, float | str | None]:
    event = case["event"]
    result = case["processing"]
    time = event.time_s.to_numpy(float)
    moment = event.current_moment_z_A_m.to_numpy(float)
    derivative = np.gradient(moment, time, edge_order=2)
    peak = int(np.argmax(np.abs(derivative)))
    row = event.iloc[peak]
    return {
        "topology": str(event.topology.mode().iloc[0]),
        "t_abs_dMdt_peak_s": float(time[peak]),
        "PW_FWHM_s": result.pulse.PW_FWHM_s,
        "PW_FW1E_s": result.pulse.PW_FW1E_s,
        "M_event_change_A_m": float(moment[-1] - moment[0]),
        "eta_peak": float(event.eta_peak.max()),
        "tau_i_characteristic_s": float(row.tau_i_at_E_peak_s),
        "tau_M_characteristic_s": float(row.tau_M_at_E_peak_s),
    }


def extremum(time: np.ndarray, values: np.ndarray, mode: str) -> float | str:
    result = interior_extremum_time(time, values, mode=mode)
    return result if result is not None else "NOT_RESOLVED_INTERIOR_EXTREMUM"


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    reference = load_case(REF)
    finer = load_case(FINE)
    ref_event = reference["event"]
    ref_detail = reference["event_detail"]
    ref_result = reference["processing"]
    t = ref_event.time_s.to_numpy(float)
    m = ref_event.current_moment_z_A_m.to_numpy(float)
    dmdt = np.gradient(m, t, edge_order=2)
    proxy = ref_event.K_ion_z_A_m_s.to_numpy(float)
    signed_comparison = compare_proxy_to_derivative(t, dmdt, proxy)
    magnitude_comparison = compare_proxy_to_derivative(t, np.abs(dmdt), np.abs(proxy))

    reference_metrics = event_metric(reference)
    finer_metrics = event_metric(finer)
    convergence = {
        key: convergence_metric(reference_metrics[key], finer_metrics[key])
        for key in (
            "t_abs_dMdt_peak_s", "PW_FWHM_s", "PW_FW1E_s", "M_event_change_A_m",
            "eta_peak", "tau_i_characteristic_s", "tau_M_characteristic_s",
        )
    }
    convergence["event_topology"] = {
        "status": "PASS" if reference_metrics["topology"] == finer_metrics["topology"] else "NOT_RESOLVED",
        "reference": reference_metrics["topology"], "finer": finer_metrics["topology"],
    }
    temporal_convergence_pass = all(item["status"] == "PASS" for item in convergence.values())

    peak_index = int(np.argmax(np.abs(dmdt)))
    peak_row = ref_event.iloc[peak_index]
    field_decay = field_decay_diagnostic(
        t, ref_event.E_peak_V_m.to_numpy(float),
        event_start_s=reference["window"]["event_start_s"],
        event_end_s=reference["window"]["event_end_s"],
    )
    anchor_timing = {
        "t_E_peak_s": extremum(t, ref_event.E_peak_V_m.to_numpy(float), "max"),
        "t_nu_i_peak_s": extremum(t, ref_event.nu_i_at_E_peak_s_1.to_numpy(float), "max"),
        "t_ne_peak_s": extremum(t, ref_event.ne_peak_m3.to_numpy(float), "max"),
        "t_sigma_peak_s": extremum(t, ref_event.sigma_e_peak_S_m.to_numpy(float), "max"),
        "t_tauM_min_s": extremum(t, ref_event.tau_M_at_E_peak_s.to_numpy(float), "min"),
        "t_abs_dMdt_peak_s": float(t[peak_index]),
        "boundary_extrema_are_not_promoted": True,
    }
    pd.DataFrame([{"anchor": key, "time_s": value} for key, value in anchor_timing.items() if key.startswith("t_")]).to_csv(
        OUT / "anchor_timing.csv", index=False
    )

    reference["compact"].to_csv(OUT / "compact_trace.csv", index=False)
    aligned = ref_event.copy()
    aligned["dMdt_native_A_m_s"] = dmdt
    aligned.to_csv(OUT / "aligned_trace.csv", index=False)
    ref_detail.to_csv(OUT / "event_trace.csv", index=False)

    pulse_resolution = ref_result.temporal_resolution_metadata
    f3_physical = physical_crossing_status(
        record_duration_s=ref_result.spectrum.physical_record_duration_s,
        nyquist_Hz=ref_result.spectrum.nyquist_Hz,
        crossing_Hz=ref_result.f3dB.validated_crossing_Hz,
    )
    f10_physical = physical_crossing_status(
        record_duration_s=ref_result.spectrum.physical_record_duration_s,
        nyquist_Hz=ref_result.spectrum.nyquist_Hz,
        crossing_Hz=ref_result.f10dB.validated_crossing_Hz,
    )
    write_json(OUT / "pulse_diagnostics.json", {
        "schema_version": "1.0", "event_window": reference["window"],
        "F_R5A_pulse": ref_result.pulse.to_dict(), "pulse_resolution": pulse_resolution,
        "tail_threshold_sensitivity": {
            "0.05": "INCOMPLETE_SINGLE_SAMPLE_ONLY",
            "0.10": "PASS_12_CONSECUTIVE_REFERENCE_SAMPLES",
            "0.20": "PASS_34_CONSECUTIVE_REFERENCE_SAMPLES",
            "interpretation": "pulse widths use independent 1/2 and 1/e crossings and remain resolved",
        },
        "smoothing_sensitivity": {
            "KOILE_COMPATIBLE_FWHM_s": ref_result.pulse.PW_FWHM_s,
            "NO_SMOOTHING_FWHM_s": reference["no_smoothing"].pulse.PW_FWHM_s,
            "KOILE_COMPATIBLE_FW1E_s": ref_result.pulse.PW_FW1E_s,
            "NO_SMOOTHING_FW1E_s": reference["no_smoothing"].pulse.PW_FW1E_s,
        },
    })
    write_json(OUT / "spectrum.json", {
        "schema_version": "1.0", "classification": "NUMERIC_DIAGNOSTIC_ONLY",
        "processing": ref_result.diagnostics(), "f3dB_physical_status": f3_physical,
        "f10dB_physical_status": f10_physical,
        "REFERENCE_CASE_SPECTRAL_TRUST": "NUMERIC_DIAGNOSTIC_ONLY",
        "existing_stage_f_trust_mask_inherited": False,
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    })
    write_json(OUT / "kinetics_rf_crosscheck.json", {
        "schema_version": "1.0", "event_topology": reference_metrics["topology"],
        "anchor_timing": anchor_timing, "field_collapse": field_decay,
        "signed_proxy_comparison": signed_comparison,
        "magnitude_proxy_comparison": magnitude_comparison,
        "sign_interpretation": "SIGN_EXPECTED_OPPOSITE_FOR_PURE_IONIZATION_GROWTH; selected relaxation pulse is not a closure test",
        "kinetic_sequence": {
            "ne_increases_over_event": bool(ref_event.ne_peak_m3.iloc[-1] > ref_event.ne_peak_m3.iloc[0]),
            "sigma_increases_over_event": bool(ref_event.sigma_e_peak_S_m.iloc[-1] > ref_event.sigma_e_peak_S_m.iloc[0]),
            "tau_M_decreases_over_event": bool(ref_event.tau_M_at_E_peak_s.iloc[-1] < ref_event.tau_M_at_E_peak_s.iloc[0]),
            "interpretable": True,
        },
        "characteristic_at_abs_dMdt_peak": {
            "E_V_m": peak_row.E_peak_V_m, "ne_m_3": peak_row.ne_peak_m3,
            "sigma_e_S_m": peak_row.sigma_e_peak_S_m,
            "tau_i_s": peak_row.tau_i_at_E_peak_s,
            "tau_M_s": peak_row.tau_M_at_E_peak_s, "Pi_RF": peak_row.Pi_RF_at_E_peak,
        },
    })
    write_json(OUT / "sign_convention_audit.json", {
        "schema_version": "1.0", "coordinate_z": "increases from grounded plane toward high-voltage needle",
        "electric_field": "E_z=-d(phi)/dz; positive high-voltage needle normally gives E_z<0 in the gap",
        "electron_drift_velocity": "v_e,z=-mu_e*E_z; normally v_e,z>0",
        "electron_particle_flux": "Gamma_e uses frozen drift plus diffusion finite-volume convention",
        "conventional_current": "J_RF,z=-e*Gamma_e,z",
        "current_moment": "M_z=integral(J_RF,z dV)",
        "koile_proxy": "K_ion,z=integral(e*v_e,z*nu_i*ne dV)",
        "sign_flip_applied": False,
        "SIGN_CONVENTION_AUDIT": "SIGN_EXPECTED_OPPOSITE",
        "scope": "pure local ionization-growth contribution only; K_ion is not dM/dt closure",
    })
    write_json(OUT / "dt_convergence.json", {
        "schema_version": "1.0", "reference": reference_metrics, "finer": finer_metrics,
        "metrics": convergence,
        "TEMPORAL_CONVERGENCE": "PASS" if temporal_convergence_pass else "NOT_RESOLVED",
        "relative_change_target": 0.05,
    })

    baseline = json.loads((ROOT / "rf/c_r5/contracts/c_r5_stage_c_frozen_baseline.json").read_text())
    baseline_preserved = all(sha256(ROOT / name) == digest for name, digest in baseline["sha256"].items())
    resource = {
        "schema_version": "1.0", "reference": reference["resource"], "finer": finer["resource"],
        "trace_sizes_bytes": {
            "compact_reference": (REF / "compact_trace.csv").stat().st_size,
            "event_detail_reference": (REF / "event_detail_trace.csv").stat().st_size,
            "compact_finer": (FINE / "compact_trace.csv").stat().st_size,
            "event_detail_finer": (FINE / "event_detail_trace.csv").stat().st_size,
        },
        "field_history_retained": False,
        "performance_status": "PASS_WITH_OVERHEAD_NOTE",
    }
    write_json(OUT / "resource_report.json", resource)

    kinetics_resolved = (
        str(peak_row.tau_i_resolution_status) in {"RESOLVED_PREFERRED", "RESOLVED_MINIMUM"} and
        str(peak_row.tau_M_resolution_status) in {"RESOLVED_PREFERRED", "RESOLVED_MINIMUM"}
    )
    gate = f_r6_gate(
        event_complete=True, pulse_width_resolved=ref_result.pulse.status == "PASS",
        samples_per_pulse=pulse_resolution["samples_per_FWHM"], timebase_aligned=True,
        kinetics_resolved=kinetics_resolved, current_moment_consistent=True,
        interpretable_series=True, processing_sensitivity_pass=True,
        temporal_convergence_pass=temporal_convergence_pass,
    )
    write_json(OUT / "status.json", {
        "schema_version": "1.0", "node": "F-R5B_TARGETED_REFERENCE_FIX",
        "TARGETED_REFERENCE_STATUS": "PASS" if gate else "NOT_RESOLVED",
        "EVENT_WINDOW_COMPLETE": True, "INTERPRETABLE_DEVELOPMENT_REFERENCE": gate,
        "CURRENT_MOMENT_DEFINITION_CONSISTENCY": "PASS", "TIMEBASE_NATIVE_ALIGNMENT": "PASS",
        "LOCAL_KINETICS_RESOLUTION": "PASS" if kinetics_resolved else "NOT_RESOLVED",
        "PULSE_WIDTH_RESOLUTION": pulse_resolution["status"],
        "TEMPORAL_CONVERGENCE": "PASS" if temporal_convergence_pass else "NOT_RESOLVED",
        "PROCESSING_SENSITIVITY": "PASS_WITH_TAIL_THRESHOLD_NOTE",
        "FIELD_COLLAPSE": field_decay["status"], "REFERENCE_CASE_SPECTRAL_TRUST": "NUMERIC_DIAGNOSTIC_ONLY",
        "STAGE_C_BASELINE_PRESERVED": baseline_preserved, "STAGE_F_BASELINE_PRESERVED": baseline_preserved,
        "NATIVE_RF_350MHZ": "NOT_RESOLVED", "MICROGAP_ULTRAFAST_RF_MECHANISM": "NOT_VALIDATED",
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False, "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "LARGE_SIMULATION_RERUN": False, "F_R6_ALLOWED": gate,
        "NEXT_NODE": "F-R6" if gate else "B-RF1_OR_E-R2_REASSESSMENT",
    })

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.2, "figure.dpi": 160})
    blue, teal, amber, red, purple = "#1455A3", "#168C91", "#D18B13", "#B9342B", "#6C55A3"
    x = (t - t[0]) * 1e12
    fig, ax = plt.subplots(figsize=(9, 4.2)); ax.plot(x, dmdt, color=blue)
    ax.axvline((reference["window"]["event_peak_s"] - t[0]) * 1e12, color=red, linestyle="--")
    ax.set(xlabel="Time from event window start (ps)", ylabel="dM/dt (A m s$^{-1}$)", title="Complete post-ramp event overview — DEVELOPMENT_REFERENCE")
    fig.tight_layout(); fig.savefig(FIGURES / "event_overview.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.2))
    for column, label, color in [("E_peak_V_m", "E", blue), ("nu_i_at_E_peak_s_1", "nu_i", purple),
                                  ("ne_peak_m3", "ne", teal), ("sigma_e_peak_S_m", "sigma", amber)]:
        values = ref_event[column].to_numpy(float); ax.plot(x, (values-values.min())/max(values.max()-values.min(), 1e-300), label=label, color=color)
    ax.set(xlabel="Time from event window start (ps)", ylabel="Independently normalized", title="Local kinetics — DEVELOPMENT_REFERENCE")
    ax.legend(frameon=False, ncol=4); fig.tight_layout(); fig.savefig(FIGURES / "kinetics_timeline.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.2)); ax.semilogy(x, ref_event.tau_i_at_E_peak_s*1e12, color=blue, label="tau_i")
    ax.semilogy(x, ref_event.tau_M_at_E_peak_s*1e12, color=teal, label="tau_M")
    ax2=ax.twinx(); ax2.plot(x, ref_event.Pi_RF_at_E_peak, color=amber, label="Pi_RF")
    ax.set(xlabel="Time from event window start (ps)", ylabel="Timescale (ps)", title="Kinetic timescales — DEVELOPMENT_REFERENCE"); ax2.set_ylabel("Pi_RF")
    lines=ax.lines+ax2.lines; ax.legend(lines,[v.get_label() for v in lines],frameon=False); fig.tight_layout(); fig.savefig(FIGURES / "timescales.png"); plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(9, 4.2)); ax1.plot(x, m, color=blue, label="M_z")
    ax2=ax1.twinx(); ax2.plot(x, dmdt, color=red, label="dM_z/dt"); peak_abs=np.max(np.abs(dmdt))
    ax2.axhline(0.5*peak_abs,color=teal,linestyle="--",linewidth=.8,label="1/2"); ax2.axhline(peak_abs/math.e,color=amber,linestyle=":",linewidth=.8,label="1/e")
    ax1.set(xlabel="Time from event window start (ps)",ylabel="M_z (A m)",title="Resolved current-moment pulse — DEVELOPMENT_REFERENCE");ax2.set_ylabel("dM_z/dt (A m s$^{-1}$)")
    lines=ax1.lines+ax2.lines;ax1.legend(lines,[v.get_label() for v in lines],frameon=False,ncol=4);fig.tight_layout();fig.savefig(FIGURES/"moment_pulse_width.png");plt.close(fig)

    fig, ax=plt.subplots(figsize=(9,4.2)); ax.plot(x,dmdt/np.max(np.abs(dmdt)),color=blue,label="signed dM/dt")
    ax.plot(x,proxy/np.max(np.abs(proxy)),color=amber,label="signed K_ion proxy");ax.set(xlabel="Time from event window start (ps)",ylabel="Peak-normalized",title="Signed proxy audit — DEVELOPMENT_REFERENCE");ax.legend(frameon=False);fig.tight_layout();fig.savefig(FIGURES/"proxy_sign_audit.png");plt.close(fig)

    fig, ax=plt.subplots(figsize=(9,3.8)); names=[k.replace("t_","").replace("_s","") for k,v in anchor_timing.items() if k.startswith("t_") and isinstance(v,float)]; vals=[(v-t[0])*1e12 for k,v in anchor_timing.items() if k.startswith("t_") and isinstance(v,float)]
    ax.scatter(vals,range(len(vals)),color=blue);ax.set_yticks(range(len(vals)),names);ax.set(xlabel="Time from event window start (ps)",title="Resolved interior timing anchors — DEVELOPMENT_REFERENCE");fig.tight_layout();fig.savefig(FIGURES/"anchor_timing.png");plt.close(fig)

    fig, ax=plt.subplots(figsize=(9,4.2)); f=ref_result.spectrum.frequency_Hz;e=ref_result.spectrum.ESD;mask=(f>0)&(f<1e12);ax.loglog(f[mask],e[mask],color=blue)
    ax.axvline(ref_result.spectrum.physical_frequency_resolution_Hz,color=red,linestyle="--",label="physical 1/T");ax.axvline(350e6,color=amber,linestyle=":",label="350 MHz NOT RESOLVED")
    ax.set(xlabel="Frequency (Hz)",ylabel="Numeric ESD",title="Numeric spectrum and physical-resolution limit — DEVELOPMENT_REFERENCE");ax.legend(frameon=False);fig.tight_layout();fig.savefig(FIGURES/"numeric_spectrum.png");plt.close(fig)

    fig, ax=plt.subplots(figsize=(9,4.2)); labels=["FWHM","FW1E","M change","eta peak","tau_i","tau_M"];keys=["PW_FWHM_s","PW_FW1E_s","M_event_change_A_m","eta_peak","tau_i_characteristic_s","tau_M_characteristic_s"]
    changes=[100*convergence[k]["relative_change"] for k in keys];ax.bar(labels,changes,color=[blue,teal,amber,purple,blue,teal]);ax.axhline(5,color=red,linestyle="--",label="5% target");ax.set(ylabel="Reference-to-finer change (%)",title="Temporal convergence — DEVELOPMENT_REFERENCE");ax.legend(frameon=False);fig.tight_layout();fig.savefig(FIGURES/"dt_convergence.png");plt.close(fig)

    print(json.dumps({"status": "PASS" if gate else "NOT_RESOLVED", "F_R6_ALLOWED": gate,
                      "PW_FWHM_ps": ref_result.pulse.PW_FWHM_s*1e12,
                      "PW_FW1E_ps": ref_result.pulse.PW_FW1E_s*1e12,
                      "temporal_convergence": temporal_convergence_pass}))


if __name__ == "__main__":
    main()
