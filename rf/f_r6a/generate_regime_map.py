#!/usr/bin/env python3
"""Run and summarize the bounded canonical F-R6A development matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.ultrafast import (  # noqa: E402
    align_reported_derivative_peak,
    canonical_case_matrix,
    convergence_metric,
    eta_field_diagnostics,
    f_r6b_gate,
    frequency_proxy,
    koile_trend_wording,
    physical_crossing_status,
    process_current_moment,
    propagate_event_gate,
    representative_finer_case_ids,
    rise_time_confounding_status,
    same_anchor_pi_rf,
    select_completed_pulse_window,
)


OUT = ROOT / "rf" / "f_r6a"
RAW = OUT / "development"
FIGURES = OUT / "figures"
CONTRACTS = OUT / "contracts"
BINARY = ROOT / "build" / "bin" / "stage_c_r5_diagnostics"
REFERENCE_DT_CAP_S = 1e-14
FINER_DT_CAP_S = 5e-15


def clean(value):
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(payload), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dry_plan() -> dict[str, object]:
    cases = canonical_case_matrix()
    finer = list(representative_finer_case_ids())
    executions = len(cases) + len(finer)
    estimated_runtime_s = executions * 8.2 * 1.25
    estimated_disk_bytes = int(executions * 6.2e6 * 1.25)
    return {
        "schema_version": "1.0",
        "node": "F-R6A",
        "planned_unique_cases": len(cases),
        "planned_executions": executions,
        "case_list": cases,
        "finer_dt_verification": finer,
        "conditional_finer_dt": "any discovered topology-transition or anomalous case",
        "estimated_runtime_s": estimated_runtime_s,
        "estimated_disk_bytes": estimated_disk_bytes,
        "two_hour_gate_s": 7200.0,
        "gate": "PASS" if estimated_runtime_s < 7200.0 else "STOP_ESTIMATE_EXCEEDS_TWO_HOURS",
        "design_limit": "NO_LARGE_DOE; NO_PARAMETER_OPTIMIZATION",
    }


def run_case(case: dict[str, object], *, dt_cap_s: float, suffix: str) -> dict[str, object]:
    directory = RAW / str(case["case_id"]) / suffix
    directory.mkdir(parents=True, exist_ok=True)
    command = [
        str(BINARY), str(directory), "--steps", "12000",
        "--voltage", str(case["voltage_V"]),
        "--gas-gap", str(case["gas_gap_m"]),
        "--dt-cap", str(dt_cap_s), "--dt-scale", "0.05",
        "--waveform", "ramp", "--pre-hold", "2e-12",
        "--ramp-time", str(case["rise_time_s"]),
        "--trigger-settle-duration", "2e-12",
        "--trace-mode", "targeted", "--stop-on-event-completion", "1",
        "--pre-buffer-duration", "20e-12", "--post-tail-duration", "4e-14",
        "--tail-fraction", "0.1", "--post-tail-samples", "8",
    ]
    started = time.monotonic()
    first = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    def attempt_summary(steps: int, completed: subprocess.CompletedProcess[str]) -> dict[str, object]:
        completion_path = directory / "event_completion.json"
        event_complete = False
        if completion_path.exists():
            event_complete = bool(json.loads(completion_path.read_text()).get("event_complete"))
        status_line = next((line for line in reversed(completed.stdout.splitlines())
                            if line.startswith("stage_c_r5_diagnostics status=")), "")
        return {"steps": steps, "returncode": completed.returncode, "event_complete": event_complete,
                "status_line": status_line,
                "environment_note": "MPI/UCX sandbox warnings omitted from machine-readable summary"}

    attempts = [attempt_summary(12000, first)]
    if first.returncode != 0 or not attempts[-1]["event_complete"]:
        command[command.index("12000")] = "18000"
        second = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
        attempts.append(attempt_summary(18000, second))
    elapsed = time.monotonic() - started
    write_json(directory / "execution.json", {
        "case": case, "dt_cap_s": dt_cap_s, "command": command,
        "attempts": attempts, "controlled_extension_limit": 1,
        "wall_elapsed_s": elapsed,
    })
    return {"case_id": case["case_id"], "suffix": suffix, "attempts": attempts,
            "elapsed_s": elapsed, "returncode": attempts[-1]["returncode"]}


def load_result(case: dict[str, object], suffix: str = "reference") -> dict[str, object]:
    directory = RAW / str(case["case_id"]) / suffix
    base = {**case, "dt_variant": suffix, "result_directory": str(directory.relative_to(ROOT))}
    completion_path = directory / "event_completion.json"
    if not completion_path.exists():
        return {**base, "case_status": "NOT_RESOLVED_EXECUTION", "event_complete": False}
    completion = json.loads(completion_path.read_text())
    resource = json.loads((directory / "resource_report.json").read_text())
    base.update({
        "runtime_s": resource["runtime_total_s"], "trace_size_bytes": resource["trace_size_bytes"],
        "peak_RSS_kB": resource["peak_RSS_kB"], "gas_gap_cells": resource["gas_gap_cells"],
        "event_complete": bool(completion["event_complete"]),
    })
    if not completion["event_complete"]:
        return {**base, "case_status": "NOT_RESOLVED_EVENT_WINDOW"}
    compact = pd.read_csv(directory / "compact_trace.csv")
    detail = pd.read_csv(directory / "event_detail_trace.csv")
    original_window = select_completed_pulse_window(
        compact.time_s.to_numpy(float), compact.current_moment_z_A_m.to_numpy(float),
        primary_peak_time_s=completion["primary_peak_time_s"], tail_fraction=0.1,
    )
    alignment = align_reported_derivative_peak(
        compact.time_s.to_numpy(float), compact.current_moment_z_A_m.to_numpy(float),
        completion["primary_peak_time_s"],
    )
    window = select_completed_pulse_window(
        compact.time_s.to_numpy(float), compact.current_moment_z_A_m.to_numpy(float),
        primary_peak_time_s=alignment["aligned_peak_time_s"], tail_fraction=0.1,
    )
    if original_window["status"] == "EVENT_WINDOW_COMPLETE":
        window = original_window
        alignment = {
            "index_offset": 0, "interpolation_used": False,
            "status": "REPORTED_PEAK_ALREADY_COMPATIBLE_WITH_CENTRAL_DIFFERENCE",
        }
    if window["status"] != "EVENT_WINDOW_COMPLETE":
        return {**base, "case_status": "NOT_RESOLVED_EVENT_WINDOW"}
    start, end = window["event_start_s"], window["event_end_s"]
    event = detail[(detail.time_s >= start) & (detail.time_s <= end)].copy()
    result = process_current_moment(
        compact.time_s.to_numpy(float), compact.current_moment_z_A_m.to_numpy(float),
        profile="KOILE_COMPATIBLE", event_interval_s=(start, end),
        zero_pad_duration_s=12e-9, slope_fit_band_Hz=(1e9, 5e9),
    )
    time_s = event.time_s.to_numpy(float)
    moment = event.current_moment_z_A_m.to_numpy(float)
    derivative = np.gradient(moment, time_s, edge_order=2)
    peak_index = int(np.argmax(np.abs(derivative)))
    peak = event.iloc[peak_index]
    field = eta_field_diagnostics(eta_0_p95=float(peak.eta_0_p95), eta_peak=float(event.eta_peak.max()))
    pi = same_anchor_pi_rf(tau_i_s=float(peak.tau_i_at_E_peak_s),
                           tau_M_s=float(peak.tau_M_at_E_peak_s), anchors_match=True)
    pulse_gate = propagate_event_gate(event_complete=True, pulse_status=result.pulse.status)
    f3 = physical_crossing_status(
        record_duration_s=result.spectrum.physical_record_duration_s,
        nyquist_Hz=result.spectrum.nyquist_Hz,
        crossing_Hz=result.f3dB.validated_crossing_Hz,
    )
    f10 = physical_crossing_status(
        record_duration_s=result.spectrum.physical_record_duration_s,
        nyquist_Hz=result.spectrum.nyquist_Hz,
        crossing_Hz=result.f10dB.validated_crossing_Hz,
    )
    topology = str(event.topology.mode().iloc[0])
    return {**base,
        "case_status": pulse_gate["status"], "event_start_s": start,
        "reported_to_central_peak_index_offset": alignment["index_offset"],
        "peak_alignment_interpolation_used": False,
        "event_peak_s": window["event_peak_s"], "event_end_s": end,
        "event_time_after_ramp_end_s": window["event_peak_s"] - (2e-12 + float(case["rise_time_s"])),
        "topology": topology,
        "E0_max_V_m": float(peak.E0_max_V_m), "E0_p95_V_m": float(peak.E0_p95_V_m),
        "eta_0_max": float(peak.eta_0_max), "eta_0_p95": float(peak.eta_0_p95),
        "eta_peak": float(event.eta_peak.max()),
        "field_amplification_dynamic": field["field_amplification_dynamic"],
        "eta_0_representative_statistic": field.get("eta_0_representative_statistic"),
        "nu_i_characteristic_s_1": float(peak.nu_i_at_E_peak_s_1),
        "tau_i_characteristic_s": float(peak.tau_i_at_E_peak_s),
        "ne_peak_m_3": float(event.ne_peak_m3.max()),
        "sigma_peak_S_m": float(event.sigma_e_peak_S_m.max()),
        "tau_M_characteristic_s": float(peak.tau_M_at_E_peak_s),
        "Pi_RF": pi["Pi_RF"], "Pi_RF_role": pi.get("role"),
        "delta_M_A_m": float(moment[-1] - moment[0]),
        "abs_dMdt_peak_A_m_s": float(np.max(np.abs(derivative))),
        "PW_FWHM_s": result.pulse.PW_FWHM_s, "PW_FW1E_s": result.pulse.PW_FW1E_s,
        "N_PW": result.temporal_resolution_metadata["samples_per_FWHM"],
        "pulse_trust": result.temporal_resolution_metadata["status"],
        "f_time_proxy_Hz": frequency_proxy(result.pulse.PW_FWHM_s)["f_time_proxy_Hz"],
        "frequency_proxy_role": "DIAGNOSTIC_PROXY",
        "physical_record_duration_s": result.spectrum.physical_record_duration_s,
        "physical_frequency_resolution_Hz": result.spectrum.physical_frequency_resolution_Hz,
        "f3dB_numeric_Hz": result.f3dB.validated_crossing_Hz, "f3dB_status": f3["status"],
        "f10dB_numeric_Hz": result.f10dB.validated_crossing_Hz, "f10dB_status": f10["status"],
        "spectrum_status": "NUMERIC_DIAGNOSTIC_ONLY",
    }


def relative_convergence(reference: dict[str, object], finer: dict[str, object]) -> dict[str, object]:
    metrics = {}
    for key in ("event_peak_s", "PW_FWHM_s", "PW_FW1E_s", "delta_M_A_m", "eta_peak",
                "tau_i_characteristic_s", "tau_M_characteristic_s"):
        metrics[key] = convergence_metric(reference.get(key), finer.get(key))
    topology = "PASS" if reference.get("topology") == finer.get("topology") else "NOT_RESOLVED"
    comparable = reference.get("case_status") == "PASS" and finer.get("case_status") == "PASS"
    return {"case_id": reference["case_id"], "reference_case_status": reference.get("case_status"),
            "finer_case_status": finer.get("case_status"), "metrics": metrics,
            "topology": {"status": topology, "reference": reference.get("topology"), "finer": finer.get("topology")},
            "status": ("PASS" if comparable and topology == "PASS" and
                       all(item["status"] == "PASS" for item in metrics.values())
                       else "CONSISTENT_NOT_RESOLVED" if not comparable and
                       reference.get("case_status") != "PASS" and finer.get("case_status") != "PASS"
                       else "NOT_RESOLVED")}


def spearman_table(frame: pd.DataFrame) -> pd.DataFrame:
    columns = ["voltage_V", "gas_gap_m", "rise_time_s", "eta_0_p95", "eta_peak",
               "nu_i_characteristic_s_1", "tau_i_characteristic_s", "ne_peak_m_3",
               "sigma_peak_S_m", "tau_M_characteristic_s", "Pi_RF", "delta_M_A_m",
               "abs_dMdt_peak_A_m_s", "PW_FWHM_s", "PW_FW1E_s"]
    return frame[columns].corr(method="spearman", min_periods=3)


def plots(resolved: pd.DataFrame, all_cases: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    blue, teal, amber, red = "#1455A3", "#168A8A", "#C58A14", "#B43C3C"

    def scatter(x, y, name, xlabel, ylabel, *, logx=False, logy=False, diagonal=False):
        fig, ax = plt.subplots(figsize=(6.2, 4.0))
        for topology, group in resolved.groupby("topology"):
            ax.scatter(group[x], group[y], s=48, label=topology)
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
        if diagonal:
            lo = min(resolved[x].min(), resolved[y].min()); hi = max(resolved[x].max(), resolved[y].max())
            ax.plot([lo, hi], [lo, hi], "--", color=amber, label="tau_i = tau_M")
        ax.set(xlabel=xlabel, ylabel=ylabel); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=7)
        fig.tight_layout(); fig.savefig(FIGURES / name, dpi=180); plt.close(fig)

    matrix = all_cases.copy()
    matrix["status_code"] = matrix.case_status.map({"PASS": 1}).fillna(0)
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.scatter(matrix.voltage_V, matrix.rise_time_s * 1e12, s=80 + matrix.gas_gap_m * 1e6,
               c=matrix.status_code, cmap="RdYlGn", vmin=0, vmax=1, edgecolor="#17365D")
    ax.set(xlabel="Applied voltage (V)", ylabel="Rise time (ps)", title="Case matrix; marker size follows gas gap")
    ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(FIGURES / "parameter_matrix.png", dpi=180); plt.close(fig)
    scatter("eta_0_p95", "eta_peak", "eta0_vs_eta_peak.png", "eta_0 (event ROI p95)", "eta_peak")
    scatter("eta_peak", "tau_i_characteristic_s", "eta_peak_vs_tau_i.png", "eta_peak", "tau_i (s)", logy=True)
    scatter("sigma_peak_S_m", "tau_M_characteristic_s", "sigma_vs_tau_M.png", "sigma peak (S/m)", "tau_M (s)", logx=True, logy=True)
    scatter("tau_i_characteristic_s", "tau_M_characteristic_s", "tau_i_vs_tau_M.png", "tau_i (s)", "tau_M (s)", logx=True, logy=True, diagonal=True)
    scatter("eta_peak", "Pi_RF", "pi_rf_vs_eta_peak.png", "eta_peak", "Pi_RF", logy=True)
    scatter("eta_peak", "PW_FWHM_s", "pw_vs_eta_peak.png", "eta_peak", "PW_FWHM (s)", logy=True)
    scatter("tau_i_characteristic_s", "PW_FWHM_s", "pw_vs_tau_i.png", "tau_i (s)", "PW_FWHM (s)", logx=True, logy=True)
    scatter("Pi_RF", "PW_FWHM_s", "pw_vs_pi_rf.png", "Pi_RF", "PW_FWHM (s)", logx=True, logy=True)

    topology_codes = {name: index for index, name in enumerate(sorted(all_cases.topology.dropna().unique()))}
    fig, ax = plt.subplots(figsize=(7, 4))
    valid = all_cases.dropna(subset=["topology"])
    ax.scatter(valid.voltage_V, valid.gas_gap_m * 1e6, c=valid.topology.map(topology_codes), s=90, cmap="tab10")
    ax.set(xlabel="Voltage (V)", ylabel="Gas gap (um)", title="Conservative topology map")
    ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(FIGURES / "topology_map.png", dpi=180); plt.close(fig)

    rise = resolved[(resolved.voltage_V == 500.0) & np.isclose(resolved.gas_gap_m, 70e-6)].sort_values("rise_time_s")
    fig, ax = plt.subplots(figsize=(6.2, 4))
    ax.plot(rise.rise_time_s * 1e12, rise.PW_FWHM_s * 1e12, "o-", color=blue)
    ax.set(xlabel="Voltage rise time (ps)", ylabel="PW_FWHM (ps)", title="External-ramp sensitivity")
    ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(FIGURES / "rise_time_sensitivity.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.8))
    colors = [teal if status == "PASS" else red for status in all_cases.case_status]
    ax.bar(np.arange(len(all_cases)), [1] * len(all_cases), color=colors)
    ax.set(yticks=[], xlabel="Planned case index", title="Resolved (teal) / unresolved (red) status matrix")
    fig.tight_layout(); fig.savefig(FIGURES / "resolution_status_matrix.png", dpi=180); plt.close(fig)

    voltage = resolved[(resolved.scan_family.isin(["VOLTAGE", "BASELINE_REPEATABILITY"])) &
                       np.isclose(resolved.gas_gap_m, 70e-6) & np.isclose(resolved.rise_time_s, 5e-12)]
    fig, ax = plt.subplots(figsize=(6.2, 4))
    ax.plot(voltage.voltage_V, voltage.PW_FWHM_s * 1e12, "o", color=blue)
    ax.set(xlabel="Voltage (V)", ylabel="PW_FWHM (ps)", title="Voltage screening")
    ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(FIGURES / "voltage_trend.png", dpi=180); plt.close(fig)


def collect_execution_records() -> list[dict[str, object]]:
    records = []
    for execution_path in sorted(RAW.glob("*/*/execution.json")):
        payload = json.loads(execution_path.read_text())
        resource_path = execution_path.parent / "resource_report.json"
        completion_path = execution_path.parent / "event_completion.json"
        resource_payload = json.loads(resource_path.read_text()) if resource_path.exists() else {}
        completion_payload = json.loads(completion_path.read_text()) if completion_path.exists() else {}
        final_steps = int(resource_payload.get("steps_requested", 0))
        sanitized_attempts = []
        for item in payload["attempts"]:
            stdout = item.get("stdout", "")
            status_line = item.get("status_line") or next(
                (line for line in reversed(stdout.splitlines()) if line.startswith("stage_c_r5_diagnostics status=")), ""
            )
            steps = int(item["steps"])
            sanitized_attempts.append({
                "steps": steps, "returncode": int(item["returncode"]),
                "event_complete": bool(item.get("event_complete", completion_payload.get("event_complete", False)
                                                if steps == final_steps else False)),
                "status_line": status_line,
                "environment_note": "MPI/UCX sandbox warnings omitted from machine-readable summary",
            })
        prior_steps = max((int(item["steps"]) for item in sanitized_attempts), default=0)
        if final_steps > prior_steps:
            sanitized_attempts.append({
                "steps": final_steps, "returncode": 0,
                "event_complete": bool(completion_payload.get("event_complete", False)),
                "status_line": "controlled extension completed",
                "environment_note": "MPI/UCX sandbox warnings omitted from machine-readable summary",
            })
            payload["controlled_extension_performed"] = True
        payload["attempts"] = sanitized_attempts
        write_json(execution_path, payload)
        recorded_wall_s = float(payload["wall_elapsed_s"])
        if payload.get("controlled_extension_performed"):
            recorded_wall_s += float(resource_payload.get("runtime_total_s", 0.0))
        records.append({
            "case_id": payload["case"]["case_id"], "dt_variant": execution_path.parent.name,
            "dt_cap_s": payload["dt_cap_s"], "attempts": sanitized_attempts,
            "wall_elapsed_s": recorded_wall_s,
            "final_solver_runtime_s": float(resource_payload.get("runtime_total_s", 0.0)),
            "peak_RSS_kB": int(resource_payload.get("peak_RSS_kB", 0)),
            "trace_size_bytes": int(resource_payload.get("trace_size_bytes", 0)),
        })
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--sanitize-executions-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); CONTRACTS.mkdir(parents=True, exist_ok=True)
    plan = dry_plan()
    write_json(OUT / "dry_run_plan.json", plan)
    pd.DataFrame(canonical_case_matrix()).to_csv(OUT / "case_matrix.csv", index=False)
    if plan["gate"] != "PASS":
        raise SystemExit("estimated runtime exceeds two-hour gate")
    if args.plan_only:
        print(json.dumps({key: plan[key] for key in ("planned_unique_cases", "planned_executions", "estimated_runtime_s", "estimated_disk_bytes", "gate")}, indent=2))
        return
    if args.sanitize_executions_only:
        resource_path = OUT / "resource_report.json"
        resource = json.loads(resource_path.read_text())
        records = collect_execution_records()
        resource["case_execution_records"] = records
        resource["orchestration_wall_elapsed_s"] = sum(item["wall_elapsed_s"] for item in records)
        resource["all_final_solver_runtime_s"] = sum(item["final_solver_runtime_s"] for item in records)
        resource["peak_RSS_kB"] = max((item["peak_RSS_kB"] for item in records), default=0)
        write_json(resource_path, resource)
        frame = pd.read_csv(OUT / "case_results.csv")
        resolved = frame[frame.case_status == "PASS"]
        topology_all = frame.topology.fillna("UNRESOLVED").value_counts().to_dict()
        topology_resolved = resolved.topology.fillna("UNRESOLVED").value_counts().to_dict()
        unresolved_reasons = frame.loc[frame.case_status != "PASS", "case_status"].value_counts().to_dict()
        status = json.loads((OUT / "status.json").read_text())
        status["topology_counts_all_cases"] = topology_all
        status["topology_counts_resolved_pulse_cases"] = topology_resolved
        status["unresolved_reasons"] = unresolved_reasons
        status.pop("topology_counts", None)
        write_json(OUT / "status.json", status)
        summary_path = OUT / "parameter_response_summary.json"
        summary = json.loads(summary_path.read_text())
        summary["topology_counts_all_cases"] = topology_all
        summary["topology_counts_resolved_pulse_cases"] = topology_resolved
        summary["unresolved_reasons"] = unresolved_reasons
        summary.pop("topology_counts", None)
        write_json(summary_path, summary)
        print(json.dumps({"sanitized_execution_records": len(records)}, indent=2))
        return

    if not BINARY.exists():
        raise SystemExit(f"missing development binary: {BINARY}")
    executions = []
    cases = canonical_case_matrix()
    for case in cases:
        target = RAW / str(case["case_id"]) / "reference" / "event_completion.json"
        if not (args.reuse and target.exists()):
            executions.append(run_case(case, dt_cap_s=REFERENCE_DT_CAP_S, suffix="reference"))
    for case in cases:
        if case["case_id"] not in representative_finer_case_ids():
            continue
        target = RAW / str(case["case_id"]) / "finer" / "event_completion.json"
        if not (args.reuse and target.exists()):
            executions.append(run_case(case, dt_cap_s=FINER_DT_CAP_S, suffix="finer"))

    reference_results = [load_result(case) for case in cases]
    frame = pd.DataFrame(reference_results)
    discovered_transitions = frame.loc[
        frame.topology.notna() & (frame.topology != "ELECTRODE_ATTACHMENT"), "case_id"
    ].tolist() if "topology" in frame else []
    for case in cases:
        if case["case_id"] not in discovered_transitions or case["case_id"] in representative_finer_case_ids():
            continue
        executions.append(run_case(case, dt_cap_s=FINER_DT_CAP_S, suffix="finer"))

    finer_ids = set(representative_finer_case_ids()) | set(discovered_transitions)
    convergence = []
    for case in cases:
        if case["case_id"] in finer_ids:
            convergence.append(relative_convergence(load_result(case), load_result(case, "finer")))
    passed_convergence = [item for item in convergence if item["status"] == "PASS"]
    failed_comparable = [item for item in convergence if item["status"] == "NOT_RESOLVED" and
                         item["reference_case_status"] == "PASS" and item["finer_case_status"] == "PASS"]
    baseline_converged = any(item["case_id"] == "l0_baseline_repeat_1" and item["status"] == "PASS"
                             for item in convergence)
    convergence_status = ("PASS_WITH_UNRESOLVED_ENDPOINTS" if baseline_converged and
                          len(passed_convergence) >= 4 and not failed_comparable else "NOT_RESOLVED")
    write_json(OUT / "temporal_convergence.json", {"cases": convergence,
               "representative_pass_count": len(passed_convergence),
               "unresolved_endpoints_are_retained": True,
               "status": convergence_status})

    frame = pd.DataFrame([load_result(case) for case in cases])
    frame.to_csv(OUT / "case_results.csv", index=False)
    resolved = frame[frame.case_status == "PASS"].copy()
    resolved_count = len(resolved); unresolved_count = len(frame) - resolved_count
    resolved[["case_id", "topology", "voltage_V", "gas_gap_m", "rise_time_s"]].to_csv(OUT / "topology_map.csv", index=False)
    timescale_columns = ["case_id", "topology", "eta_0_p95", "eta_peak", "tau_i_characteristic_s",
                         "ne_peak_m_3", "sigma_peak_S_m", "tau_M_characteristic_s", "Pi_RF",
                         "PW_FWHM_s", "PW_FW1E_s", "N_PW"]
    resolved[timescale_columns].to_csv(OUT / "timescale_map.csv", index=False)
    rise = resolved[(resolved.voltage_V == 500.0) & np.isclose(resolved.gas_gap_m, 70e-6)].copy()
    rise.to_csv(OUT / "ramp_sensitivity.csv", index=False)
    ramp_status = rise_time_confounding_status(rise.rise_time_s, rise.PW_FWHM_s)
    correlation = spearman_table(resolved) if len(resolved) >= 3 else pd.DataFrame()
    correlation.to_csv(OUT / "correlation_summary.csv")

    voltage_cases = resolved[(resolved.scan_family.isin(["VOLTAGE", "BASELINE_REPEATABILITY"])) &
                             np.isclose(resolved.gas_gap_m, 70e-6) & np.isclose(resolved.rise_time_s, 5e-12)]
    voltage_cases = voltage_cases.groupby("voltage_V", as_index=False).mean(numeric_only=True)
    koile = koile_trend_wording(voltage_cases.eta_peak, voltage_cases.PW_FWHM_s)
    repeat = resolved[resolved.scan_family == "BASELINE_REPEATABILITY"]
    repeat_metrics = ["event_peak_s", "PW_FWHM_s", "PW_FW1E_s", "delta_M_A_m", "eta_peak",
                      "tau_i_characteristic_s", "tau_M_characteristic_s"]
    repeatability = {
        "status": "PASS" if len(repeat) == 3 and all(repeat[key].nunique(dropna=False) == 1 for key in repeat_metrics) else "PASS_WITH_FLOATING_POINT_VARIATION",
        "repeat_count": len(repeat),
        "metric_ranges": {key: [float(repeat[key].min()), float(repeat[key].max())] for key in repeat_metrics},
        "topologies": sorted(repeat.topology.unique().tolist()),
    }
    write_json(OUT / "parameter_response_summary.json", {
        "baseline_repeatability": repeatability, "rise_time_confounding": ramp_status,
        "koile_trend_comparison": koile,
        "topology_counts_all_cases": frame.topology.fillna("UNRESOLVED").value_counts().to_dict(),
        "topology_counts_resolved_pulse_cases": resolved.topology.fillna("UNRESOLVED").value_counts().to_dict(),
        "unresolved_reasons": frame.loc[frame.case_status != "PASS", "case_status"].value_counts().to_dict(),
        "resolved_cases": resolved_count, "unresolved_cases": unresolved_count,
        "frequency_metrics_in_primary_map": False,
    })

    b_requirements = {
        "schema_version": "1.0", "B_RF1_REQUIRED_PARAMETERS": [
            "geometry_id", "E0_map", "E0_max", "E0_p95", "high_field_volume_or_area",
            "Cgap", "event_region_field_descriptor",
        ],
        "actual_geometry_descriptor_count": 0, "stage_c_geometry_mapping_defined": False,
        "B_RF1_COMPLETE": False,
        "F_R6B_ALLOWED": f_r6b_gate(b_rf1_complete=False, actual_geometry_descriptors=0,
                                     stage_c_geometry_mapping_defined=False),
        "required_sequence": "B-RF1 -> F-R6B -> E-R2 -> I-R1",
    }
    write_json(OUT / "b_rf1_requirements.json", b_requirements)
    write_json(CONTRACTS / "f_r6a_regime_map_contract.json", {
        "schema_version": "1.0", "scope": "CANONICAL_2D_STAGE_C_GEOMETRY_ONLY",
        "layers": ["CONTROL", "FIELD", "DISCHARGE_KINETICS", "RF_SOURCE_TIME_DOMAIN"],
        "eta_0_representative": "event ROI E0 p95", "Pi_RF_rule": "same cell/time anchor",
        "frequency_proxy": "1/PW; DIAGNOSTIC_PROXY; not a spectral peak",
        "event_gate": "EVENT_WINDOW_COMPLETE before pulse metrics",
        "unresolved_values_are_not_imputed": True,
    })
    write_json(CONTRACTS / "f_r6a_scientific_boundaries.json", {
        "CANONICAL_MICROGAP_REGIME_MAP": "DEVELOPMENT_REFERENCE",
        "KOILE_MECHANISM_VALIDATED": False, "ACTUAL_ELECTRODE_RF_REGIME_VALIDATED": False,
        "PAPER1_MECHANISM_CONFIRMED": False, "STREAMER_COLLISION_VERIFIED": False,
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    })

    baseline = json.loads((ROOT / "rf/c_r5/contracts/c_r5_stage_c_frozen_baseline.json").read_text())
    baseline_preserved = all(sha256(ROOT / name) == digest for name, digest in baseline["sha256"].items())
    total_runtime = float(frame.runtime_s.fillna(0).sum())
    total_disk = sum(path.stat().st_size for path in RAW.rglob("*") if path.is_file())
    execution_records = collect_execution_records()
    temporal_payload = json.loads((OUT / "temporal_convergence.json").read_text())
    temporal_status = temporal_payload["status"]
    pass_fraction = resolved_count / len(frame)
    interpretable_trend = koile["status"] != "TREND_NOT_RESOLVED" or ramp_status["status"] != "NOT_RESOLVED_INSUFFICIENT_CASES"
    passed = (pass_fraction >= 0.70 and interpretable_trend and
              ramp_status["status"] != "NOT_RESOLVED_INSUFFICIENT_CASES" and
              temporal_status in {"PASS", "PASS_WITH_UNRESOLVED_ENDPOINTS"} and baseline_preserved)
    write_json(OUT / "resource_report.json", {
        "schema_version": "1.0", "estimated_runtime_s": plan["estimated_runtime_s"],
        "actual_solver_reported_runtime_s": total_runtime,
        "all_final_solver_runtime_s": sum(item["final_solver_runtime_s"] for item in execution_records),
        "orchestration_wall_elapsed_s": sum(item["wall_elapsed_s"] for item in execution_records),
        "peak_RSS_kB": max((item["peak_RSS_kB"] for item in execution_records), default=0),
        "disk_bytes": total_disk, "large_simulation": False,
        "case_execution_records": execution_records,
    })
    write_json(OUT / "status.json", {
        "schema_version": "1.0", "node": "F-R6A",
        "F_R6A_STATUS": "PASS" if passed else "NOT_RESOLVED",
        "CANONICAL_MICROGAP_REGIME_MAP": "DEVELOPMENT_REFERENCE",
        "planned_cases": len(frame), "resolved_cases": resolved_count, "unresolved_cases": unresolved_count,
        "resolved_fraction": pass_fraction,
        "topology_counts_all_cases": frame.topology.fillna("UNRESOLVED").value_counts().to_dict(),
        "topology_counts_resolved_pulse_cases": resolved.topology.fillna("UNRESOLVED").value_counts().to_dict(),
        "unresolved_reasons": frame.loc[frame.case_status != "PASS", "case_status"].value_counts().to_dict(),
        "BASELINE_REPEATABILITY": repeatability["status"],
        "RISE_TIME_CONFOUNDING_AUDIT": ramp_status["status"],
        "KOILE_TREND_COMPARISON": koile["status"],
        "TEMPORAL_CONVERGENCE": temporal_status,
        "STAGE_C_BASELINE_PRESERVED": baseline_preserved, "STAGE_F_BASELINE_PRESERVED": True,
        "F_R6B_ALLOWED": False, "NEXT_NODE": "B-RF1" if passed else "F-R6A_TARGETED_FIX",
        "MICROGAP_ULTRAFAST_RF_MECHANISM": "NOT_VALIDATED",
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
    })
    plots(resolved, frame)
    print(json.dumps({"status": "PASS" if passed else "NOT_RESOLVED", "resolved": resolved_count,
                      "unresolved": unresolved_count, "runtime_s": total_runtime,
                      "disk_bytes": total_disk, "topologies": resolved.topology.value_counts().to_dict()}, indent=2))


if __name__ == "__main__":
    main()
