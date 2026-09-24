"""Deterministic governance helpers for the canonical F-R6A regime map."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def canonical_case_matrix() -> list[dict[str, object]]:
    """Return the bounded 15-case screening and interaction design."""
    cases: list[dict[str, object]] = []

    def add(case_id: str, level: str, family: str, voltage: float, gap_um: float, rise_ps: float) -> None:
        cases.append({
            "case_id": case_id,
            "design_level": level,
            "scan_family": family,
            "voltage_V": voltage,
            "gas_gap_m": gap_um * 1e-6,
            "rise_time_s": rise_ps * 1e-12,
            "geometry_id": f"canonical-2d-needle-plane-gap-{gap_um:g}um",
            "gas": "air",
            "pressure_Pa": 101325.0,
            "temperature_K": 300.0,
        })

    for repeat in range(1, 4):
        add(f"l0_baseline_repeat_{repeat}", "LEVEL_0", "BASELINE_REPEATABILITY", 500.0, 70.0, 5.0)
    for voltage in (400.0, 450.0, 550.0, 600.0):
        add(f"l1_voltage_{int(voltage)}V", "LEVEL_1", "VOLTAGE", voltage, 70.0, 5.0)
    for gap_um in (60.0, 80.0):
        add(f"l1_gap_{int(gap_um)}um", "LEVEL_1", "GAP", 500.0, gap_um, 5.0)
    for rise_ps in (2.5, 10.0):
        label = str(rise_ps).replace(".", "p")
        add(f"l1_rise_{label}ps", "LEVEL_1", "RISE_TIME", 500.0, 70.0, rise_ps)
    for voltage in (450.0, 550.0):
        for rise_ps in (2.5, 10.0):
            label = str(rise_ps).replace(".", "p")
            add(f"l2_voltage_{int(voltage)}V_rise_{label}ps", "LEVEL_2", "VOLTAGE_X_RISE", voltage, 70.0, rise_ps)
    return cases


def representative_finer_case_ids() -> tuple[str, ...]:
    return (
        "l0_baseline_repeat_1",
        "l1_voltage_400V",
        "l1_voltage_600V",
        "l1_gap_60um",
        "l1_gap_80um",
        "l1_rise_2p5ps",
        "l1_rise_10p0ps",
    )


def align_reported_derivative_peak(
    time_s: np.ndarray, moment_A_m: np.ndarray, reported_peak_time_s: float, *, search_steps: int = 3,
) -> dict[str, object]:
    """Align a backward-difference endpoint timestamp to a nearby central-difference peak.

    This selects an existing accepted-step timestamp and never interpolates the trace.
    """
    time = np.asarray(time_s, dtype=float)
    moment = np.asarray(moment_A_m, dtype=float)
    if time.ndim != 1 or moment.shape != time.shape or time.size < 7 or np.any(np.diff(time) <= 0.0):
        raise ValueError("peak-alignment inputs must be equal, finite, increasing traces")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(moment)) or not math.isfinite(reported_peak_time_s):
        raise ValueError("peak-alignment inputs must be finite")
    reported_index = int(np.argmin(np.abs(time - reported_peak_time_s)))
    derivative = np.abs(np.gradient(moment, time, edge_order=2))
    start = max(1, reported_index - search_steps)
    stop = min(time.size - 1, reported_index + search_steps + 1)
    aligned_index = start + int(np.argmax(derivative[start:stop]))
    return {
        "reported_peak_time_s": float(reported_peak_time_s),
        "reported_nearest_index": reported_index,
        "aligned_peak_time_s": float(time[aligned_index]),
        "aligned_index": aligned_index,
        "index_offset": aligned_index - reported_index,
        "time_offset_s": float(time[aligned_index] - reported_peak_time_s),
        "search_steps": search_steps,
        "interpolation_used": False,
        "status": "ALIGNED_EXISTING_ACCEPTED_STEP",
    }


def eta_field_diagnostics(*, eta_0_p95: float, eta_peak: float) -> dict[str, object]:
    if not all(math.isfinite(value) and value > 0.0 for value in (eta_0_p95, eta_peak)):
        return {"status": "NOT_RESOLVED", "field_amplification_dynamic": None}
    return {
        "status": "PASS",
        "eta_0_representative_statistic": "EVENT_ROI_E0_P95",
        "field_amplification_dynamic": eta_peak / eta_0_p95,
    }


def same_anchor_pi_rf(*, tau_i_s: float, tau_M_s: float, anchors_match: bool) -> dict[str, object]:
    if not anchors_match:
        return {"status": "NOT_RESOLVED_ANCHOR_MISMATCH", "Pi_RF": None}
    if not all(math.isfinite(value) and value > 0.0 for value in (tau_i_s, tau_M_s)):
        return {"status": "NOT_RESOLVED_INVALID_TIMESCALE", "Pi_RF": None}
    return {"status": "PASS", "Pi_RF": tau_M_s / tau_i_s, "role": "DIAGNOSTIC_CANDIDATE"}


def propagate_event_gate(*, event_complete: bool, pulse_status: str) -> dict[str, object]:
    if not event_complete:
        return {"status": "NOT_RESOLVED_EVENT_WINDOW", "PW_trusted": False}
    if pulse_status != "PASS":
        return {"status": "NOT_RESOLVED_PULSE", "PW_trusted": False}
    return {"status": "PASS", "PW_trusted": True}


def rise_time_confounding_status(rise_time_s: Iterable[float], pulse_width_s: Iterable[float]) -> dict[str, object]:
    pairs = sorted((float(r), float(p)) for r, p in zip(rise_time_s, pulse_width_s)
                   if math.isfinite(r) and math.isfinite(p) and r > 0.0 and p > 0.0)
    if len(pairs) < 3:
        return {"status": "NOT_RESOLVED_INSUFFICIENT_CASES", "relative_span": None}
    rise = [item[0] for item in pairs]
    pulse = [item[1] for item in pairs]
    rise_rank = _ranks(rise)
    pulse_rank = _ranks(pulse)
    rho = _pearson(rise_rank, pulse_rank)
    span = (max(pulse) - min(pulse)) / max(min(pulse), 1e-300)
    status = "EXTERNAL_DRIVE_SENSITIVE" if abs(rho) >= 0.8 and span >= 0.05 else "INTRINSIC_TIMESCALE_EVIDENCE_STRENGTHENED"
    return {"status": status, "spearman_rho": rho, "relative_span": span,
            "interpretation_limit": "diagnostic evidence, not proof of a fully intrinsic timescale"}


def koile_trend_wording(eta: Iterable[float], pulse_width_s: Iterable[float]) -> dict[str, object]:
    pairs = [(float(e), float(p)) for e, p in zip(eta, pulse_width_s)
             if math.isfinite(e) and math.isfinite(p) and e > 0.0 and p > 0.0]
    if len(pairs) < 3:
        return {"status": "TREND_NOT_RESOLVED", "allowed_claim": "insufficient canonical cases"}
    rho = _pearson(_ranks([p[0] for p in pairs]), _ranks([p[1] for p in pairs]))
    if rho <= -0.5:
        status = "TREND_CONSISTENT_WITH_KOILE_REFERENCE"
    elif rho >= 0.5:
        status = "TREND_DIFFERS_FROM_KOILE_REFERENCE"
    else:
        status = "TREND_NOT_RESOLVED"
    return {"status": status, "spearman_rho": rho,
            "forbidden_claim": "KOILE_MECHANISM_VALIDATED",
            "scope": "trend direction only for canonical microgap attachment development cases"}


def frequency_proxy(pulse_width_s: float | None) -> dict[str, object]:
    if pulse_width_s is None or not math.isfinite(pulse_width_s) or pulse_width_s <= 0.0:
        return {"status": "NOT_RESOLVED", "f_time_proxy_Hz": None}
    return {"status": "DIAGNOSTIC_PROXY", "f_time_proxy_Hz": 1.0 / pulse_width_s,
            "is_spectral_peak": False}


def f_r6b_gate(*, b_rf1_complete: bool, actual_geometry_descriptors: int,
               stage_c_geometry_mapping_defined: bool) -> bool:
    return bool(b_rf1_complete and actual_geometry_descriptors >= 2 and stage_c_geometry_mapping_defined)


def _ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = 0.5 * (start + end - 1) + 1.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def _pearson(left: list[float], right: list[float]) -> float:
    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    numerator = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - mean_left) ** 2 for a in left) *
                            sum((b - mean_right) ** 2 for b in right))
    return numerator / denominator if denominator > 0.0 else 0.0
