"""Same-event local-kinetics to current-moment diagnostics.

This module is a read-only bridge.  It does not define a new current source and
does not use interpolation to improve event timing claims.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


CONFIRMED = "CONFIRMED"


@dataclass(frozen=True)
class TimebaseAudit:
    status: str
    maximum_absolute_difference_s: float
    matching_tolerance_s: float
    resampling_required_for_rf: bool
    trust_upgrade_from_interpolation: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "maximum_absolute_difference_s": self.maximum_absolute_difference_s,
            "matching_tolerance_s": self.matching_tolerance_s,
            "resampling_required_for_rf": self.resampling_required_for_rf,
            "trust_upgrade_from_interpolation": self.trust_upgrade_from_interpolation,
        }


def audit_native_timebases(
    solver_time_s: np.ndarray,
    diagnostic_time_s: np.ndarray,
    moment_time_s: np.ndarray,
    *,
    matching_tolerance_s: float | None = None,
) -> TimebaseAudit:
    grids = [np.asarray(value, dtype=float) for value in (solver_time_s, diagnostic_time_s, moment_time_s)]
    if any(value.ndim != 1 or value.size < 2 or not np.all(np.isfinite(value)) for value in grids):
        raise ValueError("timebases must be finite one-dimensional arrays with at least two samples")
    if any(np.any(np.diff(value) <= 0.0) for value in grids):
        raise ValueError("native timebases must be strictly increasing")
    reference = grids[0]
    if any(value.shape != reference.shape for value in grids[1:]):
        return TimebaseAudit("NOT_ALIGNED_DIFFERENT_SAMPLE_COUNT", math.inf, 0.0, True)
    dt_min = float(np.min(np.diff(reference)))
    tolerance = matching_tolerance_s if matching_tolerance_s is not None else max(1e-30, dt_min * 1e-9)
    difference = max(float(np.max(np.abs(value - reference))) for value in grids[1:])
    if difference == 0.0:
        status = "TIMEBASE_NATIVE_ALIGNMENT_PASS"
    elif difference <= tolerance:
        status = "TIMEBASE_NEAREST_STEP_ALIGNMENT_PASS"
    else:
        status = "TIMEBASE_NOT_ALIGNED"
    return TimebaseAudit(status, difference, float(tolerance), status == "TIMEBASE_NOT_ALIGNED")


def assess_event_window(
    time_s: np.ndarray,
    topology: np.ndarray,
    topology_status: np.ndarray,
    *,
    event_topology: str,
) -> dict[str, object]:
    time = np.asarray(time_s, dtype=float)
    top = np.asarray(topology, dtype=str)
    status = np.asarray(topology_status, dtype=str)
    if time.ndim != 1 or top.shape != time.shape or status.shape != time.shape or time.size < 3:
        raise ValueError("event-window arrays must be one-dimensional and equally sized")
    selected = np.flatnonzero((top == event_topology) & (status == CONFIRMED))
    if selected.size == 0:
        return {"status": "NO_CONFIRMED_EVENT", "event_topology": event_topology}
    first, last = int(selected[0]), int(selected[-1])
    contiguous = bool(np.all(np.diff(selected) == 1))
    n_pre = first
    n_event = int(selected.size)
    n_post = int(time.size - last - 1)
    complete = contiguous and n_pre >= 1 and n_event >= 3 and n_post >= 1
    return {
        "status": "COMPLETE" if complete else "INCOMPLETE",
        "event_topology": event_topology,
        "event_start_s": float(time[first]),
        "event_peak_s": None,
        "event_end_s": float(time[last]),
        "N_pre": n_pre,
        "N_event": n_event,
        "N_post": n_post,
        "window_duration_s": float(time[last] - time[first]),
        "dt_min_s": float(np.min(np.diff(time[selected]))) if selected.size > 1 else None,
        "dt_median_s": float(np.median(np.diff(time[selected]))) if selected.size > 1 else None,
        "dt_max_s": float(np.max(np.diff(time[selected]))) if selected.size > 1 else None,
        "contiguous": contiguous,
        "requirements": "pre/rise/peak/decay/post must all be present",
    }


def interior_extremum_time(time_s: np.ndarray, values: np.ndarray, *, mode: str) -> float | None:
    time = np.asarray(time_s, dtype=float)
    data = np.asarray(values, dtype=float)
    finite = np.flatnonzero(np.isfinite(data))
    if finite.size < 3:
        return None
    local = data[finite]
    relative = int(np.argmax(local) if mode == "max" else np.argmin(local))
    index = int(finite[relative])
    if index == 0 or index == data.size - 1:
        return None
    if mode == "max" and not (data[index] > data[index - 1] and data[index] >= data[index + 1]):
        return None
    if mode == "min" and not (data[index] < data[index - 1] and data[index] <= data[index + 1]):
        return None
    return float(time[index])


def central_derivative_native(time_s: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time_s, dtype=float)
    signal = np.asarray(values, dtype=float)
    if time.ndim != 1 or signal.shape != time.shape or time.size < 3:
        raise ValueError("native derivative requires equally sized one-dimensional arrays")
    if np.any(np.diff(time) <= 0.0) or not np.all(np.isfinite(time)) or not np.all(np.isfinite(signal)):
        raise ValueError("native derivative input must be finite with increasing time")
    return time, np.gradient(signal, time, edge_order=2)


def compare_proxy_to_derivative(
    time_s: np.ndarray,
    dMdt_A_m_s: np.ndarray,
    K_ion_A_m_s: np.ndarray,
) -> dict[str, object]:
    time = np.asarray(time_s, dtype=float)
    derivative = np.asarray(dMdt_A_m_s, dtype=float)
    proxy = np.asarray(K_ion_A_m_s, dtype=float)
    if time.ndim != 1 or derivative.shape != time.shape or proxy.shape != time.shape or time.size < 3:
        raise ValueError("comparison arrays must be one-dimensional and equally sized")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(derivative)) or not np.all(np.isfinite(proxy)):
        raise ValueError("comparison arrays must be finite")
    x = derivative - np.mean(derivative)
    y = proxy - np.mean(proxy)
    xnorm, ynorm = float(np.linalg.norm(x)), float(np.linalg.norm(y))
    if xnorm == 0.0 or ynorm == 0.0:
        return {"status": "NOT_RESOLVED_CONSTANT_SIGNAL"}
    zero_lag = float(np.dot(x, y) / (xnorm * ynorm))
    correlation = np.correlate(x / xnorm, y / ynorm, mode="full")
    lags = np.arange(-time.size + 1, time.size)
    best = int(np.argmax(correlation))
    dt = float(np.median(np.diff(time)))
    x_scale = derivative / max(float(np.linalg.norm(derivative)), 1e-300)
    y_scale = proxy / max(float(np.linalg.norm(proxy)), 1e-300)
    return {
        "status": "DIAGNOSTIC_ONLY_NO_FREE_SHIFT_CLAIM",
        "zero_lag_correlation": zero_lag,
        "maximum_cross_correlation": float(correlation[best]),
        "lag_at_maximum_correlation_s": float(lags[best] * dt),
        "normalized_L2_without_time_shift": float(np.linalg.norm(x_scale - y_scale)),
        "dMdt_peak_time_s": float(time[np.argmax(np.abs(derivative))]),
        "proxy_peak_time_s": float(time[np.argmax(np.abs(proxy))]),
        "peak_time_difference_s": float(time[np.argmax(np.abs(proxy))] - time[np.argmax(np.abs(derivative))]),
        "interpretation": "K_ion is a mechanism proxy, not a closure relation for dM/dt",
    }


def pulse_resolution_status(pulse_width_s: float | None, diagnostic_trace_dt_s: float) -> dict[str, object]:
    if pulse_width_s is None or not math.isfinite(pulse_width_s) or pulse_width_s <= 0.0:
        return {"status": "NOT_EVALUABLE", "N_PW": None}
    if not math.isfinite(diagnostic_trace_dt_s) or diagnostic_trace_dt_s <= 0.0:
        raise ValueError("diagnostic trace dt must be finite and positive")
    samples = pulse_width_s / diagnostic_trace_dt_s
    if samples >= 10.0:
        status = "PULSE_RESOLVED_PREFERRED"
    elif samples >= 5.0:
        status = "PULSE_RESOLVED_MINIMUM"
    else:
        status = "NOT_RESOLVED_PULSE_TIMESCALE"
    return {"status": status, "N_PW": float(samples), "interpolation_can_upgrade_status": False}


def combined_temporal_trust(
    tau_i_status: str,
    tau_M_status: str,
    pulse_status: str,
    *,
    event_window_status: str,
) -> str:
    accepted_kinetic = {"RESOLVED_PREFERRED", "RESOLVED_MINIMUM"}
    accepted_pulse = {"PULSE_RESOLVED_PREFERRED", "PULSE_RESOLVED_MINIMUM"}
    if (tau_i_status in accepted_kinetic and tau_M_status in accepted_kinetic and
            pulse_status in accepted_pulse and event_window_status == "COMPLETE"):
        return "PASS"
    return "NOT_RESOLVED"
