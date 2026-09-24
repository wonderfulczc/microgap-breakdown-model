"""Targeted development-reference gates for the F-R5B extension."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math

import numpy as np


@dataclass
class ScalarRingBuffer:
    duration_s: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.duration_s) or self.duration_s <= 0.0:
            raise ValueError("ring-buffer duration must be finite and positive")
        self._rows: deque[tuple[float, object]] = deque()

    def append(self, time_s: float, value: object) -> None:
        if not math.isfinite(time_s):
            raise ValueError("ring-buffer time must be finite")
        if self._rows and time_s <= self._rows[-1][0]:
            raise ValueError("ring-buffer time must be strictly increasing")
        self._rows.append((time_s, value))
        while self._rows and time_s - self._rows[0][0] > self.duration_s:
            self._rows.popleft()

    def rows(self) -> tuple[tuple[float, object], ...]:
        return tuple(self._rows)


def select_completed_pulse_window(
    time_s: np.ndarray,
    moment_A_m: np.ndarray,
    *,
    primary_peak_time_s: float,
    tail_fraction: float = 0.1,
    minimum_pre_samples: int = 5,
    minimum_post_samples: int = 5,
) -> dict[str, object]:
    time = np.asarray(time_s, dtype=float)
    moment = np.asarray(moment_A_m, dtype=float)
    if time.ndim != 1 or moment.shape != time.shape or time.size < 7:
        raise ValueError("pulse-window inputs must be equally sized one-dimensional arrays")
    if np.any(np.diff(time) <= 0.0) or not np.all(np.isfinite(time)) or not np.all(np.isfinite(moment)):
        raise ValueError("pulse-window inputs must be finite with increasing time")
    if not (0.0 < tail_fraction < 0.5):
        raise ValueError("tail fraction must lie between zero and one half")
    derivative = np.gradient(moment, time, edge_order=2)
    amplitude = np.abs(derivative)
    peak_index = int(np.argmin(np.abs(time - primary_peak_time_s)))
    if peak_index <= minimum_pre_samples or peak_index >= time.size - minimum_post_samples:
        return {"status": "INCOMPLETE", "reason": "PRIMARY_PEAK_AT_WINDOW_BOUNDARY"}
    peak = float(amplitude[peak_index])
    if not (peak > 0.0 and amplitude[peak_index] >= amplitude[peak_index - 1] and
            amplitude[peak_index] >= amplitude[peak_index + 1]):
        return {"status": "INCOMPLETE", "reason": "PRIMARY_PEAK_NOT_INTERIOR_LOCAL_MAXIMUM"}
    left_tail = np.flatnonzero(amplitude[:peak_index] <= tail_fraction * peak)
    right_tail = np.flatnonzero(amplitude[peak_index + 1:] <= tail_fraction * peak)
    if left_tail.size == 0 or right_tail.size == 0:
        return {"status": "INCOMPLETE", "reason": "TAIL_THRESHOLD_NOT_CROSSED_ON_BOTH_SIDES"}
    start = int(left_tail[-1])
    end = int(peak_index + 1 + right_tail[0])
    left_half = bool(np.any(amplitude[start:peak_index] <= 0.5 * peak))
    right_half = bool(np.any(amplitude[peak_index + 1:end + 1] <= 0.5 * peak))
    left_one_e = bool(np.any(amplitude[start:peak_index] <= peak / math.e))
    right_one_e = bool(np.any(amplitude[peak_index + 1:end + 1] <= peak / math.e))
    complete = (peak_index - start >= minimum_pre_samples and end - peak_index >= minimum_post_samples and
                left_half and right_half and left_one_e and right_one_e)
    return {
        "status": "EVENT_WINDOW_COMPLETE" if complete else "INCOMPLETE",
        "start_index": start,
        "peak_index": peak_index,
        "end_index": end,
        "event_start_s": float(time[start]),
        "event_peak_s": float(time[peak_index]),
        "event_end_s": float(time[end]),
        "N_pre": int(peak_index - start),
        "N_post": int(end - peak_index),
        "peak_abs_dMdt_A_m_s": peak,
        "tail_fraction": tail_fraction,
        "left_half_crossing_available": left_half,
        "right_half_crossing_available": right_half,
        "left_one_over_e_crossing_available": left_one_e,
        "right_one_over_e_crossing_available": right_one_e,
    }


def field_decay_diagnostic(
    time_s: np.ndarray,
    field_V_m: np.ndarray,
    *,
    event_start_s: float,
    event_end_s: float,
) -> dict[str, object]:
    time = np.asarray(time_s, dtype=float)
    field = np.asarray(field_V_m, dtype=float)
    mask = (time >= event_start_s) & (time <= event_end_s) & np.isfinite(field)
    indices = np.flatnonzero(mask)
    if indices.size < 3:
        return {"status": "NOT_RESOLVED", "reason": "INSUFFICIENT_EVENT_SAMPLES"}
    local = field[indices]
    relative_peak = int(np.argmax(local))
    if relative_peak == 0 or relative_peak == local.size - 1:
        return {"status": "NOT_RESOLVED", "reason": "NO_INTERIOR_FIELD_PEAK"}
    peak_index = int(indices[relative_peak])
    threshold = field[peak_index] / math.e
    after = indices[relative_peak + 1:]
    crossings = after[field[after] <= threshold]
    if crossings.size == 0:
        return {
            "status": "NOT_RESOLVED",
            "reason": "FIELD_DOES_NOT_REACH_ONE_OVER_E_WITHIN_EVENT",
            "E_peak_V_m": float(field[peak_index]),
            "t_E_peak_s": float(time[peak_index]),
            "definition": "time from interior E peak to first E_peak/e crossing",
        }
    crossing = int(crossings[0])
    return {
        "status": "PASS",
        "E_peak_V_m": float(field[peak_index]),
        "t_E_peak_s": float(time[peak_index]),
        "t_E_one_over_e_s": float(time[crossing]),
        "field_decay_time_s": float(time[crossing] - time[peak_index]),
        "definition": "time from interior E peak to first E_peak/e crossing",
    }


def physical_crossing_status(
    *,
    record_duration_s: float,
    nyquist_Hz: float,
    crossing_Hz: float | None,
    reference_frequency_Hz: float = 1e9,
    minimum_resolution_bins: float = 3.0,
) -> dict[str, object]:
    if record_duration_s <= 0.0:
        raise ValueError("record duration must be positive")
    resolution = 1.0 / record_duration_s
    if crossing_Hz is None:
        status = "CROSSING_NOT_FOUND"
    elif reference_frequency_Hz < resolution:
        status = "NOT_RESOLVED_REFERENCE_FREQUENCY"
    elif crossing_Hz > nyquist_Hz:
        status = "OUT_OF_BAND"
    elif crossing_Hz / resolution < minimum_resolution_bins:
        status = "NOT_RESOLVED_RECORD_DURATION"
    else:
        status = "RESOLVED"
    return {
        "status": status,
        "physical_frequency_resolution_Hz": resolution,
        "display_bin_spacing_is_not_physical_resolution": True,
        "minimum_resolution_bins": minimum_resolution_bins,
    }


def convergence_metric(reference: float | None, finer: float | None, *, tolerance: float = 0.05) -> dict[str, object]:
    if reference is None or finer is None or not math.isfinite(reference) or not math.isfinite(finer):
        return {"status": "NOT_EVALUABLE", "relative_change": None, "tolerance": tolerance}
    change = abs(reference - finer) / max(abs(finer), 1e-300)
    return {"status": "PASS" if change < tolerance else "NOT_RESOLVED", "relative_change": change, "tolerance": tolerance}


def f_r6_gate(
    *, event_complete: bool, pulse_width_resolved: bool, samples_per_pulse: float | None,
    timebase_aligned: bool, kinetics_resolved: bool, current_moment_consistent: bool,
    interpretable_series: bool, processing_sensitivity_pass: bool, temporal_convergence_pass: bool,
) -> bool:
    return bool(
        event_complete and pulse_width_resolved and samples_per_pulse is not None and samples_per_pulse >= 5.0 and
        timebase_aligned and kinetics_resolved and current_moment_consistent and interpretable_series and
        processing_sensitivity_pass and temporal_convergence_pass
    )
