from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from streamer_rf.rf.ultrafast import (
    ScalarRingBuffer,
    convergence_metric,
    f_r6_gate,
    field_decay_diagnostic,
    physical_crossing_status,
    select_completed_pulse_window,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "rf/f_r5b_targeted_reference"


def load(name: str) -> dict:
    return json.loads((OUT / name).read_text())


def test_event_completion_requires_pre_peak_decay_and_post() -> None:
    time = np.linspace(-8.0, 8.0, 321) * 1e-12
    moment = 0.5 * (1.0 + np.tanh(time / 1e-12))
    complete = select_completed_pulse_window(time, moment, primary_peak_time_s=0.0)
    assert complete["status"] == "EVENT_WINDOW_COMPLETE"
    assert complete["N_pre"] >= 5 and complete["N_post"] >= 5
    truncated = select_completed_pulse_window(time[150:], moment[150:], primary_peak_time_s=0.0)
    assert truncated["status"] == "INCOMPLETE"


def test_scalar_ring_buffer_retains_only_requested_physical_duration() -> None:
    buffer = ScalarRingBuffer(2e-12)
    for index in range(6):
        buffer.append(index * 1e-12, index)
    rows = buffer.rows()
    assert [value for _, value in rows] == [3, 4, 5]
    assert np.isclose(rows[-1][0] - rows[0][0], 2e-12)


def test_compact_and_event_detail_trace_schemas_are_distinct() -> None:
    compact = pd.read_csv(OUT / "compact_trace.csv", nrows=1)
    detail = pd.read_csv(OUT / "event_trace.csv", nrows=1)
    required_compact = {"time_s", "dt_s", "topology", "E_peak_V_m", "ne_peak_m3", "sigma_e_peak_S_m", "tau_i_at_E_peak_s", "tau_M_at_E_peak_s", "Pi_RF_at_E_peak", "current_moment_z_A_m"}
    assert required_compact <= set(compact.columns)
    assert len(compact.columns) < len(detail.columns)
    assert {"roi_definition", "K_ion_z_A_m_s", "signed_proxy_status"} <= set(detail.columns)


def test_driver_uses_buffered_output_without_per_step_flush() -> None:
    source = (ROOT / "cpp/apps/stage_c_r5_diagnostics.cpp").read_text()
    assert "trace.flush()" not in source
    assert "std::deque<UltrafastEventSample> pre_event_buffer" in source
    assert "PHYSICAL_EVENT_COMPLETION_WITH_MAX_STEPS_SAFETY_CAP" in source


def test_field_decay_requires_interior_peak_and_one_over_e_crossing() -> None:
    time = np.arange(7, dtype=float)
    passed = field_decay_diagnostic(time, np.array([1.0, 2.0, 4.0, 3.0, 1.0, 0.5, 0.2]), event_start_s=0.0, event_end_s=6.0)
    assert passed["status"] == "PASS"
    unresolved = field_decay_diagnostic(time, np.arange(7, dtype=float), event_start_s=0.0, event_end_s=6.0)
    assert unresolved["status"] == "NOT_RESOLVED"


def test_physical_frequency_resolution_is_not_upgraded_by_zero_padding() -> None:
    unresolved = physical_crossing_status(record_duration_s=5e-12, nyquist_Hz=1e14, crossing_Hz=2e11)
    assert unresolved["status"] == "NOT_RESOLVED_REFERENCE_FREQUENCY"
    assert unresolved["display_bin_spacing_is_not_physical_resolution"] is True
    resolved = physical_crossing_status(record_duration_s=5e-9, nyquist_Hz=1e11, crossing_Hz=5e9)
    assert resolved["status"] == "RESOLVED"


def test_dt_convergence_and_f_r6_gate() -> None:
    assert convergence_metric(3.416e-12, 3.415e-12)["status"] == "PASS"
    assert convergence_metric(3.8e-12, 3.0e-12)["status"] == "NOT_RESOLVED"
    assert f_r6_gate(
        event_complete=True, pulse_width_resolved=True, samples_per_pulse=20.0,
        timebase_aligned=True, kinetics_resolved=True, current_moment_consistent=True,
        interpretable_series=True, processing_sensitivity_pass=True, temporal_convergence_pass=True,
    ) is True
    assert f_r6_gate(
        event_complete=True, pulse_width_resolved=True, samples_per_pulse=4.9,
        timebase_aligned=True, kinetics_resolved=True, current_moment_consistent=True,
        interpretable_series=True, processing_sensitivity_pass=True, temporal_convergence_pass=True,
    ) is False


def test_targeted_reference_status_and_scientific_boundaries() -> None:
    status = load("status.json")
    pulse = load("pulse_diagnostics.json")
    spectrum = load("spectrum.json")
    sign = load("sign_convention_audit.json")
    assert status["TARGETED_REFERENCE_STATUS"] == "PASS"
    assert status["EVENT_WINDOW_COMPLETE"] is True
    assert status["INTERPRETABLE_DEVELOPMENT_REFERENCE"] is True
    assert status["F_R6_ALLOWED"] is True
    assert pulse["F_R5A_pulse"]["status"] == "PASS"
    assert pulse["pulse_resolution"]["samples_per_FWHM"] >= 5
    assert spectrum["REFERENCE_CASE_SPECTRAL_TRUST"] == "NUMERIC_DIAGNOSTIC_ONLY"
    assert spectrum["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert sign["SIGN_CONVENTION_AUDIT"] == "SIGN_EXPECTED_OPPOSITE"
    assert sign["sign_flip_applied"] is False
    assert status["MICROGAP_ULTRAFAST_RF_MECHANISM"] == "NOT_VALIDATED"
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"


def test_temporal_convergence_and_baseline_preservation() -> None:
    convergence = load("dt_convergence.json")
    status = load("status.json")
    assert convergence["TEMPORAL_CONVERGENCE"] == "PASS"
    assert all(value["status"] == "PASS" for value in convergence["metrics"].values())
    assert status["STAGE_C_BASELINE_PRESERVED"] is True
    assert status["STAGE_F_BASELINE_PRESERVED"] is True
