from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from streamer_rf.rf.ultrafast import (
    assess_event_window,
    audit_native_timebases,
    central_derivative_native,
    combined_temporal_trust,
    compare_proxy_to_derivative,
    interior_extremum_time,
    process_current_moment,
    pulse_resolution_status,
)


ROOT = Path(__file__).resolve().parents[2]


def test_native_timebase_alignment_and_no_interpolation_upgrade() -> None:
    time = np.array([0.0, 1e-12, 2e-12, 4e-12])
    exact = audit_native_timebases(time, time.copy(), time.copy())
    assert exact.status == "TIMEBASE_NATIVE_ALIGNMENT_PASS"
    assert exact.resampling_required_for_rf is False
    assert exact.trust_upgrade_from_interpolation is False

    shifted = time.copy()
    shifted[2] += 0.2e-12
    mismatch = audit_native_timebases(time, shifted, time)
    assert mismatch.status == "TIMEBASE_NOT_ALIGNED"
    assert mismatch.resampling_required_for_rf is True
    assert mismatch.trust_upgrade_from_interpolation is False


def test_event_window_requires_pre_event_and_post_event() -> None:
    time = np.arange(7, dtype=float) * 1e-12
    topology = np.array(["UNRESOLVED", "PROPAGATION", "PROPAGATION", "PROPAGATION", "UNRESOLVED", "UNRESOLVED", "UNRESOLVED"])
    status = np.array(["UNRESOLVED", "CONFIRMED", "CONFIRMED", "CONFIRMED", "UNRESOLVED", "UNRESOLVED", "UNRESOLVED"])
    complete = assess_event_window(time, topology, status, event_topology="PROPAGATION")
    assert complete["status"] == "COMPLETE"
    assert complete["N_pre"] == 1 and complete["N_event"] == 3 and complete["N_post"] == 3

    incomplete = assess_event_window(time, np.full(7, "PROPAGATION"), np.full(7, "CONFIRMED"), event_topology="PROPAGATION")
    assert incomplete["status"] == "INCOMPLETE"


def test_anchor_separation_and_boundary_extrema_remain_unresolved() -> None:
    time = np.arange(5, dtype=float)
    assert interior_extremum_time(time, np.array([0.0, 1.0, 4.0, 2.0, 1.0]), mode="max") == 2.0
    assert interior_extremum_time(time, np.array([5.0, 4.0, 3.0, 2.0, 1.0]), mode="max") is None
    e = np.array([1.0, 5.0, 2.0])
    ne = np.array([1.0, 2.0, 8.0])
    sigma = np.array([7.0, 1.0, 2.0])
    assert (int(np.argmax(e)), int(np.argmax(ne)), int(np.argmax(sigma))) == (1, 2, 0)


def test_native_derivative_and_known_lag() -> None:
    time = np.arange(9, dtype=float) * 1e-12
    _, derivative = central_derivative_native(time, time**2)
    assert np.allclose(derivative, 2.0 * time, rtol=1e-12, atol=1e-30)

    reference = np.array([0.0, 0.0, 1.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    delayed = np.roll(reference, 2)
    comparison = compare_proxy_to_derivative(time, reference, delayed)
    assert comparison["lag_at_maximum_correlation_s"] == pytest.approx(-2e-12)
    assert comparison["peak_time_difference_s"] == pytest.approx(2e-12)
    assert comparison["status"] == "DIAGNOSTIC_ONLY_NO_FREE_SHIFT_CLAIM"


def test_ionization_proxy_value_unit_and_sign_convention() -> None:
    elementary_charge = 1.602176634e-19
    ne, nu_i, mu_e, ez, volume = 2e17, 3e10, 0.04, 5e6, 2e-18
    electron_drift_z = -mu_e * ez
    proxy = elementary_charge * electron_drift_z * nu_i * ne * volume
    assert proxy == pytest.approx(-3.8452239216e-4)
    assert abs(proxy) == pytest.approx(3.8452239216e-4)


def test_f_r5a_is_reused_and_pulse_resolution_uses_native_dt() -> None:
    time = np.linspace(-40e-12, 40e-12, 801)
    moment = 0.5 * (1.0 + np.tanh(time / 5e-12))
    result = process_current_moment(time, moment, profile="NO_SMOOTHING")
    assert result.processing_provenance["derivative_scheme"] == "PROJECT_OPERATIONAL_CHOICE"
    resolution = pulse_resolution_status(result.pulse.PW_FWHM_s, float(np.median(np.diff(time))))
    assert resolution["status"] == "PULSE_RESOLVED_PREFERRED"
    assert resolution["interpolation_can_upgrade_status"] is False


@pytest.mark.parametrize(
    ("width", "expected"),
    [(10e-12, "PULSE_RESOLVED_PREFERRED"), (7e-12, "PULSE_RESOLVED_MINIMUM"), (4e-12, "NOT_RESOLVED_PULSE_TIMESCALE")],
)
def test_pulse_resolution_thresholds(width: float, expected: str) -> None:
    assert pulse_resolution_status(width, 1e-12)["status"] == expected
    assert pulse_resolution_status(None, 1e-12)["status"] == "NOT_EVALUABLE"


def test_combined_temporal_trust_requires_complete_window_and_all_timescales() -> None:
    assert combined_temporal_trust(
        "RESOLVED_MINIMUM", "RESOLVED_PREFERRED", "PULSE_RESOLVED_MINIMUM", event_window_status="COMPLETE"
    ) == "PASS"
    assert combined_temporal_trust(
        "RESOLVED_MINIMUM", "RESOLVED_PREFERRED", "PULSE_RESOLVED_MINIMUM", event_window_status="INCOMPLETE"
    ) == "NOT_RESOLVED"


def test_contracts_preserve_stage_and_scientific_boundaries() -> None:
    bridge = json.loads((ROOT / "rf/f_r5b/contracts/same_event_bridge_contract.json").read_text())
    proxy = json.loads((ROOT / "rf/f_r5b/contracts/mechanism_proxy_contract.json").read_text())
    status = json.loads((ROOT / "rf/f_r5b/development/status.json").read_text())
    assert bridge["stage"] == "F"
    assert bridge["current_moment_definition"] == "M(t)=integral(J_RF dV)"
    assert bridge["f_r5a_reuse"] is True
    assert proxy["role"] == "MECHANISM_DIAGNOSTIC_PROXY"
    assert proxy["closure_relation"] is False
    assert status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert status["STAGE_F_BASELINE_PRESERVED"] is True
    assert status["KOILE_MECHANISM_VALIDATED"] is False
    assert status["F_R6_ALLOWED"] is False
