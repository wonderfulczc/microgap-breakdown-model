from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from streamer_rf.rf.ultrafast import (
    align_reported_derivative_peak,
    canonical_case_matrix,
    eta_field_diagnostics,
    f_r6b_gate,
    frequency_proxy,
    koile_trend_wording,
    propagate_event_gate,
    representative_finer_case_ids,
    rise_time_confounding_status,
    same_anchor_pi_rf,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "rf/f_r6a"


def load(name: str) -> dict:
    return json.loads((OUT / name).read_text())


def test_case_matrix_is_bounded_hierarchical_and_provenanced() -> None:
    cases = canonical_case_matrix()
    assert len(cases) == 15
    assert {case["design_level"] for case in cases} == {"LEVEL_0", "LEVEL_1", "LEVEL_2"}
    assert sum(case["scan_family"] == "BASELINE_REPEATABILITY" for case in cases) == 3
    assert {case["voltage_V"] for case in cases} >= {400.0, 500.0, 600.0}
    gaps = np.array([case["gas_gap_m"] for case in cases])
    rises = np.array([case["rise_time_s"] for case in cases])
    assert all(np.any(np.isclose(gaps, value, rtol=0.0, atol=1e-18)) for value in (60e-6, 70e-6, 80e-6))
    assert all(np.any(np.isclose(rises, value, rtol=0.0, atol=1e-24)) for value in (2.5e-12, 5e-12, 10e-12))
    assert all(case["gas"] == "air" and case["pressure_Pa"] == 101325.0 for case in cases)
    assert len(representative_finer_case_ids()) == 7


def test_reported_backward_difference_peak_aligns_without_interpolation() -> None:
    time = np.linspace(-5.0, 5.0, 101) * 1e-12
    moment = np.tanh(time / 1e-12)
    result = align_reported_derivative_peak(time, moment, time[51])
    assert result["aligned_index"] == 50
    assert result["index_offset"] == -1
    assert result["interpolation_used"] is False


def test_unresolved_event_and_pulse_states_propagate() -> None:
    assert propagate_event_gate(event_complete=False, pulse_status="PASS")["status"] == "NOT_RESOLVED_EVENT_WINDOW"
    assert propagate_event_gate(event_complete=True, pulse_status="NOT_RESOLVED_MULTI_PEAK")["status"] == "NOT_RESOLVED_PULSE"
    assert propagate_event_gate(event_complete=True, pulse_status="PASS")["PW_trusted"] is True


def test_eta_definition_and_same_anchor_pi_rule() -> None:
    field = eta_field_diagnostics(eta_0_p95=2.0, eta_peak=3.0)
    assert field["eta_0_representative_statistic"] == "EVENT_ROI_E0_P95"
    assert field["field_amplification_dynamic"] == 1.5
    assert np.isclose(same_anchor_pi_rf(tau_i_s=2e-12, tau_M_s=4e-9, anchors_match=True)["Pi_RF"], 2000.0)
    assert same_anchor_pi_rf(tau_i_s=2e-12, tau_M_s=4e-9, anchors_match=False)["status"] == "NOT_RESOLVED_ANCHOR_MISMATCH"


def test_rise_time_confounding_and_frequency_proxy_labels() -> None:
    driven = rise_time_confounding_status([2.5, 5.0, 10.0], [1.0, 2.0, 4.0])
    assert driven["status"] == "EXTERNAL_DRIVE_SENSITIVE"
    stable = rise_time_confounding_status([2.5, 5.0, 10.0], [1.0, 1.01, 1.02])
    assert stable["status"] == "INTRINSIC_TIMESCALE_EVIDENCE_STRENGTHENED"
    proxy = frequency_proxy(2e-12)
    assert proxy["status"] == "DIAGNOSTIC_PROXY"
    assert proxy["f_time_proxy_Hz"] == 5e11
    assert proxy["is_spectral_peak"] is False


def test_koile_wording_and_f_r6b_gate_boundaries() -> None:
    trend = koile_trend_wording([1.0, 2.0, 3.0], [3.0, 2.0, 1.0])
    assert trend["status"] == "TREND_CONSISTENT_WITH_KOILE_REFERENCE"
    assert trend["forbidden_claim"] == "KOILE_MECHANISM_VALIDATED"
    assert f_r6b_gate(b_rf1_complete=False, actual_geometry_descriptors=2, stage_c_geometry_mapping_defined=True) is False
    assert f_r6b_gate(b_rf1_complete=True, actual_geometry_descriptors=1, stage_c_geometry_mapping_defined=True) is False
    assert f_r6b_gate(b_rf1_complete=True, actual_geometry_descriptors=2, stage_c_geometry_mapping_defined=True) is True


def test_generated_map_retains_unresolved_cases_and_conservative_topology() -> None:
    frame = pd.read_csv(OUT / "case_results.csv")
    status = load("status.json")
    assert len(frame) == 15
    assert status["resolved_cases"] + status["unresolved_cases"] == 15
    assert status["resolved_fraction"] >= 0.70
    assert (frame.case_status != "PASS").any()
    assert "STREAMER_COLLISION" not in set(frame.topology.dropna())
    assert set(frame.frequency_proxy_role.dropna()) == {"DIAGNOSTIC_PROXY"}
    assert set(frame.spectrum_status.dropna()) == {"NUMERIC_DIAGNOSTIC_ONLY"}


def test_rise_audit_b_rf1_gate_and_scientific_boundaries() -> None:
    summary = load("parameter_response_summary.json")
    requirements = load("b_rf1_requirements.json")
    boundaries = load("contracts/f_r6a_scientific_boundaries.json")
    assert summary["rise_time_confounding"]["status"] in {
        "EXTERNAL_DRIVE_SENSITIVE", "INTRINSIC_TIMESCALE_EVIDENCE_STRENGTHENED"
    }
    assert requirements["F_R6B_ALLOWED"] is False
    assert boundaries["KOILE_MECHANISM_VALIDATED"] is False
    assert boundaries["ACTUAL_ELECTRODE_RF_REGIME_VALIDATED"] is False
    assert boundaries["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"


def test_frozen_stage_c_hashes_and_driver_scope() -> None:
    baseline = load("../c_r5/contracts/c_r5_stage_c_frozen_baseline.json")
    for name, expected in baseline["sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    source = (ROOT / "cpp/apps/stage_c_r5_diagnostics.cpp").read_text()
    assert "--gas-gap" in source
    assert "gas_gap_cells" in source
    assert "electron drift-diffusion equation" not in source
