import json
from pathlib import Path

import numpy as np
import pytest

from streamer_rf.validation.wp_i_b import (
    ANALYSIS_BANDS_HZ,
    SYNTHETIC_INPUT,
    band_metrics,
    enforce_synthetic_status,
    event_metrics,
    normalized_magnitude_comparison,
    read_scope_waveform,
    scope_quality_gate,
    synthetic_loading_diagnostic,
    validate_vna_input,
    verify_synthetic_bundle,
)


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "validation/stage_i/system_350mhz/synthetic_input"
OUTPUT = ROOT / "validation/stage_i/system_350mhz"


def load_output(name):
    return json.loads((OUTPUT / name).read_text())


def test_synthetic_bundle_provenance_and_manifest_integrity():
    result = verify_synthetic_bundle(INPUT)
    assert result["verified_file_count"] == 40
    assert result["metadata"]["data_origin"] == SYNTHETIC_INPUT
    assert result["metadata"]["scientific_validation_allowed"] is False


@pytest.mark.parametrize("status", ["MEASURED", "VALIDATED", "VALIDATED_WITH_VNA", "SYSTEM_350MHZ_VALIDATED"])
def test_synthetic_input_cannot_satisfy_scientific_status(status):
    with pytest.raises(ValueError, match="CANNOT_VALIDATE_SCIENCE"):
        enforce_synthetic_status(status)


def test_vna_files_reuse_touchstone_parser_and_exact_grid():
    for name in ("TX_S11.s1p", "RX_S11.s1p", "TX_RX_S21.s2p"):
        data = validate_vna_input(INPUT / "vna" / name)
        assert len(data.frequency_Hz) == 901
        assert data.frequency_Hz[0] == 50e6
        assert data.frequency_Hz[-1] == 500e6
        assert np.all(np.diff(data.frequency_Hz) == 0.5e6)
        assert data.Z0_ohm == 50.0


def test_scope_quality_and_frequency_bands():
    metadata = json.loads((INPUT / "measurement_metadata.json").read_text())["oscilloscope"]
    path = INPUT / "scope_raw/0.6m/CURVED_NEEDLE_0.6m_rep01.csv"
    time_s, voltage = read_scope_waveform(path)
    quality = scope_quality_gate(time_s, voltage, metadata, baseline_available=True)
    assert quality["passed"]
    metrics, frequency, spectrum = event_metrics(time_s, voltage, metadata)
    assert metrics["FULL_50_500MHZ"]["frequency_low_Hz"] == pytest.approx(50e6)
    assert metrics["LOW_50_100MHZ"]["frequency_high_Hz"] == pytest.approx(100e6)
    assert metrics["H3_OVERLAP_200_500MHZ"]["frequency_low_Hz"] == pytest.approx(200e6)
    assert np.isfinite(spectrum).all()
    assert frequency[-1] == pytest.approx(2.5e9)


def test_scope_quality_rejects_nonmonotonic_time():
    metadata = json.loads((INPUT / "measurement_metadata.json").read_text())["oscilloscope"]
    time_s = np.arange(10000) * 2e-10
    time_s[20] = time_s[19]
    quality = scope_quality_gate(time_s, np.sin(np.arange(10000)), metadata, baseline_available=True)
    assert not quality["passed"]


def test_band_metric_keeps_56mhz_component_outside_h3_band():
    frequency = np.arange(0.0, 501e6, 0.5e6)
    spectrum = np.zeros(frequency.shape, dtype=complex)
    spectrum[np.argmin(abs(frequency - 56e6))] = 1.0
    low = band_metrics(frequency, spectrum, *ANALYSIS_BANDS_HZ["LOW_50_100MHZ"])
    high = band_metrics(frequency, spectrum, *ANALYSIS_BANDS_HZ["H3_OVERLAP_200_500MHZ"])
    assert low["peak_frequency_Hz"] == 56e6
    assert high["band_integrated_V2_Hz"] == 0.0


def test_normalized_spectral_comparison_is_scale_independent():
    frequency = np.linspace(200e6, 500e6, 31)
    spectrum = np.exp(-((frequency - 350e6) / 50e6) ** 2) * np.exp(1j * frequency / 1e8)
    metrics = normalized_magnitude_comparison(frequency, spectrum, 3.0 * spectrum)
    assert metrics["normalized_shape_correlation"] == pytest.approx(1.0)
    assert metrics["normalized_shape_L2"] == pytest.approx(0.0)


def test_loading_diagnostic_uses_vna_admittance_without_feedback():
    frequency = np.array([200e6, 350e6, 500e6])
    s11 = np.zeros(3, dtype=complex)
    voltage = np.array([1.0, 2.0, 1.0], dtype=complex)
    frozen_current = voltage / 50.0
    impedance, _, implied, metrics = synthetic_loading_diagnostic(
        frequency, s11, 50.0, frequency, voltage, frozen_current
    )
    assert np.allclose(impedance, 50.0)
    assert np.allclose(implied, frozen_current)
    assert metrics["normalized_L2_current_mismatch"] == pytest.approx(0.0)


def test_wp_i_b_outputs_preserve_dry_run_statuses():
    summary = load_output("wp_i_b_summary.json")
    assert summary["WP_I_B_INPUT_TYPE"] == SYNTHETIC_INPUT
    assert summary["SCIENTIFIC_VALIDATION_ALLOWED"] is False
    assert summary["WP_I_B_PIPELINE_DRY_RUN"] == "PASS"
    assert summary["STAGE_I_END_TO_END_DRY_RUN"] == "PASS"
    assert summary["STAGE_I_EXPERIMENTAL_DATA"] == "NOT_PROVIDED"
    assert summary["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"
    assert summary["H2_SCIENTIFIC_VALIDATION"] == "VNA_MEASUREMENT_PENDING"
    assert summary["output_data_layer"] == "DERIVED"
    assert summary["derived_parent_hash"] == summary["raw_input_manifest_hash"]


def test_quality_gate_retains_per_event_results():
    quality = load_output("wp_i_b_quality_gate.json")
    assert quality["status"] == "PASS"
    assert len(quality["event_results"]) == 30
    assert all(event["passed"] for event in quality["event_results"])
    assert quality["reference_plane_status"] == "SYNTHETIC_IDEAL_REFERENCE_PLANE"


def test_h3_comparison_never_extends_below_200mhz():
    summary = load_output("wp_i_b_summary.json")
    assert summary["frequency_policy"]["H3_COMPARISON_BAND_Hz"] == [200e6, 500e6]
    assert summary["frequency_policy"]["50_TO_200MHZ_H3_STATUS"] == "OUTSIDE_H3_COMPARISON_SUPPORT"
    data = np.genfromtxt(OUTPUT / "wp_i_b_h3_comparison.csv", delimiter=",", names=True)
    assert data["frequency_Hz"].min() >= 200e6
    assert data["frequency_Hz"].max() <= 500e6


def test_repetitions_remain_individually_addressable():
    data = np.genfromtxt(OUTPUT / "wp_i_b_event_metrics.csv", delimiter=",", names=True, dtype=None, encoding=None)
    assert len(data) == 30
    for distance in (0.6, 1.0):
        selected = data[np.isclose(data["distance_m"], distance)]
        assert len(selected) == 15
        assert set(selected["repetition_index"]) == set(range(1, 16))
        assert len(set(selected["raw_data_hash"])) == 15


def test_discrepancy_ledger_is_synthetic_only():
    ledger = load_output("wp_i_b_discrepancy_ledger.json")
    assert ledger["scientific_interpretation_allowed"] is False
    assert ledger["entries"]
    assert all(entry["evidence_type"] == "SYNTHETIC_DRY_RUN" for entry in ledger["entries"])


def test_future_real_data_replacement_does_not_require_api_redesign():
    replacement = load_output("wp_i_b_summary.json")["future_real_data_replacement"]
    assert replacement["status"] == "ARCHITECTURE_READY"
    assert replacement["core_parser_redesign_required"] is False
    assert replacement["core_comparison_API_redesign_required"] is False
