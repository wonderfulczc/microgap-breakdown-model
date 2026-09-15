import json
from pathlib import Path

import numpy as np
import pytest

from streamer_rf.validation.wp_i_b import event_metrics, read_scope_waveform, scope_quality_gate
from streamer_rf.validation.wp_i_d import (
    SYNTHETIC_INPUT,
    SYNTHETIC_UNCERTAINTY,
    condition_statistics,
    cosine_floor_fit,
    distance_frequency_stability,
    enforce_synthetic_development_status,
    merge_condition_records,
    monte_carlo_power_law_uncertainty,
    polarization_contrast_db,
    power_law_fit,
    verify_supplement_bundle,
)


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "validation/stage_i/system_350mhz/wp_i_d_synthetic"
OUTPUT = ROOT / "validation/stage_i/wp_i_d"


def load_json(name):
    return json.loads((OUTPUT / name).read_text())


def load_csv(name):
    return np.genfromtxt(OUTPUT / name, delimiter=",", names=True, dtype=None, encoding=None)


def test_supplement_provenance_manifest_and_status_protection():
    verified = verify_supplement_bundle(INPUT)
    assert verified["verified_file_count"] == 64
    assert verified["metadata"]["data_origin"] == SYNTHETIC_INPUT
    assert verified["uncertainty"]["data_origin"] == SYNTHETIC_UNCERTAINTY
    for status in ("MEASURED", "VALIDATED", "EXPERIMENTAL", "VALIDATED_WITH_VNA"):
        with pytest.raises(ValueError, match="CANNOT_VALIDATE"):
            enforce_synthetic_development_status(status)


def test_condition_index_has_unique_events_and_retains_all_repetitions():
    index = load_csv("wp_i_d_condition_index.csv")
    assert len(index) == 90
    identities = {(row["distance_m"], row["orientation_deg"], row["repetition_index"]) for row in index}
    assert len(identities) == 90
    existing = index[np.isclose(index["distance_m"], 0.6) & np.isclose(index["orientation_deg"], 0.0)]
    assert len(existing) == 15
    assert len(set(existing["raw_data_hash"])) == 15


def test_condition_merging_rejects_duplicate_event():
    event = {"distance_m": 0.6, "orientation_deg": 0.0, "repetition_index": 1}
    with pytest.raises(ValueError, match="DUPLICATE"):
        merge_condition_records([event], [event])


def test_new_raw_event_passes_quality_gate_and_three_band_processing():
    metadata = json.loads((INPUT / "supplement_metadata.json").read_text())
    record = metadata["records"][0]
    time_s, voltage = read_scope_waveform(INPUT / record["file"])
    scope = {
        key: metadata[key]
        for key in ("sample_interval_s", "record_length_samples", "record_duration_s", "trigger_time_s")
    }
    quality = scope_quality_gate(time_s, voltage, scope, baseline_available=True)
    metrics, _, _ = event_metrics(time_s, voltage, scope)
    assert quality["passed"]
    assert metrics["FULL_50_500MHZ"]["frequency_low_Hz"] == pytest.approx(50e6)
    assert metrics["LOW_50_100MHZ"]["frequency_high_Hz"] == pytest.approx(100e6)
    assert metrics["H3_OVERLAP_200_500MHZ"]["frequency_low_Hz"] == pytest.approx(200e6)


def test_condition_statistics_and_repeatability_matrix():
    rows = load_csv("wp_i_d_event_metrics.csv")
    stats = condition_statistics(rows, 0.3, 0.0)
    assert stats["event_count"] == 10
    assert stats["peak_voltage_V_mean"] > 0
    matrix = load_csv("wp_i_d_repeatability_matrix.csv")
    assert len(matrix) == 8
    assert set(matrix["event_count"]) == {10, 15}


def test_power_law_fit_recovers_known_exponent():
    distance = np.array([0.3, 0.5, 0.6, 0.8, 1.0])
    amplitude = 0.04 * distance ** -1.25
    fit = power_law_fit(distance, amplitude)
    assert fit["A"] == pytest.approx(0.04)
    assert fit["n"] == pytest.approx(1.25)
    assert fit["R_squared_amplitude_space"] == pytest.approx(1.0)


def test_distance_frequency_stability_uses_explicit_reference():
    result = distance_frequency_stability(
        [0.3, 0.6, 1.0], [351e6, 350e6, 349e6], [345e6, 344e6, 343e6]
    )
    assert result["peak_frequency_shift_Hz"] == [1e6, 0.0, -1e6]
    assert result["maximum_absolute_centroid_shift_Hz"] == 1e6
    assert result["scientific_interpretation_allowed"] is False


def test_cosine_floor_fit_recovers_finite_floor_and_contrast():
    angle = np.array([0.0, 30.0, 60.0, 90.0])
    amplitude = np.sqrt(2.0**2 * np.cos(np.deg2rad(angle)) ** 2 + 0.25**2)
    fit = cosine_floor_fit(angle, amplitude)
    assert fit["A_parallel"] == pytest.approx(2.0, rel=2e-3)
    assert fit["A_floor"] == pytest.approx(0.25, rel=2e-3)
    assert fit["R_squared"] > 0.99999
    assert polarization_contrast_db(amplitude[0], amplitude[-1]) > 18.0


def test_fixed_seed_monte_carlo_is_reproducible():
    assumptions = json.loads((INPUT / "synthetic_uncertainty_assumptions.json").read_text())
    args = ([0.3, 0.5, 0.8, 1.0], [0.12, 0.08, 0.05, 0.04], [0.002] * 4, assumptions)
    first = monte_carlo_power_law_uncertainty(*args, samples=200, seed=17)
    second = monte_carlo_power_law_uncertainty(*args, samples=200, seed=17)
    assert first == second
    assert first["convergence_sanity"]["relative_change"] < 0.25


def test_h3_distance_and_orientation_boundaries_are_not_fabricated():
    distance = load_json("wp_i_d_distance_fit.json")["H3_distance_boundary"]
    orientation = load_json("wp_i_d_orientation_fit.json")
    assert distance["direct_reference_distance_m"] == 1.0
    assert distance["other_distances"] == "NO_DIRECT_H3_FULLWAVE_REFERENCE"
    assert distance["inverse_distance_H3_extrapolation_used"] is False
    assert orientation["H3_orientation_comparison_status"] == "NO_DIRECT_H3_ORIENTATION_REFERENCE"


def test_uncertainty_and_discrepancy_outputs_remain_synthetic_only():
    uncertainty = load_json("wp_i_d_uncertainty_summary.json")
    ledger = load_json("wp_i_d_discrepancy_ledger.json")
    assert uncertainty["data_origin"] == SYNTHETIC_UNCERTAINTY
    assert uncertainty["scientific_use_allowed"] is False
    assert uncertainty["distance_exponent_monte_carlo"]["seed"] == 20260915
    assert len(ledger["entries"]) == 5
    assert all(entry["evidence_type"] == "SYNTHETIC_DRY_RUN" for entry in ledger["entries"])
    assert all(entry["scientific_interpretation_allowed"] is False for entry in ledger["entries"])


def test_wp_i_d_status_and_future_replacement_contract():
    status = load_json("wp_i_d_status.json")
    assert status["WP_I_D_PIPELINE_DRY_RUN"] == "PASS"
    assert status["DISTANCE_ANALYSIS_PIPELINE_READY"] is True
    assert status["ORIENTATION_ANALYSIS_PIPELINE_READY"] is True
    assert status["REPEATABILITY_PIPELINE_READY"] is True
    assert status["UNCERTAINTY_PIPELINE_READY"] is True
    assert status["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"
    assert status["STAGE_I_EXPERIMENTAL_DATA"] == "NOT_PROVIDED"
    replacement = status["future_real_data_replacement"]
    assert replacement["status"] == "ARCHITECTURE_READY"
    assert replacement["core_analysis_code_change_required"] is False
