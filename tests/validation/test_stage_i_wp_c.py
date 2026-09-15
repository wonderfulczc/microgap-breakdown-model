import json
from pathlib import Path

import numpy as np
import pytest

from streamer_rf.validation.native_ghz import (
    NATIVE_DISCREPANCY_CATEGORIES,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE4_TRUSTED_BINS_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BINS_HZ,
    comparison_mask,
    frequency_coverage_status,
    sparse_trusted_bin_metrics,
    trusted_bin_support,
    validate_native_discrepancy_entries,
    validate_native_measurement_contract,
    validate_synthetic_native_fixture,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "validation/stage_i/native_ghz"


def load(name):
    return json.loads((OUTPUT / name).read_text())


def test_stage4_trust_mask_is_exactly_preserved():
    profile = load("native_ghz_stage4_profile.json")
    assert tuple(profile["trusted_frequency_band_Hz"]) == STAGE4_TRUSTED_BAND_HZ
    assert np.array_equal(profile["trusted_frequency_bins_Hz"], STAGE4_TRUSTED_BINS_HZ)
    assert profile["trusted_bin_count"] == 5


def test_stage5_trust_mask_is_exactly_preserved():
    profile = load("native_ghz_stage5_profile.json")
    assert tuple(profile["trusted_frequency_band_Hz"]) == STAGE5_TRUSTED_BAND_HZ
    assert np.array_equal(profile["trusted_frequency_bins_Hz"], STAGE5_TRUSTED_BINS_HZ)
    assert profile["trusted_bin_count"] == 7


def test_measurement_trust_receiver_mask_intersection():
    measurement = np.array([True, True, False, True])
    trust = np.array([False, True, True, True])
    receiver = np.array([True, True, True, False])
    assert np.array_equal(comparison_mask(measurement, trust, receiver), [False, True, False, False])


def test_partial_instrument_bandwidth_handling():
    assert frequency_coverage_status([0.0, 8e9], STAGE4_TRUSTED_BAND_HZ) == "FULL_TRUST_BAND_COVERAGE"
    assert frequency_coverage_status([0.0, 8e9], STAGE5_TRUSTED_BAND_HZ) == "PARTIAL_TRUST_BAND_COVERAGE"


def test_8ghz_scope_cannot_claim_complete_stage5_coverage():
    policy = load("native_ghz_frequency_coverage_policy.json")
    scope = policy["instruments"]["OSCILLOSCOPE_WIDEBAND_WAVEFORM"]
    assert scope["Stage5_coverage"] == "PARTIAL_TRUST_BAND_COVERAGE"
    assert scope["Stage5_trusted_bins_covered"] == 5
    assert np.count_nonzero(trusted_bin_support(STAGE5_TRUSTED_BINS_HZ, [0.0, 8e9], [2.8e9, 10.5e9])) == 5


def test_spectrum_analyzer_capability_covers_stage5():
    analyzer = load("native_ghz_frequency_coverage_policy.json")["instruments"]["SPECTRUM_ANALYZER"]
    assert analyzer["capability_frequency_Hz"] == [2.0, 67e9]
    assert analyzer["Stage4_coverage"] == "FULL_TRUST_BAND_COVERAGE"
    assert analyzer["Stage5_coverage"] == "FULL_TRUST_BAND_COVERAGE"
    assert analyzer["Stage5_trusted_bins_covered"] == 7


def test_sparse_trusted_bin_metrics_do_not_require_dense_interpolation():
    frequency = STAGE4_TRUSTED_BINS_HZ
    simulation = np.array([1, 2, 4, 2, 1], dtype=complex)
    measurement = 3 * simulation * np.exp(0.2j)
    result = sparse_trusted_bin_metrics(frequency, simulation, measurement)
    assert result["trusted_bin_count"] == 5
    assert result["normalized_trusted_bin_correlation"] == pytest.approx(1.0)
    assert result["relative_amplitude_pattern_L2"] == pytest.approx(0.0)


def test_stage4_and_stage5_profiles_are_independent():
    stage4 = load("native_ghz_stage4_profile.json")
    stage5 = load("native_ghz_stage5_profile.json")
    assert stage4["VALIDATION_PROFILE"] == "NATIVE_GHZ_STAGE4"
    assert stage5["VALIDATION_PROFILE"] == "NATIVE_GHZ_STAGE5"
    assert stage4["source_case_id"] != stage5["source_case_id"]
    assert stage4["received_native_spectrum_hash"] != stage5["received_native_spectrum_hash"]


def test_full_maxwell_pending_is_propagated_only_to_stage5():
    assert load("native_ghz_stage4_profile.json")["full_maxwell_status"] is None
    assert load("native_ghz_stage5_profile.json")["full_maxwell_status"] == "FULL_MAXWELL_REFERENCE_PENDING"
    assert load("native_ghz_validation_status.json")["FULL_MAXWELL_REFERENCE_PENDING"] is True


def test_synthetic_native_provenance_protection():
    rows = [{"frequency_Hz": 4e9, "value_real": 1.0, "value_imag": 0.0, "trust_valid": 1, "repetition_index": 1}]
    valid = {
        "data_origin": "SYNTHETIC_NATIVE_GHZ_DRY_RUN",
        "scientific_validation_allowed": False,
        "frequency_unit": "Hz",
        "quantity_unit": "V_s",
        "trust_mask_semantics": "STAGE_F_EXACT_TRUSTED_BINS_ONLY",
    }
    assert validate_synthetic_native_fixture(valid, rows)
    invalid = dict(valid, data_origin="MEASURED")
    with pytest.raises(ValueError, match="PROVENANCE_REQUIRED"):
        validate_synthetic_native_fixture(invalid, rows)
    with pytest.raises(ValueError, match="CANNOT_VALIDATE"):
        validate_synthetic_native_fixture(dict(valid, scientific_validation_allowed=True), rows)
    with pytest.raises(ValueError, match="CANNOT_VALIDATE"):
        validate_synthetic_native_fixture(dict(valid, validation_status="VALIDATED"), rows)


def test_system_350mhz_profile_is_isolated_and_unmodified():
    system = json.loads((ROOT / "validation/stage_i/stage_i_350mhz_profile.json").read_text())
    status = load("native_ghz_validation_status.json")
    assert system["VALIDATION_PROFILE"] == "SYSTEM_350MHZ"
    assert system["analysis_band_Hz"] == [200e6, 500e6]
    assert system["target_center_Hz"] == 350e6
    assert status["SYSTEM_350MHZ_PROFILE_MODIFIED"] is False
    assert status["NATIVE_350MHZ_RF_TRUST"] == "NOT_RESOLVED"


def test_native_discrepancy_ledger_categories_and_paths():
    ledger = load("native_ghz_discrepancy_ledger.json")
    assert set(ledger["allowed_possible_sources"]) == set(NATIVE_DISCREPANCY_CATEGORIES)
    assert ledger["entries"] == []
    sample = [{
        "pathway": "NATIVE_STAGE5",
        "quantity": "TRUSTED_BIN_PATTERN",
        "possible_source": "FULL_MAXWELL_REFERENCE",
        "evidence_type": "FUTURE_MEASUREMENT",
        "action_required": "REVIEW",
    }]
    assert validate_native_discrepancy_entries(sample)


def test_measurement_templates_remain_not_measured():
    contract = load("native_ghz_measurement_contract.json")
    for profile, template in contract["templates"].items():
        assert template["validation_profile"] == profile
        assert template["measurement_status"] == "NOT_MEASURED"
        assert validate_native_measurement_contract(template)


def test_h4_frozen_hashes_and_observer_metadata_are_recorded():
    provenance = load("native_ghz_h4_input_provenance.json")
    assert provenance["receiver_transfer_hash"] == "dbfe52d59367877d5b5d9575bb431b6a1744ccaabe5c276ea7f02bc9c7fa3db6"
    for stage in ("Stage4", "Stage5"):
        source = provenance["stages"][stage]
        assert source["observer_position_m"] == [0.2, 0.0, 0.005]
        assert source["NATIVE_PROPAGATION_ALREADY_INCLUDED"] is True
