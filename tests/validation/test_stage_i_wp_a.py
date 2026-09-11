import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.receiver import read_touchstone  # noqa: E402
from streamer_rf.validation.stage_i import (  # noqa: E402
    NOT_PROVIDED,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    amplitude_comparison_class,
    assert_raw_not_overwritten,
    compare_spectra,
    compare_waveforms,
    frequency_to_hz,
    reference_plane_compatible,
    summarize_repetitions,
    validate_data_layer,
    validate_discrepancy_ledger,
    validate_measurement_contract,
    validate_repetition_records,
    validate_scope_contract,
    validate_status_transition,
    validate_vna_contract,
    validate_waveform_metadata,
)


OUT = ROOT / "validation/stage_i"


def load(name):
    return json.loads((OUT / name).read_text())


def test_measurement_contract_accepts_explicit_missing_values():
    record = load("stage_i_measurement_contract.json")
    assert validate_measurement_contract(record)
    assert record["raw_data_path"] == NOT_PROVIDED
    assert record["raw_data_hash"] == NOT_PROVIDED


def test_vna_contract_reuses_h2_touchstone_parser_without_raw_smoothing():
    record = load("stage_i_vna_contract.json")
    assert validate_vna_contract(record)
    assert record["expected_files"] == ["RX_S11.s1p", "TX_S11.s1p", "TX_RX_S21.s2p"]
    assert record["raw_smoothing"] is False


def test_oscilloscope_contract_preserves_raw_voltage():
    record = load("stage_i_scope_contract.json")
    assert validate_scope_contract(record)
    assert record["raw_waveform_columns"] == ["time_s", "voltage_V", "channel_id"]
    assert record["primary_values_normalized"] is False


def test_raw_calibrated_derived_layer_separation():
    assert validate_data_layer({"layer": "RAW", "immutable": True, "parent_hash": None, "operations": []})
    parent = "a" * 64
    assert validate_data_layer({"layer": "CALIBRATED", "parent_hash": parent, "operations": ["SOLT"]})
    assert validate_data_layer({"layer": "DERIVED", "parent_hash": parent, "operations": ["FFT"]})
    with pytest.raises(ValueError, match="PARENT_HASH"):
        validate_data_layer({"layer": "DERIVED", "parent_hash": NOT_PROVIDED, "operations": ["FFT"]})


def test_raw_files_cannot_be_overwritten():
    assert assert_raw_not_overwritten("a" * 64, "a" * 64)
    with pytest.raises(ValueError, match="OVERWRITE"):
        assert_raw_not_overwritten("a" * 64, "b" * 64)


def test_reference_plane_compatibility_is_explicit():
    assert reference_plane_compatible("RX_FEED", "RX_FEED")
    assert not reference_plane_compatible("RX_FEED", "INSTRUMENT_INPUT")
    with pytest.raises(ValueError, match="UNKNOWN"):
        reference_plane_compatible("RX_FEED", "UNKNOWN")


def test_incompatible_reference_planes_reject_absolute_comparison():
    simulation = {"reference_plane": "RX_FEED"}
    measurement = {"reference_plane": "INSTRUMENT_INPUT"}
    assert amplitude_comparison_class(simulation, measurement) == "NOT_COMPARABLE"


def test_absolute_comparison_requires_matching_load_geometry_and_calibration():
    common = {
        "reference_plane": "RX_FEED",
        "receiver_id": "RX1",
        "load_impedance_ohm": 50.0,
        "geometry_id": "G1",
        "instrument_impedance_ohm": 50.0,
    }
    measurement = dict(common, cable_status="CALIBRATED_OUT", calibration_traceable=True)
    assert amplitude_comparison_class(common, measurement) == "ABSOLUTE_AMPLITUDE_VALID"
    measurement["geometry_id"] = "G2"
    assert amplitude_comparison_class(common, measurement) == "NORMALIZED_SHAPE_ONLY"


def test_existing_touchstone_parser_handles_wp_a_dry_run_files():
    root = ROOT / "fullwave/h2/theory_inputs"
    for name in ("H2_THEORY_350MHz_RX_S11.s1p", "H2_THEORY_350MHz_TX_S11.s1p", "H2_THEORY_350MHz_TX_RX_S21.s2p"):
        data = read_touchstone(root / name)
        assert data.frequency_Hz.size == 601
        assert (data.frequency_Hz[0], data.frequency_Hz[-1]) == pytest.approx((200e6, 500e6))
        assert data.Z0_ohm == 50.0


def test_waveform_metadata_validation_and_repetition_index():
    record = {
        "sample_interval_s": 1.25e-11,
        "sample_rate_Hz": 80e9,
        "input_impedance_ohm": 50.0,
        "record_length": 1000,
        "pretrigger_fraction": 0.2,
        "experiment_repetition_index": 3,
    }
    assert validate_waveform_metadata(record)
    record["sample_rate_Hz"] = 40e9
    with pytest.raises(ValueError, match="MISMATCH"):
        validate_waveform_metadata(record)


def test_frequency_units_are_unambiguous():
    assert np.array_equal(frequency_to_hz([200, 350, 500], "MHz"), [200e6, 350e6, 500e6])
    assert np.array_equal(frequency_to_hz([3, 10], "GHz"), [3e9, 10e9])
    with pytest.raises(ValueError, match="UNIT"):
        frequency_to_hz([350], "M")


def test_repeated_measurements_keep_unique_raw_events():
    records = [
        {"repetition_index": 0, "raw_data_hash": "a" * 64},
        {"repetition_index": 1, "raw_data_hash": "b" * 64},
    ]
    assert validate_repetition_records(records)
    with pytest.raises(ValueError, match="UNIQUE"):
        validate_repetition_records([records[0], records[0]])


def test_uncertainty_placeholders_and_repeat_statistics():
    record = load("stage_i_uncertainty_contract.json")
    assert all(value == NOT_PROVIDED for value in record["components"].values())
    result = summarize_repetitions([1.0, 2.0, 3.0])
    assert result["mean"] == pytest.approx(2.0)
    assert result["standard_deviation"] == pytest.approx(1.0)
    assert len(result["confidence_interval"]) == 2


def test_validation_status_cannot_jump_from_data_available_to_validated():
    assert validate_status_transition("NOT_MEASURED", "DATA_AVAILABLE")
    assert validate_status_transition("COMPARISON_READY", "VALIDATED")
    with pytest.raises(ValueError, match="TRANSITION"):
        validate_status_transition("DATA_AVAILABLE", "VALIDATED")


def test_h3_h4_mapping_profiles_are_separate():
    rows = np.genfromtxt(OUT / "stage_i_simulation_measurement_mapping.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    mapping = {row["measured_quantity"]: row for row in rows}
    assert mapping["RECEIVED_SPECTRUM_200_500MHZ"]["validation_profile"] == "SYSTEM_350MHZ"
    assert mapping["MEASURED_GHZ_NATIVE_SPECTRUM"]["validation_profile"] == "NATIVE_GHZ"
    assert mapping["RECEIVED_SPECTRUM_200_500MHZ"]["simulation_quantity"].startswith("H3")
    assert mapping["MEASURED_GHZ_NATIVE_SPECTRUM"]["simulation_quantity"].startswith("H4")


def test_stage_f_native_trust_masks_are_preserved():
    record = load("stage_i_native_ghz_profile.json")
    assert tuple(record["Stage4_trusted_frequency_mask_Hz"]) == STAGE4_TRUSTED_BAND_HZ
    assert tuple(record["Stage5_trusted_frequency_mask_Hz"]) == STAGE5_TRUSTED_BAND_HZ
    assert record["Stage5_status"] == "FULL_MAXWELL_REFERENCE_PENDING"


def test_350mhz_native_path_remains_not_resolved():
    system = load("stage_i_350mhz_profile.json")
    native = load("stage_i_native_ghz_profile.json")
    assert system["native_350MHz_status"] == "NOT_RESOLVED"
    assert native["native_350MHz_status"] == "NOT_RESOLVED"


def test_discrepancy_ledger_schema_is_empty_and_neutral():
    ledger = load("stage_i_discrepancy_ledger.json")
    assert validate_discrepancy_ledger(ledger)
    assert ledger["entries"] == []
    assert ledger["automatic_blame_assignment"] is False


def test_metric_apis_have_no_embedded_acceptance_thresholds():
    frequency = np.array([1.0, 2.0, 3.0])
    spectrum = compare_spectra(frequency, [1, 2, 1], [1, 1.5, 1])
    waveform = compare_waveforms(frequency, [0, 1, 0], [0, 0.8, 0])
    assert "status" not in spectrum and "status" not in waveform
    assert np.isfinite(spectrum["normalized_spectral_correlation"])
    assert waveform["normalized_RMSE"] == pytest.approx(0.25)

