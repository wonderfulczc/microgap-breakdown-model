import csv
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.evidence import (  # noqa: E402
    H3_DEVELOPMENT_BAND_HZ,
    SIGNAL_LEVEL_DEFINITIONS,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    STAGE_I_REQUIRED_INPUT_IDS,
    TRUST_HIERARCHY,
    coherent_summation_permitted,
    common_spectral_metrics,
    exact_frequency_support,
    native_propagation_action,
    validate_signal_levels,
    validate_stage_h_transfer_contract,
    validate_stage_i_requirements,
)


H5 = ROOT / "fullwave/h5"


def contract():
    return json.loads((H5 / "h5_stage_h_transfer_contract.json").read_text())


def test_signal_0_1_2_semantics_are_frozen_and_sources_differ():
    assert validate_signal_levels(SIGNAL_LEVEL_DEFINITIONS)
    assert SIGNAL_LEVEL_DEFINITIONS["SIGNAL_0"] == "SOURCE_LEVEL_PHYSICS"
    assert SIGNAL_LEVEL_DEFINITIONS["SIGNAL_1"] == "PROPAGATION_STRUCTURE_LEVEL"
    assert SIGNAL_LEVEL_DEFINITIONS["SIGNAL_2"] == "RECEIVER_TERMINAL_LEVEL"
    assert SIGNAL_LEVEL_DEFINITIONS["SIGNAL_0_NATIVE"] != SIGNAL_LEVEL_DEFINITIONS["SIGNAL_0_PORT"]


def test_h3_and_h4_pathways_remain_separate():
    with (H5 / "h5_pathway_matrix.csv").open(newline="") as handle:
        rows = {row["pathway"]: row for row in csv.DictReader(handle)}
    assert set(rows) == {"NATIVE_STAGE4", "NATIVE_STAGE5", "POST_BREAKDOWN_G3"}
    assert rows["NATIVE_STAGE4"]["propagation_model"].startswith("STAGE_F")
    assert rows["POST_BREAKDOWN_G3"]["propagation_model"].startswith("H3")
    assert rows["NATIVE_STAGE4"]["receiver_model"] != rows["POST_BREAKDOWN_G3"]["receiver_model"]


def test_current_h3_h4_coherent_summation_is_not_permitted():
    paths = [
        {"physical_tx_rx_geometry": "h3", "receiver_id": "rx350"},
        {"physical_tx_rx_geometry": "h4", "receiver_id": "rxGHz"},
    ]
    assert not coherent_summation_permitted(paths)
    assert contract()["CROSS_PATH_COHERENT_SUMMATION"] == "NOT_PERMITTED_CURRENT_CONFIGURATION"


def test_native_propagation_cannot_be_applied_twice():
    assert (
        native_propagation_action(True, "APPLY_RECEIVER_RESPONSE_ONLY")
        == "NATIVE_PROPAGATION_OWNED_BY_STAGE_F"
    )
    with pytest.raises(ValueError, match="DOUBLE"):
        native_propagation_action(True, "APPLY_FRIIS_PROPAGATION")


def test_350mhz_native_status_remains_not_resolved():
    record = contract()
    assert record["native_350MHz_status"] == "NOT_RESOLVED"
    assert record["frequency_support"]["NATIVE_200_500MHz"] == "NOT_RESOLVED"


def test_stage_f_trust_masks_are_exactly_unchanged():
    assert exact_frequency_support("NATIVE_STAGE4") == STAGE4_TRUSTED_BAND_HZ
    assert exact_frequency_support("NATIVE_STAGE5") == STAGE5_TRUSTED_BAND_HZ
    record = contract()["frequency_support"]
    assert tuple(record["NATIVE_STAGE4_Hz"]) == STAGE4_TRUSTED_BAND_HZ
    assert tuple(record["NATIVE_STAGE5_Hz"]) == STAGE5_TRUSTED_BAND_HZ


def test_h3_practical_band_and_frequency_gap_are_unchanged():
    assert exact_frequency_support("POST_BREAKDOWN_G3") == H3_DEVELOPMENT_BAND_HZ
    record = contract()["frequency_support"]
    assert tuple(record["POST_BREAKDOWN_G3_Hz"]) == (200e6, 500e6)
    assert record["UNSUPPORTED_500MHz_TO_STAGE4_LOW"] == "NO_TRUST_INTERPOLATION"


def test_sparse_h4_bins_make_time_waveform_secondary_only():
    native = contract()["native_paths"]
    assert native["Stage4"]["trusted_bin_count"] == 5
    assert native["Stage5"]["trusted_bin_count"] == 7
    assert all(path["frequency_domain_role"] == "PRIMARY_NATIVE_RF_OUTPUT" for path in native.values())
    assert all(path["time_domain_role"] == "SECONDARY_RECONSTRUCTION_ONLY" for path in native.values())


def test_h3_loading_debts_are_propagated():
    debts = contract()["known_debts"]
    assert "FULL_WAVE_LOADING_MISMATCH_HIGH" in debts
    assert "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED" in debts


def test_h4_mesh_sensitivity_and_absolute_status_are_propagated():
    record = contract()
    assert "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT" in record["known_debts"]
    assert all(
        path["absolute_amplitude_status"] == "NUMERICAL_REFERENCE_ONLY"
        for path in record["native_paths"].values()
    )


def test_stage5_full_maxwell_debt_is_propagated():
    assert "FULL_MAXWELL_REFERENCE_PENDING" in contract()["known_debts"]


def test_vna_and_production_geometry_debts_are_propagated():
    record = contract()
    assert record["H2_SCIENTIFIC_VALIDATION"] == "VNA_MEASUREMENT_PENDING"
    assert record["STAGE_H_PRODUCTION_GEOMETRY_PENDING"] is True
    assert "PRODUCTION_TX_RX_GEOMETRY_PENDING" in record["known_debts"]


def test_common_trust_hierarchy_and_metrics_preserve_units():
    assert TRUST_HIERARCHY == (
        "TRUSTED_PHYSICS",
        "DEVELOPMENT_VERIFIED",
        "NUMERICAL_REFERENCE_ONLY",
        "NOT_RESOLVED",
        "EXPERIMENTAL_VALIDATION_PENDING",
    )
    frequency = np.array([1.0, 2.0, 3.0])
    result = common_spectral_metrics(
        frequency, [1, 2, 1], [1, 1, 3], source_unit="V_s_per_m", receiver_unit="V_s"
    )
    assert result["source_unit"] == "V_s_per_m"
    assert result["receiver_unit"] == "V_s"
    assert result["receiver_induced_peak_shift_Hz"] == 1.0
    assert np.max(result["source_normalized_shape"]) == pytest.approx(1.0)


def test_stage_i_required_input_contract_is_complete_and_unfabricated():
    record = json.loads((H5 / "h5_stage_i_required_inputs.json").read_text())
    assert validate_stage_i_requirements(record)
    assert tuple(item["input_id"] for item in record["required_inputs"]) == STAGE_I_REQUIRED_INPUT_IDS
    assert not any(item["available"] for item in record["required_inputs"])


def test_unified_stage_h_contract_validates_without_scientific_promotion():
    record = contract()
    assert validate_stage_h_transfer_contract(record)
    assert record["STAGE_H_TOOL_DEVELOPMENT"] == "PASS"
    assert record["STAGE_H_SCIENTIFIC_VALIDATION"] == "PENDING_STAGE_I"
    assert record["CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION"] == "NOT_RESOLVED"

