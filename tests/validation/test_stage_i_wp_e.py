import csv
import json
from pathlib import Path

import pytest

from streamer_rf.validation.final_audit import (
    REENTRY_STEPS,
    determine_stage_i_status,
    enforce_scientific_validation_evidence,
    stage_j_packaging_allowed,
    validate_debt_ledger,
    validate_export_inventory,
    validate_reentry_contract,
    validate_synthetic_discrepancy_separation,
    validate_validation_matrix,
    verify_real_data_replacement,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "validation/stage_i/final"


def load_json(name):
    return json.loads((OUTPUT / name).read_text())


def load_matrix():
    with (OUTPUT / "stage_i_validation_matrix.csv").open() as stream:
        return list(csv.DictReader(stream))


def test_final_validation_matrix_schema_and_required_rows():
    rows = load_matrix()
    assert validate_validation_matrix(rows)
    quantities = {row["quantity"] for row in rows}
    assert len(rows) == 16
    assert {
        "350MHZ_TX_S11",
        "DISTANCE_DEPENDENCE",
        "NATIVE_STAGE5_GHZ_SPECTRUM",
        "CROSS_PATH_ABSOLUTE_CONTRIBUTION",
    } <= quantities


def test_stage_i_dual_axis_status_requires_all_development_pipelines():
    passing = {name: True for name in ("WP_I_A", "WP_I_B", "WP_I_C", "WP_I_D")}
    status = determine_stage_i_status(passing, False)
    assert status["STAGE_I_TOOL_DEVELOPMENT"] == "PASS"
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    passing["WP_I_D"] = False
    assert determine_stage_i_status(passing, False)["STAGE_I_TOOL_DEVELOPMENT"] == "NOT_RESOLVED"


def test_scientific_validation_cannot_pass_without_real_data():
    for status in ("PASS", "VALIDATED"):
        with pytest.raises(ValueError, match="REAL_EXPERIMENT_REQUIRED"):
            enforce_scientific_validation_evidence(status, False)


def test_real_data_replacement_readiness_preserves_core_apis():
    contract = load_json("stage_i_final_contract.json")
    assert contract["REAL_DATA_REPLACEMENT_READY"] is True
    replacement = {
        "replace_only": [
            "RAW_MEASUREMENT_FILES", "ACTUAL_GEOMETRY_METADATA", "INSTRUMENT_METADATA",
            "CALIBRATION_METADATA", "UNCERTAINTY_METADATA", "PROVENANCE_STATUS",
        ],
        "redesign_not_required": [
            "PARSER_APIS", "COMPARISON_APIS", "SPECTRAL_METRICS",
            "UNCERTAINTY_ENGINE", "DISCREPANCY_LEDGER", "STAGE_H_CONTRACTS",
        ],
    }
    assert verify_real_data_replacement(replacement)


def test_all_required_scientific_debts_are_propagated():
    ledger = load_json("stage_i_final_debt_ledger.json")
    assert validate_debt_ledger(ledger["debts"])
    identifiers = {entry["debt_id"] for entry in ledger["debts"]}
    assert {
        "H2_VNA_MEASUREMENT_PENDING",
        "H3_FULL_WAVE_LOADING_MISMATCH_HIGH",
        "H3_FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED",
        "H4_RECEIVER_TRANSFER_MESH_SENSITIVITY",
        "H4_ABSOLUTE_AMPLITUDE_NUMERICAL_REFERENCE_ONLY",
        "STAGE5_FULL_MAXWELL_REFERENCE_PENDING",
        "STAGE_F_SPARSE_TRUSTED_FFT_BINS",
        "STAGE_I_EXPERIMENTAL_DATA_NOT_PROVIDED",
        "PRODUCTION_TX_RX_GEOMETRY_NOT_PROVIDED",
    } == identifiers
    assert ledger["blocks_tool_packaging"] is False


def test_stage_j_permission_is_independent_of_scientific_validation():
    assert stage_j_packaging_allowed("PASS", "PASS") is True
    status = load_json("stage_i_final_status.json")
    assert status["STAGE_J_TOOL_PACKAGING_ALLOWED"] is True
    assert status["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False


def test_frequency_and_cross_path_boundaries_are_frozen():
    contract = load_json("stage_i_final_contract.json")
    bands = contract["frequency_boundaries"]
    assert bands["system_experimental_context_Hz"] == [50e6, 500e6]
    assert bands["H3_formal_comparison_Hz"] == [200e6, 500e6]
    assert bands["system_target_center_Hz"] == 350e6
    assert bands["native_350MHz"] == "NOT_RESOLVED"
    assert bands["Stage4_trusted_Hz"] == [2941408508.9091916, 7966314711.629051]
    assert bands["Stage5_trusted_Hz"] == [3047273105.1868486, 10233758844.919172]
    assert contract["CROSS_PATH_COHERENT_SUMMATION"] == "NOT_PERMITTED_CURRENT_CONFIGURATION"
    assert contract["CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION"] == "NOT_RESOLVED"


def test_stage5_full_maxwell_debt_remains_visible():
    matrix = {row["quantity"]: row for row in load_matrix()}
    assert matrix["NATIVE_STAGE5_GHZ_SPECTRUM"]["main_limitation"] == "FULL_MAXWELL_REFERENCE_PENDING"


def test_synthetic_discrepancies_remain_separate_from_future_experiment():
    ledger = load_json("stage_i_final_discrepancy_ledger.json")
    assert validate_synthetic_discrepancy_separation(ledger)
    assert ledger["development_discrepancies"]
    assert ledger["future_experimental_discrepancies"] == []
    assert all(entry["evidence_type"] == "SYNTHETIC_DRY_RUN" for entry in ledger["development_discrepancies"])


def test_stage_j_export_inventory_schema_and_external_boundaries():
    inventory = load_json("stage_j_export_inventory.json")
    assert validate_export_inventory(inventory["items"])
    assert inventory["external_source_policy"] == {
        "Afivo": "MUST_REMAIN_EXTERNAL",
        "openEMS": "MUST_REMAIN_EXTERNAL",
    }
    assert any(item["classification"] == "DO_NOT_PACKAGE_RAW_TEMP" for item in inventory["items"])


def test_real_experiment_reentry_contract_is_ordered_and_no_replay():
    contract = load_json("stage_i_real_data_reentry_contract.json")
    assert contract["workflow"] == list(REENTRY_STEPS)
    assert validate_reentry_contract(contract)
    assert contract["rerun_tool_development_dry_runs"] is False


def test_final_contract_hashes_and_statuses_are_complete():
    contract = load_json("stage_i_final_contract.json")
    assert contract["STAGE_I_TOOL_DEVELOPMENT"] == "PASS"
    assert contract["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert contract["STAGE_I_EXPERIMENTAL_DATA"] == "NOT_PROVIDED"
    assert len(contract["frozen_upstream_hashes"]) == 5
    assert all(len(value) == 64 for value in contract["frozen_upstream_hashes"].values())
