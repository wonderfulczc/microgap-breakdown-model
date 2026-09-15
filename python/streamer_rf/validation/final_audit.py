"""Stage-I final audit contracts and status guards."""
from __future__ import annotations


MATRIX_COLUMNS = (
    "quantity",
    "pathway",
    "simulation_tool_status",
    "synthetic_dry_run_status",
    "real_measurement_status",
    "scientific_validation_status",
    "main_limitation",
    "required_future_input",
)

EXPORT_CLASSES = {
    "CORE_SOURCE",
    "EXTERNAL_BACKEND_REFERENCE",
    "CONFIG",
    "TEST",
    "DOCUMENTATION",
    "REFERENCE_DATA",
    "SYNTHETIC_FIXTURE",
    "DERIVED_RESULT",
    "DO_NOT_PACKAGE_RAW_TEMP",
}

REENTRY_STEPS = (
    "REAL_DATA_INGESTION",
    "PROVENANCE_GATE",
    "QUALITY_GATE",
    "CALIBRATION_REFERENCE_PLANE_CHECK",
    "STAGE_I_COMPARISON",
    "DISCREPANCY_LEDGER_UPDATE",
    "VALIDATION_STATUS_UPDATE",
)


def validate_validation_matrix(rows):
    if not rows:
        raise ValueError("EMPTY_STAGE_I_VALIDATION_MATRIX")
    quantities = set()
    for row in rows:
        missing = set(MATRIX_COLUMNS) - set(row)
        if missing:
            raise ValueError(f"VALIDATION_MATRIX_FIELDS_MISSING:{sorted(missing)}")
        if row["quantity"] in quantities:
            raise ValueError("DUPLICATE_VALIDATION_MATRIX_QUANTITY")
        quantities.add(row["quantity"])
    return True


def determine_stage_i_status(upstream_pass, real_experimental_data_available):
    required = {"WP_I_A", "WP_I_B", "WP_I_C", "WP_I_D"}
    if set(upstream_pass) != required or not all(upstream_pass.values()):
        return {
            "STAGE_I_TOOL_DEVELOPMENT": "NOT_RESOLVED",
            "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        }
    return {
        "STAGE_I_TOOL_DEVELOPMENT": "PASS",
        "STAGE_I_SCIENTIFIC_VALIDATION": (
            "COMPARISON_READY" if real_experimental_data_available else "PENDING_REAL_EXPERIMENT"
        ),
    }


def verify_real_data_replacement(replacement_contract):
    required_replace = {
        "RAW_MEASUREMENT_FILES",
        "ACTUAL_GEOMETRY_METADATA",
        "INSTRUMENT_METADATA",
        "CALIBRATION_METADATA",
        "UNCERTAINTY_METADATA",
        "PROVENANCE_STATUS",
    }
    protected_apis = {
        "PARSER_APIS",
        "COMPARISON_APIS",
        "SPECTRAL_METRICS",
        "UNCERTAINTY_ENGINE",
        "DISCREPANCY_LEDGER",
        "STAGE_H_CONTRACTS",
    }
    return (
        set(replacement_contract.get("replace_only", ())) == required_replace
        and set(replacement_contract.get("redesign_not_required", ())) == protected_apis
    )


def stage_j_packaging_allowed(tool_status, reproducibility_status):
    return tool_status == "PASS" and reproducibility_status == "PASS"


def validate_debt_ledger(ledger):
    allowed = {
        "BLOCKS_TOOL_PACKAGING",
        "BLOCKS_SCIENTIFIC_VALIDATION",
        "LIMITS_ABSOLUTE_AMPLITUDE",
        "LIMITS_MECHANISM_INTERPRETATION",
        "LIMITS_PRODUCTION_PREDICTION",
        "NON_BLOCKING_DOCUMENTED_DEBT",
    }
    if not ledger:
        raise ValueError("EMPTY_SCIENTIFIC_DEBT_LEDGER")
    for entry in ledger:
        categories = set(entry.get("severity_categories", ()))
        if not entry.get("debt_id") or not categories or not categories <= allowed:
            raise ValueError("INVALID_SCIENTIFIC_DEBT_ENTRY")
    return True


def validate_synthetic_discrepancy_separation(contract):
    for entry in contract.get("development_discrepancies", ()):
        if entry.get("evidence_type") != "SYNTHETIC_DRY_RUN":
            raise ValueError("SYNTHETIC_DISCREPANCY_PROVENANCE_LOST")
        if entry.get("scientific_interpretation_allowed") is not False:
            raise ValueError("SYNTHETIC_DISCREPANCY_CANNOT_BE_SCIENTIFIC")
    if any(entry.get("evidence_type") == "SYNTHETIC_DRY_RUN" for entry in contract.get("future_experimental_discrepancies", ())):
        raise ValueError("SYNTHETIC_ENTRY_IN_EXPERIMENTAL_LEDGER")
    return True


def validate_export_inventory(items):
    if not items:
        raise ValueError("EMPTY_STAGE_J_EXPORT_INVENTORY")
    for item in items:
        if not item.get("path") or item.get("classification") not in EXPORT_CLASSES:
            raise ValueError("INVALID_STAGE_J_EXPORT_ITEM")
    return True


def validate_reentry_contract(contract):
    if tuple(contract.get("workflow", ())) != REENTRY_STEPS:
        raise ValueError("INVALID_STAGE_I_REENTRY_WORKFLOW")
    if contract.get("rerun_tool_development_dry_runs") is not False:
        raise ValueError("REENTRY_MUST_NOT_REQUIRE_DRY_RUN_REPLAY")
    return True


def enforce_scientific_validation_evidence(status, real_data_available):
    if status in {"VALIDATED", "PASS"} and not real_data_available:
        raise ValueError("REAL_EXPERIMENT_REQUIRED_FOR_SCIENTIFIC_VALIDATION")
    return status
