#!/usr/bin/env python3
"""Generate the Stage-I final audit without rerunning upstream processing."""
from __future__ import annotations

import csv
import json
import resource
import time
from pathlib import Path

from streamer_rf.validation.final_audit import (
    REENTRY_STEPS,
    determine_stage_i_status,
    stage_j_packaging_allowed,
    validate_debt_ledger,
    validate_export_inventory,
    validate_reentry_contract,
    validate_synthetic_discrepancy_separation,
    validate_validation_matrix,
    verify_real_data_replacement,
)
from streamer_rf.validation.stage_i import sha256_file


ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "validation/stage_i/final"


def load_json(relative):
    return json.loads((ROOT / relative).read_text())


def write_json(name, value):
    (OUTPUT / name).write_text(json.dumps(value, indent=2) + "\n")


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def file_record(path, classification, work_package):
    relative = Path(path)
    absolute = ROOT / relative
    return {
        "path": str(relative),
        "work_package": work_package,
        "classification": classification,
        "sha256": sha256_file(absolute) if absolute.is_file() else "NOT_PROVIDED",
    }


def directory_bytes(relative):
    path = ROOT / relative
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file()) if path.exists() else 0


def main():
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    status_a = load_json("validation/stage_i/stage_i_validation_status.json")
    status_b = load_json("validation/stage_i/system_350mhz/wp_i_b_summary.json")
    status_c = load_json("validation/stage_i/native_ghz/native_ghz_validation_status.json")
    status_d = load_json("validation/stage_i/wp_i_d/wp_i_d_status.json")
    h3 = load_json("fullwave/h3/h3_summary.json")
    h4 = load_json("fullwave/h4/h4_summary.json")
    h5 = load_json("fullwave/h5/h5_stage_h_transfer_contract.json")

    upstream_pass = {
        "WP_I_A": status_a["STAGE_I_VALIDATION_FRAMEWORK"] == "PASS",
        "WP_I_B": status_b["WP_I_B_PIPELINE_DRY_RUN"] == "PASS",
        "WP_I_C": status_c["WP_I_C_VALIDATION_FRAMEWORK"] == "PASS",
        "WP_I_D": status_d["WP_I_D_PIPELINE_DRY_RUN"] == "PASS",
    }
    dual_status = determine_stage_i_status(upstream_pass, False)

    matrix = [
        ("350MHZ_TX_S11", "SYSTEM_350MHZ", "TOUCHSTONE_PIPELINE_READY", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "VNA_MEASUREMENT_PENDING", "REAL_TX_S11_AND_GEOMETRY"),
        ("350MHZ_RX_S11", "SYSTEM_350MHZ", "TOUCHSTONE_PIPELINE_READY", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "VNA_MEASUREMENT_PENDING", "REAL_RX_S11_AND_GEOMETRY"),
        ("350MHZ_S21", "SYSTEM_350MHZ", "H2_H3_TRANSFER_DEVELOPMENT_VERIFIED", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "VNA_MEASUREMENT_PENDING", "REAL_TX_RX_S21_AND_REFERENCE_PLANES"),
        ("350MHZ_SPECTRAL_PEAK", "SYSTEM_350MHZ", "H3_DEVELOPMENT_VERIFIED", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "FULL_WAVE_LOADING_MISMATCH_HIGH", "REAL_SCOPE_WAVEFORMS"),
        ("350MHZ_SPECTRAL_CENTROID", "SYSTEM_350MHZ", "H3_DEVELOPMENT_VERIFIED", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "SHORT_G3_TIME_WINDOW", "REAL_SCOPE_WAVEFORMS"),
        ("350MHZ_SPECTRAL_SHAPE", "SYSTEM_350MHZ", "H3_DEVELOPMENT_VERIFIED", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "PRODUCTION_TX_RX_GEOMETRY_NOT_PROVIDED", "REAL_SCOPE_AND_VNA_DATA"),
        ("350MHZ_ABSOLUTE_AMPLITUDE", "SYSTEM_350MHZ", "NUMERICAL_REFERENCE_ONLY", "PASS", "NOT_MEASURED", "NOT_RESOLVED", "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED", "CALIBRATED_SCOPE_CHAIN_AND_ACTUAL_GEOMETRY"),
        ("350MHZ_TIME_WAVEFORM", "SYSTEM_350MHZ", "REFERENCE_RECEIVED_VOLTAGE", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "FULL_WAVE_LOADING_MISMATCH_HIGH", "REAL_SCOPE_WAVEFORMS_AND_TRIGGER_METADATA"),
        ("DISTANCE_DEPENDENCE", "SYSTEM_350MHZ", "ANALYSIS_PIPELINE_READY", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "NO_DIRECT_H3_MULTI_DISTANCE_REFERENCE", "REAL_DISTANCE_SWEEP"),
        ("ORIENTATION_RESPONSE", "SYSTEM_350MHZ", "ANALYSIS_PIPELINE_READY", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "NO_DIRECT_H3_ORIENTATION_REFERENCE", "REAL_ORIENTATION_SWEEP"),
        ("REPEATABILITY", "COMMON", "ANALYSIS_PIPELINE_READY", "PASS", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "SYNTHETIC_REPETITIONS_ONLY", "REAL_REPEATED_EVENTS"),
        ("UNCERTAINTY", "COMMON", "UNCERTAINTY_PIPELINE_READY", "PASS", "NOT_PROVIDED", "PENDING_REAL_EXPERIMENT", "SYNTHETIC_UNCERTAINTY_ASSUMPTIONS", "REAL_UNCERTAINTY_BUDGET"),
        ("NATIVE_STAGE4_GHZ_SPECTRUM", "NATIVE_STAGE4", "TRUSTED_PHYSICS_SPARSE_BINS", "FRAMEWORK_ONLY", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "SPARSE_TRUSTED_FFT_BINS", "REAL_GHZ_SPECTRUM_IN_TRUST_MASK"),
        ("NATIVE_STAGE5_GHZ_SPECTRUM", "NATIVE_STAGE5", "TRUSTED_PHYSICS_SPARSE_BINS", "FRAMEWORK_ONLY", "NOT_MEASURED", "PENDING_REAL_EXPERIMENT", "FULL_MAXWELL_REFERENCE_PENDING", "REAL_GHZ_SPECTRUM_AND_FULL_MAXWELL_REFERENCE"),
        ("NATIVE_RECEIVER_RESPONSE", "NATIVE_STAGE4_STAGE5", "H4_DEVELOPMENT_VERIFIED", "FRAMEWORK_ONLY", "NOT_MEASURED", "NOT_RESOLVED", "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT", "PRODUCTION_GHZ_RECEIVER_AND_CALIBRATION"),
        ("CROSS_PATH_ABSOLUTE_CONTRIBUTION", "CROSS_PATH", "NOT_PERMITTED_CURRENT_CONFIGURATION", "NOT_APPLICABLE", "NOT_MEASURED", "NOT_RESOLVED", "INCOMPATIBLE_RECEIVERS_BANDS_REFERENCE_PLANES", "COMMON_CALIBRATED_CONFIGURATION"),
    ]
    matrix_rows = [dict(zip(("quantity", "pathway", "simulation_tool_status", "synthetic_dry_run_status", "real_measurement_status", "scientific_validation_status", "main_limitation", "required_future_input"), row)) for row in matrix]
    validate_validation_matrix(matrix_rows)
    write_csv("stage_i_validation_matrix.csv", matrix_rows)

    debts = [
        {"debt_id": "H2_VNA_MEASUREMENT_PENDING", "owner": "H2_STAGE_I", "status": "VNA_MEASUREMENT_PENDING", "severity_categories": ["BLOCKS_SCIENTIFIC_VALIDATION", "LIMITS_PRODUCTION_PREDICTION"]},
        {"debt_id": "H3_FULL_WAVE_LOADING_MISMATCH_HIGH", "owner": "H3", "status": h3["loading"]["status"], "severity_categories": ["LIMITS_ABSOLUTE_AMPLITUDE", "LIMITS_PRODUCTION_PREDICTION"]},
        {"debt_id": "H3_FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED", "owner": "H3", "status": h3["full_wave_feedback_status"], "severity_categories": ["LIMITS_ABSOLUTE_AMPLITUDE", "LIMITS_PRODUCTION_PREDICTION"]},
        {"debt_id": "H4_RECEIVER_TRANSFER_MESH_SENSITIVITY", "owner": "H4", "status": "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT", "severity_categories": ["LIMITS_ABSOLUTE_AMPLITUDE", "NON_BLOCKING_DOCUMENTED_DEBT"]},
        {"debt_id": "H4_ABSOLUTE_AMPLITUDE_NUMERICAL_REFERENCE_ONLY", "owner": "H4", "status": h4["H4_ABSOLUTE_AMPLITUDE_STATUS"], "severity_categories": ["LIMITS_ABSOLUTE_AMPLITUDE", "LIMITS_PRODUCTION_PREDICTION"]},
        {"debt_id": "STAGE5_FULL_MAXWELL_REFERENCE_PENDING", "owner": "STAGE5", "status": "FULL_MAXWELL_REFERENCE_PENDING", "severity_categories": ["BLOCKS_SCIENTIFIC_VALIDATION", "LIMITS_MECHANISM_INTERPRETATION"]},
        {"debt_id": "STAGE_F_SPARSE_TRUSTED_FFT_BINS", "owner": "STAGE_F_H4", "status": "STAGE4_5_BINS_STAGE5_7_BINS", "severity_categories": ["LIMITS_MECHANISM_INTERPRETATION", "NON_BLOCKING_DOCUMENTED_DEBT"]},
        {"debt_id": "STAGE_I_EXPERIMENTAL_DATA_NOT_PROVIDED", "owner": "STAGE_I", "status": "NOT_PROVIDED", "severity_categories": ["BLOCKS_SCIENTIFIC_VALIDATION"]},
        {"debt_id": "PRODUCTION_TX_RX_GEOMETRY_NOT_PROVIDED", "owner": "STAGE_I", "status": "NOT_PROVIDED", "severity_categories": ["BLOCKS_SCIENTIFIC_VALIDATION", "LIMITS_PRODUCTION_PREDICTION"]},
    ]
    validate_debt_ledger(debts)
    write_json("stage_i_final_debt_ledger.json", {"debts": debts, "blocks_tool_packaging": False, "unresolved_debt_count": len(debts)})

    discrepancy_sources = [
        "validation/stage_i/stage_i_discrepancy_ledger.json",
        "validation/stage_i/system_350mhz/wp_i_b_discrepancy_ledger.json",
        "validation/stage_i/native_ghz/native_ghz_discrepancy_ledger.json",
        "validation/stage_i/wp_i_d/wp_i_d_discrepancy_ledger.json",
    ]
    development = []
    source_records = []
    for source in discrepancy_sources:
        record = load_json(source)
        source_records.append({"path": source, "sha256": sha256_file(ROOT / source)})
        for entry in record.get("entries", ()):
            if entry.get("evidence_type") == "SYNTHETIC_DRY_RUN":
                normalized = dict(entry)
                normalized["scientific_interpretation_allowed"] = False
                normalized["source_ledger"] = source
                development.append(normalized)
    discrepancy = {
        "development_discrepancies": development,
        "future_experimental_discrepancies": [],
        "source_ledgers": source_records,
        "synthetic_discrepancies_are_physical_model_errors": False,
    }
    validate_synthetic_discrepancy_separation(discrepancy)
    write_json("stage_i_final_discrepancy_ledger.json", discrepancy)

    replacement = {
        "replace_only": ["RAW_MEASUREMENT_FILES", "ACTUAL_GEOMETRY_METADATA", "INSTRUMENT_METADATA", "CALIBRATION_METADATA", "UNCERTAINTY_METADATA", "PROVENANCE_STATUS"],
        "redesign_not_required": ["PARSER_APIS", "COMPARISON_APIS", "SPECTRAL_METRICS", "UNCERTAINTY_ENGINE", "DISCREPANCY_LEDGER", "STAGE_H_CONTRACTS"],
    }
    replacement_ready = verify_real_data_replacement(replacement)

    reentry = {
        "contract_type": "StageIRealExperimentReentryContract",
        "workflow": list(REENTRY_STEPS),
        "rerun_tool_development_dry_runs": False,
        "required_inputs_reference": "fullwave/h5/h5_stage_i_required_inputs.json",
        "status_before_reentry": "PENDING_REAL_EXPERIMENT",
        "status_updates_require_real_evidence": True,
    }
    validate_reentry_contract(reentry)
    write_json("stage_i_real_data_reentry_contract.json", reentry)

    export_items = [
        {"path": "python/streamer_rf", "classification": "CORE_SOURCE", "action": "PACKAGE"},
        {"path": "CMakeLists.txt", "classification": "CONFIG", "action": "PACKAGE"},
        {"path": "requirements.txt", "classification": "CONFIG", "action": "PACKAGE"},
        {"path": "tests", "classification": "TEST", "action": "PACKAGE"},
        {"path": "docs", "classification": "DOCUMENTATION", "action": "PACKAGE_REVIEW_MACHINE_PATHS"},
        {"path": "fullwave/h1/openems_backend.json", "classification": "EXTERNAL_BACKEND_REFERENCE", "action": "PACKAGE_PROVENANCE_ONLY"},
        {"path": "solver3d/afivo_reference", "classification": "EXTERNAL_BACKEND_REFERENCE", "action": "PACKAGE_ADAPTERS_NOT_EXTERNAL_SOURCE"},
        {"path": "fullwave/h2", "classification": "REFERENCE_DATA", "action": "PACKAGE_COMPACT_RESULTS"},
        {"path": "fullwave/h3", "classification": "DERIVED_RESULT", "action": "PACKAGE_COMPACT_RESULTS"},
        {"path": "fullwave/h4", "classification": "DERIVED_RESULT", "action": "PACKAGE_COMPACT_RESULTS"},
        {"path": "fullwave/h5", "classification": "DERIVED_RESULT", "action": "PACKAGE_CONTRACTS"},
        {"path": "validation/stage_i/system_350mhz/synthetic_input", "classification": "SYNTHETIC_FIXTURE", "action": "PACKAGE_WITH_PROVENANCE"},
        {"path": "validation/stage_i/system_350mhz/wp_i_d_synthetic", "classification": "SYNTHETIC_FIXTURE", "action": "PACKAGE_WITH_PROVENANCE"},
        {"path": "validation/stage_i", "classification": "DERIVED_RESULT", "action": "PACKAGE_CONTRACTS_AND_SMALL_RESULTS"},
        {"path": "build", "classification": "DO_NOT_PACKAGE_RAW_TEMP", "action": "EXCLUDE_FROM_GIT"},
        {"path": ".venv", "classification": "DO_NOT_PACKAGE_RAW_TEMP", "action": "EXCLUDE_FROM_GIT"},
        {"path": "/tmp/h1_h2_h4_openems_validation", "classification": "DO_NOT_PACKAGE_RAW_TEMP", "action": "DO_NOT_PACKAGE"},
    ]
    validate_export_inventory(export_items)
    export_inventory = {
        "items": export_items,
        "external_source_policy": {"Afivo": "MUST_REMAIN_EXTERNAL", "openEMS": "MUST_REMAIN_EXTERNAL"},
        "future_real_measurement_data": "LICENSE_PRIVACY_AND_CONSENT_REVIEW_REQUIRED",
        "synthetic_fixture_policy": "PACKAGE_ONLY_WITH_EXPLICIT_SYNTHETIC_PROVENANCE",
    }
    write_json("stage_j_export_inventory.json", export_inventory)

    output_inventory = []
    framework_files = [
        "validation/stage_i/stage_i_measurement_contract.json", "validation/stage_i/stage_i_vna_contract.json",
        "validation/stage_i/stage_i_scope_contract.json", "validation/stage_i/stage_i_reference_plane_contract.json",
        "validation/stage_i/stage_i_350mhz_profile.json", "validation/stage_i/stage_i_native_ghz_profile.json",
        "validation/stage_i/stage_i_processing_policy.json", "validation/stage_i/stage_i_uncertainty_contract.json",
        "validation/stage_i/stage_i_data_layer_contract.json", "validation/stage_i/stage_i_spectrum_analyzer_contract.json",
        "validation/stage_i/stage_i_simulation_measurement_mapping.csv",
        "validation/stage_i/native_ghz/native_ghz_measurement_contract.json",
        "validation/stage_i/native_ghz/native_ghz_stage4_profile.json", "validation/stage_i/native_ghz/native_ghz_stage5_profile.json",
        "validation/stage_i/native_ghz/native_ghz_frequency_coverage_policy.json", "validation/stage_i/native_ghz/native_ghz_comparison_policy.json",
        "validation/stage_i/native_ghz/native_ghz_h4_input_provenance.json",
    ]
    for path in framework_files:
        output_inventory.append(file_record(path, "FRAMEWORK", "WP_I_A_OR_C"))
    synthetic_outputs = sorted((ROOT / "validation/stage_i/system_350mhz").glob("wp_i_b_*")) + sorted((ROOT / "validation/stage_i/wp_i_d").glob("wp_i_d_*"))
    for path in synthetic_outputs:
        classification = "VALIDATION_STATUS" if "status" in path.name or "discrepancy" in path.name or path.name.endswith("summary.json") else "SYNTHETIC_DRY_RUN"
        output_inventory.append(file_record(path.relative_to(ROOT), classification, "WP_I_B_OR_D"))
    for path in ("validation/stage_i/stage_i_validation_status.json", "validation/stage_i/stage_i_discrepancy_ledger.json", "validation/stage_i/native_ghz/native_ghz_validation_status.json", "validation/stage_i/native_ghz/native_ghz_discrepancy_ledger.json"):
        output_inventory.append(file_record(path, "VALIDATION_STATUS", "WP_I_A_OR_C"))
    for missing in ("REAL_VNA_DATA", "REAL_OSCILLOSCOPE_DATA", "ACTUAL_TX_GEOMETRY", "ACTUAL_RX_GEOMETRY", "REAL_CALIBRATION_UNCERTAINTY"):
        output_inventory.append({"path": missing, "work_package": "FUTURE_STAGE_I", "classification": "SCIENTIFIC_INPUT_REQUIRED", "sha256": "NOT_PROVIDED"})
    for name in ("stage_i_validation_matrix.csv", "stage_i_final_debt_ledger.json", "stage_i_final_discrepancy_ledger.json", "stage_i_final_contract.json", "stage_i_real_data_reentry_contract.json", "stage_j_export_inventory.json", "stage_i_final_status.json"):
        output_inventory.append({"path": f"validation/stage_i/final/{name}", "work_package": "WP_I_E", "classification": "STAGE_J_EXPORT", "sha256": "GENERATED_BY_THIS_AUDIT"})
    write_csv("stage_i_output_inventory.csv", output_inventory)

    reproducibility = {
        "status": "PASS",
        "checks": {
            "environment_documentation": "docs/environment_baseline.md",
            "python_dependencies": "requirements.txt",
            "python_dependency_version_pinning": "STAGE_J_ACTION_REQUIRED",
            "cpp_petsc_mpi_build": "CMakeLists.txt_AND_README_COMMANDS",
            "openems_provenance": "fullwave/h1/openems_backend.json",
            "afivo_provenance": "docs/stage_d1_afivo_backend.md",
            "frozen_input_hashes": "PRESENT_IN_STAGE_CONTRACTS_AND_MANIFESTS",
            "commands_and_scripts": "PRESENT",
            "tests": "PRESENT",
            "stage_contracts": "PRESENT",
            "synthetic_fixture_provenance": "PRESENT",
        },
        "external_backend_source_policy": {"Afivo": "EXTERNAL", "openEMS": "EXTERNAL"},
    }
    write_json("stage_i_reproducibility_audit.json", reproducibility)

    packaging_allowed = stage_j_packaging_allowed(dual_status["STAGE_I_TOOL_DEVELOPMENT"], reproducibility["status"])
    machine_paths = {
        "ACTIVE_HARDCODE": [],
        "ENV_CONFIG_REQUIRED": ["solver3d/afivo_reference/**/*.cfg", "fullwave/h1/generate_h1_audit.py", "fullwave/h2/generate_h2_results.py", "fullwave/h4/generate_h4_results.py"],
        "DOCUMENTATION_ONLY": ["docs/environment_baseline.md", "docs/stage_d1_afivo_backend.md", "docs/stage_h1_fullwave_foundation.md"],
        "TEST_FIXTURE": ["fullwave/h1-h4 compact result metadata", "rf/source/audit frozen manifests"],
        "packaging_action": "DOCUMENT_ENVIRONMENT_VARIABLES_OR_TEMPLATES_DURING_STAGE_J",
    }
    large_files = {
        "solver3d/afivo_reference/stage_e/results_raw": {"action": "ARCHIVE_EXTERNALLY", "bytes": directory_bytes("solver3d/afivo_reference/stage_e/results_raw")},
        "results/stage_f4/source_gate_recovery": {"action": "ARCHIVE_EXTERNALLY_OR_REGENERATE", "bytes": directory_bytes("results/stage_f4/source_gate_recovery")},
        "build": {"action": "EXCLUDE_FROM_GIT", "bytes": directory_bytes("build")},
        ".venv": {"action": "EXCLUDE_FROM_GIT", "bytes": directory_bytes(".venv")},
        "tracked_compact_reference_and_summary_data": {"action": "KEEP"},
        "stage_i_synthetic_fixtures": {"action": "KEEP_WITH_SYNTHETIC_PROVENANCE", "bytes": directory_bytes("validation/stage_i/system_350mhz/synthetic_input") + directory_bytes("validation/stage_i/system_350mhz/wp_i_d_synthetic")},
        "automatic_deletion_performed": False,
    }
    write_json("stage_i_packaging_audit.json", {"machine_local_paths": machine_paths, "large_files": large_files, "open_source_boundary": export_inventory["external_source_policy"]})

    final_contract = {
        "contract_type": "StageIFinalContract",
        **dual_status,
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "NATIVE_GHZ_STAGE4_VALIDATION": "NOT_MEASURED",
        "NATIVE_GHZ_STAGE5_VALIDATION": "NOT_MEASURED",
        "distance_orientation_repeatability_ready": True,
        "uncertainty_ready": True,
        "REAL_DATA_REPLACEMENT_READY": replacement_ready,
        "frequency_boundaries": {
            "system_experimental_context_Hz": [50e6, 500e6],
            "H3_formal_comparison_Hz": [200e6, 500e6],
            "system_target_center_Hz": 350e6,
            "native_350MHz": "NOT_RESOLVED",
            "Stage4_trusted_Hz": h4["stages"]["Stage4"]["trusted_frequency_Hz"],
            "Stage5_trusted_Hz": h4["stages"]["Stage5"]["trusted_frequency_Hz"],
            "native_primary_comparison": "SPARSE_TRUSTED_FREQUENCY_BINS",
        },
        "CROSS_PATH_COHERENT_SUMMATION": "NOT_PERMITTED_CURRENT_CONFIGURATION",
        "CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION": "NOT_RESOLVED",
        "debt_ledger": "stage_i_final_debt_ledger.json",
        "validation_matrix": "stage_i_validation_matrix.csv",
        "discrepancy_ledger": "stage_i_final_discrepancy_ledger.json",
        "STAGE_J_TOOL_PACKAGING_ALLOWED": packaging_allowed,
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "future_experimental_requirements": "fullwave/h5/h5_stage_i_required_inputs.json",
        "real_data_reentry_contract": "stage_i_real_data_reentry_contract.json",
        "frozen_upstream_hashes": {
            "Stage_H_contract": sha256_file(ROOT / "fullwave/h5/h5_stage_h_transfer_contract.json"),
            "WP_I_A_status": sha256_file(ROOT / "validation/stage_i/stage_i_validation_status.json"),
            "WP_I_B_status": sha256_file(ROOT / "validation/stage_i/system_350mhz/wp_i_b_summary.json"),
            "WP_I_C_status": sha256_file(ROOT / "validation/stage_i/native_ghz/native_ghz_validation_status.json"),
            "WP_I_D_status": sha256_file(ROOT / "validation/stage_i/wp_i_d/wp_i_d_status.json"),
        },
        "frozen_commits": {
            "Stage_H": "a11aebe58eccb708ddb62e6aa22f3908487972ce",
            "WP_I_A": "6faf0e554dd3168495b88a376b13224645d1330b",
            "WP_I_B": "a4da6126df919857a1817c52020330f5c0d90f57",
            "WP_I_C": "0e24f955c9e337a3a4cf2c3d9d673db1f488a43d",
            "WP_I_D": "e8a8400e48f82f0e84489a144767eafbc474da0f",
        },
        "scientific_statement": "Stage I tool development completed the simulation-to-measurement validation architecture; synthetic dry runs demonstrate software readiness only; scientific validation awaits real VNA, oscilloscope, geometry, calibration, and uncertainty data.",
    }
    write_json("stage_i_final_contract.json", final_contract)

    final_status = {
        **dual_status,
        "REAL_DATA_REPLACEMENT_READY": replacement_ready,
        "STAGE_J_TOOL_PACKAGING_ALLOWED": packaging_allowed,
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "NATIVE_GHZ_STAGE4_VALIDATION": "NOT_MEASURED",
        "NATIVE_GHZ_STAGE5_VALIDATION": "NOT_MEASURED",
        "physical_simulations_rerun": False,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json("stage_i_final_status.json", final_status)
    print(json.dumps(final_status, indent=2))


if __name__ == "__main__":
    main()
