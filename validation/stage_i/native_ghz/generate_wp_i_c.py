#!/usr/bin/env python3
"""Generate WP-I-C native-GHz validation contracts from frozen H4 outputs."""
from __future__ import annotations

import json
from pathlib import Path
import resource
import time

import numpy as np

from streamer_rf.validation.native_ghz import (
    NATIVE_DISCREPANCY_CATEGORIES,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE4_TRUSTED_BINS_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BINS_HZ,
    frequency_coverage_status,
    trusted_bin_support,
    validate_native_measurement_contract,
)
from streamer_rf.validation.stage_i import NOT_PROVIDED, sha256_file


ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "validation/stage_i/native_ghz"
H4_CONTRACT_PATH = ROOT / "fullwave/h4/h4_result_contract.json"
H4_SUMMARY_PATH = ROOT / "fullwave/h4/h4_summary.json"


def load_json(path):
    return json.loads(Path(path).read_text())


def write_json(name, record):
    (OUTPUT / name).write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def verify_h4_inputs():
    contract = load_json(H4_CONTRACT_PATH)
    summary = load_json(H4_SUMMARY_PATH)
    receiver_path = ROOT / contract["receiver_transfer_path"]
    if sha256_file(receiver_path) != contract["receiver_transfer_hash"]:
        raise ValueError("H4_RECEIVER_TRANSFER_HASH_MISMATCH")
    receiver = np.genfromtxt(receiver_path, delimiter=",", names=True)
    receiver_support = [float(receiver["frequency_Hz"].min()), float(receiver["frequency_Hz"].max())]
    stages = {}
    expected = {
        "Stage4": (STAGE4_TRUSTED_BAND_HZ, STAGE4_TRUSTED_BINS_HZ),
        "Stage5": (STAGE5_TRUSTED_BAND_HZ, STAGE5_TRUSTED_BINS_HZ),
    }
    for stage, (band, bins) in expected.items():
        source = contract["source_contracts"][stage]
        if tuple(source["trusted_frequency_mask_Hz"]) != band:
            raise ValueError(f"{stage.upper()}_TRUST_BAND_CHANGED")
        if sha256_file(ROOT / source["source_path"]) != source["source_hash"]:
            raise ValueError(f"{stage.upper()}_SOURCE_HASH_MISMATCH")
        if sha256_file(ROOT / source["rf_trust_report_path"]) != source["rf_trust_report_hash"]:
            raise ValueError(f"{stage.upper()}_RF_TRUST_HASH_MISMATCH")
        spectrum_path = ROOT / contract["native_received_spectrum_paths"][stage]
        spectrum = np.genfromtxt(spectrum_path, delimiter=",", names=True)
        trusted = spectrum["frequency_Hz"][spectrum["trust_valid"].astype(bool)]
        if not np.array_equal(trusted, bins):
            raise ValueError(f"{stage.upper()}_TRUSTED_BINS_CHANGED")
        if summary["stages"][stage]["trusted_bin_count"] != len(bins):
            raise ValueError(f"{stage.upper()}_TRUSTED_BIN_COUNT_CHANGED")
        stages[stage] = {
            "source_case_id": source["source_case_id"],
            "source_path": source["source_path"],
            "source_hash": source["source_hash"],
            "rf_trust_report_path": source["rf_trust_report_path"],
            "rf_trust_report_hash": source["rf_trust_report_hash"],
            "observer_id": source["observer_id"],
            "observer_position_m": source["observer_position_m"],
            "observer_distance_m": source["observer_distance_m"],
            "field_components": source["field_components"],
            "receiver_axis": source["receiver_axis"],
            "NATIVE_PROPAGATION_ALREADY_INCLUDED": source["NATIVE_PROPAGATION_ALREADY_INCLUDED"],
            "trusted_frequency_band_Hz": list(band),
            "trusted_frequency_bins_Hz": bins.tolist(),
            "trusted_bin_count": len(bins),
            "received_native_spectrum_path": contract["native_received_spectrum_paths"][stage],
            "received_native_spectrum_hash": sha256_file(spectrum_path),
            "received_time_semantics": summary["stages"][stage]["received_time_semantics"],
        }
    if contract["Stage5_full_maxwell_status"] != "FULL_MAXWELL_REFERENCE_PENDING":
        raise ValueError("STAGE5_FULL_MAXWELL_DEBT_MISSING")
    return {
        "H4_contract_path": str(H4_CONTRACT_PATH.relative_to(ROOT)),
        "H4_contract_hash": sha256_file(H4_CONTRACT_PATH),
        "H4_summary_path": str(H4_SUMMARY_PATH.relative_to(ROOT)),
        "H4_summary_hash": sha256_file(H4_SUMMARY_PATH),
        "receiver_id": contract["receiver_id"],
        "receiver_transfer_path": contract["receiver_transfer_path"],
        "receiver_transfer_hash": contract["receiver_transfer_hash"],
        "receiver_transfer_support_Hz": receiver_support,
        "receiver_absolute_amplitude_status": contract["H4_ABSOLUTE_AMPLITUDE_STATUS"],
        "receiver_mesh_status": "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
        "stages": stages,
        "Stage5_full_maxwell_status": contract["Stage5_full_maxwell_status"],
        "inputs_recomputed": False,
    }


def measurement_template(profile):
    record = {
        "contract_type": "StageINativeGHzMeasurementContract",
        "experiment_id": NOT_PROVIDED,
        "validation_profile": profile,
        "receiver_id": NOT_PROVIDED,
        "receiver_geometry": NOT_PROVIDED,
        "receiver_orientation": NOT_PROVIDED,
        "distance_m": NOT_PROVIDED,
        "instrument": NOT_PROVIDED,
        "frequency_range_Hz": NOT_PROVIDED,
        "frequency_resolution_Hz": NOT_PROVIDED,
        "reference_plane": NOT_PROVIDED,
        "receiver_load_impedance_ohm": NOT_PROVIDED,
        "cable_metadata": NOT_PROVIDED,
        "calibration_metadata": NOT_PROVIDED,
        "raw_spectrum_path": NOT_PROVIDED,
        "raw_waveform_path": NOT_PROVIDED,
        "repetition_index": NOT_PROVIDED,
        "uncertainty_metadata": NOT_PROVIDED,
        "measurement_status": "NOT_MEASURED",
    }
    validate_native_measurement_contract(record)
    return record


def main():
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    h4 = verify_h4_inputs()
    write_json("native_ghz_h4_input_provenance.json", h4)

    measurement_contract = {
        "contract_type": "StageINativeGHzMeasurementContractSet",
        "profile_separation": "SYSTEM_350MHZ_AND_NATIVE_GHZ_INDEPENDENT",
        "templates": {
            "NATIVE_GHZ_STAGE4": measurement_template("NATIVE_GHZ_STAGE4"),
            "NATIVE_GHZ_STAGE5": measurement_template("NATIVE_GHZ_STAGE5"),
        },
        "raw_values_normalized": False,
        "missing_values": "NOT_PROVIDED",
        "scientific_data_present": False,
    }
    write_json("native_ghz_measurement_contract.json", measurement_contract)

    instruments = {
        "OSCILLOSCOPE_WIDEBAND_WAVEFORM": {
            "capability_bandwidth_Hz": 8e9,
            "capability_sample_rate_Hz": 80e9,
            "capability_vertical_bits": 12,
            "capability_only": True,
            "measurement_settings": NOT_PROVIDED,
            "Stage4_coverage": frequency_coverage_status([0.0, 8e9], STAGE4_TRUSTED_BAND_HZ),
            "Stage5_coverage": frequency_coverage_status([0.0, 8e9], STAGE5_TRUSTED_BAND_HZ),
            "Stage4_trusted_bins_covered": int(
                np.count_nonzero(trusted_bin_support(STAGE4_TRUSTED_BINS_HZ, [0.0, 8e9], h4["receiver_transfer_support_Hz"]))
            ),
            "Stage5_trusted_bins_covered": int(
                np.count_nonzero(trusted_bin_support(STAGE5_TRUSTED_BINS_HZ, [0.0, 8e9], h4["receiver_transfer_support_Hz"]))
            ),
        },
        "SPECTRUM_ANALYZER": {
            "capability_frequency_Hz": [2.0, 67e9],
            "capability_only": True,
            "measurement_settings": NOT_PROVIDED,
            "Stage4_coverage": frequency_coverage_status([2.0, 67e9], STAGE4_TRUSTED_BAND_HZ),
            "Stage5_coverage": frequency_coverage_status([2.0, 67e9], STAGE5_TRUSTED_BAND_HZ),
            "Stage4_trusted_bins_covered": 5,
            "Stage5_trusted_bins_covered": 7,
        },
    }
    coverage_policy = {
        "contract_type": "NativeGHzFrequencyCoveragePolicy",
        "coverage_statuses": [
            "FULL_TRUST_BAND_COVERAGE",
            "PARTIAL_TRUST_BAND_COVERAGE",
            "NO_TRUST_BAND_COVERAGE",
        ],
        "classification_basis": "MEASUREMENT_SUPPORT_VERSUS_EXACT_RFTRUSTREPORT_BAND",
        "comparison_bins": "EXACT_FROZEN_STAGE_F_TRUSTED_FFT_BINS",
        "instruments": instruments,
    }
    write_json("native_ghz_frequency_coverage_policy.json", coverage_policy)

    common = {
        "receiver_id": h4["receiver_id"],
        "receiver_transfer_hash": h4["receiver_transfer_hash"],
        "receiver_transfer_support_Hz": h4["receiver_transfer_support_Hz"],
        "primary_comparison": "FREQUENCY_DOMAIN_TRUSTED_BINS",
        "time_waveform_status": "SECONDARY_RECONSTRUCTION_ONLY",
        "absolute_amplitude_status": "NUMERICAL_REFERENCE_ONLY_AND_CALIBRATION_REQUIRED",
        "validation_status": "NOT_MEASURED",
        "production_receiver_status": "NOT_RESOLVED",
    }
    stage4_profile = {
        "VALIDATION_PROFILE": "NATIVE_GHZ_STAGE4",
        **common,
        **h4["stages"]["Stage4"],
        "full_maxwell_status": None,
    }
    stage5_profile = {
        "VALIDATION_PROFILE": "NATIVE_GHZ_STAGE5",
        **common,
        **h4["stages"]["Stage5"],
        "full_maxwell_status": "FULL_MAXWELL_REFERENCE_PENDING",
    }
    write_json("native_ghz_stage4_profile.json", stage4_profile)
    write_json("native_ghz_stage5_profile.json", stage5_profile)

    comparison_policy = {
        "contract_type": "NativeGHzComparisonPolicy",
        "comparison_mask": "MEASUREMENT_SUPPORT_INTERSECT_RFTRUSTREPORT_MASK_INTERSECT_RECEIVER_TRANSFER_SUPPORT",
        "interpolation_across_missing_trusted_bins": "FORBIDDEN",
        "validity_extrapolation": "FORBIDDEN",
        "NATIVE_GHZ_PRIMARY_COMPARISON": "FREQUENCY_DOMAIN_TRUSTED_BINS",
        "H4_time_waveform": "SECONDARY_RECONSTRUCTION_ONLY",
        "metrics": [
            "TRUSTED_BIN_PEAK_FREQUENCY",
            "TRUSTED_BIN_SPECTRAL_CENTROID",
            "NORMALIZED_TRUSTED_BIN_CORRELATION",
            "RELATIVE_TRUSTED_BIN_AMPLITUDE_PATTERN",
            "ABSOLUTE_VOLTAGE_ONLY_WITH_COMPATIBLE_CALIBRATION",
            "POLARIZATION_CONTRAST_WHERE_AVAILABLE",
        ],
        "dense_spectrum_assumption": False,
        "synthetic_hook": {
            "status": "READY_NO_FIXTURE_PROVIDED",
            "required_provenance": "SYNTHETIC_NATIVE_GHZ_DRY_RUN",
            "scientific_validation_allowed": False,
            "frequency_unit": "Hz",
            "allowed_quantity_units": ["V", "V_s", "V_per_m", "V_s_per_m"],
            "trust_mask_semantics": "STAGE_F_EXACT_TRUSTED_BINS_ONLY",
            "repetition_index_required": True,
        },
    }
    write_json("native_ghz_comparison_policy.json", comparison_policy)

    ledger = {
        "ledger_type": "StageINativeGHzDiscrepancyLedger",
        "pathways": ["NATIVE_STAGE4", "NATIVE_STAGE5"],
        "allowed_possible_sources": list(NATIVE_DISCREPANCY_CATEGORIES),
        "entries": [],
        "empty_reason": "STAGE_I_EXPERIMENTAL_DATA_NOT_PROVIDED",
        "automatic_cause_assignment": False,
    }
    write_json("native_ghz_discrepancy_ledger.json", ledger)

    status = {
        "WP_I_C_VALIDATION_FRAMEWORK": "PASS",
        "NATIVE_GHZ_PIPELINE_READY": True,
        "NATIVE_GHZ_STAGE4_VALIDATION": "NOT_MEASURED",
        "NATIVE_GHZ_STAGE5_VALIDATION": "NOT_MEASURED",
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "SYSTEM_350MHZ_PROFILE_MODIFIED": False,
        "H3_OUTPUTS_MODIFIED": False,
        "NATIVE_350MHZ_RF_TRUST": "NOT_RESOLVED",
        "PRODUCTION_NATIVE_RF_RECEIVER": "NOT_RESOLVED",
        "H4_ABSOLUTE_AMPLITUDE_STATUS": "NUMERICAL_REFERENCE_ONLY",
        "FULL_MAXWELL_REFERENCE_PENDING": True,
        "synthetic_fixture_status": "READY_NO_FIXTURE_PROVIDED",
        "physical_simulations_rerun": False,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json("native_ghz_validation_status.json", status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
