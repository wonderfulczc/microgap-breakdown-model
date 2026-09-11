#!/usr/bin/env python3
"""Generate WP-I-A contracts without creating experimental data."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.receiver import read_touchstone  # noqa: E402
from streamer_rf.validation.stage_i import (  # noqa: E402
    DISCREPANCY_SOURCES,
    NOT_PROVIDED,
    REFERENCE_PLANES,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    amplitude_comparison_class,
    compare_spectra,
    sha256_file,
    validate_data_layer,
    validate_discrepancy_ledger,
    validate_measurement_contract,
    validate_scope_contract,
    validate_vna_contract,
)


OUT = ROOT / "validation/stage_i"


def write_json(name, record):
    (OUT / name).write_text(json.dumps(record, indent=2) + "\n")


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    upstream_path = ROOT / "fullwave/h5/h5_stage_i_required_inputs.json"
    upstream = json.loads(upstream_path.read_text())

    inventory = {
        "streamer_rf.fullwave.receiver.read_touchstone": "REUSABLE",
        "fullwave/h2/h2_350mhz_future_vna_contract.json": "PARTIAL",
        "fullwave/h2/h2_vna_measurement_protocol.json": "PARTIAL",
        "fullwave/h2/theory_inputs": "PARTIAL",
        "fullwave/h2/h2_reference_sparameters.csv": "LEGACY",
        "fullwave/h3_h4_simulation_waveforms": "REUSABLE",
        "physical_VNA_measurements": "NOT_APPLICABLE",
        "physical_oscilloscope_waveforms": "NOT_APPLICABLE",
        "calibration_records": "NOT_APPLICABLE",
    }

    measurement = {
        "contract_type": "StageIMeasurementContract",
        "contract_status": "MEASUREMENT_TEMPLATE_NO_DATA",
        "upstream_requirements_path": "fullwave/h5/h5_stage_i_required_inputs.json",
        "upstream_requirements_hash": sha256_file(upstream_path),
        "experiment_id": "STAGE_I_PENDING_EXPERIMENT",
        "timestamp": NOT_PROVIDED,
        "hardware_configuration": NOT_PROVIDED,
        "transmitter_id": NOT_PROVIDED,
        "receiver_id": NOT_PROVIDED,
        "actual_tx_geometry": NOT_PROVIDED,
        "actual_rx_geometry": NOT_PROVIDED,
        "distance_m": NOT_PROVIDED,
        "tx_orientation": NOT_PROVIDED,
        "rx_orientation": NOT_PROVIDED,
        "polarization": NOT_PROVIDED,
        "environment": NOT_PROVIDED,
        "reference_planes": NOT_PROVIDED,
        "cable_metadata": NOT_PROVIDED,
        "connector_metadata": NOT_PROVIDED,
        "instrument_metadata": NOT_PROVIDED,
        "calibration_metadata": NOT_PROVIDED,
        "measurement_bandwidth_Hz": NOT_PROVIDED,
        "sample_rate_Hz": NOT_PROVIDED,
        "number_of_samples": NOT_PROVIDED,
        "averaging": NOT_PROVIDED,
        "repetition_index": NOT_PROVIDED,
        "raw_data_path": NOT_PROVIDED,
        "raw_data_hash": NOT_PROVIDED,
        "uncertainty_metadata": NOT_PROVIDED,
    }
    validate_measurement_contract(measurement)
    write_json("stage_i_measurement_contract.json", measurement)

    vna = {
        "contract_type": "StageIVNAContract",
        "instrument_type": "VNA",
        "measurement_status": "NOT_MEASURED",
        "capability_frequency_Hz": [500.0, 67e9],
        "capability_only_not_measurement": True,
        "expected_files": ["RX_S11.s1p", "TX_S11.s1p", "TX_RX_S21.s2p"],
        "parser": "streamer_rf.fullwave.receiver.read_touchstone",
        "frequency_unit": "READ_FROM_TOUCHSTONE_AND_CONVERTED_TO_Hz",
        "reference_impedance_ohm": NOT_PROVIDED,
        "complex_parameters": ["S11", "S21"],
        "calibration_type": NOT_PROVIDED,
        "calibration_reference_plane": NOT_PROVIDED,
        "IF_bandwidth_Hz": NOT_PROVIDED,
        "source_power_dBm": NOT_PROVIDED,
        "averaging": NOT_PROVIDED,
        "cable_adapter_identifiers": NOT_PROVIDED,
        "raw_smoothing": False,
    }
    validate_vna_contract(vna)
    write_json("stage_i_vna_contract.json", vna)

    scope = {
        "contract_type": "StageIOscilloscopeContract",
        "instrument_type": "OSCILLOSCOPE",
        "measurement_status": "NOT_MEASURED",
        "capability_bandwidth_Hz": 8e9,
        "capability_sample_rate_Hz": 80e9,
        "capability_vertical_bits": 12,
        "capability_only_not_measurement": True,
        "raw_waveform_columns": ["time_s", "voltage_V", "channel_id"],
        "required_waveform_metadata": [
            "sample_interval_s",
            "sample_rate_Hz",
            "vertical_scale_V_per_div",
            "input_impedance_ohm",
            "bandwidth_limit_Hz",
            "trigger_mode",
            "trigger_source",
            "pretrigger_fraction",
            "probe_attenuator_metadata",
            "record_length",
            "experiment_repetition_index",
        ],
        "primary_values_normalized": False,
        "raw_waveform_path": NOT_PROVIDED,
    }
    validate_scope_contract(scope)
    write_json("stage_i_scope_contract.json", scope)

    analyzer = {
        "contract_type": "StageISpectrumAnalyzerContract",
        "instrument_type": "SPECTRUM_ANALYZER",
        "measurement_status": "NOT_MEASURED",
        "capability_frequency_Hz": [2.0, 67e9],
        "capability_only_not_measurement": True,
        "trace_format": NOT_PROVIDED,
        "resolution_bandwidth_Hz": NOT_PROVIDED,
        "video_bandwidth_Hz": NOT_PROVIDED,
        "detector": NOT_PROVIDED,
        "reference_level": NOT_PROVIDED,
    }
    write_json("stage_i_spectrum_analyzer_contract.json", analyzer)

    reference_planes = {
        "contract_type": "StageIReferencePlaneContract",
        "planes": {
            "SOURCE_PORT": "PHYSICS_DERIVED_SOURCE_TERMINAL",
            "TX_FEED": "TRANSMITTER_FEED_OR_CALIBRATED_CABLE_END",
            "FREE_SPACE_REFERENCE": "DECLARED_FIELD_LOCATION_AND_POLARIZATION",
            "RX_FEED": "RECEIVER_FEED_OR_CALIBRATED_CABLE_END",
            "INSTRUMENT_INPUT": "VNA_SCOPE_OR_ANALYZER_CONNECTOR_PLANE",
        },
        "ordered_plane_ids": list(REFERENCE_PLANES),
        "direct_absolute_comparison_requires_same_plane": True,
        "cross_plane_comparison_requires_traceable_transfer": True,
        "incompatible_absolute_comparison": "REJECT",
        "cable_connector_states": [
            "CALIBRATED_OUT",
            "MEASURED_TRANSFER_AVAILABLE",
            "MODELED",
            "UNRESOLVED",
        ],
        "automatic_cable_loss_removal": False,
        "deembedding_requirements": [
            "PARENT_RAW_HASH",
            "COMPLEX_TRANSFER_FILE_HASH",
            "OPERATION_PARAMETERS",
            "INPUT_OUTPUT_REFERENCE_PLANES",
            "REVERSIBLE_OPERATION_RECORD",
        ],
    }
    write_json("stage_i_reference_plane_contract.json", reference_planes)

    profile_350 = {
        "VALIDATION_PROFILE": "SYSTEM_350MHZ",
        "status": "NOT_MEASURED",
        "analysis_band_Hz": [200e6, 500e6],
        "target_center_Hz": 350e6,
        "simulation_target": "H3_REFERENCE_RECEIVED_SPECTRUM",
        "native_350MHz_status": "NOT_RESOLVED",
        "required_future_measurements": [
            "ACTUAL_TX_GEOMETRY",
            "ACTUAL_RX_GEOMETRY",
            "VNA_S11_S21",
            "RECEIVED_OSCILLOSCOPE_WAVEFORM",
            "DISTANCE",
            "ORIENTATION",
            "REPETITION_COUNT",
            "MEASUREMENT_METADATA",
        ],
        "provided_values": NOT_PROVIDED,
    }
    write_json("stage_i_350mhz_profile.json", profile_350)

    profile_native = {
        "VALIDATION_PROFILE": "NATIVE_GHZ",
        "status": "NOT_MEASURED",
        "Stage4_trusted_frequency_mask_Hz": list(STAGE4_TRUSTED_BAND_HZ),
        "Stage5_trusted_frequency_mask_Hz": list(STAGE5_TRUSTED_BAND_HZ),
        "Stage5_status": "FULL_MAXWELL_REFERENCE_PENDING",
        "native_350MHz_status": "NOT_RESOLVED",
        "trust_interpolation": "FORBIDDEN",
        "simulation_targets": [
            "fullwave/h4/h4_stage4_native_received_spectrum.csv",
            "fullwave/h4/h4_stage5_native_received_spectrum.csv",
        ],
        "provided_values": NOT_PROVIDED,
    }
    write_json("stage_i_native_ghz_profile.json", profile_native)

    layers = {
        "contract_type": "StageIDataLayerContract",
        "layers": {
            "RAW": {
                "immutable": True,
                "parent_hash": None,
                "operations": [],
                "meaning": "INSTRUMENT_EXPORT_EXACTLY_AS_RECORDED",
            },
            "CALIBRATED": {
                "immutable_parent_required": True,
                "parent_hash": "REQUIRED_SHA256",
                "operations": ["TRACEABLE_CALIBRATION_OR_DEEMBEDDING_ONLY"],
            },
            "DERIVED": {
                "parent_hash": "REQUIRED_SHA256",
                "operations": ["FFT_PSD_ESD_BAND_ENERGY_FEATURE_METRICS"],
            },
        },
        "raw_overwrite_policy": "FORBIDDEN",
        "derived_parent_hash_required": True,
    }
    validate_data_layer({"layer": "RAW", **layers["layers"]["RAW"]})
    write_json("stage_i_data_layer_contract.json", layers)

    processing = {
        "contract_type": "StageISpectralProcessingPolicy",
        "raw_data_modified": False,
        "time_window": "PRESERVE_ACQUIRED_INTERVAL_AND_RECORD_ANY_DERIVED_SUBWINDOW",
        "DC_policy": "RAW_RETAINED;MEAN_REMOVAL_ONLY_AS_RECORDED_DERIVED_OPERATION",
        "pairwise_policy": "SIMULATION_AND_MEASUREMENT_USE_IDENTICAL_PROCESSING_FOR_COMPARISON",
        "frequency_unit": "Hz",
        "profiles": {
            "SYSTEM_350MHZ": {
                "normalization": "H3_COMPLEX_SINGLE_SIDED_PEAK_AMPLITUDE_2_OVER_N",
                "window": "RECTANGULAR_UNLESS_A_RECORDED_COMMON_ALTERNATIVE_IS_USED",
                "zero_padding": "LINEAR_CONVOLUTION_OR_INTERPOLATION_ONLY_NO_NEW_INFORMATION",
            },
            "NATIVE_GHZ": {
                "normalization": "STAGE_F_DT_TIMES_RFFT_WITH_ABSOLUTE_TIME_PHASE",
                "window": "HANN",
                "sidedness": "ONE_SIDED",
                "ESD_policy": "STAGE_F_ONE_SIDED_ESD_WITH_WINDOW_ENERGY_CORRECTION",
            },
        },
        "arbitrary_filtering": "FORBIDDEN",
        "bandpass_metadata_required": ["filter_type", "pass_band_Hz", "order", "phase_behavior"],
        "alignment_modes": {
            "ABSOLUTE_TRIGGER_TIME": "ABSOLUTE_TIME_ELIGIBLE",
            "PROPAGATION_CORRECTED_TIME": "ABSOLUTE_TIME_ELIGIBLE_WITH_RECORDED_DELAY",
            "FEATURE_ALIGNED_FOR_SHAPE_ONLY": "SHAPE_COMPARISON_ONLY",
        },
        "arbitrary_time_shift": "FORBIDDEN",
        "amplitude_classes": {
            "ABSOLUTE_AMPLITUDE_VALID": "MATCHED_PLANE_LOAD_GEOMETRY_AND_TRACEABLE_CABLE_CALIBRATION",
            "RELATIVE_AMPLITUDE_ONLY": "REFERENCE_PLANE_MATCHED_BUT_CABLE_OR_SCALE_NOT_FULLY_CALIBRATED",
            "NORMALIZED_SHAPE_ONLY": "HARDWARE_OR_LOAD_CONFIGURATION_DIFFERS",
            "NOT_COMPARABLE": "REFERENCE_PLANES_INCOMPATIBLE_OR_PROVENANCE_INSUFFICIENT",
        },
        "universal_acceptance_thresholds": None,
    }
    write_json("stage_i_processing_policy.json", processing)

    uncertainty = {
        "contract_type": "StageIUncertaintyContract",
        "components": {
            key: NOT_PROVIDED
            for key in (
                "repeatability",
                "instrument_uncertainty",
                "geometry_tolerance",
                "distance_uncertainty",
                "orientation_uncertainty",
                "calibration_uncertainty",
                "sampling_uncertainty",
            )
        },
        "supported_ensemble_statistics": [
            "MEAN",
            "MEDIAN",
            "STANDARD_DEVIATION",
            "CONFIDENCE_INTERVAL",
            "COEFFICIENT_OF_VARIATION",
        ],
        "raw_repetitions_averaged_before_storage": False,
        "individual_event_hash_required": True,
    }
    write_json("stage_i_uncertainty_contract.json", uncertainty)

    mapping_rows = [
        ["VNA_RX_TX_S11", "H2_H3_TX_RX_INPUT_RESPONSE", "SYSTEM_350MHZ", "COMPLEX_S_PARAMETER"],
        ["VNA_TX_RX_S21", "H2_H3_REFERENCE_TRANSFER", "SYSTEM_350MHZ", "COMPLEX_S_PARAMETER"],
        ["RECEIVED_SPECTRUM_200_500MHZ", "H3_REFERENCE_RECEIVED_SPECTRUM", "SYSTEM_350MHZ", "RECEIVER_VOLTAGE"],
        ["MEASURED_GHZ_NATIVE_SPECTRUM", "H4_TRUSTED_NATIVE_RECEIVED_SPECTRUM", "NATIVE_GHZ", "RECEIVER_VOLTAGE_WITHIN_TRUST_MASK"],
        ["ORIENTATION_RESPONSE", "H2_H4_POLARIZATION_RESPONSE", "PROFILE_SPECIFIC", "NORMALIZED_OR_CALIBRATED_CONTRAST"],
        ["DISTANCE_DEPENDENCE", "APPROPRIATE_FULLWAVE_OR_STAGE_F_FIELD_MODEL", "PROFILE_SPECIFIC", "DISTANCE_LAW"],
    ]
    with (OUT / "stage_i_simulation_measurement_mapping.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["measured_quantity", "simulation_quantity", "validation_profile", "comparison_semantics"])
        writer.writerows(mapping_rows)

    ledger = {
        "ledger_type": "StageIModelDiscrepancyLedger",
        "status": "EMPTY_NO_EXPERIMENTAL_COMPARISON",
        "possible_source_categories": list(DISCREPANCY_SOURCES),
        "automatic_blame_assignment": False,
        "entries": [],
    }
    validate_discrepancy_ledger(ledger)
    write_json("stage_i_discrepancy_ledger.json", ledger)

    theory_root = ROOT / "fullwave/h2/theory_inputs"
    theory_meta = json.loads((theory_root / "H2_THEORY_350MHz_metadata.json").read_text())
    rx = read_touchstone(theory_root / "H2_THEORY_350MHz_RX_S11.s1p")
    tx = read_touchstone(theory_root / "H2_THEORY_350MHz_TX_S11.s1p")
    pair = read_touchstone(theory_root / "H2_THEORY_350MHz_TX_RX_S21.s2p")
    openems = np.genfromtxt(ROOT / "fullwave/h2/h2_350mhz_openems_sparameters.csv", delimiter=",", names=True)
    openems_s11 = openems["S11_real"] + 1j * openems["S11_imag"]
    openems_s21 = openems["S21_real"] + 1j * openems["S21_imag"]
    dry_run = {
        "status": "SOFTWARE_DRY_RUN_PASS",
        "input_provenance": theory_meta["status"],
        "experimental_evidence": False,
        "VNA_gate_satisfied": False,
        "parser_reused": "streamer_rf.fullwave.receiver.read_touchstone",
        "files": {
            path.name: {"sha256": sha256_file(path), "samples": int(data.frequency_Hz.size), "Z0_ohm": data.Z0_ohm}
            for path, data in (
                (theory_root / "H2_THEORY_350MHz_RX_S11.s1p", rx),
                (theory_root / "H2_THEORY_350MHz_TX_S11.s1p", tx),
                (theory_root / "H2_THEORY_350MHz_TX_RX_S21.s2p", pair),
            )
        },
        "frequency_support_Hz": [float(pair.frequency_Hz[0]), float(pair.frequency_Hz[-1])],
        "S11_metric_API": compare_spectra(rx.frequency_Hz, openems_s11, rx.parameters["S11"]),
        "S21_metric_API": compare_spectra(pair.frequency_Hz, openems_s21, pair.parameters["S21"]),
        "reference_plane_logic": {
            "same_plane_example": "TX_FEED",
            "compatible": True,
            "provenance_limited_comparison_class": "NORMALIZED_SHAPE_ONLY",
        },
    }

    status = {
        "STAGE_I_VALIDATION_FRAMEWORK": "PASS",
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "SYSTEM_350MHZ": "NOT_MEASURED",
        "NATIVE_GHZ": "NOT_MEASURED",
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
        "STAGE_H_SCIENTIFIC_VALIDATION": "PENDING_STAGE_I",
        "inherited_debts": [
            "FULL_WAVE_LOADING_MISMATCH_HIGH",
            "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED",
            "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
            "H4_ABSOLUTE_AMPLITUDE_STATUS_NUMERICAL_REFERENCE_ONLY",
            "FULL_MAXWELL_REFERENCE_PENDING",
        ],
        "validation_status_vocabulary": [
            "NOT_MEASURED",
            "DATA_AVAILABLE",
            "CALIBRATION_PENDING",
            "COMPARISON_READY",
            "PARTIALLY_VALIDATED",
            "VALIDATED",
            "NOT_RESOLVED",
            "MODEL_DISCREPANCY",
        ],
        "infrastructure_inventory": inventory,
        "upstream_all_inputs_available": all(item["available"] for item in upstream["required_inputs"]),
        "dry_run": dry_run,
        "physical_comparison_campaign_started": False,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json("stage_i_validation_status.json", status)
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
