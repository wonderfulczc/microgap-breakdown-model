#!/usr/bin/env python3
"""Integrate frozen H3/H4 evidence without combining their physical signals."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.evidence import (  # noqa: E402
    H3_DEVELOPMENT_BAND_HZ,
    SIGNAL_LEVEL_DEFINITIONS,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    STAGE_I_REQUIRED_INPUT_IDS,
    TRUST_HIERARCHY,
    common_spectral_metrics,
    validate_stage_h_transfer_contract,
    validate_stage_i_requirements,
)


OUT = ROOT / "fullwave/h5"


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text())


def complex_columns(data, real_name, imag_name):
    return np.asarray(data[real_name]) + 1j * np.asarray(data[imag_name])


def write_dict_rows(path: Path, rows: list[dict]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def calculate_path_metrics(pathway: str, path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if pathway == "POST_BREAKDOWN_G3":
        frequency = data["frequency_Hz"]
        source = complex_columns(data, "V_G3_real_V", "V_G3_imag_V")
        receiver = complex_columns(data, "V_rx_real_V", "V_rx_imag_V")
        source_unit, receiver_unit = "V", "V"
    else:
        valid = data["trust_valid"].astype(bool)
        frequency = data["frequency_Hz"][valid]
        source = complex_columns(data, "E_parallel_real_V_s_m", "E_parallel_imag_V_s_m")[valid]
        receiver = complex_columns(data, "V_rx_real_V_s", "V_rx_imag_V_s")[valid]
        source_unit, receiver_unit = "V_s_per_m", "V_s"
    metrics = common_spectral_metrics(
        frequency, source, receiver, source_unit=source_unit, receiver_unit=receiver_unit
    )
    shape_rows = [
        {
            "pathway": pathway,
            "frequency_Hz": float(frequency[index]),
            "source_normalized_shape": float(metrics["source_normalized_shape"][index]),
            "receiver_normalized_shape": float(metrics["receiver_normalized_shape"][index]),
            "source_unit": source_unit,
            "receiver_unit": receiver_unit,
            "absolute_cross_path_comparison_permitted": False,
        }
        for index in range(frequency.size)
    ]
    del metrics["source_normalized_shape"], metrics["receiver_normalized_shape"]
    return metrics, shape_rows


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    h3_contract_path = ROOT / "fullwave/h3/h3_result_contract.json"
    h3_summary_path = ROOT / "fullwave/h3/h3_summary.json"
    h4_contract_path = ROOT / "fullwave/h4/h4_result_contract.json"
    h4_summary_path = ROOT / "fullwave/h4/h4_summary.json"
    h3 = load_json(h3_contract_path)
    h3_summary = load_json(h3_summary_path)
    h4 = load_json(h4_contract_path)
    h4_summary = load_json(h4_summary_path)

    pathways = [
        {
            "pathway": "NATIVE_STAGE4",
            "source_quantity": "SIGNAL_0_NATIVE:TOTAL_STAGE_F_E_B_OBSERVER_FIELD",
            "source_physical_mechanism": "COLD_DISCHARGE_NATIVE_CHARGE_CURRENT_EVOLUTION",
            "source_frequency_support_Hz": f"{STAGE4_TRUSTED_BAND_HZ[0]}:{STAGE4_TRUSTED_BAND_HZ[1]}",
            "propagation_model": "STAGE_F_RETARDED_JEFIMENKO_TO_OBSERVER",
            "structure_model": "NO_ADDITIONAL_TX_PROPAGATION",
            "receiver_model": h4["receiver_id"],
            "receiver_output": "PRIMARY_NATIVE_RF_OUTPUT:STAGE4_RECEIVED_SPECTRUM",
            "trust_status": "TRUSTED_PHYSICS_SOURCE__NUMERICAL_REFERENCE_ONLY_SIGNAL_2",
            "numerical_limitation": "SPARSE_5_TRUSTED_FFT_BINS;RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
            "experimental_validation_status": "EXPERIMENTAL_VALIDATION_PENDING",
        },
        {
            "pathway": "NATIVE_STAGE5",
            "source_quantity": "SIGNAL_0_NATIVE:TOTAL_STAGE_F_E_B_OBSERVER_FIELD",
            "source_physical_mechanism": "COLLISION_NATIVE_CHARGE_CURRENT_EVOLUTION",
            "source_frequency_support_Hz": f"{STAGE5_TRUSTED_BAND_HZ[0]}:{STAGE5_TRUSTED_BAND_HZ[1]}",
            "propagation_model": "STAGE_F_RETARDED_JEFIMENKO_TO_OBSERVER",
            "structure_model": "NO_ADDITIONAL_TX_PROPAGATION",
            "receiver_model": h4["receiver_id"],
            "receiver_output": "PRIMARY_NATIVE_RF_OUTPUT:STAGE5_RECEIVED_SPECTRUM",
            "trust_status": "TRUSTED_PHYSICS_SOURCE__NUMERICAL_REFERENCE_ONLY_SIGNAL_2",
            "numerical_limitation": "SPARSE_7_TRUSTED_FFT_BINS;FULL_MAXWELL_REFERENCE_PENDING;RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
            "experimental_validation_status": "EXPERIMENTAL_VALIDATION_PENDING",
        },
        {
            "pathway": "POST_BREAKDOWN_G3",
            "source_quantity": "SIGNAL_0_PORT:G3_V_PORT_I_PORT_TRANSIENT",
            "source_physical_mechanism": "THERMAL_CHANNEL_RLC_PORT_TRANSIENT",
            "source_frequency_support_Hz": f"{H3_DEVELOPMENT_BAND_HZ[0]}:{H3_DEVELOPMENT_BAND_HZ[1]}",
            "propagation_model": "H3_PASSIVE_FULLWAVE_REFERENCE_STRUCTURE_AND_FREE_SPACE",
            "structure_model": h3["reference_structure_id"],
            "receiver_model": h3["receiver_id"],
            "receiver_output": "REFERENCE_RECEIVED_VOLTAGE",
            "trust_status": "DEVELOPMENT_VERIFIED",
            "numerical_limitation": "FULL_WAVE_LOADING_MISMATCH_HIGH;FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED;G3_10NS_WINDOW",
            "experimental_validation_status": "EXPERIMENTAL_VALIDATION_PENDING",
        },
    ]
    write_dict_rows(OUT / "h5_pathway_matrix.csv", pathways)

    frequency_rows = [
        {
            "pathway": "POST_BREAKDOWN_G3",
            "f_low_Hz": 200e6,
            "f_high_Hz": 500e6,
            "center_Hz": 350e6,
            "status": "DEVELOPMENT_VERIFIED",
            "trust_interpolation": "FORBIDDEN",
        },
        {
            "pathway": "NATIVE_STAGE4",
            "f_low_Hz": STAGE4_TRUSTED_BAND_HZ[0],
            "f_high_Hz": STAGE4_TRUSTED_BAND_HZ[1],
            "center_Hz": "",
            "status": "TRUSTED_PHYSICS_WITHIN_RFTRUSTREPORT_MASK",
            "trust_interpolation": "FORBIDDEN",
        },
        {
            "pathway": "NATIVE_STAGE5",
            "f_low_Hz": STAGE5_TRUSTED_BAND_HZ[0],
            "f_high_Hz": STAGE5_TRUSTED_BAND_HZ[1],
            "center_Hz": "",
            "status": "TRUSTED_PHYSICS_WITHIN_RFTRUSTREPORT_MASK",
            "trust_interpolation": "FORBIDDEN",
        },
        {
            "pathway": "NATIVE_200_500MHZ",
            "f_low_Hz": 200e6,
            "f_high_Hz": 500e6,
            "center_Hz": 350e6,
            "status": "NOT_RESOLVED",
            "trust_interpolation": "FORBIDDEN",
        },
        {
            "pathway": "UNSUPPORTED_FREQUENCY_GAP",
            "f_low_Hz": 500e6,
            "f_high_Hz": STAGE4_TRUSTED_BAND_HZ[0],
            "center_Hz": "",
            "status": "NO_CROSS_PATH_TRUST_CLAIM",
            "trust_interpolation": "FORBIDDEN",
        },
    ]
    write_dict_rows(OUT / "h5_frequency_support_matrix.csv", frequency_rows)

    metric_sources = {
        "POST_BREAKDOWN_G3": ROOT / "fullwave/h3/h3_received_spectrum.csv",
        "NATIVE_STAGE4": ROOT / "fullwave/h4/h4_stage4_native_received_spectrum.csv",
        "NATIVE_STAGE5": ROOT / "fullwave/h4/h4_stage5_native_received_spectrum.csv",
    }
    metrics = {}
    shape_rows = []
    for pathway, path in metric_sources.items():
        metrics[pathway], rows = calculate_path_metrics(pathway, path)
        shape_rows.extend(rows)
    metric_rows = [{"pathway": pathway, **values} for pathway, values in metrics.items()]
    write_dict_rows(OUT / "h5_common_spectral_metrics.csv", metric_rows)
    write_dict_rows(OUT / "h5_figure_spectral_shapes.csv", shape_rows)

    stage_i = {
        "status": "PENDING_STAGE_I",
        "fabricated_inputs": False,
        "required_inputs": [
            {
                "input_id": input_id,
                "available": False,
                "required_provenance": "MEASURED_OR_PHYSICALLY_SURVEYED_WITH_TRACEABLE_METADATA",
            }
            for input_id in STAGE_I_REQUIRED_INPUT_IDS
        ],
        "comparison_targets": {
            "200_500MHz": "MEASURED_RECEIVER_VOLTAGE_SPECTRUM_VS_H3_REFERENCE_SPECTRUM",
            "trusted_native_GHz": "MEASURED_GHZ_RECEIVER_SPECTRUM_VS_H4_TRUSTED_BAND_PREDICTION_WHERE_HARDWARE_EXISTS",
            "metrics": [
                "FEATURE_FREQUENCY",
                "SPECTRAL_SHAPE",
                "CALIBRATED_AMPLITUDE",
                "REFERENCE_PLANE_COMPATIBLE_PHASE",
                "DISTANCE_DEPENDENCE",
                "ORIENTATION_POLARIZATION_DEPENDENCE",
            ],
        },
    }
    validate_stage_i_requirements(stage_i)
    (OUT / "h5_stage_i_required_inputs.json").write_text(json.dumps(stage_i, indent=2) + "\n")

    contract = {
        "contract_type": "StageHTransferContract",
        "signal_level_definition": SIGNAL_LEVEL_DEFINITIONS,
        "native_paths": {
            stage: {
                "source_contract": h4["source_contracts"][stage],
                "source_file_hash": h4["source_contracts"][stage]["source_hash"],
                "propagation_model": "STAGE_F_RETARDED_JEFIMENKO_TO_OBSERVER",
                "NATIVE_PROPAGATION_ALREADY_INCLUDED": True,
                "structure_action": "NO_SECOND_PROPAGATION",
                "receiver_model": h4["receiver_id"],
                "transfer_function": "H_rx_E_COMPLEX",
                "receiver_output": h4["native_received_spectrum_paths"][stage],
                "frequency_domain_role": "PRIMARY_NATIVE_RF_OUTPUT",
                "time_domain_role": "SECONDARY_RECONSTRUCTION_ONLY",
                "trusted_bin_count": h4_summary["stages"][stage]["trusted_bin_count"],
                "absolute_amplitude_status": h4["H4_ABSOLUTE_AMPLITUDE_STATUS"],
            }
            for stage in ("Stage4", "Stage5")
        },
        "post_breakdown_path": {
            "source_case_id": h3["source_case_id"],
            "G3_input_hash": h3["G3_input_hash"],
            "source_file": "thermal/g3_port/g3_port_uniform.csv",
            "source_quantity": "V_G3_COMPLEX_WITH_I_G3_LOADING_DIAGNOSTIC_ONLY",
            "propagation_structure_model": h3["reference_structure_id"],
            "receiver_model": h3["receiver_id"],
            "transfer_function": "H_V_COMPLEX_TOTAL_PORT_VOLTAGE_TO_50OHM_RECEIVER",
            "receiver_output": f"fullwave/h3/{h3['V_received_spectrum_file']}",
            "trust_status": "DEVELOPMENT_VERIFIED",
            "time_window_s": h3_summary["time_transform"]["original_duration_s"],
            "intrinsic_Fourier_scale_Hz": h3_summary["time_transform"]["intrinsic_Fourier_scale_Hz"],
            "zero_padding_status": h3_summary["time_transform"]["zero_padding_status"],
        },
        "frequency_support": {
            "POST_BREAKDOWN_G3_Hz": list(H3_DEVELOPMENT_BAND_HZ),
            "SYSTEM_TARGET_CENTER_Hz": 350e6,
            "SYSTEM_MAXIMUM_Hz": 500e6,
            "NATIVE_STAGE4_Hz": list(STAGE4_TRUSTED_BAND_HZ),
            "NATIVE_STAGE5_Hz": list(STAGE5_TRUSTED_BAND_HZ),
            "NATIVE_200_500MHz": "NOT_RESOLVED",
            "UNSUPPORTED_500MHz_TO_STAGE4_LOW": "NO_TRUST_INTERPOLATION",
        },
        "trust_hierarchy": list(TRUST_HIERARCHY),
        "trust_mapping": {
            "Stage_F_native_within_mask": "TRUSTED_PHYSICS",
            "H3_reference_pathway": "DEVELOPMENT_VERIFIED",
            "H4_absolute_receiver_voltage": "NUMERICAL_REFERENCE_ONLY",
            "native_350MHz": "NOT_RESOLVED",
            "experimental_receiver_output": "EXPERIMENTAL_VALIDATION_PENDING",
        },
        "CROSS_PATH_COHERENT_SUMMATION": "NOT_PERMITTED_CURRENT_CONFIGURATION",
        "coherent_summation_requirements": [
            "SAME_PHYSICAL_TX_RX_GEOMETRY",
            "SAME_RECEIVER",
            "SAME_REFERENCE_PLANE",
            "COMPATIBLE_TIME_ORIGIN",
            "COMPATIBLE_FREQUENCY_COVERAGE",
            "COMPATIBLE_PHASE_CONVENTION",
            "VALIDATED_TRANSFER_FUNCTIONS",
        ],
        "CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION": "NOT_RESOLVED",
        "native_350MHz_status": "NOT_RESOLVED",
        "known_debts": [
            "FULL_WAVE_LOADING_MISMATCH_HIGH",
            "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED",
            "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
            "FULL_MAXWELL_REFERENCE_PENDING",
            "SPARSE_STAGE4_5_AND_STAGE5_7_TRUSTED_FFT_BINS",
            "VNA_MEASUREMENT_PENDING",
            "PRODUCTION_TX_RX_GEOMETRY_PENDING",
        ],
        "upstream_hashes": {
            "H3_contract": sha(h3_contract_path),
            "H3_summary": sha(h3_summary_path),
            "H3_received_spectrum": sha(ROOT / "fullwave/h3/h3_received_spectrum.csv"),
            "H4_contract": sha(h4_contract_path),
            "H4_summary": sha(h4_summary_path),
            "H4_receiver_transfer": sha(ROOT / h4["receiver_transfer_path"]),
        },
        "stage_i_required_inputs_file": "fullwave/h5/h5_stage_i_required_inputs.json",
        "STAGE_H_TOOL_DEVELOPMENT": "PASS",
        "STAGE_H_SCIENTIFIC_VALIDATION": "PENDING_STAGE_I",
        "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING": True,
        "STAGE_H_PRODUCTION_GEOMETRY_PENDING": True,
        "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED": True,
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
    }
    validate_stage_h_transfer_contract(contract)
    (OUT / "h5_stage_h_transfer_contract.json").write_text(json.dumps(contract, indent=2) + "\n")

    summary = {
        "STAGE_H_TOOL_DEVELOPMENT": "PASS",
        "STAGE_H_SCIENTIFIC_VALIDATION": "PENDING_STAGE_I",
        "CROSS_PATH_COHERENT_SUMMATION": contract["CROSS_PATH_COHERENT_SUMMATION"],
        "CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION": "NOT_RESOLVED",
        "system_350MHz_interpretation": (
            "SUPPORTED_BY_G3_H3_POST_BREAKDOWN_REFERENCE_PATH_NATIVE_STAGE_F_CONTRIBUTION_NOT_RESOLVED"
        ),
        "native_GHz_interpretation": (
            "RECEIVER_SELECTION_ALTERS_SHAPE_AND_CENTROID_WITHOUT_CHANGING_STAGE_F_TRUST_MASKS"
        ),
        "common_spectral_metrics": metrics,
        "simulation_rerun": False,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUT / "h5_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
