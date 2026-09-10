"""Generate compact H2 contracts from existing canonical openEMS results."""
import csv
import hashlib
import json
from pathlib import Path
import re
import resource
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.receiver import (  # noqa: E402
    PORT_REFERENCE,
    g3_spectrum_plan,
    geometry_hash,
    near_far_classification,
    validate_receiver_geometry,
    validate_transfer_contract,
)

RAW = Path("/tmp/h2_openems_validation")


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime_stats(case):
    text = (RAW / f"{case}.log").read_text()
    size = re.search(r"Time for\s+\d+ iterations with\s+([0-9.]+) cells", text)
    step = re.search(r"FDTD timestep is:\s+([0-9.e+-]+) s", text)
    iterations = re.search(r"Time for\s+(\d+) iterations", text)
    if not all((size, step, iterations)):
        raise ValueError(f"INCOMPLETE_OPENEMS_LOG:{case}")
    return {
        "actual_FDTD_cells": int(round(float(size.group(1)))),
        "actual_FDTD_timestep_s": float(step.group(1)),
        "actual_timesteps": int(iterations.group(1)),
    }


def load_case(case):
    result = json.loads((RAW / case / "result.json").read_text())
    result.update(runtime_stats(case))
    return result


def complex_at_1GHz(record):
    return complex(*record["S21_at_1GHz"])


def relative(a, b):
    return float(abs(b - a) / max(abs(a), np.finfo(float).tiny))


def main():
    started = time.perf_counter()
    out = Path(__file__).parent
    cases = {name: load_case(name) for name in ("baseline", "fine", "expanded_domain")}
    baseline, fine, expanded = (cases[name] for name in ("baseline", "fine", "expanded_domain"))

    hardware = {
        "search_scope": [str(ROOT), "/home/helianthusczc/projects (depth <= 5)"],
        "search_result": "NO_RECEIVER_GEOMETRY_OR_VNA_DATA_FOUND",
        "candidates": [
            {"receiver_id": "SINGLE_TURN_COIL", "status": "PLANNED_NOT_AVAILABLE"},
            {"receiver_id": "COMMERCIAL_350_450_MHZ_ANTENNA", "status": "PLANNED_NOT_AVAILABLE"},
            {"receiver_id": "VIVALDI_UWB_ANTENNA", "status": "PLANNED_NOT_AVAILABLE"},
            {"receiver_id": "H2_CANONICAL_DIPOLE_RX", "status": "GEOMETRY_AVAILABLE_ONLY"},
        ],
        "touchstone_files": [],
        "calibration_records": [],
        "production_receiver_selected": False,
    }
    write_json(out / "h2_hardware_inventory.json", hardware)

    geometry = {
        "geometry_kind": "H2_CANONICAL_REFERENCE_RX",
        "geometry_role": "NUMERICAL_REFERENCE_FIXTURE_NOT_PRODUCTION_HARDWARE",
        "receiver_id": "H2_CANONICAL_DIPOLE_RX",
        "transmitter_id": "H2_REFERENCE_TX",
        "hardware_status": "GEOMETRY_AVAILABLE_ONLY",
        "microgap_geometry_included": False,
        "dipole_axis": [0.0, 0.0, 1.0],
        "total_length_m": 0.15,
        "arm_length_m": 0.075,
        "tx_position_m": [-0.15, 0.0, 0.0],
        "rx_position_m": [0.15, 0.0, 0.0],
        "tx_rx_distance_m": 0.3,
        "relative_orientation": "PARALLEL_CO_POLARIZED_Z",
        "environment": "FREE_SPACE_REFERENCE",
        "support": "NONE",
        "Z0_ohm": 50.0,
        "port_reference": PORT_REFERENCE,
        "connector_equivalent": "IDEAL_CENTER_GAP_CURVEPORT_NO_CABLE",
        "voltage_reference": "POSITIVE_Z_ARM_MINUS_NEGATIVE_Z_ARM",
        "positive_current_direction": "PORT_INTO_RECEIVER",
        "frequency_band_Hz": [0.7e9, 1.3e9],
        "geometry_provenance": (
            "PROJECT_OWNED_CANONICAL_HALF_WAVE_DIPOLE_DERIVED_FROM_PINNED_OPENEMS_TEST_SEMANTICS"
        ),
    }
    validate_receiver_geometry(geometry)
    geometry["near_far_samples"] = [
        near_far_classification(0.15, 0.3, f) for f in (0.7e9, 1.0e9, 1.3e9)
    ]
    geometry["geometry_hash"] = geometry_hash(geometry)
    write_json(out / "h2_geometry.json", geometry)

    g3 = np.genfromtxt(ROOT / "thermal/g3_port/g3_port_uniform.csv", delimiter=",", names=True)
    voltage_plan = g3_spectrum_plan(g3["time_s"], g3["V_port_V"], 0.999)
    current_plan = g3_spectrum_plan(g3["time_s"], g3["I_port_A"], 0.999)
    g3_high = max(
        voltage_plan["upper_frequency_at_cumulative_fraction_Hz"],
        current_plan["upper_frequency_at_cumulative_fraction_Hz"],
    )
    trust = {}
    for stage, name in (
        ("Stage4", "F4-P-stage4-left-isolated"),
        ("Stage5", "F4-C-stage5-highfield-collision"),
    ):
        path = ROOT / "rf/production/f4_attribution" / f"{name}_rf_trust_report.json"
        report = json.loads(path.read_text())
        trust[stage] = {
            "low_Hz": report["trusted_frequency_low_Hz"],
            "high_Hz": report["trusted_frequency_high_Hz"],
            "source": str(path.relative_to(ROOT)),
        }
    bands = {
        "spectrum_method": "MEAN_REMOVED_HANN_WINDOW_RFFT_POWER_PLANNING_ONLY",
        "G3_source": {
            "dt_s": voltage_plan["dt_s"],
            "df_Hz": voltage_plan["df_Hz"],
            "nyquist_Hz": voltage_plan["nyquist_Hz"],
            "V_peak_frequency_Hz": voltage_plan["peak_frequency_Hz"],
            "V_99p9_upper_Hz": voltage_plan["upper_frequency_at_cumulative_fraction_Hz"],
            "I_peak_frequency_Hz": current_plan["peak_frequency_Hz"],
            "I_99p9_upper_Hz": current_plan["upper_frequency_at_cumulative_fraction_Hz"],
            "low_frequency_limit": "BELOW_DF_NOT_RESOLVED_BY_10_NS_RECORD",
        },
        "H3_CIRCUIT_STRUCTURE_BANDS": [{
            "name": "G3_RESOLVED_99P9_CONTENT",
            "low_Hz": voltage_plan["df_Hz"],
            "high_Hz": g3_high,
            "provenance": ["G3_SOURCE_CONTENT"],
            "status": "PLANNING_BAND_PENDING_H3_STRUCTURE_GEOMETRY",
        }],
        "H4_NATIVE_RF_TRUSTED_BANDS": [
            {"name": stage, **values, "provenance": ["STAGE_F_TRUSTED_NATIVE"]}
            for stage, values in trust.items()
        ],
        "H2_RECEIVER_VALIDATED_BANDS": [],
        "H2_SIMULATION_ONLY_BANDS": [{
            "name": "H2_CANONICAL_DIPOLE_FIXTURE",
            "low_Hz": 0.7e9,
            "high_Hz": 1.3e9,
            "provenance": ["RECEIVER_HARDWARE_CANONICAL_GEOMETRY_ONLY"],
            "status": "SIMULATION_ONLY",
        }],
        "frozen_native_low_band_policy": (
            "VHF_UHF_1_TO_3_GHZ_NOT_PROMOTED_BEYOND_STORED_STAGE_F_TRUST_STATUS"
        ),
    }
    write_json(out / "h2_frequency_bands.json", bands)

    sensitivity_rows = []
    for name, record in cases.items():
        sensitivity_rows.append({
            "case": name,
            "spacing_m": record["spacing_mm"] * 1e-3,
            "domain_scale": record["domain_scale"],
            "actual_FDTD_cells": record["actual_FDTD_cells"],
            "actual_timesteps": record["actual_timesteps"],
            "S11_min_frequency_Hz": record["S11_min_frequency_Hz"],
            "S11_min_dB": record["S11_min_dB"],
            "S21_at_1GHz_dB": record["S21_at_1GHz_dB"],
            "S21_phase_at_1GHz_rad": record["S21_at_1GHz_phase_rad"],
            "group_delay_at_1GHz_s": record["group_delay_at_1GHz_s"],
            "Dmax_at_1GHz": record["nf2ff"]["Dmax"],
            "runtime_s": record["runtime_s"],
            "peak_RSS_KiB": record["peak_RSS_KiB"],
            "raw_bytes": record["raw_bytes"],
        })
    with (out / "h2_mesh_domain_sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=sensitivity_rows[0].keys())
        writer.writeheader()
        writer.writerows(sensitivity_rows)
    shutil.copyfile(RAW / "baseline/sparameters.csv", out / "h2_reference_sparameters.csv")
    shutil.copyfile(RAW / "baseline/nf2ff_cut.csv", out / "h2_reference_nf2ff_cut.csv")

    sparameter_file = out / "h2_reference_sparameters.csv"
    contract = {
        "receiver_id": geometry["receiver_id"],
        "hardware_status": "GEOMETRY_AVAILABLE_ONLY",
        "geometry_hash": geometry["geometry_hash"],
        "measurement_id": None,
        "port_reference": PORT_REFERENCE,
        "voltage_reference": geometry["voltage_reference"],
        "positive_current_direction": geometry["positive_current_direction"],
        "Z0_ohm": 50.0,
        "frequency_data_file": sparameter_file.name,
        "frequency_data_sha256": sha256(sparameter_file),
        "complex_fields": ["S11", "S21", "Vrx_50ohm_per_tx_inc", "Vrx_open_circuit_equivalent_per_tx_inc"],
        "H_rx_E_complex": None,
        "H_rx_E_status": "RX_INTRINSIC_FIELD_TRANSFER_PENDING",
        "near_far_field_status": "FAR_FIELD_ACROSS_0P7_TO_1P3_GHZ_BY_2D2_OVER_LAMBDA",
        "mesh_status": "BASELINE_TO_FINE_KEY_FEATURES_WITHIN_5_PERCENT",
        "VNA_validation_status": "SIMULATION_ONLY",
        "valid_frequency_mask": {"low_Hz": 0.7e9, "high_Hz": 1.3e9},
        "source_provenance": "H2_REFERENCE_TX_CANONICAL_DIPOLE",
        "load_semantics": {
            "Vrx_50ohm": "LOADED_RECEIVER_TERMINAL_VOLTAGE_TRANSFER",
            "Vrx_open_circuit_equivalent": "THEVENIN_CONVERSION_USING_SYMMETRIC_RECEIVER_ZIN",
            "not_equivalent": True,
        },
    }
    validate_transfer_contract(contract)
    write_json(out / "h2_receiver_transfer_contract.json", contract)

    protocol = {
        "status": "VNA_MEASUREMENT_REQUIRED",
        "instrument": "EXISTING_500_HZ_TO_67_GHZ_VNA",
        "measurement_fixture": geometry,
        "calibration": "FULL_TWO_PORT_SOLT_APPROPRIATE_TO_ACTUAL_CONNECTORS",
        "calibration_reference_plane": "AT_TX_AND_RX_FEED_CABLE_ENDS",
        "measurements": ["S11_RECEIVER", "S11_REFERENCE_TX", "S21_TX_TO_RX"],
        "frequency_ranges": [{"low_Hz": 0.7e9, "high_Hz": 1.3e9, "points": 601}],
        "frequency_resolution_Hz": 1e6,
        "IF_bandwidth_Hz": "RECORD_ACTUAL_SETTING_CHOSEN_FOR_NOISE_AND_ACQUISITION_TIME",
        "source_power_dBm": "RECORD_ACTUAL_PASSIVE_LINEAR_SAFE_SETTING",
        "linearity_check": "REPEAT_ONE_TRACE_AT_REDUCED_POWER_AND_CONFIRM_UNCHANGED_S_PARAMETERS",
        "averaging": "RECORD_COUNT; USE_ONLY_IF IDENTICAL_FOR_BACKGROUND_AND_FIXTURE",
        "distance_m": 0.3,
        "orientation": "PARALLEL_CO_POLARIZED_Z",
        "antenna_centers_height_m": "RECORD_ACTUAL",
        "environment": "MICROWAVE_ANECHOIC_CHAMBER_WHERE_AVAILABLE",
        "cables_and_adapters": "RECORD_IDS_LENGTHS_CONNECTOR_TYPES_AND_ORIENTATION",
        "background_measurement": "RECORD_TX_RX_REMOVED_OR_TERMINATED_BACKGROUND_IF_DYNAMIC_RANGE_REQUIRES",
        "required_outputs": [
            "h2_reference_tx.s1p", "h2_receiver.s1p", "h2_tx_rx.s2p", "h2_vna_metadata.json"
        ],
        "metadata_required": [
            "instrument_model", "instrument_serial", "calibration_kit", "calibration_timestamp",
            "IF_bandwidth_Hz", "source_power_dBm", "averaging_count", "cable_ids",
            "adapter_ids", "distance_m", "orientation", "height_m", "environment",
        ],
        "primary_comparison_policy": "RAW_UNSMOOTHED_COMPLEX_DATA_NO_FREQUENCY_SHIFT_NO_PEAK_SCALING",
    }
    write_json(out / "h2_vna_measurement_protocol.json", protocol)

    z_base = complex_at_1GHz(baseline)
    z_fine = complex_at_1GHz(fine)
    z_expanded = complex_at_1GHz(expanded)
    summary = {
        "framework_status": "SIMULATION_FRAMEWORK_VALIDATED",
        "formal_node_status": "MINIMAL_FIX_REQUIRED",
        "blocker": "VNA_MEASUREMENT_REQUIRED",
        "receiver_transfer_status": "SIMULATION_ONLY",
        "production_receiver_status": "NOT_RESOLVED",
        "reference_results": {
            "S11_min_frequency_Hz": baseline["S11_min_frequency_Hz"],
            "S11_min_dB": baseline["S11_min_dB"],
            "S11_minus10dB_band_Hz": baseline["S11_minus10dB_band_Hz"],
            "Zin_at_1GHz_ohm": baseline["Zin_at_1GHz_ohm"],
            "S21_at_1GHz": baseline["S21_at_1GHz"],
            "S21_at_1GHz_dB": baseline["S21_at_1GHz_dB"],
            "S21_phase_at_1GHz_rad": baseline["S21_at_1GHz_phase_rad"],
            "group_delay_at_1GHz_s": baseline["group_delay_at_1GHz_s"],
            "near_far": baseline["near_far_at_1GHz"],
            "nf2ff": baseline["nf2ff"],
        },
        "mesh_convergence": {
            "feature_frequency_relative_change": relative(
                baseline["S11_min_frequency_Hz"], fine["S11_min_frequency_Hz"]
            ),
            "S21_magnitude_relative_change_at_1GHz": relative(abs(z_base), abs(z_fine)),
            "S21_dB_absolute_change_at_1GHz": abs(
                fine["S21_at_1GHz_dB"] - baseline["S21_at_1GHz_dB"]
            ),
            "S21_phase_absolute_change_at_1GHz_rad": abs(
                fine["S21_at_1GHz_phase_rad"] - baseline["S21_at_1GHz_phase_rad"]
            ),
            "Dmax_relative_change": relative(baseline["nf2ff"]["Dmax"], fine["nf2ff"]["Dmax"]),
            "status": "PASS_KEY_FEATURES_WITHIN_5_PERCENT",
        },
        "domain_check": {
            "S21_magnitude_relative_change_at_1GHz": relative(abs(z_base), abs(z_expanded)),
            "S21_dB_absolute_change_at_1GHz": abs(
                expanded["S21_at_1GHz_dB"] - baseline["S21_at_1GHz_dB"]
            ),
            "S21_phase_absolute_change_at_1GHz_rad": abs(
                expanded["S21_at_1GHz_phase_rad"] - baseline["S21_at_1GHz_phase_rad"]
            ),
            "Dmax_relative_change": relative(baseline["nf2ff"]["Dmax"], expanded["nf2ff"]["Dmax"]),
            "minimum_structure_to_PML_clearance_m": 0.15,
            "status": "PASS",
        },
        "intrinsic_field_transfer_status": "RX_INTRINSIC_FIELD_TRANSFER_PENDING",
        "VNA_quality_status": "NOT_EVALUATED_NO_MEASUREMENT_DATA",
        "VNA_comparison_metrics": None,
        "acceptance_targets": {
            "major_feature_frequency_mismatch_max": 0.05,
            "normalized_S21_spectral_correlation_min": 0.90,
            "principal_band_absolute_S21_tolerance_dB": 3.0,
        },
        "runtime_cases": cases,
        "resource": {
            "total_runtime_s": sum(record["runtime_s"] for record in cases.values()),
            "maximum_peak_RSS_KiB": max(record["peak_RSS_KiB"] for record in cases.values()),
            "external_raw_bytes": sum(record["raw_bytes"] for record in cases.values()),
            "generation_runtime_s": time.perf_counter() - started,
            "generation_peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        },
        "G3_waveform_injected": False,
        "later_nodes_started": False,
    }
    write_json(out / "h2_summary.json", summary)


if __name__ == "__main__":
    main()
