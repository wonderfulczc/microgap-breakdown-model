#!/usr/bin/env python3
"""Generate the WP-I-D synthetic multi-condition pipeline dry run."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import time

import numpy as np

from streamer_rf.validation.stage_i import sha256_file
from streamer_rf.validation.wp_i_b import event_metrics, read_scope_waveform, scope_quality_gate
from streamer_rf.validation.wp_i_d import (
    SYNTHETIC_INPUT,
    condition_statistics,
    cosine_floor_fit,
    distance_frequency_stability,
    merge_condition_records,
    monte_carlo_power_law_uncertainty,
    polarization_contrast_db,
    power_law_fit,
    verify_supplement_bundle,
)


ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "validation/stage_i/wp_i_d"
SUPPLEMENT = ROOT / "validation/stage_i/system_350mhz/wp_i_d_synthetic"
WP_B_INPUT = ROOT / "validation/stage_i/system_350mhz/synthetic_input"
WP_B_METRICS = ROOT / "validation/stage_i/system_350mhz/wp_i_b_event_metrics.csv"
WP_B_STATUS = ROOT / "validation/stage_i/system_350mhz/wp_i_b_summary.json"


def write_json(name, record):
    (OUTPUT / name).write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def write_csv(name, rows):
    if not rows:
        raise ValueError("EMPTY_WP_I_D_OUTPUT")
    with (OUTPUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_wp_b_rows():
    with WP_B_METRICS.open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    rows = []
    numeric_fields = [name for name in raw[0] if name not in ("raw_data_path", "raw_data_hash")]
    for record in raw:
        converted = {name: float(record[name]) for name in numeric_fields}
        converted["repetition_index"] = int(converted["repetition_index"])
        converted.update(
            {
                "orientation_deg": 0.0,
                "source_bundle": "WP_I_B_SYNTHETIC_FROZEN",
                "group": "distance_sweep",
                "raw_data_path": str(Path("validation/stage_i/system_350mhz/synthetic_input") / record["raw_data_path"]),
                "raw_data_hash": record["raw_data_hash"],
            }
        )
        rows.append(converted)
    return rows


def flatten_supplement_event(record, metrics):
    row = {
        "distance_m": float(record["distance_m"]),
        "orientation_deg": float(record["orientation_deg"]),
        "repetition_index": int(record["repetition_index"]),
        "source_bundle": "WP_I_D_SYNTHETIC_SUPPLEMENT",
        "group": record["group"],
        "raw_data_path": str(Path("validation/stage_i/system_350mhz/wp_i_d_synthetic") / record["file"]),
        "raw_data_hash": sha256_file(SUPPLEMENT / record["file"]),
        "peak_abs_voltage_V": metrics["peak_abs_voltage_V"],
        "rms_voltage_V": metrics["rms_voltage_V"],
        "noise_rms_V": metrics["noise_rms_V"],
        "peak_to_noise_SNR_dB": metrics["peak_to_noise_SNR_dB"],
    }
    mapping = {
        "FULL_50_500MHZ": "full_50_500mhz",
        "LOW_50_100MHZ": "low_50_100mhz",
        "H3_OVERLAP_200_500MHZ": "h3_overlap_200_500mhz",
    }
    for source, prefix in mapping.items():
        row[f"{prefix}_peak_frequency_Hz"] = metrics[source]["peak_frequency_Hz"]
        row[f"{prefix}_spectral_centroid_Hz"] = metrics[source]["spectral_centroid_Hz"]
        row[f"{prefix}_band_integrated_V2_Hz"] = metrics[source]["band_integrated_V2_Hz"]
    return row


def stats_output_row(stats, comparison_status):
    return {**stats, "H3_comparison_status": comparison_status, "evidence_type": "SYNTHETIC_DRY_RUN"}


def main():
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    supplement = verify_supplement_bundle(SUPPLEMENT)
    metadata = supplement["metadata"]
    uncertainty = supplement["uncertainty"]
    wp_b_status = json.loads(WP_B_STATUS.read_text())
    if wp_b_status["WP_I_B_PIPELINE_DRY_RUN"] != "PASS":
        raise ValueError("WP_I_B_NOT_FROZEN_PASS")

    wp_b_rows = read_wp_b_rows()
    if len(wp_b_rows) != 30:
        raise ValueError("WP_I_B_EVENT_COUNT_CHANGED")
    quality_rows = []
    supplemental_rows = []
    quality_metadata = {
        "sample_interval_s": metadata["sample_interval_s"],
        "record_length_samples": metadata["record_length_samples"],
        "record_duration_s": metadata["record_duration_s"],
        "trigger_time_s": metadata["trigger_time_s"],
    }
    background_available = (WP_B_INPUT / "scope_raw/BACKGROUND_NO_TRIGGER.csv").is_file()
    for record in metadata["records"]:
        time_s, voltage = read_scope_waveform(SUPPLEMENT / record["file"])
        quality = scope_quality_gate(time_s, voltage, quality_metadata, baseline_available=background_available)
        if not quality["passed"]:
            raise ValueError(f"WP_I_D_RAW_QUALITY_FAILURE:{record['file']}")
        metrics, _, _ = event_metrics(time_s, voltage, quality_metadata)
        supplemental_rows.append(flatten_supplement_event(record, metrics))
        quality_rows.append({"raw_data_path": record["file"], **quality})

    all_rows = merge_condition_records(wp_b_rows, supplemental_rows)
    if len(all_rows) != 90:
        raise ValueError("WP_I_D_MERGED_EVENT_COUNT_MISMATCH")
    all_rows.sort(key=lambda row: (row["distance_m"], row["orientation_deg"], row["repetition_index"]))
    condition_index = [
        {
            "condition_id": f"D{row['distance_m']:.1f}m_A{int(row['orientation_deg']):02d}",
            "distance_m": row["distance_m"],
            "orientation_deg": row["orientation_deg"],
            "repetition_index": row["repetition_index"],
            "source_bundle": row["source_bundle"],
            "group": row["group"],
            "raw_data_path": row["raw_data_path"],
            "raw_data_hash": row["raw_data_hash"],
            "data_origin": SYNTHETIC_INPUT,
            "scientific_validation_allowed": False,
        }
        for row in all_rows
    ]
    write_csv("wp_i_d_condition_index.csv", condition_index)
    write_csv("wp_i_d_event_metrics.csv", all_rows)

    distances = [0.3, 0.5, 0.6, 0.8, 1.0]
    distance_stats = [
        stats_output_row(
            condition_statistics(all_rows, distance, 0.0),
            "FROZEN_H3_REFERENCE_AVAILABLE" if np.isclose(distance, 1.0) else "NO_DIRECT_H3_FULLWAVE_REFERENCE",
        )
        for distance in distances
    ]
    write_csv("wp_i_d_distance_statistics.csv", distance_stats)

    peak_amplitude = np.asarray([row["peak_voltage_V_mean"] for row in distance_stats])
    rms_amplitude = np.asarray([row["RMS_voltage_V_mean"] for row in distance_stats])
    peak_fit = power_law_fit(distances, peak_amplitude)
    rms_fit = power_law_fit(distances, rms_amplitude)
    stability = distance_frequency_stability(
        distances,
        [row["peak_frequency_Hz_mean"] for row in distance_stats],
        [row["spectral_centroid_Hz_mean"] for row in distance_stats],
    )
    distance_fit = {
        "model": "V(d)=A*d^(-n)",
        "evidence_type": "SYNTHETIC_DRY_RUN",
        "scientific_interpretation_allowed": False,
        "distance_m": distances,
        "peak_voltage_fit": peak_fit,
        "RMS_voltage_fit": rms_fit,
        "frequency_stability": stability,
        "H3_distance_boundary": {
            "direct_reference_distance_m": 1.0,
            "other_distances": "NO_DIRECT_H3_FULLWAVE_REFERENCE",
            "inverse_distance_H3_extrapolation_used": False,
        },
    }
    write_json("wp_i_d_distance_fit.json", distance_fit)

    orientations = [0.0, 30.0, 60.0, 90.0]
    orientation_stats = [
        stats_output_row(
            condition_statistics(all_rows, 0.6, angle), "NO_DIRECT_H3_ORIENTATION_REFERENCE"
        )
        for angle in orientations
    ]
    write_csv("wp_i_d_orientation_statistics.csv", orientation_stats)
    peak_orientation = np.asarray([row["peak_voltage_V_mean"] for row in orientation_stats])
    rms_orientation = np.asarray([row["RMS_voltage_V_mean"] for row in orientation_stats])
    orientation_fit = {
        "model": "sqrt(A_parallel^2*cos(theta)^2+A_floor^2)",
        "orientation_deg": orientations,
        "peak_voltage_fit": cosine_floor_fit(orientations, peak_orientation),
        "RMS_voltage_fit": cosine_floor_fit(orientations, rms_orientation),
        "SYNTHETIC_POLARIZATION_CONTRAST_peak_dB": polarization_contrast_db(
            peak_orientation[0], peak_orientation[-1]
        ),
        "SYNTHETIC_POLARIZATION_CONTRAST_RMS_dB": polarization_contrast_db(
            rms_orientation[0], rms_orientation[-1]
        ),
        "H3_orientation_comparison_status": "NO_DIRECT_H3_ORIENTATION_REFERENCE",
        "evidence_type": "SYNTHETIC_DRY_RUN",
        "scientific_interpretation_allowed": False,
    }
    write_json("wp_i_d_orientation_fit.json", orientation_fit)

    repeatability_matrix = []
    condition_pairs = [(distance, 0.0) for distance in distances] + [(0.6, angle) for angle in orientations[1:]]
    for distance, angle in condition_pairs:
        stats = condition_statistics(all_rows, distance, angle)
        repeatability_matrix.append(
            {
                "condition_id": f"D{distance:.1f}m_A{int(angle):02d}",
                "distance_m": distance,
                "orientation_deg": angle,
                "event_count": stats["event_count"],
                "peak_frequency_mean_Hz": stats["peak_frequency_Hz_mean"],
                "peak_frequency_SD_Hz": stats["peak_frequency_Hz_standard_deviation"],
                "peak_frequency_CV": stats["peak_frequency_Hz_coefficient_of_variation"],
                "peak_voltage_mean_V": stats["peak_voltage_V_mean"],
                "peak_voltage_CV": stats["peak_voltage_V_coefficient_of_variation"],
                "RMS_voltage_mean_V": stats["RMS_voltage_V_mean"],
                "RMS_voltage_CV": stats["RMS_voltage_V_coefficient_of_variation"],
                "spectral_centroid_mean_Hz": stats["spectral_centroid_Hz_mean"],
                "spectral_centroid_SD_Hz": stats["spectral_centroid_Hz_standard_deviation"],
            }
        )
    write_csv("wp_i_d_repeatability_matrix.csv", repeatability_matrix)

    peak_standard_error = np.asarray(
        [row["peak_voltage_V_standard_deviation"] / np.sqrt(row["event_count"]) for row in distance_stats]
    )
    exponent_uncertainty = monte_carlo_power_law_uncertainty(
        distances, peak_amplitude, peak_standard_error, uncertainty, samples=2000, seed=20260915
    )
    ref = next(row for row in distance_stats if np.isclose(row["distance_m"], 0.6))
    orient0, orient90 = orientation_stats[0], orientation_stats[-1]
    factor_db = 20.0 / np.log(10.0)
    contrast_repeatability = factor_db * np.sqrt(
        (orient0["peak_voltage_V_standard_deviation"] / np.sqrt(orient0["event_count"]) / orient0["peak_voltage_V_mean"]) ** 2
        + (orient90["peak_voltage_V_standard_deviation"] / np.sqrt(orient90["event_count"]) / orient90["peak_voltage_V_mean"]) ** 2
    )
    rng = np.random.default_rng(20260915)
    angle_sigma = uncertainty["orientation_standard_uncertainty_deg"]
    fit = orientation_fit["peak_voltage_fit"]
    perturbed0 = rng.normal(0.0, angle_sigma, 2000)
    perturbed90 = rng.normal(90.0, angle_sigma, 2000)
    model0 = np.sqrt(fit["A_parallel"] ** 2 * np.cos(np.deg2rad(perturbed0)) ** 2 + fit["A_floor"] ** 2)
    model90 = np.sqrt(fit["A_parallel"] ** 2 * np.cos(np.deg2rad(perturbed90)) ** 2 + fit["A_floor"] ** 2)
    contrast_orientation = float(np.std(20 * np.log10(model0 / model90), ddof=1))
    contrast_scale = factor_db * np.sqrt(2) * uncertainty["voltage_scale_relative_standard_uncertainty"]

    budgets = []

    def add_budget(quantity, components):
        combined = float(np.sqrt(sum(value**2 for value in components.values())))
        for component, value in components.items():
            budgets.append(
                {
                    "quantity": quantity,
                    "component": component,
                    "standard_uncertainty": value,
                    "unit": "dB" if quantity == "peak_orientation_contrast" else ("1" if quantity == "distance_exponent_n" else ("V" if quantity == "peak_voltage_0p6m" else "Hz")),
                    "provenance": "SYNTHETIC_UNCERTAINTY_ASSUMPTION_OR_SYNTHETIC_REPEATABILITY",
                }
            )
        budgets.append(
            {
                "quantity": quantity,
                "component": "COMBINED_RSS",
                "standard_uncertainty": combined,
                "unit": budgets[-1]["unit"],
                "provenance": "DERIVED_SYNTHETIC_DRY_RUN",
            }
        )
        return combined

    n_components = {
        name: values["standard_uncertainty"]
        for name, values in exponent_uncertainty["components"].items()
        if name != "combined"
    }
    combined_summary = {
        "distance_exponent_n": add_budget("distance_exponent_n", n_components),
        "peak_voltage_0p6m": add_budget(
            "peak_voltage_0p6m",
            {
                "repeatability": ref["peak_voltage_V_standard_deviation"] / np.sqrt(ref["event_count"]),
                "voltage_scale": ref["peak_voltage_V_mean"] * uncertainty["voltage_scale_relative_standard_uncertainty"],
            },
        ),
        "peak_frequency_0p6m": add_budget(
            "peak_frequency_0p6m",
            {
                "repeatability": ref["peak_frequency_Hz_standard_deviation"] / np.sqrt(ref["event_count"]),
                "frequency_reference": ref["peak_frequency_Hz_mean"] * uncertainty["frequency_reference_relative_standard_uncertainty"],
                "sampling": ref["peak_frequency_Hz_mean"] * uncertainty["sample_interval_relative_standard_uncertainty"],
            },
        ),
        "spectral_centroid_0p6m": add_budget(
            "spectral_centroid_0p6m",
            {
                "repeatability": ref["spectral_centroid_Hz_standard_deviation"] / np.sqrt(ref["event_count"]),
                "frequency_reference": ref["spectral_centroid_Hz_mean"] * uncertainty["frequency_reference_relative_standard_uncertainty"],
                "sampling": ref["spectral_centroid_Hz_mean"] * uncertainty["sample_interval_relative_standard_uncertainty"],
            },
        ),
        "peak_orientation_contrast": add_budget(
            "peak_orientation_contrast",
            {
                "repeatability": contrast_repeatability,
                "orientation": contrast_orientation,
                "voltage_scale": contrast_scale,
            },
        ),
    }
    write_csv("wp_i_d_uncertainty_budget.csv", budgets)
    dominant = {}
    for quantity in combined_summary:
        candidates = [row for row in budgets if row["quantity"] == quantity and row["component"] != "COMBINED_RSS"]
        dominant[quantity] = max(candidates, key=lambda row: row["standard_uncertainty"])["component"]
    uncertainty_summary = {
        "input_path": "validation/stage_i/system_350mhz/wp_i_d_synthetic/synthetic_uncertainty_assumptions.json",
        "input_hash": sha256_file(SUPPLEMENT / "synthetic_uncertainty_assumptions.json"),
        "data_origin": uncertainty["data_origin"],
        "scientific_use_allowed": False,
        "assumptions": uncertainty,
        "distance_exponent_monte_carlo": exponent_uncertainty,
        "combined_standard_uncertainty": combined_summary,
        "dominant_terms": dominant,
        "orientation_monte_carlo_seed": 20260915,
        "orientation_monte_carlo_samples": 2000,
    }
    write_json("wp_i_d_uncertainty_summary.json", uncertainty_summary)

    discrepancy_entries = [
        ("DISTANCE_TREND", "UNRESOLVED", "REPLACE_WITH_REAL_DISTANCE_SWEEP"),
        ("DISTANCE_FREQUENCY_STABILITY", "UNRESOLVED", "REPLACE_WITH_REAL_DISTANCE_SWEEP"),
        ("ORIENTATION_RESPONSE", "UNRESOLVED", "REPLACE_WITH_REAL_ORIENTATION_SWEEP"),
        ("REPEATABILITY", "INSTRUMENT", "REPLACE_WITH_REAL_REPETITIONS"),
        ("UNCERTAINTY_PROPAGATION", "CALIBRATION", "REPLACE_SYNTHETIC_UNCERTAINTY_ASSUMPTIONS"),
    ]
    ledger = {
        "ledger_type": "StageIModelDiscrepancyLedger",
        "entries": [
            {
                "quantity": quantity,
                "simulation_value": "NOT_APPLICABLE_PIPELINE_DRY_RUN",
                "measurement_value": "SYNTHETIC_DEVELOPMENT_RESULT",
                "difference": "NOT_INTERPRETED",
                "uncertainty": "SYNTHETIC_UNCERTAINTY_ASSUMPTION",
                "possible_source": source,
                "classification": "DEVELOPMENT_DIAGNOSTIC",
                "action_required": action,
                "evidence_type": "SYNTHETIC_DRY_RUN",
                "scientific_interpretation_allowed": False,
            }
            for quantity, source, action in discrepancy_entries
        ],
    }
    write_json("wp_i_d_discrepancy_ledger.json", ledger)

    status = {
        "WP_I_D_PIPELINE_DRY_RUN": "PASS",
        "DISTANCE_ANALYSIS_PIPELINE_READY": True,
        "ORIENTATION_ANALYSIS_PIPELINE_READY": True,
        "REPEATABILITY_PIPELINE_READY": True,
        "UNCERTAINTY_PIPELINE_READY": True,
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "input_provenance": "SYNTHETIC_DEVELOPMENT_INPUT",
        "scientific_validation_allowed": False,
        "condition_count": 8,
        "event_count": 90,
        "supplemental_event_count": 60,
        "new_raw_quality_pass_count": sum(row["passed"] for row in quality_rows),
        "H3_comparison_band_Hz": [200e6, 500e6],
        "H3_below_200MHz_extrapolation": False,
        "H3_distance_reference_boundary": "ONLY_FROZEN_1P0M_COPOLARIZED_REFERENCE",
        "H3_orientation_reference_boundary": "NO_DIRECT_H3_ORIENTATION_REFERENCE",
        "future_real_data_replacement": {
            "status": "ARCHITECTURE_READY",
            "core_analysis_code_change_required": False,
            "replace_only": [
                "RAW_DATA",
                "DISTANCE_ORIENTATION_METADATA",
                "INSTRUMENT_CALIBRATION_METADATA",
                "UNCERTAINTY_VALUES",
                "PROVENANCE",
                "VALIDATION_STATUS",
            ],
        },
        "supplement_manifest_hash": sha256_file(SUPPLEMENT / "sha256_manifest.json"),
        "WP_I_B_metrics_hash": sha256_file(WP_B_METRICS),
        "physical_simulations_rerun": False,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json("wp_i_d_status.json", status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
