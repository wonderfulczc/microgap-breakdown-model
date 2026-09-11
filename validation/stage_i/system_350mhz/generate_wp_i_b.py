#!/usr/bin/env python3
"""Run the WP-I-B synthetic end-to-end validation pipeline dry run."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import resource
import time

import numpy as np

from streamer_rf.validation.stage_i import sha256_file
from streamer_rf.validation.wp_i_b import (
    ANALYSIS_BANDS_HZ,
    SYNTHETIC_CALIBRATION,
    SYNTHETIC_INPUT,
    event_metrics,
    normalized_magnitude_comparison,
    read_scope_waveform,
    repeatability_summary,
    scope_quality_gate,
    synthetic_loading_diagnostic,
    validate_vna_input,
    verify_synthetic_bundle,
)


ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "validation/stage_i/system_350mhz/synthetic_input"
OUTPUT = ROOT / "validation/stage_i/system_350mhz"
H3_RECEIVED = ROOT / "fullwave/h3/h3_received_spectrum.csv"
H3_SOURCE = ROOT / "fullwave/h3/h3_g3_source_spectrum.csv"


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def flatten_event(event, metrics):
    row = {
        "distance_m": event["distance_m"],
        "repetition_index": event["repetition_index"],
        "raw_data_path": event["file"],
        "raw_data_hash": sha256_file(INPUT / event["file"]),
        "peak_abs_voltage_V": metrics["peak_abs_voltage_V"],
        "rms_voltage_V": metrics["rms_voltage_V"],
        "noise_rms_V": metrics["noise_rms_V"],
        "peak_to_noise_SNR_dB": metrics["peak_to_noise_SNR_dB"],
    }
    for band_name, values in metrics.items():
        if band_name not in ANALYSIS_BANDS_HZ:
            continue
        prefix = band_name.lower()
        row[f"{prefix}_peak_frequency_Hz"] = values["peak_frequency_Hz"]
        row[f"{prefix}_spectral_centroid_Hz"] = values["spectral_centroid_Hz"]
        row[f"{prefix}_band_integrated_V2_Hz"] = values["band_integrated_V2_Hz"]
    return row


def main():
    started = time.perf_counter()
    bundle = verify_synthetic_bundle(INPUT)
    metadata = bundle["metadata"]

    vna_files = metadata["vna"]["files"]
    vna = {name: validate_vna_input(INPUT / name) for name in vna_files}
    vna_summary = {
        "parser": "streamer_rf.fullwave.receiver.read_touchstone",
        "raw_smoothing": False,
        "files": {
            name: {
                "sha256": sha256_file(INPUT / name),
                "frequency_min_Hz": float(data.frequency_Hz[0]),
                "frequency_max_Hz": float(data.frequency_Hz[-1]),
                "points": int(data.frequency_Hz.size),
                "frequency_step_Hz": float(np.median(np.diff(data.frequency_Hz))),
                "Z0_ohm": data.Z0_ohm,
                "parameters": sorted(data.parameters),
            }
            for name, data in vna.items()
        },
    }

    scope_metadata = metadata["oscilloscope"]
    background_path = INPUT / scope_metadata["background_file"]
    background_time, background_voltage = read_scope_waveform(background_path)
    background_quality = scope_quality_gate(
        background_time, background_voltage, scope_metadata, baseline_available=True
    )
    background_noise_rms = float(np.std(background_voltage, ddof=1))

    event_rows = []
    quality_rows = []
    spectra = {distance: [] for distance in scope_metadata["distances_m"]}
    frequency_axis = None
    for event in metadata["event_metadata"]:
        time_s, voltage = read_scope_waveform(INPUT / event["file"])
        quality = scope_quality_gate(time_s, voltage, scope_metadata, baseline_available=background_path.is_file())
        metrics, frequency, spectrum = event_metrics(time_s, voltage, scope_metadata)
        if frequency_axis is None:
            frequency_axis = frequency
        elif not np.array_equal(frequency_axis, frequency):
            raise ValueError("SCOPE_FREQUENCY_GRID_MISMATCH")
        spectra[event["distance_m"]].append(spectrum)
        event_rows.append(flatten_event(event, metrics))
        quality_rows.append({"file": event["file"], **quality})

    if len(event_rows) != 30 or any(len(items) != 15 for items in spectra.values()):
        raise ValueError("SCOPE_REPETITION_COUNT_MISMATCH")
    if not background_quality["passed"] or not all(row["passed"] for row in quality_rows):
        raise ValueError("SCOPE_QUALITY_GATE_FAILED")

    event_fields = list(event_rows[0])
    write_csv(OUTPUT / "wp_i_b_event_metrics.csv", event_fields, event_rows)

    repeat_rows = []
    repeat_quantities = (
        "full_50_500mhz_peak_frequency_Hz",
        "full_50_500mhz_spectral_centroid_Hz",
        "peak_abs_voltage_V",
        "rms_voltage_V",
    )
    for distance in scope_metadata["distances_m"]:
        for quantity in repeat_quantities:
            stats = repeatability_summary(event_rows, distance, quantity)
            repeat_rows.append(
                {
                    "distance_m": distance,
                    "quantity": quantity,
                    "count": stats["count"],
                    "mean": stats["mean"],
                    "standard_deviation": stats["standard_deviation"],
                    "coefficient_of_variation": stats["coefficient_of_variation"],
                    "confidence_95_low": stats["confidence_interval"][0],
                    "confidence_95_high": stats["confidence_interval"][1],
                }
            )
    write_csv(OUTPUT / "wp_i_b_repeatability.csv", list(repeat_rows[0]), repeat_rows)

    by_distance = {
        distance: {row["quantity"]: row for row in repeat_rows if np.isclose(row["distance_m"], distance)}
        for distance in scope_metadata["distances_m"]
    }
    d06, d10 = by_distance[0.6], by_distance[1.0]
    ensemble_magnitude = {
        distance: np.mean(np.abs(np.asarray(items)), axis=0) for distance, items in spectra.items()
    }
    full_mask = (frequency_axis >= 50e6) & (frequency_axis <= 500e6)
    distance_shape_correlation = float(
        np.corrcoef(ensemble_magnitude[0.6][full_mask], ensemble_magnitude[1.0][full_mask])[0, 1]
    )
    distance_trend = {
        "evidence_type": "SYNTHETIC_DRY_RUN",
        "distances_m": [0.6, 1.0],
        "peak_voltage_ratio_1p0m_over_0p6m": d10["peak_abs_voltage_V"]["mean"]
        / d06["peak_abs_voltage_V"]["mean"],
        "RMS_voltage_ratio_1p0m_over_0p6m": d10["rms_voltage_V"]["mean"]
        / d06["rms_voltage_V"]["mean"],
        "peak_frequency_shift_1p0m_minus_0p6m_Hz": d10[
            "full_50_500mhz_peak_frequency_Hz"
        ]["mean"]
        - d06["full_50_500mhz_peak_frequency_Hz"]["mean"],
        "spectral_centroid_shift_1p0m_minus_0p6m_Hz": d10[
            "full_50_500mhz_spectral_centroid_Hz"
        ]["mean"]
        - d06["full_50_500mhz_spectral_centroid_Hz"]["mean"],
        "ensemble_spectral_shape_correlation_50_500MHz": distance_shape_correlation,
        "development_trend": {
            "frequency_relatively_stable": bool(
                abs(
                    d10["full_50_500mhz_peak_frequency_Hz"]["mean"]
                    - d06["full_50_500mhz_peak_frequency_Hz"]["mean"]
                )
                <= 0.02 * 350e6
            ),
            "amplitude_reduced_at_1p0m": bool(
                d10["peak_abs_voltage_V"]["mean"] < d06["peak_abs_voltage_V"]["mean"]
            ),
            "scientific_evidence": False,
        },
    }
    write_json(OUTPUT / "wp_i_b_distance_trend.json", distance_trend)

    h3_received = np.genfromtxt(H3_RECEIVED, delimiter=",", names=True)
    h3_frequency = h3_received["frequency_Hz"]
    if h3_frequency.min() < 200e6 or h3_frequency.max() > 500e6:
        raise ValueError("H3_COMPARISON_SUPPORT_VIOLATION")
    h3_voltage = h3_received["V_rx_real_V"] + 1j * h3_received["V_rx_imag_V"]
    synthetic_1m = np.interp(h3_frequency, frequency_axis, ensemble_magnitude[1.0])
    comparison = normalized_magnitude_comparison(h3_frequency, h3_voltage, synthetic_1m)
    comparison_rows = [
        {
            "frequency_Hz": frequency,
            "H3_received_real_V": h3.real,
            "H3_received_imag_V": h3.imag,
            "H3_received_magnitude_V": abs(h3),
            "synthetic_1p0m_ensemble_magnitude_V": synthetic,
            "comparison_support": "H3_OVERLAP_200_500MHZ",
        }
        for frequency, h3, synthetic in zip(h3_frequency, h3_voltage, synthetic_1m)
    ]
    write_csv(OUTPUT / "wp_i_b_h3_comparison.csv", list(comparison_rows[0]), comparison_rows)
    amplitude_ratio = float(np.max(synthetic_1m) / np.max(np.abs(h3_voltage)))
    comparison.update(
        {
            "H3_COMPARISON_BAND_Hz": [200e6, 500e6],
            "evaluated_frequency_support_Hz": [float(h3_frequency[0]), float(h3_frequency[-1])],
            "50_TO_200MHZ_H3_STATUS": "OUTSIDE_H3_COMPARISON_SUPPORT",
            "amplitude_comparison_status": "DEVELOPMENT_ABSOLUTE_DIAGNOSTIC_ONLY",
            "synthetic_to_H3_peak_amplitude_ratio": amplitude_ratio,
            "scientific_amplitude_validation_allowed": False,
        }
    )

    h3_source = np.genfromtxt(H3_SOURCE, delimiter=",", names=True)
    if not np.array_equal(h3_source["frequency_Hz"], h3_frequency):
        raise ValueError("H3_SOURCE_RECEIVER_FREQUENCY_MISMATCH")
    v_g3 = h3_source["V_G3_real_V"] + 1j * h3_source["V_G3_imag_V"]
    i_g3 = h3_source["I_G3_real_A"] + 1j * h3_source["I_G3_imag_A"]
    tx_data = vna["vna/TX_S11.s1p"]
    impedance, admittance, implied, loading = synthetic_loading_diagnostic(
        tx_data.frequency_Hz, tx_data.parameters["S11"], tx_data.Z0_ohm, h3_frequency, v_g3, i_g3
    )
    loading_rows = [
        {
            "frequency_Hz": frequency,
            "Zin_real_ohm": z.real,
            "Zin_imag_ohm": z.imag,
            "Yin_real_S": y.real,
            "Yin_imag_S": y.imag,
            "I_implied_real_A": current.real,
            "I_implied_imag_A": current.imag,
            "I_G3_real_A": frozen.real,
            "I_G3_imag_A": frozen.imag,
        }
        for frequency, z, y, current, frozen in zip(h3_frequency, impedance, admittance, implied, i_g3)
    ]
    write_csv(OUTPUT / "wp_i_b_loading_diagnostic.csv", list(loading_rows[0]), loading_rows)
    loading.update(
        {
            "input_provenance": SYNTHETIC_INPUT,
            "evidence_type": "SYNTHETIC_DRY_RUN",
            "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED": True,
            "G2_G3_modified": False,
        }
    )

    discrepancies = [
        {
            "quantity": "H3_OVERLAP_SPECTRAL_PEAK_FREQUENCY_Hz",
            "simulation_value": comparison["simulation_peak_frequency_Hz"],
            "measurement_value": comparison["synthetic_peak_frequency_Hz"],
            "difference": comparison["simulation_peak_frequency_Hz"]
            - comparison["synthetic_peak_frequency_Hz"],
            "uncertainty": "SYNTHETIC_REPEATABILITY_ONLY",
            "possible_source": "SOURCE_MODEL",
            "classification": "DEVELOPMENT_DIAGNOSTIC",
            "action_required": "REPEAT_WITH_REAL_MEASUREMENT",
            "evidence_type": "SYNTHETIC_DRY_RUN",
        },
        {
            "quantity": "H3_OVERLAP_SPECTRAL_CENTROID_Hz",
            "simulation_value": comparison["simulation_centroid_Hz"],
            "measurement_value": comparison["synthetic_centroid_Hz"],
            "difference": comparison["simulation_centroid_Hz"] - comparison["synthetic_centroid_Hz"],
            "uncertainty": "SYNTHETIC_REPEATABILITY_ONLY",
            "possible_source": "SOURCE_MODEL",
            "classification": "DEVELOPMENT_DIAGNOSTIC",
            "action_required": "REPEAT_WITH_REAL_MEASUREMENT",
            "evidence_type": "SYNTHETIC_DRY_RUN",
        },
        {
            "quantity": "H3_OVERLAP_NORMALIZED_SHAPE_CORRELATION",
            "simulation_value": 1.0,
            "measurement_value": comparison["normalized_shape_correlation"],
            "difference": 1.0 - comparison["normalized_shape_correlation"],
            "uncertainty": "SYNTHETIC_REPEATABILITY_ONLY",
            "possible_source": "RX_GEOMETRY",
            "classification": "DEVELOPMENT_DIAGNOSTIC",
            "action_required": "REPEAT_WITH_REAL_MEASUREMENT",
            "evidence_type": "SYNTHETIC_DRY_RUN",
        },
        {
            "quantity": "FULL_WAVE_LOADING_CURRENT_NORMALIZED_L2",
            "simulation_value": 0.0,
            "measurement_value": loading["normalized_L2_current_mismatch"],
            "difference": loading["normalized_L2_current_mismatch"],
            "uncertainty": "SYNTHETIC_VNA_ONLY",
            "possible_source": "FULL_WAVE_LOADING",
            "classification": loading["status"],
            "action_required": "RETAIN_FEEDBACK_NOT_COUPLED_DEBT",
            "evidence_type": "SYNTHETIC_DRY_RUN",
        },
    ]
    ledger = {
        "ledger_type": "StageIModelDiscrepancyLedger",
        "dataset_provenance": SYNTHETIC_INPUT,
        "scientific_interpretation_allowed": False,
        "entries": discrepancies,
    }
    write_json(OUTPUT / "wp_i_b_discrepancy_ledger.json", ledger)

    quality_summary = {
        "status": "PASS",
        "input_provenance": SYNTHETIC_INPUT,
        "verified_manifest_files": bundle["verified_file_count"],
        "event_count": len(event_rows),
        "events_per_distance": 15,
        "sample_rate_Hz": scope_metadata["sample_rate_Hz"],
        "sample_interval_s": scope_metadata["sample_interval_s"],
        "record_length_samples": scope_metadata["record_length_samples"],
        "record_duration_s": scope_metadata["record_duration_s"],
        "background_noise_rms_V": background_noise_rms,
        "all_scope_events_passed": True,
        "clipping_test": "NO_REPEATED_EXTREMA_DETECTED",
        "frequency_coverage_Hz": metadata["frequency_coverage_Hz"],
        "reference_plane_status": SYNTHETIC_CALIBRATION,
        "event_results": quality_rows,
        "vna": vna_summary,
    }
    write_json(OUTPUT / "wp_i_b_quality_gate.json", quality_summary)

    summary = {
        "WP_I_B_INPUT_TYPE": SYNTHETIC_INPUT,
        "SCIENTIFIC_VALIDATION_ALLOWED": False,
        "allowed_use": "SOFTWARE_PIPELINE_DRY_RUN_ONLY",
        "raw_data_preserved": True,
        "raw_input_manifest_hash": sha256_file(INPUT / "sha256_manifest.json"),
        "output_data_layer": "DERIVED",
        "derived_parent_hash": sha256_file(INPUT / "sha256_manifest.json"),
        "processing_policy_path": "validation/stage_i/stage_i_processing_policy.json",
        "processing_policy_hash": sha256_file(ROOT / "validation/stage_i/stage_i_processing_policy.json"),
        "frequency_policy": {
            "raw_and_descriptive_support_Hz": [50e6, 500e6],
            "separate_descriptive_bands_Hz": {
                name: list(bounds) for name, bounds in ANALYSIS_BANDS_HZ.items()
            },
            "H3_COMPARISON_BAND_Hz": [200e6, 500e6],
            "50_TO_200MHZ_H3_STATUS": "OUTSIDE_H3_COMPARISON_SUPPORT",
            "weak_56MHz_component_preserved": True,
        },
        "vna_ingestion": vna_summary,
        "scope_ingestion": {
            "distances_m": scope_metadata["distances_m"],
            "events_per_distance": 15,
            "total_events": 30,
            "individual_raw_events_preserved": True,
        },
        "quality_gate": quality_summary,
        "distance_trend": distance_trend,
        "H3_comparison": comparison,
        "loading_diagnostic": loading,
        "discrepancy_ledger": "wp_i_b_discrepancy_ledger.json",
        "future_real_data_replacement": {
            "status": "ARCHITECTURE_READY",
            "core_parser_redesign_required": False,
            "core_comparison_API_redesign_required": False,
            "replace_only": [
                "DATA_FILES",
                "GEOMETRY_METADATA",
                "CALIBRATION_METADATA",
                "PROVENANCE",
                "UNCERTAINTY",
                "VALIDATION_STATUS",
            ],
        },
        "WP_I_B_PIPELINE_DRY_RUN": "PASS",
        "STAGE_I_END_TO_END_DRY_RUN": "PASS",
        "STAGE_I_EXPERIMENTAL_DATA": "NOT_PROVIDED",
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
        "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING": True,
        "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED": True,
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json(OUTPUT / "wp_i_b_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
