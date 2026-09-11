"""WP-I-B synthetic ingestion and comparison helpers.

The provenance gate in this module is deliberately stricter than the generic
Stage-I contracts: synthetic development inputs can exercise software but can
never satisfy a scientific validation status.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from streamer_rf.fullwave.receiver import bounded_complex_interpolate, read_touchstone
from streamer_rf.fullwave.transient import implied_loading_current, loading_metrics, loading_status

from .stage_i import sha256_file, spectral_centroid, summarize_repetitions


SYNTHETIC_INPUT = "SYNTHETIC_DEVELOPMENT_INPUT"
SYNTHETIC_CALIBRATION = "SYNTHETIC_IDEAL_REFERENCE_PLANE"
FORBIDDEN_SCIENTIFIC_STATUSES = {
    "MEASURED",
    "VALIDATED",
    "VALIDATED_WITH_VNA",
    "SYSTEM_350MHZ_VALIDATED",
}
ANALYSIS_BANDS_HZ = {
    "FULL_50_500MHZ": (50e6, 500e6),
    "LOW_50_100MHZ": (50e6, 100e6),
    "H3_OVERLAP_200_500MHZ": (200e6, 500e6),
}


def _load_json(path):
    return json.loads(Path(path).read_text())


def verify_synthetic_bundle(input_dir):
    root = Path(input_dir)
    metadata = _load_json(root / "measurement_metadata.json")
    geometry = _load_json(root / "geometry_contract.json")
    calibration = _load_json(root / "calibration_metadata.json")
    manifest = _load_json(root / "sha256_manifest.json")
    origins = (
        metadata.get("data_origin"),
        geometry.get("data_origin"),
        calibration.get("data_origin"),
        manifest.get("bundle_status"),
    )
    if any(origin != SYNTHETIC_INPUT for origin in origins):
        raise ValueError("SYNTHETIC_PROVENANCE_REQUIRED")
    if metadata.get("scientific_validation_allowed") is not False or manifest.get(
        "scientific_validation_allowed"
    ) is not False:
        raise ValueError("SYNTHETIC_SCIENTIFIC_VALIDATION_FORBIDDEN")
    if not FORBIDDEN_SCIENTIFIC_STATUSES.issubset(metadata.get("forbidden_statuses", ())):
        raise ValueError("SYNTHETIC_FORBIDDEN_STATUS_GUARD_INCOMPLETE")
    if calibration.get("vna_calibration_status") != SYNTHETIC_CALIBRATION:
        raise ValueError("INVALID_SYNTHETIC_CALIBRATION")
    if metadata.get("vna", {}).get("reference_plane") != "TX_FEED/RX_FEED":
        raise ValueError("INVALID_SYNTHETIC_VNA_REFERENCE_PLANE")
    if calibration.get("reference_planes") != {
        "vna": "TX_FEED_to_RX_FEED",
        "scope": "INSTRUMENT_INPUT",
    }:
        raise ValueError("INVALID_SYNTHETIC_REFERENCE_PLANES")

    entries = manifest.get("files", ())
    if manifest.get("file_count") != len(entries):
        raise ValueError("MANIFEST_FILE_COUNT_MISMATCH")
    verified = []
    for entry in entries:
        path = root / entry["path"]
        if not path.is_file():
            raise ValueError(f"MANIFEST_FILE_MISSING:{entry['path']}")
        if path.stat().st_size != entry["bytes"] or sha256_file(path) != entry["sha256"]:
            raise ValueError(f"MANIFEST_INTEGRITY_FAILURE:{entry['path']}")
        verified.append(entry["path"])
    expected_events = {event["file"] for event in metadata.get("event_metadata", ())}
    if not expected_events or not expected_events.issubset(verified):
        raise ValueError("EVENT_MANIFEST_COVERAGE_FAILURE")
    return {
        "metadata": metadata,
        "geometry": geometry,
        "calibration": calibration,
        "manifest": manifest,
        "verified_file_count": len(verified),
    }


def enforce_synthetic_status(status):
    if status in FORBIDDEN_SCIENTIFIC_STATUSES:
        raise ValueError("SYNTHETIC_INPUT_CANNOT_VALIDATE_SCIENCE")
    return status


def validate_vna_input(path, *, low_Hz=50e6, high_Hz=500e6, points=901, step_Hz=0.5e6):
    data = read_touchstone(path)
    frequency = data.frequency_Hz
    if frequency.shape != (points,) or not np.all(np.isfinite(frequency)):
        raise ValueError("VNA_FREQUENCY_GRID_MISMATCH")
    if not np.isclose(frequency[0], low_Hz) or not np.isclose(frequency[-1], high_Hz):
        raise ValueError("VNA_FREQUENCY_SUPPORT_MISMATCH")
    if not np.allclose(np.diff(frequency), step_Hz, rtol=0, atol=1e-6):
        raise ValueError("VNA_FREQUENCY_STEP_MISMATCH")
    if not np.isclose(data.Z0_ohm, 50.0):
        raise ValueError("VNA_REFERENCE_IMPEDANCE_MISMATCH")
    if any(np.any(~np.isfinite(values)) for values in data.parameters.values()):
        raise ValueError("NONFINITE_VNA_DATA")
    return data


def read_scope_waveform(path):
    path = Path(path)
    with path.open(newline="") as stream:
        header = next(csv.reader(stream))
    if header != ["time_s", "voltage_V"]:
        raise ValueError("INVALID_SCOPE_COLUMNS")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("INVALID_SCOPE_DATA")
    return data[:, 0], data[:, 1]


def scope_quality_gate(time_s, voltage_V, metadata, *, baseline_available):
    time = np.asarray(time_s, dtype=float)
    voltage = np.asarray(voltage_V, dtype=float)
    if time.ndim != 1 or voltage.shape != time.shape or time.size < 2:
        raise ValueError("INVALID_SCOPE_SHAPE")
    finite = bool(np.all(np.isfinite(time)) and np.all(np.isfinite(voltage)))
    monotonic = bool(np.all(np.diff(time) > 0))
    dt = float(np.median(np.diff(time)))
    dt_match = bool(np.allclose(np.diff(time), metadata["sample_interval_s"], rtol=1e-10, atol=1e-18))
    length_match = time.size == metadata["record_length_samples"]
    duration_match = bool(np.isclose(time.size * dt, metadata["record_duration_s"], rtol=1e-10))
    extrema_repeats = max(np.count_nonzero(voltage == np.max(voltage)), np.count_nonzero(voltage == np.min(voltage)))
    clipping_status = "NO_REPEATED_EXTREMA_DETECTED" if extrema_repeats == 1 else "POSSIBLE_CLIPPING"
    return {
        "finite": finite,
        "monotonic_time": monotonic,
        "sample_interval_s": dt,
        "sample_interval_valid": dt_match,
        "record_length_samples": int(time.size),
        "record_length_valid": bool(length_match and duration_match),
        "clipping_status": clipping_status,
        "baseline_available": bool(baseline_available),
        "passed": bool(
            finite
            and monotonic
            and dt_match
            and length_match
            and duration_match
            and extrema_repeats == 1
            and baseline_available
        ),
    }


def scope_spectrum(time_s, voltage_V, *, baseline_end_s):
    time = np.asarray(time_s, dtype=float)
    voltage = np.asarray(voltage_V, dtype=float)
    if time.shape != voltage.shape or time.size < 4 or np.any(np.diff(time) <= 0):
        raise ValueError("INVALID_SCOPE_SPECTRUM_INPUT")
    baseline_mask = time < baseline_end_s
    if np.count_nonzero(baseline_mask) < 2:
        raise ValueError("BASELINE_INTERVAL_UNAVAILABLE")
    centered = voltage - np.mean(voltage[baseline_mask])
    dt = float(np.median(np.diff(time)))
    spectrum = 2.0 / time.size * np.fft.rfft(centered)
    spectrum[0] *= 0.5
    if time.size % 2 == 0:
        spectrum[-1] *= 0.5
    return np.fft.rfftfreq(time.size, dt), spectrum, centered


def band_metrics(frequency_Hz, spectrum, low_Hz, high_Hz):
    frequency = np.asarray(frequency_Hz, dtype=float)
    values = np.asarray(spectrum, dtype=complex)
    grid_step = float(np.median(np.diff(frequency)))
    tolerance = max(
        64 * np.finfo(float).eps * max(abs(low_Hz), abs(high_Hz), 1.0),
        abs(grid_step) * 1e-9,
    )
    mask = (frequency >= low_Hz - tolerance) & (frequency <= high_Hz + tolerance)
    if np.count_nonzero(mask) < 2:
        raise ValueError("INSUFFICIENT_BAND_SUPPORT")
    selected_f = frequency[mask]
    selected = values[mask]
    magnitude = np.abs(selected)
    return {
        "frequency_low_Hz": float(selected_f[0]),
        "frequency_high_Hz": float(selected_f[-1]),
        "frequency_bins": int(selected_f.size),
        "peak_frequency_Hz": float(selected_f[np.argmax(magnitude)]),
        "spectral_centroid_Hz": spectral_centroid(selected_f, selected),
        "band_integrated_V2_Hz": float(np.trapezoid(magnitude**2, selected_f)),
    }


def event_metrics(time_s, voltage_V, scope_metadata):
    frequency, spectrum, centered = scope_spectrum(
        time_s, voltage_V, baseline_end_s=scope_metadata["trigger_time_s"]
    )
    baseline = np.asarray(time_s) < scope_metadata["trigger_time_s"]
    noise_rms = float(np.sqrt(np.mean(centered[baseline] ** 2)))
    peak = float(np.max(np.abs(centered)))
    rms = float(np.sqrt(np.mean(centered**2)))
    result = {
        "peak_abs_voltage_V": peak,
        "rms_voltage_V": rms,
        "noise_rms_V": noise_rms,
        "peak_to_noise_SNR_dB": float(20 * np.log10(peak / noise_rms)) if noise_rms > 0 else np.nan,
    }
    for name, (low, high) in ANALYSIS_BANDS_HZ.items():
        result[name] = band_metrics(frequency, spectrum, low, high)
    return result, frequency, spectrum


def repeatability_summary(event_rows, distance_m, quantity):
    values = [row[quantity] for row in event_rows if np.isclose(row["distance_m"], distance_m)]
    return summarize_repetitions(values)


def normalized_magnitude_comparison(frequency_Hz, simulation, synthetic):
    frequency = np.asarray(frequency_Hz, dtype=float)
    sim = np.abs(np.asarray(simulation, dtype=complex))
    syn = np.abs(np.asarray(synthetic, dtype=complex))
    if frequency.shape != sim.shape or sim.shape != syn.shape or sim.size < 2:
        raise ValueError("INVALID_SYNTHETIC_SPECTRAL_COMPARISON")
    sim_norm = sim / np.linalg.norm(sim)
    syn_norm = syn / np.linalg.norm(syn)
    return {
        "simulation_peak_frequency_Hz": float(frequency[np.argmax(sim)]),
        "synthetic_peak_frequency_Hz": float(frequency[np.argmax(syn)]),
        "simulation_centroid_Hz": spectral_centroid(frequency, sim),
        "synthetic_centroid_Hz": spectral_centroid(frequency, syn),
        "normalized_shape_correlation": float(np.corrcoef(sim_norm, syn_norm)[0, 1]),
        "normalized_shape_L2": float(np.linalg.norm(sim_norm - syn_norm)),
    }


def input_impedance_from_s11(s11, z0_ohm):
    reflection = np.asarray(s11, dtype=complex)
    if np.any(~np.isfinite(reflection)) or not np.isfinite(z0_ohm) or z0_ohm <= 0:
        raise ValueError("INVALID_S11_IMPEDANCE_INPUT")
    if np.any(np.isclose(1.0 - reflection, 0.0)):
        raise ValueError("S11_IMPEDANCE_SINGULAR")
    return z0_ohm * (1.0 + reflection) / (1.0 - reflection)


def synthetic_loading_diagnostic(vna_frequency, tx_s11, z0_ohm, h3_frequency, v_g3, i_g3):
    s11 = bounded_complex_interpolate(vna_frequency, tx_s11, h3_frequency)
    impedance = input_impedance_from_s11(s11, z0_ohm)
    admittance = 1.0 / impedance
    implied = implied_loading_current(admittance, v_g3)
    metrics = loading_metrics(implied, i_g3)
    return impedance, admittance, implied, {**metrics, "status": loading_status(metrics["normalized_L2_current_mismatch"])}
