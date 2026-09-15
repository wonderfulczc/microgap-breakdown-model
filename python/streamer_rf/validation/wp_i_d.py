"""WP-I-D distance, orientation, repeatability, and uncertainty helpers."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .stage_i import sha256_file, summarize_repetitions


SYNTHETIC_INPUT = "SYNTHETIC_DEVELOPMENT_INPUT"
SYNTHETIC_UNCERTAINTY = "SYNTHETIC_UNCERTAINTY_ASSUMPTION"
FORBIDDEN_SYNTHETIC_STATUSES = {"MEASURED", "VALIDATED", "EXPERIMENTAL", "VALIDATED_WITH_VNA"}


def _load_json(path):
    return json.loads(Path(path).read_text())


def verify_supplement_bundle(root):
    root = Path(root)
    metadata = _load_json(root / "supplement_metadata.json")
    uncertainty = _load_json(root / "synthetic_uncertainty_assumptions.json")
    manifest = _load_json(root / "sha256_manifest.json")
    if metadata.get("data_origin") != SYNTHETIC_INPUT or manifest.get("bundle_status") != SYNTHETIC_INPUT:
        raise ValueError("WP_I_D_SYNTHETIC_PROVENANCE_REQUIRED")
    if metadata.get("scientific_validation_allowed") is not False or manifest.get(
        "scientific_validation_allowed"
    ) is not False:
        raise ValueError("WP_I_D_SCIENTIFIC_VALIDATION_FORBIDDEN")
    if uncertainty.get("data_origin") != SYNTHETIC_UNCERTAINTY or uncertainty.get(
        "scientific_use_allowed"
    ) is not False:
        raise ValueError("WP_I_D_UNCERTAINTY_PROVENANCE_REQUIRED")
    entries = manifest.get("files", ())
    if manifest.get("file_count") != len(entries):
        raise ValueError("WP_I_D_MANIFEST_COUNT_MISMATCH")
    verified = set()
    for entry in entries:
        path = root / entry["path"]
        if not path.is_file():
            raise ValueError(f"WP_I_D_MANIFEST_FILE_MISSING:{entry['path']}")
        if path.stat().st_size != entry["bytes"] or sha256_file(path) != entry["sha256"]:
            raise ValueError(f"WP_I_D_MANIFEST_INTEGRITY_FAILURE:{entry['path']}")
        verified.add(entry["path"])
    record_paths = {record["file"] for record in metadata.get("records", ())}
    if len(record_paths) != 60 or not record_paths.issubset(verified):
        raise ValueError("WP_I_D_RECORD_MANIFEST_COVERAGE_FAILURE")
    return {
        "metadata": metadata,
        "uncertainty": uncertainty,
        "manifest": manifest,
        "verified_file_count": len(verified),
    }


def enforce_synthetic_development_status(status):
    if status in FORBIDDEN_SYNTHETIC_STATUSES:
        raise ValueError("SYNTHETIC_DEVELOPMENT_INPUT_CANNOT_VALIDATE")
    return status


def merge_condition_records(existing, supplemental):
    merged = []
    identities = set()
    for record in (*existing, *supplemental):
        identity = (
            float(record["distance_m"]),
            float(record["orientation_deg"]),
            int(record["repetition_index"]),
        )
        if identity in identities:
            raise ValueError("DUPLICATE_DEVELOPMENT_EVENT")
        identities.add(identity)
        merged.append(dict(record))
    return merged


def condition_statistics(rows, distance_m, orientation_deg):
    selected = [
        row
        for row in rows
        if np.isclose(row["distance_m"], distance_m) and np.isclose(row["orientation_deg"], orientation_deg)
    ]
    if len(selected) < 2:
        raise ValueError("INSUFFICIENT_CONDITION_REPETITIONS")
    quantities = {
        "peak_frequency_Hz": "full_50_500mhz_peak_frequency_Hz",
        "spectral_centroid_Hz": "full_50_500mhz_spectral_centroid_Hz",
        "peak_voltage_V": "peak_abs_voltage_V",
        "RMS_voltage_V": "rms_voltage_V",
        "SNR_dB": "peak_to_noise_SNR_dB",
    }
    result = {
        "distance_m": float(distance_m),
        "orientation_deg": float(orientation_deg),
        "event_count": len(selected),
    }
    for output_name, field in quantities.items():
        stats = summarize_repetitions([row[field] for row in selected])
        result[f"{output_name}_mean"] = stats["mean"]
        result[f"{output_name}_standard_deviation"] = stats["standard_deviation"]
        result[f"{output_name}_coefficient_of_variation"] = stats["coefficient_of_variation"]
        result[f"{output_name}_confidence_95_low"] = stats["confidence_interval"][0]
        result[f"{output_name}_confidence_95_high"] = stats["confidence_interval"][1]
    return result


def power_law_fit(distance_m, amplitude):
    distance = np.asarray(distance_m, dtype=float)
    values = np.asarray(amplitude, dtype=float)
    if distance.ndim != 1 or values.shape != distance.shape or distance.size < 3:
        raise ValueError("INVALID_POWER_LAW_DATA")
    if np.any(~np.isfinite(distance)) or np.any(~np.isfinite(values)) or np.any(distance <= 0) or np.any(values <= 0):
        raise ValueError("INVALID_POWER_LAW_DATA")
    slope_intercept, covariance = np.polyfit(np.log(distance), np.log(values), 1, cov=True)
    slope, intercept = slope_intercept
    exponent = -float(slope)
    prefactor = float(np.exp(intercept))
    prediction = prefactor * distance ** (-exponent)
    residual = values - prediction
    denominator = float(np.sum((values - np.mean(values)) ** 2))
    r_squared = 1.0 - float(np.sum(residual**2)) / denominator if denominator else np.nan
    n_se = float(np.sqrt(covariance[0, 0]))
    log_a_se = float(np.sqrt(covariance[1, 1]))
    return {
        "A": prefactor,
        "n": exponent,
        "A_confidence_95": [float(np.exp(intercept - 1.96 * log_a_se)), float(np.exp(intercept + 1.96 * log_a_se))],
        "n_confidence_95": [exponent - 1.96 * n_se, exponent + 1.96 * n_se],
        "A_standard_error_log_space": log_a_se,
        "n_standard_error": n_se,
        "R_squared_amplitude_space": r_squared,
        "residuals": residual.tolist(),
        "predicted": prediction.tolist(),
        "interpretation": "SYNTHETIC_DEVELOPMENT_EMPIRICAL_FIT_ONLY",
    }


def distance_frequency_stability(distance_m, peak_frequency_Hz, centroid_Hz, reference_distance_m=0.6):
    distance = np.asarray(distance_m, dtype=float)
    peak = np.asarray(peak_frequency_Hz, dtype=float)
    centroid = np.asarray(centroid_Hz, dtype=float)
    if distance.shape != peak.shape or peak.shape != centroid.shape or distance.size < 3:
        raise ValueError("INVALID_DISTANCE_FREQUENCY_DATA")
    reference = np.flatnonzero(np.isclose(distance, reference_distance_m))
    if reference.size != 1:
        raise ValueError("DISTANCE_FREQUENCY_REFERENCE_MISSING")
    peak_shift = peak - peak[reference[0]]
    centroid_shift = centroid - centroid[reference[0]]
    return {
        "reference_distance_m": reference_distance_m,
        "peak_frequency_shift_Hz": peak_shift.tolist(),
        "spectral_centroid_shift_Hz": centroid_shift.tolist(),
        "maximum_absolute_peak_shift_Hz": float(np.max(np.abs(peak_shift))),
        "maximum_absolute_centroid_shift_Hz": float(np.max(np.abs(centroid_shift))),
        "peak_frequency_across_distance_CV": float(np.std(peak, ddof=1) / np.mean(peak)),
        "spectral_centroid_across_distance_CV": float(np.std(centroid, ddof=1) / np.mean(centroid)),
        "distance_peak_frequency_correlation": float(np.corrcoef(distance, peak)[0, 1]),
        "distance_centroid_correlation": float(np.corrcoef(distance, centroid)[0, 1]),
        "scientific_interpretation_allowed": False,
    }


def cosine_floor_fit(orientation_deg, amplitude):
    angle = np.asarray(orientation_deg, dtype=float)
    values = np.asarray(amplitude, dtype=float)
    if angle.ndim != 1 or values.shape != angle.shape or angle.size < 3:
        raise ValueError("INVALID_ORIENTATION_DATA")
    if np.any(~np.isfinite(angle)) or np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("INVALID_ORIENTATION_DATA")
    cosine_squared = np.cos(np.deg2rad(angle)) ** 2
    parallel = float(np.max(values))
    floor = float(np.min(values))
    parallel_span = 1.5 * parallel
    floor_span = parallel_span
    for _ in range(5):
        parallel_grid = np.linspace(max(0.0, parallel - parallel_span), parallel + parallel_span, 121)
        floor_grid = np.linspace(max(0.0, floor - floor_span), floor + floor_span, 121)
        predicted_grid = np.sqrt(
            parallel_grid[:, None, None] ** 2 * cosine_squared[None, None, :]
            + floor_grid[None, :, None] ** 2
        )
        objective = np.sum((predicted_grid - values[None, None, :]) ** 2, axis=2)
        best_parallel, best_floor = np.unravel_index(np.argmin(objective), objective.shape)
        parallel = float(parallel_grid[best_parallel])
        floor = float(floor_grid[best_floor])
        parallel_span /= 10.0
        floor_span /= 10.0
    predicted = np.sqrt(parallel**2 * cosine_squared + floor**2)
    residual = values - predicted
    denominator = float(np.sum((values - np.mean(values)) ** 2))
    r_squared = 1.0 - float(np.sum(residual**2)) / denominator if denominator else np.nan
    observed_normalized = values / values[0]
    cosine_reference = np.abs(np.cos(np.deg2rad(angle)))
    return {
        "A_parallel": parallel,
        "A_floor": floor,
        "cross_polarization_floor_ratio": floor / parallel if parallel else np.nan,
        "predicted": predicted.tolist(),
        "residuals": residual.tolist(),
        "R_squared": r_squared,
        "observed_normalized": observed_normalized.tolist(),
        "absolute_cosine_reference": cosine_reference.tolist(),
        "fit_method": "BOUNDED_NONNEGATIVE_GRID_REFINEMENT",
        "interpretation": "SYNTHETIC_POLARIZATION_MODEL_DRY_RUN_ONLY",
    }


def polarization_contrast_db(amplitude_0deg, amplitude_90deg):
    values = np.asarray([amplitude_0deg, amplitude_90deg], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("INVALID_POLARIZATION_CONTRAST")
    return float(20 * np.log10(values[0] / values[1]))


def monte_carlo_power_law_uncertainty(
    distance_m,
    amplitude,
    amplitude_standard_error,
    assumptions,
    *,
    samples=2000,
    seed=20260915,
):
    distance = np.asarray(distance_m, dtype=float)
    values = np.asarray(amplitude, dtype=float)
    standard_error = np.asarray(amplitude_standard_error, dtype=float)
    if distance.shape != values.shape or values.shape != standard_error.shape or samples < 100:
        raise ValueError("INVALID_MONTE_CARLO_INPUT")
    rng = np.random.default_rng(seed)
    z_repeat = rng.normal(size=(samples, distance.size))
    z_distance = rng.normal(size=(samples, distance.size))
    z_scale = rng.normal(size=(samples, distance.size))

    def exponents(use_repeat, use_distance, use_scale):
        perturbed_distance = np.broadcast_to(distance, (samples, distance.size)).copy()
        perturbed_values = np.broadcast_to(values, (samples, values.size)).copy()
        if use_distance:
            perturbed_distance += z_distance * assumptions["distance_standard_uncertainty_m"]
        if use_repeat:
            perturbed_values += z_repeat * standard_error
        if use_scale:
            perturbed_values = perturbed_values * (
                1 + z_scale * assumptions["voltage_scale_relative_standard_uncertainty"]
            )
        if np.any(perturbed_distance <= 0) or np.any(perturbed_values <= 0):
            raise ValueError("NONPOSITIVE_MONTE_CARLO_SAMPLE")
        x = np.log(perturbed_distance)
        y = np.log(perturbed_values)
        x_centered = x - np.mean(x, axis=1, keepdims=True)
        y_centered = y - np.mean(y, axis=1, keepdims=True)
        return -np.sum(x_centered * y_centered, axis=1) / np.sum(x_centered**2, axis=1)

    component_samples = {
        "repeatability": exponents(True, False, False),
        "distance_geometry": exponents(False, True, False),
        "voltage_scale": exponents(False, False, True),
        "combined": exponents(True, True, True),
    }
    result = {
        "method": "FIXED_SEED_MONTE_CARLO",
        "seed": seed,
        "samples": samples,
        "components": {},
    }
    for name, draws in component_samples.items():
        result["components"][name] = {
            "mean": float(np.mean(draws)),
            "standard_uncertainty": float(np.std(draws, ddof=1)),
            "confidence_95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        }
    combined = component_samples["combined"]
    first_half_std = float(np.std(combined[: samples // 2], ddof=1))
    full_std = float(np.std(combined, ddof=1))
    result["convergence_sanity"] = {
        "first_half_standard_uncertainty": first_half_std,
        "full_standard_uncertainty": full_std,
        "relative_change": abs(first_half_std - full_std) / full_std if full_std else 0.0,
    }
    return result
