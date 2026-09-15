"""Stage-I native-GHz validation contracts and sparse-bin comparisons."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .stage_i import spectral_centroid


STAGE4_TRUSTED_BAND_HZ = (2941408508.9091916, 7966314711.629051)
STAGE5_TRUSTED_BAND_HZ = (3047273105.1868486, 10233758844.919172)
STAGE4_TRUSTED_BINS_HZ = np.asarray(
    [
        3771036549.883574,
        4713795687.354467,
        5656554824.825361,
        6599313962.296254,
        7542073099.767148,
    ]
)
STAGE5_TRUSTED_BINS_HZ = np.asarray(
    [
        3936061094.1996818,
        4920076367.749602,
        5904091641.299522,
        6888106914.849443,
        7872122188.3993635,
        8856137461.949284,
        9840152735.499205,
    ]
)
NATIVE_VALIDATION_STATUSES = (
    "NOT_MEASURED",
    "DATA_AVAILABLE",
    "COMPARISON_READY",
    "PARTIALLY_VALIDATED",
    "VALIDATED",
    "MODEL_DISCREPANCY",
    "NOT_RESOLVED",
)
SYNTHETIC_NATIVE_PROVENANCE = "SYNTHETIC_NATIVE_GHZ_DRY_RUN"
NATIVE_DISCREPANCY_CATEGORIES = (
    "SOURCE_MODEL",
    "RECEIVER_TRANSFER",
    "MESH_NUMERICS",
    "FULL_MAXWELL_REFERENCE",
    "POLARIZATION",
    "INSTRUMENT_BANDWIDTH",
    "CALIBRATION",
    "UNRESOLVED",
)


def frequency_coverage_status(measurement_support_Hz, trusted_band_Hz):
    measurement = np.asarray(measurement_support_Hz, dtype=float)
    trusted = np.asarray(trusted_band_Hz, dtype=float)
    if (
        measurement.shape != (2,)
        or trusted.shape != (2,)
        or np.any(~np.isfinite(measurement))
        or np.any(~np.isfinite(trusted))
        or measurement[0] > measurement[1]
        or trusted[0] >= trusted[1]
    ):
        raise ValueError("INVALID_FREQUENCY_COVERAGE")
    overlap_low = max(measurement[0], trusted[0])
    overlap_high = min(measurement[1], trusted[1])
    if overlap_high < overlap_low:
        return "NO_TRUST_BAND_COVERAGE"
    if measurement[0] <= trusted[0] and measurement[1] >= trusted[1]:
        return "FULL_TRUST_BAND_COVERAGE"
    return "PARTIAL_TRUST_BAND_COVERAGE"


def trusted_bin_support(trusted_bins_Hz, measurement_support_Hz, receiver_support_Hz):
    bins = np.asarray(trusted_bins_Hz, dtype=float)
    measurement = np.asarray(measurement_support_Hz, dtype=float)
    receiver = np.asarray(receiver_support_Hz, dtype=float)
    if (
        bins.ndim != 1
        or not bins.size
        or np.any(~np.isfinite(bins))
        or np.any(np.diff(bins) <= 0)
        or measurement.shape != (2,)
        or receiver.shape != (2,)
        or np.any(~np.isfinite(measurement))
        or np.any(~np.isfinite(receiver))
    ):
        raise ValueError("INVALID_TRUST_INTERSECTION")
    return (
        (bins >= measurement[0])
        & (bins <= measurement[1])
        & (bins >= receiver[0])
        & (bins <= receiver[1])
    )


def comparison_mask(measurement_valid, rf_trust_valid, receiver_valid):
    masks = tuple(np.asarray(mask, dtype=bool) for mask in (measurement_valid, rf_trust_valid, receiver_valid))
    if any(mask.ndim != 1 for mask in masks) or len({mask.shape for mask in masks}) != 1:
        raise ValueError("INVALID_COMPARISON_MASK")
    return masks[0] & masks[1] & masks[2]


def sparse_trusted_bin_metrics(frequency_Hz, simulation, measurement, valid_mask=None):
    frequency = np.asarray(frequency_Hz, dtype=float)
    simulation = np.asarray(simulation, dtype=complex)
    measurement = np.asarray(measurement, dtype=complex)
    if frequency.ndim != 1 or simulation.shape != frequency.shape or measurement.shape != frequency.shape:
        raise ValueError("INVALID_SPARSE_NATIVE_METRICS")
    valid = np.ones(frequency.shape, dtype=bool) if valid_mask is None else np.asarray(valid_mask, dtype=bool)
    if valid.shape != frequency.shape:
        raise ValueError("INVALID_SPARSE_NATIVE_METRICS")
    f, sim, measured = frequency[valid], simulation[valid], measurement[valid]
    if f.size < 2 or np.any(~np.isfinite(f)) or np.any(~np.isfinite(sim)) or np.any(~np.isfinite(measured)):
        raise ValueError("INSUFFICIENT_SPARSE_NATIVE_BINS")
    sim_mag, measured_mag = np.abs(sim), np.abs(measured)
    sim_pattern = sim_mag / np.linalg.norm(sim_mag)
    measured_pattern = measured_mag / np.linalg.norm(measured_mag)
    correlation = np.corrcoef(sim_pattern, measured_pattern)[0, 1]
    return {
        "trusted_bin_count": int(f.size),
        "simulation_peak_frequency_Hz": float(f[np.argmax(sim_mag)]),
        "measurement_peak_frequency_Hz": float(f[np.argmax(measured_mag)]),
        "simulation_spectral_centroid_Hz": spectral_centroid(f, sim),
        "measurement_spectral_centroid_Hz": spectral_centroid(f, measured),
        "normalized_trusted_bin_correlation": float(correlation),
        "relative_amplitude_pattern_L2": float(np.linalg.norm(sim_pattern - measured_pattern)),
    }


def validate_native_measurement_contract(record):
    required = {
        "experiment_id",
        "validation_profile",
        "receiver_id",
        "receiver_geometry",
        "receiver_orientation",
        "distance_m",
        "instrument",
        "frequency_range_Hz",
        "frequency_resolution_Hz",
        "reference_plane",
        "receiver_load_impedance_ohm",
        "cable_metadata",
        "calibration_metadata",
        "raw_spectrum_path",
        "raw_waveform_path",
        "repetition_index",
        "uncertainty_metadata",
        "measurement_status",
    }
    if not required.issubset(record) or record.get("validation_profile") not in (
        "NATIVE_GHZ_STAGE4",
        "NATIVE_GHZ_STAGE5",
    ):
        raise ValueError("INVALID_NATIVE_GHZ_MEASUREMENT_CONTRACT")
    if record.get("measurement_status") not in NATIVE_VALIDATION_STATUSES:
        raise ValueError("INVALID_NATIVE_GHZ_VALIDATION_STATUS")
    return True


def validate_synthetic_native_fixture(metadata, rows):
    if metadata.get("data_origin") != SYNTHETIC_NATIVE_PROVENANCE:
        raise ValueError("SYNTHETIC_NATIVE_GHZ_PROVENANCE_REQUIRED")
    if metadata.get("scientific_validation_allowed") is not False:
        raise ValueError("SYNTHETIC_NATIVE_GHZ_CANNOT_VALIDATE")
    forbidden_statuses = {"MEASURED", "VALIDATED", "VALIDATED_WITH_VNA", "EXPERIMENTAL"}
    if any(value in forbidden_statuses for value in metadata.values() if isinstance(value, str)):
        raise ValueError("SYNTHETIC_NATIVE_GHZ_CANNOT_VALIDATE")
    if metadata.get("frequency_unit") != "Hz" or metadata.get("quantity_unit") not in (
        "V",
        "V_s",
        "V_per_m",
        "V_s_per_m",
    ):
        raise ValueError("INVALID_SYNTHETIC_NATIVE_UNITS")
    if metadata.get("trust_mask_semantics") != "STAGE_F_EXACT_TRUSTED_BINS_ONLY":
        raise ValueError("INVALID_SYNTHETIC_NATIVE_TRUST_SEMANTICS")
    required = {"frequency_Hz", "value_real", "value_imag", "trust_valid", "repetition_index"}
    if not rows or any(not required.issubset(row) for row in rows):
        raise ValueError("INVALID_SYNTHETIC_NATIVE_ROWS")
    for row in rows:
        values = np.asarray([row["frequency_Hz"], row["value_real"], row["value_imag"]], dtype=float)
        if np.any(~np.isfinite(values)) or row["trust_valid"] not in (0, 1, False, True):
            raise ValueError("INVALID_SYNTHETIC_NATIVE_ROWS")
        if not isinstance(row["repetition_index"], int) or row["repetition_index"] < 0:
            raise ValueError("INVALID_SYNTHETIC_NATIVE_REPETITION")
    return True


def ingest_synthetic_native_fixture(metadata_path, spectrum_path):
    metadata = json.loads(Path(metadata_path).read_text())
    with Path(spectrum_path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    converted = [
        {
            "frequency_Hz": float(row["frequency_Hz"]),
            "value_real": float(row["value_real"]),
            "value_imag": float(row["value_imag"]),
            "trust_valid": int(row["trust_valid"]),
            "repetition_index": int(row["repetition_index"]),
        }
        for row in rows
    ]
    validate_synthetic_native_fixture(metadata, converted)
    return metadata, converted


def validate_native_discrepancy_entries(entries):
    required = {"pathway", "quantity", "possible_source", "evidence_type", "action_required"}
    for entry in entries:
        if not required.issubset(entry):
            raise ValueError("INVALID_NATIVE_DISCREPANCY_ENTRY")
        if entry["pathway"] not in ("NATIVE_STAGE4", "NATIVE_STAGE5"):
            raise ValueError("INVALID_NATIVE_DISCREPANCY_PATHWAY")
        if entry["possible_source"] not in NATIVE_DISCREPANCY_CATEGORIES:
            raise ValueError("INVALID_NATIVE_DISCREPANCY_SOURCE")
    return True
