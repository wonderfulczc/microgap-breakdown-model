"""H5 source-to-receiver evidence contracts and pathway separation rules."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SIGNAL_LEVEL_DEFINITIONS = {
    "SIGNAL_0": "SOURCE_LEVEL_PHYSICS",
    "SIGNAL_0_NATIVE": "STAGE_F_TRUSTED_NATIVE_E_B_FIELD_OR_SPECTRUM",
    "SIGNAL_0_PORT": "G3_PHYSICS_DERIVED_V_PORT_I_PORT_TRANSIENT",
    "SIGNAL_1": "PROPAGATION_STRUCTURE_LEVEL",
    "SIGNAL_2": "RECEIVER_TERMINAL_LEVEL",
}

TRUST_HIERARCHY = (
    "TRUSTED_PHYSICS",
    "DEVELOPMENT_VERIFIED",
    "NUMERICAL_REFERENCE_ONLY",
    "NOT_RESOLVED",
    "EXPERIMENTAL_VALIDATION_PENDING",
)

STAGE4_TRUSTED_BAND_HZ = (2941408508.9091916, 7966314711.629051)
STAGE5_TRUSTED_BAND_HZ = (3047273105.1868486, 10233758844.919172)
H3_DEVELOPMENT_BAND_HZ = (200e6, 500e6)

STAGE_I_REQUIRED_INPUT_IDS = (
    "ACTUAL_TRANSMITTER_GEOMETRY",
    "ACTUAL_RECEIVER_GEOMETRY",
    "VNA_S11_S21_TOUCHSTONE",
    "CABLE_CONNECTOR_REFERENCE_PLANE_METADATA",
    "SOURCE_RECEIVER_DISTANCE_ORIENTATION",
    "OSCILLOSCOPE_RECEIVED_WAVEFORM",
    "MEASUREMENT_BANDWIDTH_SAMPLE_RATE",
    "CALIBRATION_ENVIRONMENT_METADATA",
    "REPEATED_MEASUREMENTS_UNCERTAINTY",
)


def validate_signal_levels(levels):
    if levels != SIGNAL_LEVEL_DEFINITIONS:
        raise ValueError("SIGNAL_LEVEL_SEMANTICS_CHANGED")
    if levels["SIGNAL_0_NATIVE"] == levels["SIGNAL_0_PORT"]:
        raise ValueError("NATIVE_AND_PORT_SOURCE_QUANTITIES_CONFLATED")
    return True


def native_propagation_action(native_propagation_already_included: bool, requested_action: str):
    if native_propagation_already_included:
        if requested_action != "APPLY_RECEIVER_RESPONSE_ONLY":
            raise ValueError("NATIVE_PROPAGATION_DOUBLE_APPLICATION_FORBIDDEN")
        return "NATIVE_PROPAGATION_OWNED_BY_STAGE_F"
    if requested_action == "APPLY_RECEIVER_RESPONSE_ONLY":
        raise ValueError("MISSING_NATIVE_PROPAGATION_MODEL")
    return "PROPAGATION_REQUIRED_UPSTREAM"


def exact_frequency_support(pathway: str):
    supports = {
        "POST_BREAKDOWN_G3": H3_DEVELOPMENT_BAND_HZ,
        "NATIVE_STAGE4": STAGE4_TRUSTED_BAND_HZ,
        "NATIVE_STAGE5": STAGE5_TRUSTED_BAND_HZ,
    }
    if pathway not in supports:
        raise ValueError("UNKNOWN_STAGE_H_PATHWAY")
    return supports[pathway]


def coherent_summation_permitted(paths):
    paths = tuple(paths)
    required = (
        "physical_tx_rx_geometry",
        "receiver_id",
        "reference_plane",
        "time_origin",
        "frequency_support",
        "phase_convention",
        "transfer_validation",
    )
    if len(paths) < 2 or any(any(not path.get(key) for key in required) for path in paths):
        return False
    return all(all(path[key] == paths[0][key] for path in paths[1:]) for key in required[:-1]) and all(
        path["transfer_validation"] == "VALIDATED" for path in paths
    )


def common_spectral_metrics(
    frequency_Hz,
    source_spectrum,
    receiver_spectrum,
    *,
    source_unit: str,
    receiver_unit: str,
):
    frequency = np.asarray(frequency_Hz, dtype=float)
    source = np.asarray(source_spectrum, dtype=complex)
    receiver = np.asarray(receiver_spectrum, dtype=complex)
    if (
        frequency.ndim != 1
        or source.shape != frequency.shape
        or receiver.shape != frequency.shape
        or frequency.size < 2
        or np.any(~np.isfinite(frequency))
        or np.any(np.diff(frequency) <= 0)
        or np.any(~np.isfinite(source))
        or np.any(~np.isfinite(receiver))
        or not source_unit
        or not receiver_unit
    ):
        raise ValueError("INVALID_COMMON_SPECTRAL_METRIC_INPUT")
    source_abs = np.abs(source)
    receiver_abs = np.abs(receiver)
    source_norm = source_abs / np.max(source_abs) if np.max(source_abs) > 0 else np.zeros_like(source_abs)
    receiver_norm = receiver_abs / np.max(receiver_abs) if np.max(receiver_abs) > 0 else np.zeros_like(receiver_abs)
    source_weight = source_abs**2
    receiver_weight = receiver_abs**2
    source_centroid = (
        float(np.sum(frequency * source_weight) / np.sum(source_weight)) if np.sum(source_weight) else np.nan
    )
    receiver_centroid = (
        float(np.sum(frequency * receiver_weight) / np.sum(receiver_weight)) if np.sum(receiver_weight) else np.nan
    )
    source_peak = float(frequency[np.argmax(source_abs)])
    receiver_peak = float(frequency[np.argmax(receiver_abs)])
    correlation = (
        float(np.corrcoef(source_abs, receiver_abs)[0, 1])
        if np.std(source_abs) > 0 and np.std(receiver_abs) > 0
        else np.nan
    )
    return {
        "valid_frequency_band_Hz": [float(frequency[0]), float(frequency[-1])],
        "frequency_samples": int(frequency.size),
        "source_unit": source_unit,
        "receiver_unit": receiver_unit,
        "source_peak_frequency_Hz": source_peak,
        "receiver_peak_frequency_Hz": receiver_peak,
        "receiver_induced_peak_shift_Hz": receiver_peak - source_peak,
        "source_centroid_frequency_Hz": source_centroid,
        "receiver_centroid_frequency_Hz": receiver_centroid,
        "receiver_induced_centroid_shift_Hz": receiver_centroid - source_centroid,
        "source_band_integrated_metric": float(np.trapezoid(source_weight, frequency)),
        "receiver_band_integrated_metric": float(np.trapezoid(receiver_weight, frequency)),
        "normalized_shape_correlation": correlation,
        "source_normalized_shape": source_norm,
        "receiver_normalized_shape": receiver_norm,
    }


def validate_stage_i_requirements(record):
    items = record.get("required_inputs")
    if not isinstance(items, list) or tuple(item.get("input_id") for item in items) != STAGE_I_REQUIRED_INPUT_IDS:
        raise ValueError("INCOMPLETE_STAGE_I_REQUIRED_INPUTS")
    if any(item.get("available") is not False for item in items):
        raise ValueError("FABRICATED_STAGE_I_INPUT")
    return True


def validate_stage_h_transfer_contract(record):
    if record.get("signal_level_definition") != SIGNAL_LEVEL_DEFINITIONS:
        raise ValueError("INVALID_STAGE_H_SIGNAL_LEVELS")
    if record.get("CROSS_PATH_COHERENT_SUMMATION") != "NOT_PERMITTED_CURRENT_CONFIGURATION":
        raise ValueError("CROSS_PATH_COHERENT_SUMMATION_FORBIDDEN")
    if record.get("CROSS_MECHANISM_ABSOLUTE_CONTRIBUTION") != "NOT_RESOLVED":
        raise ValueError("CROSS_MECHANISM_CONTRIBUTION_PROMOTION_FORBIDDEN")
    support = record.get("frequency_support", {})
    if tuple(support.get("POST_BREAKDOWN_G3_Hz", ())) != H3_DEVELOPMENT_BAND_HZ:
        raise ValueError("H3_FREQUENCY_SUPPORT_CHANGED")
    if tuple(support.get("NATIVE_STAGE4_Hz", ())) != STAGE4_TRUSTED_BAND_HZ:
        raise ValueError("STAGE4_TRUST_SUPPORT_CHANGED")
    if tuple(support.get("NATIVE_STAGE5_Hz", ())) != STAGE5_TRUSTED_BAND_HZ:
        raise ValueError("STAGE5_TRUST_SUPPORT_CHANGED")
    if record.get("native_350MHz_status") != "NOT_RESOLVED":
        raise ValueError("NATIVE_350MHZ_TRUST_PROMOTION_FORBIDDEN")
    if record.get("STAGE_H_TOOL_DEVELOPMENT") != "PASS":
        raise ValueError("STAGE_H_TOOL_DEVELOPMENT_INCOMPLETE")
    if record.get("STAGE_H_SCIENTIFIC_VALIDATION") != "PENDING_STAGE_I":
        raise ValueError("STAGE_H_SCIENTIFIC_VALIDATION_PROMOTION_FORBIDDEN")
    for path in record.get("native_paths", {}).values():
        if path.get("NATIVE_PROPAGATION_ALREADY_INCLUDED") is not True:
            raise ValueError("NATIVE_PROPAGATION_OWNERSHIP_MISSING")
        if path.get("time_domain_role") != "SECONDARY_RECONSTRUCTION_ONLY":
            raise ValueError("SPARSE_NATIVE_TIME_WAVEFORM_PROMOTED")
    debts = set(record.get("known_debts", ()))
    required_debts = {
        "FULL_WAVE_LOADING_MISMATCH_HIGH",
        "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED",
        "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
        "FULL_MAXWELL_REFERENCE_PENDING",
        "VNA_MEASUREMENT_PENDING",
        "PRODUCTION_TX_RX_GEOMETRY_PENDING",
    }
    if not required_debts.issubset(debts):
        raise ValueError("INHERITED_STAGE_H_DEBT_MISSING")
    return True
