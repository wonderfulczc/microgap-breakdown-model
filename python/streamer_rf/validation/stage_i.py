"""Stage-I measurement contracts, traceability gates, and neutral metrics."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np


NOT_PROVIDED = "NOT_PROVIDED"
REFERENCE_PLANES = (
    "SOURCE_PORT",
    "TX_FEED",
    "FREE_SPACE_REFERENCE",
    "RX_FEED",
    "INSTRUMENT_INPUT",
)
DATA_LAYERS = ("RAW", "CALIBRATED", "DERIVED")
AMPLITUDE_CLASSES = (
    "ABSOLUTE_AMPLITUDE_VALID",
    "RELATIVE_AMPLITUDE_ONLY",
    "NORMALIZED_SHAPE_ONLY",
    "NOT_COMPARABLE",
)
VALIDATION_STATUSES = (
    "NOT_MEASURED",
    "DATA_AVAILABLE",
    "CALIBRATION_PENDING",
    "COMPARISON_READY",
    "PARTIALLY_VALIDATED",
    "VALIDATED",
    "NOT_RESOLVED",
    "MODEL_DISCREPANCY",
)
DISCREPANCY_SOURCES = (
    "SOURCE_MODEL",
    "FULL_WAVE_LOADING",
    "TX_GEOMETRY",
    "RX_GEOMETRY",
    "CABLE_CONNECTOR",
    "REFERENCE_PLANE",
    "MESH_NUMERICS",
    "RECEIVER_TRANSFER",
    "FULL_MAXWELL_REFERENCE",
    "POLARIZATION",
    "INSTRUMENT_BANDWIDTH",
    "CALIBRATION",
    "INSTRUMENT",
    "ENVIRONMENT",
    "UNRESOLVED",
)
STAGE4_TRUSTED_BAND_HZ = (2941408508.9091916, 7966314711.629051)
STAGE5_TRUSTED_BAND_HZ = (3047273105.1868486, 10233758844.919172)

MEASUREMENT_REQUIRED_FIELDS = (
    "experiment_id",
    "timestamp",
    "hardware_configuration",
    "transmitter_id",
    "receiver_id",
    "actual_tx_geometry",
    "actual_rx_geometry",
    "distance_m",
    "tx_orientation",
    "rx_orientation",
    "polarization",
    "environment",
    "reference_planes",
    "cable_metadata",
    "connector_metadata",
    "instrument_metadata",
    "calibration_metadata",
    "measurement_bandwidth_Hz",
    "sample_rate_Hz",
    "number_of_samples",
    "averaging",
    "repetition_index",
    "raw_data_path",
    "raw_data_hash",
    "uncertainty_metadata",
)


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _missing_or_finite(value, *, positive=False):
    if value == NOT_PROVIDED:
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(number) and (number > 0 if positive else True)


def validate_measurement_contract(record):
    if any(field not in record for field in MEASUREMENT_REQUIRED_FIELDS):
        raise ValueError("INCOMPLETE_STAGE_I_MEASUREMENT_CONTRACT")
    if record.get("contract_status") != "MEASUREMENT_TEMPLATE_NO_DATA":
        raise ValueError("MEASUREMENT_DATA_STATUS_CONFLICT")
    for field in ("distance_m", "measurement_bandwidth_Hz", "sample_rate_Hz", "number_of_samples"):
        if not _missing_or_finite(record[field], positive=True):
            raise ValueError(f"INVALID_MEASUREMENT_FIELD:{field}")
    repetition = record["repetition_index"]
    if repetition != NOT_PROVIDED and (not isinstance(repetition, int) or repetition < 0):
        raise ValueError("INVALID_REPETITION_INDEX")
    raw_path, raw_hash = record["raw_data_path"], record["raw_data_hash"]
    if (raw_path == NOT_PROVIDED) != (raw_hash == NOT_PROVIDED):
        raise ValueError("RAW_DATA_PATH_HASH_MISMATCH")
    if raw_hash != NOT_PROVIDED and (len(raw_hash) != 64 or any(c not in "0123456789abcdef" for c in raw_hash)):
        raise ValueError("INVALID_RAW_DATA_HASH")
    return True


def validate_vna_contract(record):
    if record.get("instrument_type") != "VNA" or record.get("measurement_status") != "NOT_MEASURED":
        raise ValueError("INVALID_VNA_CONTRACT")
    if record.get("parser") != "streamer_rf.fullwave.receiver.read_touchstone":
        raise ValueError("TOUCHSTONE_PARSER_DUPLICATION_FORBIDDEN")
    if record.get("raw_smoothing") is not False:
        raise ValueError("RAW_TOUCHSTONE_SMOOTHING_FORBIDDEN")
    if record.get("expected_files") != ["RX_S11.s1p", "TX_S11.s1p", "TX_RX_S21.s2p"]:
        raise ValueError("INVALID_VNA_FILE_CONTRACT")
    if record.get("capability_frequency_Hz") != [500.0, 67e9]:
        raise ValueError("INVALID_VNA_CAPABILITY")
    return True


def validate_scope_contract(record):
    if record.get("instrument_type") != "OSCILLOSCOPE" or record.get("measurement_status") != "NOT_MEASURED":
        raise ValueError("INVALID_OSCILLOSCOPE_CONTRACT")
    if record.get("capability_bandwidth_Hz") != 8e9 or record.get("capability_sample_rate_Hz") != 80e9:
        raise ValueError("INVALID_OSCILLOSCOPE_CAPABILITY")
    if record.get("capability_vertical_bits") != 12:
        raise ValueError("INVALID_OSCILLOSCOPE_CAPABILITY")
    required_columns = {"time_s", "voltage_V", "channel_id"}
    if not required_columns.issubset(record.get("raw_waveform_columns", ())):
        raise ValueError("INVALID_SCOPE_WAVEFORM_SCHEMA")
    required_metadata = {
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
    }
    if not required_metadata.issubset(record.get("required_waveform_metadata", ())):
        raise ValueError("INCOMPLETE_SCOPE_WAVEFORM_METADATA")
    if record.get("primary_values_normalized") is not False:
        raise ValueError("RAW_SCOPE_NORMALIZATION_FORBIDDEN")
    return True


def validate_waveform_metadata(record):
    numeric_positive = ("sample_interval_s", "sample_rate_Hz", "input_impedance_ohm", "record_length")
    if any(not _missing_or_finite(record.get(key, NOT_PROVIDED), positive=True) for key in numeric_positive):
        raise ValueError("INVALID_WAVEFORM_METADATA")
    dt, rate = record.get("sample_interval_s"), record.get("sample_rate_Hz")
    if dt != NOT_PROVIDED and rate != NOT_PROVIDED and not np.isclose(float(dt) * float(rate), 1.0, rtol=1e-9):
        raise ValueError("WAVEFORM_SAMPLE_INTERVAL_RATE_MISMATCH")
    pretrigger = record.get("pretrigger_fraction", NOT_PROVIDED)
    if pretrigger != NOT_PROVIDED and not 0 <= float(pretrigger) <= 1:
        raise ValueError("INVALID_PRETRIGGER_FRACTION")
    repetition = record.get("experiment_repetition_index", NOT_PROVIDED)
    if repetition != NOT_PROVIDED and (not isinstance(repetition, int) or repetition < 0):
        raise ValueError("INVALID_REPETITION_INDEX")
    return True


def frequency_to_hz(values, unit):
    scale = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}.get(unit)
    data = np.asarray(values, dtype=float)
    if scale is None or np.any(~np.isfinite(data)):
        raise ValueError("INVALID_FREQUENCY_UNIT")
    return data * scale


def validate_data_layer(record):
    layer = record.get("layer")
    if layer not in DATA_LAYERS:
        raise ValueError("INVALID_MEASUREMENT_DATA_LAYER")
    if layer == "RAW":
        if record.get("parent_hash") not in (None, NOT_PROVIDED) or record.get("immutable") is not True:
            raise ValueError("RAW_LAYER_MUST_BE_IMMUTABLE_AND_PARENTLESS")
        if record.get("operations") not in ([], None):
            raise ValueError("RAW_LAYER_OPERATION_FORBIDDEN")
    else:
        parent_hash = record.get("parent_hash")
        if not isinstance(parent_hash, str) or len(parent_hash) != 64:
            raise ValueError("TRACEABLE_PARENT_HASH_REQUIRED")
        if not record.get("operations"):
            raise ValueError("TRACEABLE_PROCESSING_OPERATION_REQUIRED")
    return True


def assert_raw_not_overwritten(original_hash, current_hash):
    if original_hash != current_hash:
        raise ValueError("RAW_MEASUREMENT_OVERWRITE_FORBIDDEN")
    return True


def reference_plane_compatible(simulation, measurement):
    if simulation not in REFERENCE_PLANES or measurement not in REFERENCE_PLANES:
        raise ValueError("UNKNOWN_REFERENCE_PLANE")
    return simulation == measurement


def amplitude_comparison_class(simulation, measurement):
    sim_plane = simulation.get("reference_plane")
    meas_plane = measurement.get("reference_plane")
    if not reference_plane_compatible(sim_plane, meas_plane):
        return "NOT_COMPARABLE"
    if not all(
        (
            simulation.get("receiver_id") == measurement.get("receiver_id"),
            simulation.get("load_impedance_ohm") == measurement.get("load_impedance_ohm"),
            simulation.get("geometry_id") == measurement.get("geometry_id"),
            simulation.get("instrument_impedance_ohm") == measurement.get("instrument_impedance_ohm"),
        )
    ):
        return "NORMALIZED_SHAPE_ONLY"
    cable_status = measurement.get("cable_status")
    if cable_status in ("CALIBRATED_OUT", "MEASURED_TRANSFER_AVAILABLE") and measurement.get(
        "calibration_traceable"
    ):
        return "ABSOLUTE_AMPLITUDE_VALID"
    if cable_status in ("MODELED", "UNRESOLVED"):
        return "RELATIVE_AMPLITUDE_ONLY"
    return "NOT_COMPARABLE"


def validate_time_alignment(mode, *, applied_shift_s=0.0, absolute_validation=False):
    modes = ("ABSOLUTE_TRIGGER_TIME", "PROPAGATION_CORRECTED_TIME", "FEATURE_ALIGNED_FOR_SHAPE_ONLY")
    if mode not in modes or not np.isfinite(applied_shift_s):
        raise ValueError("INVALID_TIME_ALIGNMENT")
    if absolute_validation and mode == "FEATURE_ALIGNED_FOR_SHAPE_ONLY":
        raise ValueError("FEATURE_ALIGNMENT_CANNOT_VALIDATE_ABSOLUTE_TIME")
    return "SHAPE_COMPARISON_ONLY" if mode == "FEATURE_ALIGNED_FOR_SHAPE_ONLY" else "ABSOLUTE_TIME_ELIGIBLE"


def feature_frequency_error(simulation_Hz, measurement_Hz):
    values = np.asarray([simulation_Hz, measurement_Hz], dtype=float)
    if np.any(~np.isfinite(values)) or measurement_Hz <= 0:
        raise ValueError("INVALID_FEATURE_FREQUENCY")
    return {
        "signed_error_Hz": float(simulation_Hz - measurement_Hz),
        "absolute_error_Hz": float(abs(simulation_Hz - measurement_Hz)),
        "relative_error": float(abs(simulation_Hz - measurement_Hz) / measurement_Hz),
    }


def spectral_centroid(frequency_Hz, spectrum):
    frequency = np.asarray(frequency_Hz, dtype=float)
    values = np.asarray(spectrum, dtype=complex)
    if frequency.ndim != 1 or values.shape != frequency.shape or np.any(~np.isfinite(frequency)) or np.any(~np.isfinite(values)):
        raise ValueError("INVALID_SPECTRUM")
    weight = np.abs(values) ** 2
    return float(np.sum(frequency * weight) / np.sum(weight)) if np.sum(weight) else np.nan


def compare_spectra(frequency_Hz, simulation, measurement):
    frequency = np.asarray(frequency_Hz, dtype=float)
    simulation = np.asarray(simulation, dtype=complex)
    measurement = np.asarray(measurement, dtype=complex)
    if frequency.ndim != 1 or simulation.shape != frequency.shape or measurement.shape != frequency.shape:
        raise ValueError("INVALID_SPECTRAL_COMPARISON")
    if np.any(~np.isfinite(frequency)) or np.any(~np.isfinite(simulation)) or np.any(~np.isfinite(measurement)):
        raise ValueError("INVALID_SPECTRAL_COMPARISON")
    sim_mag, meas_mag = np.abs(simulation), np.abs(measurement)
    correlation = np.corrcoef(sim_mag, meas_mag)[0, 1] if np.std(sim_mag) and np.std(meas_mag) else np.nan
    complex_error = np.linalg.norm(simulation - measurement) / np.linalg.norm(measurement)
    peak_sim, peak_meas = float(frequency[np.argmax(sim_mag)]), float(frequency[np.argmax(meas_mag)])
    amplitude_ratio = float(np.max(sim_mag) / np.max(meas_mag)) if np.max(meas_mag) > 0 else np.nan
    return {
        "feature_frequency": feature_frequency_error(peak_sim, peak_meas),
        "spectral_centroid_signed_error_Hz": spectral_centroid(frequency, simulation)
        - spectral_centroid(frequency, measurement),
        "normalized_spectral_correlation": float(correlation),
        "complex_normalized_RMSE": float(complex_error),
        "peak_amplitude_ratio": amplitude_ratio,
        "peak_amplitude_error_dB": float(20 * np.log10(amplitude_ratio)) if amplitude_ratio > 0 else np.nan,
    }


def compare_waveforms(time_s, simulation, measurement):
    time = np.asarray(time_s, dtype=float)
    simulation = np.asarray(simulation, dtype=float)
    measurement = np.asarray(measurement, dtype=float)
    if time.ndim != 1 or simulation.shape != time.shape or measurement.shape != time.shape or time.size < 2:
        raise ValueError("INVALID_WAVEFORM_COMPARISON")
    if np.any(~np.isfinite(time)) or np.any(~np.isfinite(simulation)) or np.any(~np.isfinite(measurement)):
        raise ValueError("INVALID_WAVEFORM_COMPARISON")
    denominator = np.sqrt(np.mean(measurement**2))
    return {
        "normalized_RMSE": float(np.sqrt(np.mean((simulation - measurement) ** 2)) / denominator)
        if denominator > 0
        else np.nan,
        "waveform_correlation": float(np.corrcoef(simulation, measurement)[0, 1])
        if np.std(simulation) and np.std(measurement)
        else np.nan,
        "arrival_time_error_s": float(time[np.argmax(np.abs(simulation))] - time[np.argmax(np.abs(measurement))]),
    }


def distance_law_fit(distance_m, amplitude):
    distance = np.asarray(distance_m, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    if distance.shape != amplitude.shape or distance.ndim != 1 or distance.size < 2:
        raise ValueError("INVALID_DISTANCE_SERIES")
    if np.any(~np.isfinite(distance)) or np.any(~np.isfinite(amplitude)) or np.any(distance <= 0) or np.any(amplitude <= 0):
        raise ValueError("INVALID_DISTANCE_SERIES")
    slope, intercept = np.polyfit(np.log(distance), np.log(amplitude), 1)
    return {"log_log_slope": float(slope), "log_intercept": float(intercept)}


def polarization_contrast(parallel_amplitude, orthogonal_amplitude):
    values = np.asarray([parallel_amplitude, orthogonal_amplitude], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("INVALID_POLARIZATION_AMPLITUDE")
    return {
        "amplitude_ratio_parallel_over_orthogonal": float(parallel_amplitude / orthogonal_amplitude),
        "contrast_dB": float(20 * np.log10(parallel_amplitude / orthogonal_amplitude)),
    }


def summarize_repetitions(values, confidence_z=1.96):
    data = np.asarray(values, dtype=float)
    if data.ndim != 1 or data.size < 2 or np.any(~np.isfinite(data)) or confidence_z <= 0:
        raise ValueError("INSUFFICIENT_REPETITION_DATA")
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=1))
    half = float(confidence_z * std / np.sqrt(data.size))
    return {
        "count": int(data.size),
        "mean": mean,
        "median": float(np.median(data)),
        "standard_deviation": std,
        "confidence_interval": [mean - half, mean + half],
        "coefficient_of_variation": float(std / abs(mean)) if mean != 0 else np.nan,
    }


def validate_repetition_records(records):
    indices = [record.get("repetition_index") for record in records]
    if any(not isinstance(index, int) or index < 0 for index in indices) or len(indices) != len(set(indices)):
        raise ValueError("REPETITION_INDEX_NOT_UNIQUE")
    if any(not record.get("raw_data_hash") for record in records):
        raise ValueError("REPETITION_RAW_EVENT_LINK_REQUIRED")
    return True


def validate_status_transition(before, after):
    allowed = {
        "NOT_MEASURED": {"DATA_AVAILABLE"},
        "DATA_AVAILABLE": {"CALIBRATION_PENDING", "COMPARISON_READY", "NOT_RESOLVED"},
        "CALIBRATION_PENDING": {"COMPARISON_READY", "NOT_RESOLVED"},
        "COMPARISON_READY": {"PARTIALLY_VALIDATED", "VALIDATED", "MODEL_DISCREPANCY", "NOT_RESOLVED"},
        "PARTIALLY_VALIDATED": {"VALIDATED", "MODEL_DISCREPANCY", "NOT_RESOLVED"},
        "VALIDATED": {"MODEL_DISCREPANCY"},
        "NOT_RESOLVED": {"DATA_AVAILABLE"},
        "MODEL_DISCREPANCY": {"COMPARISON_READY", "PARTIALLY_VALIDATED"},
    }
    if before not in VALIDATION_STATUSES or after not in VALIDATION_STATUSES or after not in allowed[before]:
        raise ValueError("INVALID_VALIDATION_STATUS_TRANSITION")
    return True


def validate_discrepancy_ledger(record):
    if record.get("ledger_type") != "StageIModelDiscrepancyLedger" or not isinstance(record.get("entries"), list):
        raise ValueError("INVALID_STAGE_I_DISCREPANCY_LEDGER")
    required = {
        "quantity",
        "simulation_value",
        "measurement_value",
        "difference",
        "uncertainty",
        "possible_source",
        "classification",
        "action_required",
    }
    for entry in record["entries"]:
        if not required.issubset(entry) or entry["possible_source"] not in DISCREPANCY_SOURCES:
            raise ValueError("INVALID_STAGE_I_DISCREPANCY_ENTRY")
    return True
