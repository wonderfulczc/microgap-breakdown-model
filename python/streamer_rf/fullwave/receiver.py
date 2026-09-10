"""Backend-independent H2 receiver contracts, Touchstone parsing and audits."""
from dataclasses import dataclass
from pathlib import Path
import hashlib

import numpy as np

from .foundation import C0


PORT_REFERENCE = "CURVEPORT_CENTER_GAP_CONNECTOR_EQUIVALENT"


def validate_receiver_geometry(record):
    if record.get("geometry_kind") != "H2_CANONICAL_REFERENCE_RX":
        raise ValueError("UNDOCUMENTED_RECEIVER_GEOMETRY")
    if record.get("hardware_status") != "GEOMETRY_AVAILABLE_ONLY":
        raise ValueError("HARDWARE_STATUS_CONFLICT")
    if record.get("microgap_geometry_included", True):
        raise ValueError("STAGE_H_MICROGAP_REFINEMENT_FORBIDDEN")
    values = np.asarray(
        [record.get("total_length_m"), record.get("tx_rx_distance_m"), record.get("Z0_ohm")],
        dtype=float,
    )
    if np.any(~np.isfinite(values) | (values <= 0)):
        raise ValueError("INVALID_RECEIVER_GEOMETRY")
    if record.get("port_reference") != PORT_REFERENCE:
        raise ValueError("RECEIVER_REFERENCE_PLANE_MISMATCH")
    if record.get("voltage_reference") != "POSITIVE_Z_ARM_MINUS_NEGATIVE_Z_ARM":
        raise ValueError("RECEIVER_VOLTAGE_POLARITY_MISMATCH")
    if record.get("positive_current_direction") != "PORT_INTO_RECEIVER":
        raise ValueError("RECEIVER_CURRENT_DIRECTION_MISMATCH")
    return True


def geometry_hash(record):
    import json

    encoded = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def band_mask(frequency_Hz, low_Hz, high_Hz):
    f = np.asarray(frequency_Hz, dtype=float)
    bounds = np.asarray([low_Hz, high_Hz], dtype=float)
    if f.ndim != 1 or np.any(~np.isfinite(f)) or np.any(np.diff(f) <= 0):
        raise ValueError("INVALID_FREQUENCY_AXIS")
    if np.any(~np.isfinite(bounds)) or not 0 <= low_Hz < high_Hz:
        raise ValueError("INVALID_FREQUENCY_BAND")
    return (f >= low_Hz) & (f <= high_Hz)


def near_far_classification(maximum_dimension_m, distance_m, frequency_Hz):
    values = np.asarray([maximum_dimension_m, distance_m, frequency_Hz], dtype=float)
    if np.any(~np.isfinite(values) | (values <= 0)):
        raise ValueError("INVALID_NEAR_FAR_INPUT")
    wavelength = C0 / frequency_Hz
    reactive_limit = 0.62 * np.sqrt(maximum_dimension_m**3 / wavelength)
    far_limit = 2 * maximum_dimension_m**2 / wavelength
    if distance_m < reactive_limit:
        status = "REACTIVE_NEAR_FIELD"
    elif distance_m < far_limit:
        status = "RADIATING_NEAR_FIELD"
    else:
        status = "FAR_FIELD"
    return {
        "frequency_Hz": frequency_Hz,
        "wavelength_m": wavelength,
        "maximum_dimension_m": maximum_dimension_m,
        "distance_m": distance_m,
        "reactive_near_field_limit_m": reactive_limit,
        "far_field_limit_2D2_over_lambda_m": far_limit,
        "status": status,
    }


def loaded_to_open_circuit(V_loaded, Z_receiver, Z_load):
    voltage = np.asarray(V_loaded, dtype=complex)
    impedance = np.asarray(Z_receiver, dtype=complex)
    load = complex(Z_load)
    if not np.isfinite(load.real) or not np.isfinite(load.imag) or abs(load) == 0:
        raise ValueError("INVALID_RECEIVER_LOAD")
    if np.any(~np.isfinite(voltage)) or np.any(~np.isfinite(impedance)):
        raise ValueError("NONFINITE_RECEIVER_QUANTITY")
    return voltage * (impedance + load) / load


def bounded_complex_interpolate(source_frequency, source_values, target_frequency):
    source_frequency = np.asarray(source_frequency, dtype=float)
    source_values = np.asarray(source_values, dtype=complex)
    target_frequency = np.asarray(target_frequency, dtype=float)
    if (
        source_frequency.ndim != 1
        or source_values.shape != source_frequency.shape
        or np.any(~np.isfinite(source_frequency))
        or np.any(~np.isfinite(source_values))
        or np.any(np.diff(source_frequency) <= 0)
        or target_frequency.ndim != 1
        or np.any(~np.isfinite(target_frequency))
    ):
        raise ValueError("INVALID_INTERPOLATION_INPUT")
    if target_frequency.size and (
        target_frequency.min() < source_frequency[0]
        or target_frequency.max() > source_frequency[-1]
    ):
        raise ValueError("FREQUENCY_EXTRAPOLATION_FORBIDDEN")
    return np.interp(target_frequency, source_frequency, source_values.real) + 1j * np.interp(
        target_frequency, source_frequency, source_values.imag
    )


def spectral_metrics(simulated, measured):
    a = np.asarray(simulated, dtype=complex)
    b = np.asarray(measured, dtype=complex)
    if a.shape != b.shape or a.ndim != 1 or a.size < 2 or np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("INVALID_SPECTRAL_METRIC_INPUT")
    an = np.abs(a)
    bn = np.abs(b)
    denom = np.linalg.norm(bn)
    l2 = np.linalg.norm(an - bn) / denom if denom else np.nan
    correlation = np.corrcoef(an, bn)[0, 1] if np.std(an) and np.std(bn) else np.nan
    db_error = 20 * np.log10(np.maximum(an, np.finfo(float).tiny) / np.maximum(bn, np.finfo(float).tiny))
    return {
        "normalized_magnitude_L2": float(l2),
        "magnitude_correlation": float(correlation),
        "median_absolute_dB_error": float(np.median(np.abs(db_error))),
    }


def g3_spectrum_plan(time_s, waveform, cumulative_fraction=0.999):
    time_s = np.asarray(time_s, dtype=float)
    waveform = np.asarray(waveform, dtype=float)
    if (
        time_s.ndim != 1
        or waveform.shape != time_s.shape
        or time_s.size < 4
        or np.any(~np.isfinite(time_s))
        or np.any(~np.isfinite(waveform))
        or np.any(np.diff(time_s) <= 0)
        or not 0 < cumulative_fraction < 1
    ):
        raise ValueError("INVALID_G3_SPECTRUM_INPUT")
    dt = np.diff(time_s)
    if not np.allclose(dt, np.median(dt), rtol=1e-8, atol=0):
        raise ValueError("NONUNIFORM_G3_WAVEFORM")
    signal = (waveform - np.mean(waveform)) * np.hanning(waveform.size)
    frequency = np.fft.rfftfreq(waveform.size, float(np.median(dt)))
    energy = np.abs(np.fft.rfft(signal)) ** 2
    energy[0] = 0
    if not np.sum(energy) > 0:
        raise ValueError("ZERO_G3_SPECTRUM")
    cumulative = np.cumsum(energy) / np.sum(energy)
    return {
        "dt_s": float(np.median(dt)),
        "df_Hz": float(frequency[1]),
        "nyquist_Hz": float(frequency[-1]),
        "peak_frequency_Hz": float(frequency[np.argmax(energy)]),
        "upper_frequency_at_cumulative_fraction_Hz": float(
            frequency[np.searchsorted(cumulative, cumulative_fraction)]
        ),
        "cumulative_fraction": cumulative_fraction,
        "frequency_Hz": frequency,
        "spectral_energy": energy,
    }


@dataclass(frozen=True)
class TouchstoneData:
    frequency_Hz: np.ndarray
    parameters: dict
    Z0_ohm: float
    source_format: str


def _to_complex(a, b, data_format):
    if data_format == "RI":
        return a + 1j * b
    angle = np.deg2rad(b)
    magnitude = a if data_format == "MA" else 10 ** (a / 20)
    return magnitude * np.exp(1j * angle)


def read_touchstone(path):
    path = Path(path)
    suffix = path.suffix.lower()
    ports = {".s1p": 1, ".s2p": 2}.get(suffix)
    if ports is None:
        raise ValueError("UNSUPPORTED_TOUCHSTONE_PORT_COUNT")
    option = None
    rows = []
    for raw in path.read_text().splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            option = line.upper().split()
        else:
            rows.extend(float(token) for token in line.split())
    if option is None or len(option) < 6 or option[2] != "S" or option[3] not in ("RI", "MA", "DB"):
        raise ValueError("INVALID_TOUCHSTONE_OPTION")
    unit_scale = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}.get(option[1])
    if unit_scale is None or "R" not in option:
        raise ValueError("INVALID_TOUCHSTONE_OPTION")
    z0 = float(option[option.index("R") + 1])
    columns = 1 + 2 * ports * ports
    if len(rows) == 0 or len(rows) % columns:
        raise ValueError("INVALID_TOUCHSTONE_DATA")
    data = np.asarray(rows, dtype=float).reshape(-1, columns)
    frequency = data[:, 0] * unit_scale
    if np.any(~np.isfinite(data)) or np.any(np.diff(frequency) <= 0) or not np.isfinite(z0) or z0 <= 0:
        raise ValueError("INVALID_TOUCHSTONE_DATA")
    names = ("S11",) if ports == 1 else ("S11", "S21", "S12", "S22")
    parameters = {
        name: _to_complex(data[:, 1 + 2 * index], data[:, 2 + 2 * index], option[3])
        for index, name in enumerate(names)
    }
    return TouchstoneData(frequency, parameters, z0, option[3])


def validate_transfer_contract(record):
    required = (
        "receiver_id", "hardware_status", "geometry_hash", "port_reference",
        "Z0_ohm", "VNA_validation_status", "H_rx_E_status",
    )
    if any(key not in record for key in required):
        raise ValueError("INCOMPLETE_RECEIVER_TRANSFER_CONTRACT")
    if record["VNA_validation_status"] not in ("VALIDATED_WITH_VNA", "SIMULATION_ONLY", "NOT_RESOLVED"):
        raise ValueError("INVALID_VNA_STATUS")
    if record["VNA_validation_status"] == "VALIDATED_WITH_VNA" and not record.get("measurement_id"):
        raise ValueError("MISSING_VNA_MEASUREMENT_ID")
    if not np.isfinite(record["Z0_ohm"]) or record["Z0_ohm"] <= 0:
        raise ValueError("INVALID_REFERENCE_IMPEDANCE")
    return True


def validate_exact_frequency_grid(frequency_Hz, low_Hz, high_Hz, points):
    frequency = np.asarray(frequency_Hz, dtype=float)
    if frequency.shape != (points,) or np.any(~np.isfinite(frequency)):
        raise ValueError("FREQUENCY_GRID_SHAPE_MISMATCH")
    expected = np.linspace(low_Hz, high_Hz, points)
    if not np.allclose(frequency, expected, rtol=0, atol=max(1e-6, abs(high_Hz) * 1e-12)):
        raise ValueError("FREQUENCY_GRID_VALUE_MISMATCH")
    return True


def receiver_validation_status(provenance, calibration_metadata=None):
    if provenance == "THEORETICAL_SURROGATE_ONLY":
        return "SIMULATION_ONLY"
    if provenance != "VNA_MEASUREMENT":
        raise ValueError("INVALID_RECEIVER_DATA_PROVENANCE")
    required = ("instrument_model", "calibration_type", "calibration_reference_plane", "Z0_ohm")
    if calibration_metadata is None or any(not calibration_metadata.get(key) for key in required):
        return "NOT_RESOLVED"
    if not np.isfinite(calibration_metadata["Z0_ohm"]) or calibration_metadata["Z0_ohm"] <= 0:
        raise ValueError("INVALID_REFERENCE_IMPEDANCE")
    return "VALIDATED_WITH_VNA"
