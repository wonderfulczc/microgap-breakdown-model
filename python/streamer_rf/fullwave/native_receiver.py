"""H4 trusted native-field to local receiver utilities."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from streamer_rf.rf.spectral.fft import Spectrum, compute_one_sided_spectrum

from .receiver import bounded_complex_interpolate


STAGE4_TRUSTED_BAND_HZ = (2941408508.9091916, 7966314711.629051)
STAGE5_TRUSTED_BAND_HZ = (3047273105.1868486, 10233758844.919172)
SYSTEM_BAND_HZ = (200e6, 500e6)
RECEIVER_AXIS = np.array([0.0, 0.0, 1.0])


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass(frozen=True)
class NativeRFSourceContract:
    source_case_id: str
    source_path: str
    source_hash: str
    observer_id: str
    observer_position_m: tuple[float, float, float]
    observer_distance_m: float
    time_start_s: float
    time_end_s: float
    dt_s: float
    samples: int
    field_components: tuple[str, ...]
    dominant_trusted_component: str
    receiver_axis: tuple[float, float, float]
    trusted_frequency_Hz: tuple[float, float]
    rf_trust_report_path: str
    rf_trust_report_hash: str
    native_propagation_already_included: bool = True

    def validate(self) -> None:
        values = np.asarray(
            [
                *self.observer_position_m,
                self.observer_distance_m,
                self.time_start_s,
                self.time_end_s,
                self.dt_s,
                *self.receiver_axis,
                *self.trusted_frequency_Hz,
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(values)) or self.samples < 4 or self.dt_s <= 0:
            raise ValueError("INVALID_NATIVE_RF_SOURCE_CONTRACT")
        if self.time_end_s <= self.time_start_s or self.observer_distance_m <= 0:
            raise ValueError("INVALID_NATIVE_RF_SOURCE_CONTRACT")
        if not np.isclose(np.linalg.norm(self.receiver_axis), 1.0):
            raise ValueError("INVALID_RECEIVER_AXIS")
        if self.trusted_frequency_Hz[0] >= self.trusted_frequency_Hz[1]:
            raise ValueError("INVALID_NATIVE_TRUST_BAND")
        if self.field_components != ("Ex_total", "Ey_total", "Ez_total"):
            raise ValueError("TOTAL_NATIVE_FIELD_REQUIRED")
        if self.dominant_trusted_component not in self.field_components:
            raise ValueError("INVALID_DOMINANT_FIELD_COMPONENT")
        if not self.native_propagation_already_included:
            raise ValueError("NATIVE_PROPAGATION_CONTRACT_REQUIRED")

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return {
            "source_case_id": self.source_case_id,
            "source_path": self.source_path,
            "source_hash": self.source_hash,
            "observer_id": self.observer_id,
            "observer_position_m": list(self.observer_position_m),
            "observer_distance_m": self.observer_distance_m,
            "time_grid": {
                "start_s": self.time_start_s,
                "end_s": self.time_end_s,
                "dt_s": self.dt_s,
                "samples": self.samples,
            },
            "field_components": list(self.field_components),
            "polarization_metadata": {
                "basis": "CARTESIAN_XYZ",
                "dominant_trusted_component": self.dominant_trusted_component,
                "receiver_axis_policy": "FIXED_ALONG_DOMINANT_TRUSTED_E_COMPONENT",
            },
            "receiver_axis": list(self.receiver_axis),
            "trusted_frequency_mask_Hz": list(self.trusted_frequency_Hz),
            "rf_trust_report_path": self.rf_trust_report_path,
            "rf_trust_report_hash": self.rf_trust_report_hash,
            "NATIVE_PROPAGATION_ALREADY_INCLUDED": True,
        }


def load_trusted_band(trust_path: str | Path, expected: tuple[float, float]) -> tuple[float, float]:
    record = json.loads(Path(trust_path).read_text())
    if not record.get("scientific_rf_valid") or not record.get("source_valid"):
        raise ValueError("NATIVE_RF_SOURCE_NOT_TRUSTED")
    actual = (float(record["trusted_frequency_low_Hz"]), float(record["trusted_frequency_high_Hz"]))
    if actual != tuple(expected):
        raise ValueError("FROZEN_NATIVE_TRUST_BAND_CHANGED")
    return actual


def exact_trusted_mask(frequency_Hz, trusted_band_Hz):
    frequency = np.asarray(frequency_Hz, dtype=float)
    low, high = np.asarray(trusted_band_Hz, dtype=float)
    if frequency.ndim != 1 or np.any(~np.isfinite(frequency)) or np.any(np.diff(frequency) <= 0):
        raise ValueError("INVALID_NATIVE_FREQUENCY_AXIS")
    if not np.isfinite(low + high) or not 0 <= low < high:
        raise ValueError("INVALID_NATIVE_TRUST_BAND")
    return (frequency >= low) & (frequency <= high)


def project_field(field_xyz, receiver_axis=RECEIVER_AXIS):
    field = np.asarray(field_xyz)
    axis = np.asarray(receiver_axis, dtype=float)
    if field.ndim != 2 or field.shape[1] != 3 or axis.shape != (3,):
        raise ValueError("INVALID_FIELD_PROJECTION_INPUT")
    if np.any(~np.isfinite(field)) or np.any(~np.isfinite(axis)) or not np.isclose(np.linalg.norm(axis), 1.0):
        raise ValueError("INVALID_FIELD_PROJECTION_INPUT")
    return field @ axis


def native_field_spectrum(time_s, field_xyz) -> Spectrum:
    return compute_one_sided_spectrum(
        time_s,
        field_xyz,
        window="hann",
        component_names=("Ex_total", "Ey_total", "Ez_total"),
    )


def apply_trusted_receiver_transfer(
    native_spectrum: Spectrum,
    receiver_frequency_Hz,
    H_rx_E,
    trusted_band_Hz,
    receiver_axis=RECEIVER_AXIS,
):
    trusted = exact_trusted_mask(native_spectrum.frequency_Hz, trusted_band_Hz)
    projected = project_field(native_spectrum.complex_spectrum, receiver_axis)
    voltage = np.full(projected.shape, np.nan + 1j * np.nan)
    transfer = np.full(projected.shape, np.nan + 1j * np.nan)
    if np.any(trusted):
        transfer[trusted] = bounded_complex_interpolate(
            receiver_frequency_Hz, H_rx_E, native_spectrum.frequency_Hz[trusted]
        )
        voltage[trusted] = transfer[trusted] * projected[trusted]
    return projected, transfer, voltage, trusted


def trusted_bandlimited_timeseries(
    voltage_spectrum, trusted_mask, n_samples: int, dt_s: float, time_start_s: float = 0.0
):
    voltage = np.asarray(voltage_spectrum, dtype=complex)
    trusted = np.asarray(trusted_mask, dtype=bool)
    expected = n_samples // 2 + 1
    if (
        voltage.shape != (expected,)
        or trusted.shape != voltage.shape
        or dt_s <= 0
        or not np.isfinite(time_start_s)
    ):
        raise ValueError("INVALID_BANDLIMITED_TIMESERIES_INPUT")
    finite = np.zeros_like(voltage)
    finite[trusted] = voltage[trusted]
    if np.any(~np.isfinite(finite)):
        raise ValueError("NONFINITE_TRUSTED_VOLTAGE")
    frequency = np.fft.rfftfreq(n_samples, dt_s)
    physical_phase_inverse = np.exp(2j * np.pi * frequency * time_start_s)
    return np.fft.irfft(finite * physical_phase_inverse / dt_s, n=n_samples)


def spectral_centroid(frequency_Hz, values):
    frequency = np.asarray(frequency_Hz, dtype=float)
    magnitude = np.abs(np.asarray(values, dtype=complex))
    weight = magnitude * magnitude
    return float(np.sum(frequency * weight) / np.sum(weight)) if np.sum(weight) > 0 else np.nan


def normalized_shape_metrics(frequency_Hz, incident, received):
    frequency = np.asarray(frequency_Hz, dtype=float)
    a = np.abs(np.asarray(incident, dtype=complex))
    b = np.abs(np.asarray(received, dtype=complex))
    if frequency.shape != a.shape or a.shape != b.shape or a.size < 2:
        raise ValueError("INVALID_NATIVE_SPECTRAL_METRIC_INPUT")
    if np.any(~np.isfinite(frequency)) or np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("INVALID_NATIVE_SPECTRAL_METRIC_INPUT")
    correlation = np.corrcoef(a, b)[0, 1] if np.std(a) > 0 and np.std(b) > 0 else np.nan
    return {
        "normalized_shape_correlation": float(correlation),
        "incident_peak_frequency_Hz": float(frequency[np.argmax(a)]),
        "received_peak_frequency_Hz": float(frequency[np.argmax(b)]),
        "peak_shift_Hz": float(frequency[np.argmax(b)] - frequency[np.argmax(a)]),
        "incident_centroid_frequency_Hz": spectral_centroid(frequency, a),
        "received_centroid_frequency_Hz": spectral_centroid(frequency, b),
        "centroid_shift_Hz": spectral_centroid(frequency, b) - spectral_centroid(frequency, a),
        "incident_band_integrated_metric": float(np.trapezoid(a * a, frequency)),
        "received_band_integrated_metric": float(np.trapezoid(b * b, frequency)),
    }


def validate_h4_result_contract(record):
    required = (
        "source_contracts",
        "receiver_id",
        "receiver_axis",
        "receiver_transfer_hash",
        "native_path_status",
        "native_350MHz_status",
        "experimental_validation_status",
        "H4_ABSOLUTE_AMPLITUDE_STATUS",
    )
    if any(key not in record for key in required):
        raise ValueError("INCOMPLETE_H4_RESULT_CONTRACT")
    if record["native_350MHz_status"] != "NOT_RESOLVED":
        raise ValueError("NATIVE_350MHZ_TRUST_PROMOTION_FORBIDDEN")
    if record["native_path_status"] != "TRUSTED_BAND_DEVELOPMENT_VERIFIED":
        raise ValueError("INVALID_H4_NATIVE_PATH_STATUS")
    if record.get("mechanism_attribution_consumed", True):
        raise ValueError("F_R4_MECHANISM_ATTRIBUTION_FORBIDDEN")
    if record.get("h2_h3_numerical_inputs_consumed", True):
        raise ValueError("H2_H3_NATIVE_PATH_MIXING_FORBIDDEN")
    if record.get("RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT") is not True:
        raise ValueError("H4_MESH_SENSITIVITY_MUST_PROPAGATE")
    if record["H4_ABSOLUTE_AMPLITUDE_STATUS"] != "NUMERICAL_REFERENCE_ONLY":
        raise ValueError("H4_ABSOLUTE_AMPLITUDE_TRUST_PROMOTION_FORBIDDEN")
    return True
