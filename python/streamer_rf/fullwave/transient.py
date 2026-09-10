"""H3 transient-to-passive-transfer utilities."""
from pathlib import Path

import numpy as np

from .receiver import bounded_complex_interpolate


H3_BAND_HZ = (200e6, 500e6)
LOADING_THRESHOLDS = (0.25, 0.75)
BASELINE_PROVENANCE = "NUMERICAL_SOURCE_SHAPE_BASELINE_ONLY"


def read_g3_waveform(path):
    data = np.genfromtxt(Path(path), delimiter=",", names=True)
    required = {"time_s", "V_port_V", "I_port_A"}
    if not required.issubset(data.dtype.names or ()):
        raise ValueError("INCOMPLETE_G3_WAVEFORM")
    time = np.asarray(data["time_s"], dtype=float)
    voltage = np.asarray(data["V_port_V"], dtype=float)
    current = np.asarray(data["I_port_A"], dtype=float)
    if (
        time.ndim != 1 or time.size < 4 or np.any(~np.isfinite(time))
        or np.any(~np.isfinite(voltage)) or np.any(~np.isfinite(current))
        or np.any(np.diff(time) <= 0)
    ):
        raise ValueError("INVALID_G3_WAVEFORM")
    dt = np.diff(time)
    if not np.allclose(dt, np.median(dt), rtol=1e-8, atol=0):
        raise ValueError("NONUNIFORM_G3_WAVEFORM")
    return time, voltage, current


def total_port_voltage_transfer(s11, s21):
    s11 = np.asarray(s11, dtype=complex)
    s21 = np.asarray(s21, dtype=complex)
    if s11.shape != s21.shape or np.any(~np.isfinite(s11)) or np.any(~np.isfinite(s21)):
        raise ValueError("INVALID_S_PARAMETERS")
    denominator = 1.0 + s11
    if np.any(np.abs(denominator) <= np.finfo(float).eps):
        raise ValueError("ZERO_TOTAL_TX_PORT_VOLTAGE")
    return s21 / denominator


def input_admittance(impedance):
    impedance = np.asarray(impedance, dtype=complex)
    if np.any(~np.isfinite(impedance)) or np.any(np.abs(impedance) <= np.finfo(float).eps):
        raise ValueError("INVALID_INPUT_IMPEDANCE")
    return 1.0 / impedance


def zero_padded_rfft(time_s, waveform, multiplier=16):
    time_s = np.asarray(time_s, dtype=float)
    waveform = np.asarray(waveform, dtype=float)
    if time_s.shape != waveform.shape or time_s.ndim != 1 or multiplier < 1:
        raise ValueError("INVALID_TRANSFORM_INPUT")
    dt = float(np.median(np.diff(time_s)))
    if not np.allclose(np.diff(time_s), dt, rtol=1e-8, atol=0):
        raise ValueError("NONUNIFORM_TRANSFORM_INPUT")
    nfft = int(multiplier) * time_s.size
    frequency = np.fft.rfftfreq(nfft, dt)
    # Complex single-sided peak-amplitude convention; zero padding does not
    # alter amplitudes at frequencies that coincide with the original DFT bins.
    spectrum = (2.0 / waveform.size) * np.fft.rfft(waveform, n=nfft)
    spectrum[0] *= 0.5
    if nfft % 2 == 0:
        spectrum[-1] *= 0.5
    return frequency, spectrum, dt, nfft


def exact_band_mask(frequency_Hz, low_Hz=H3_BAND_HZ[0], high_Hz=H3_BAND_HZ[1]):
    frequency = np.asarray(frequency_Hz, dtype=float)
    if frequency.ndim != 1 or np.any(~np.isfinite(frequency)) or np.any(np.diff(frequency) <= 0):
        raise ValueError("INVALID_FREQUENCY_AXIS")
    if not 0 <= low_Hz < high_Hz:
        raise ValueError("INVALID_H3_BAND")
    return (frequency >= low_Hz) & (frequency <= high_Hz)


def interpolate_transfer(source_frequency, source_transfer, target_frequency):
    return bounded_complex_interpolate(source_frequency, source_transfer, target_frequency)


def apply_transfer(source_spectrum, transfer):
    source = np.asarray(source_spectrum, dtype=complex)
    response = np.asarray(transfer, dtype=complex)
    if source.shape != response.shape or np.any(~np.isfinite(source)) or np.any(~np.isfinite(response)):
        raise ValueError("INVALID_TRANSFER_INPUT")
    return source * response


def retained_spectral_fraction(full_spectrum, mask):
    spectrum = np.asarray(full_spectrum, dtype=complex)
    mask = np.asarray(mask, dtype=bool)
    if spectrum.shape != mask.shape or np.any(~np.isfinite(spectrum)):
        raise ValueError("INVALID_SPECTRAL_FRACTION_INPUT")
    total = float(np.sum(np.abs(spectrum[1:]) ** 2))
    return float(np.sum(np.abs(spectrum[mask]) ** 2) / total) if total > 0 else np.nan


def causal_band_impulse(frequency, transfer_on_grid, nfft, dt, minimum_delay_s, energy_fraction=0.999):
    frequency = np.asarray(frequency, dtype=float)
    transfer = np.asarray(transfer_on_grid, dtype=complex)
    if transfer.shape != frequency.shape or frequency.size != nfft // 2 + 1:
        raise ValueError("INVALID_IMPULSE_TRANSFER")
    if np.any(~np.isfinite(transfer)) or minimum_delay_s < 0 or not 0 < energy_fraction <= 1:
        raise ValueError("INVALID_IMPULSE_CONFIGURATION")
    impulse = np.fft.irfft(transfer, n=nfft)
    start = int(np.ceil(minimum_delay_s / dt))
    impulse[:start] = 0.0
    # The second FFT half represents wrapped negative-time content.
    impulse[nfft // 2:] = 0.0
    energy = impulse * impulse
    if not np.sum(energy) > 0:
        raise ValueError("ZERO_CAUSAL_IMPULSE")
    end = int(np.searchsorted(np.cumsum(energy) / np.sum(energy), energy_fraction)) + 1
    end = max(end, start + 1)
    return impulse[:end]


def causal_linear_convolution(source, impulse):
    source = np.asarray(source, dtype=float)
    impulse = np.asarray(impulse, dtype=float)
    if source.ndim != 1 or impulse.ndim != 1 or not source.size or not impulse.size:
        raise ValueError("INVALID_CONVOLUTION_INPUT")
    if np.any(~np.isfinite(source)) or np.any(~np.isfinite(impulse)):
        raise ValueError("NONFINITE_CONVOLUTION_INPUT")
    return np.convolve(source, impulse, mode="full")


def analytic_envelope(waveform):
    waveform = np.asarray(waveform, dtype=float)
    if waveform.ndim != 1 or waveform.size < 2 or np.any(~np.isfinite(waveform)):
        raise ValueError("INVALID_ENVELOPE_INPUT")
    weights = np.zeros(waveform.size)
    weights[0] = 1.0
    if waveform.size % 2 == 0:
        weights[1 : waveform.size // 2] = 2.0
        weights[waveform.size // 2] = 1.0
    else:
        weights[1 : (waveform.size + 1) // 2] = 2.0
    return np.abs(np.fft.ifft(np.fft.fft(waveform) * weights))


def implied_loading_current(admittance, source_voltage_spectrum):
    return apply_transfer(source_voltage_spectrum, admittance)


def loading_metrics(implied_current, frozen_current, relative_gate=1e-6):
    implied = np.asarray(implied_current, dtype=complex)
    frozen = np.asarray(frozen_current, dtype=complex)
    if implied.shape != frozen.shape or implied.size < 2 or np.any(~np.isfinite(implied)) or np.any(~np.isfinite(frozen)):
        raise ValueError("INVALID_LOADING_DIAGNOSTIC")
    denominator = np.linalg.norm(frozen)
    mismatch = float(np.linalg.norm(implied - frozen) / denominator) if denominator > 0 else np.nan
    corr = np.vdot(frozen, implied) / (np.linalg.norm(frozen) * np.linalg.norm(implied))
    valid = np.abs(frozen) >= relative_gate * np.max(np.abs(frozen))
    ratio = np.abs(implied[valid] / frozen[valid])
    return {
        "normalized_L2_current_mismatch": mismatch,
        "complex_correlation": [float(corr.real), float(corr.imag)],
        "complex_correlation_magnitude": float(abs(corr)),
        "complex_correlation_phase_rad": float(np.angle(corr)),
        "magnitude_ratio_median": float(np.median(ratio)),
        "magnitude_ratio_p10": float(np.percentile(ratio, 10)),
        "magnitude_ratio_p90": float(np.percentile(ratio, 90)),
        "stable_ratio_fraction": float(np.mean(valid)),
    }


def loading_status(normalized_l2):
    if not np.isfinite(normalized_l2) or normalized_l2 < 0:
        raise ValueError("INVALID_LOADING_MISMATCH")
    if normalized_l2 <= LOADING_THRESHOLDS[0]:
        return "FULL_WAVE_LOADING_MISMATCH_LOW"
    if normalized_l2 <= LOADING_THRESHOLDS[1]:
        return "FULL_WAVE_LOADING_MISMATCH_MODERATE"
    return "FULL_WAVE_LOADING_MISMATCH_HIGH"


def waveform_impedance(voltage, current, relative_gate=1e-6):
    voltage = np.asarray(voltage, dtype=complex)
    current = np.asarray(current, dtype=complex)
    if voltage.shape != current.shape or np.any(~np.isfinite(voltage)) or np.any(~np.isfinite(current)):
        raise ValueError("INVALID_WAVEFORM_IMPEDANCE_INPUT")
    peak = float(np.max(np.abs(current)))
    valid = np.zeros(current.shape, dtype=bool) if peak == 0 else np.abs(current) >= relative_gate * peak
    impedance = np.full(current.shape, np.nan + 1j * np.nan)
    impedance[valid] = voltage[valid] / current[valid]
    return impedance, valid


def matched_rise_smooth_step(time_s, reference_voltage):
    time = np.asarray(time_s, dtype=float)
    voltage = np.asarray(reference_voltage, dtype=float)
    if time.shape != voltage.shape or np.any(~np.isfinite(time)) or np.any(~np.isfinite(voltage)):
        raise ValueError("INVALID_BASELINE_INPUT")
    peak = float(np.max(np.abs(voltage)))
    initial = abs(voltage[0])
    crossings = []
    for fraction in (0.9, 0.1):
        indices = np.flatnonzero(np.abs(voltage) <= fraction * initial)
        if not indices.size:
            raise ValueError("G3_TRANSITION_TIMESCALE_NOT_RESOLVED")
        crossings.append(float(time[indices[0]]))
    rise_duration = crossings[1] - crossings[0]
    if rise_duration <= 0:
        raise ValueError("G3_TRANSITION_TIMESCALE_NOT_RESOLVED")
    x = np.clip((time - time[0]) / rise_duration, 0.0, 1.0)
    baseline = peak * (3.0 * x**2 - 2.0 * x**3)
    return baseline, rise_duration
