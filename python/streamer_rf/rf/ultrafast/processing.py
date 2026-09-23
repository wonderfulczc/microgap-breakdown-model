from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy.constants import c, epsilon_0
from scipy.signal import find_peaks
from scipy.signal.windows import tukey


PAPER_REPORTED = "PAPER_REPORTED"
PROJECT_OPERATIONAL_CHOICE = "PROJECT_OPERATIONAL_CHOICE"
PULSE_WIDTH_AMBIGUOUS = "AMBIGUOUS_PENDING_WAVEFORM_REPRODUCTION"


@dataclass(frozen=True)
class ResamplingResult:
    time_s: np.ndarray
    values: np.ndarray
    strategy: str
    original_dt_min_s: float
    original_dt_median_s: float
    original_dt_max_s: float
    uniform_dt_s: float
    upsampling_factor: float
    resolution_status: str

    def metadata(self) -> dict[str, object]:
        data = asdict(self)
        data.pop("time_s")
        data.pop("values")
        return data


@dataclass(frozen=True)
class PulseDiagnostics:
    status: str
    peaks: tuple[dict[str, float], ...]
    primary_peak_index: int | None
    primary_peak_time_s: float | None
    primary_peak_abs_A_m_s: float | None
    PW_FWHM_s: float | None
    PW_FW1E_s: float | None
    interpretation: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CrossingResult:
    level_dB: float
    reference_frequency_Hz: float
    reference_ESD: float
    threshold_ESD: float
    raw_crossing_candidates_Hz: tuple[float, ...]
    validated_crossing_Hz: float | None
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SpectrumResult:
    frequency_Hz: np.ndarray
    transform_A_m: np.ndarray
    ESD: np.ndarray
    dt_s: float
    n_physical_samples: int
    n_fft: int
    physical_record_duration_s: float
    zero_padded_duration_s: float
    physical_frequency_resolution_Hz: float
    display_fft_bin_spacing_Hz: float
    nyquist_Hz: float
    window_name: str
    tukey_alpha: float
    parseval_time_integral: float
    parseval_spectral_integral: float
    parseval_relative_error: float
    normalization: str
    one_sided_convention: str
    dc_handling: str
    nyquist_handling: str

    def metadata(self) -> dict[str, object]:
        data = asdict(self)
        for key in ("frequency_Hz", "transform_A_m", "ESD"):
            data.pop(key)
        return data


@dataclass(frozen=True)
class ProcessingResult:
    original_time_s: np.ndarray
    original_M_A_m: np.ndarray
    uniform_time_s: np.ndarray
    M_uniform_A_m: np.ndarray
    processed_time_s: np.ndarray
    M_processed_A_m: np.ndarray
    derivative_time_s: np.ndarray
    dMdt_A_m_s: np.ndarray
    windowed_dMdt_A_m_s: np.ndarray
    resampling: ResamplingResult
    smoothing_metadata: dict[str, object]
    derivative_metadata: dict[str, object]
    temporal_resolution_metadata: dict[str, object]
    pulse: PulseDiagnostics
    spectrum: SpectrumResult
    f3dB: CrossingResult
    f10dB: CrossingResult
    spectral_slope: dict[str, object]
    processing_provenance: dict[str, object]

    def diagnostics(self) -> dict[str, object]:
        return {
            "resampling": self.resampling.metadata(),
            "smoothing": self.smoothing_metadata,
            "derivative": self.derivative_metadata,
            "temporal_resolution": self.temporal_resolution_metadata,
            "pulse": self.pulse.to_dict(),
            "spectrum": self.spectrum.metadata(),
            "f3dB": self.f3dB.to_dict(),
            "f10dB": self.f10dB.to_dict(),
            "spectral_slope": self.spectral_slope,
            "processing_provenance": self.processing_provenance,
        }


def _validate_input(time_s: np.ndarray, values: np.ndarray, *, minimum_samples: int = 5) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time_s, dtype=float)
    moment = np.asarray(values, dtype=float)
    if time.ndim != 1 or moment.ndim != 1 or time.shape != moment.shape:
        raise ValueError("time_s and M_A_m must be one-dimensional arrays of equal length")
    if time.size < minimum_samples:
        raise ValueError(f"at least {minimum_samples} physical samples are required")
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(moment)):
        raise ValueError("time_s and M_A_m must contain only finite values")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("time_s must be strictly monotonically increasing")
    return time, moment


def resample_uniform(
    time_s: np.ndarray,
    values: np.ndarray,
    *,
    strategy: str = "MEDIAN_INPUT_DT",
    user_dt_s: float | None = None,
    reference_dt_s: float | None = None,
) -> ResamplingResult:
    time, moment = _validate_input(time_s, values)
    intervals = np.diff(time)
    dt_min = float(np.min(intervals))
    dt_median = float(np.median(intervals))
    dt_max = float(np.max(intervals))
    if strategy == "MEDIAN_INPUT_DT":
        target_dt = dt_median
    elif strategy == "USER_SPECIFIED_DT":
        target_dt = user_dt_s
    elif strategy == "REFERENCE_PROFILE_DT":
        target_dt = reference_dt_s
    else:
        raise ValueError("unknown uniform resampling strategy")
    if target_dt is None or not math.isfinite(target_dt) or target_dt <= 0.0:
        raise ValueError("resampling strategy requires a finite positive target dt")
    span = float(time[-1] - time[0])
    count = int(math.floor(span / target_dt + 1e-12)) + 1
    if count < 5:
        raise ValueError("target dt leaves too few uniform samples")
    uniform_time = time[0] + np.arange(count, dtype=float) * target_dt
    uniform_values = np.interp(uniform_time, time, moment)
    factor = dt_median / target_dt
    resolution_status = (
        "INTERPOLATED_NOT_PHYSICALLY_RESOLVED"
        if target_dt < dt_min * (1.0 - 1e-9)
        else "ORIGINAL_PHYSICAL_TIME_RESOLUTION_PRESERVED"
    )
    return ResamplingResult(
        time_s=uniform_time,
        values=uniform_values,
        strategy=strategy,
        original_dt_min_s=dt_min,
        original_dt_median_s=dt_median,
        original_dt_max_s=dt_max,
        uniform_dt_s=float(target_dt),
        upsampling_factor=float(factor),
        resolution_status=resolution_status,
    )


def _moving_average_valid(time_s: np.ndarray, values: np.ndarray, samples: int) -> tuple[np.ndarray, np.ndarray]:
    if samples < 1 or samples > values.size:
        raise ValueError("invalid moving-average window length")
    if samples == 1:
        return time_s.copy(), values.copy()
    kernel = np.full(samples, 1.0 / samples)
    averaged = np.convolve(values, kernel, mode="valid")
    averaged_time = np.convolve(time_s, kernel, mode="valid")
    return averaged_time, averaged


def central_difference_second_order(time_s: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    time, signal = _validate_input(time_s, values)
    intervals = np.diff(time)
    if not np.allclose(intervals, intervals[0], rtol=1e-9, atol=max(1e-30, abs(intervals[0]) * 1e-12)):
        raise ValueError("central difference requires a uniform time grid")
    dt = float(intervals[0])
    derivative = (signal[2:] - signal[:-2]) / (2.0 * dt)
    return time[1:-1], derivative


def _threshold_width(time_s: np.ndarray, amplitude: np.ndarray, peak_index: int, fraction: float) -> float | None:
    threshold = amplitude[peak_index] * fraction
    left_candidates = np.flatnonzero(amplitude[:peak_index] <= threshold)
    right_candidates = np.flatnonzero(amplitude[peak_index + 1 :] <= threshold)
    if left_candidates.size == 0 or right_candidates.size == 0:
        return None
    left_low = int(left_candidates[-1])
    left_high = left_low + 1
    right_high = int(peak_index + 1 + right_candidates[0])
    right_low = right_high - 1

    def crossing(i0: int, i1: int) -> float:
        y0, y1 = amplitude[i0], amplitude[i1]
        if y1 == y0:
            return float(0.5 * (time_s[i0] + time_s[i1]))
        weight = (threshold - y0) / (y1 - y0)
        return float(time_s[i0] + weight * (time_s[i1] - time_s[i0]))

    return crossing(right_low, right_high) - crossing(left_low, left_high)


def analyze_pulses(
    time_s: np.ndarray,
    dMdt_A_m_s: np.ndarray,
    *,
    primary_peak_index: int | None = None,
    minimum_relative_height: float = 0.2,
    minimum_prominence_fraction: float = 0.1,
) -> PulseDiagnostics:
    time, derivative = _validate_input(time_s, dMdt_A_m_s)
    amplitude = np.abs(derivative)
    peak = float(np.max(amplitude))
    if peak <= 0.0:
        return PulseDiagnostics("NOT_RESOLVED_NO_PULSE", (), None, None, None, None, None, PULSE_WIDTH_AMBIGUOUS)
    indices, properties = find_peaks(
        amplitude,
        height=minimum_relative_height * peak,
        prominence=minimum_prominence_fraction * peak,
    )
    peaks = tuple(
        {
            "index": int(idx),
            "time_s": float(time[idx]),
            "abs_dMdt_A_m_s": float(amplitude[idx]),
            "prominence_A_m_s": float(properties["prominences"][j]),
        }
        for j, idx in enumerate(indices)
    )
    if primary_peak_index is None:
        if len(indices) == 0:
            return PulseDiagnostics("NOT_RESOLVED_NO_INTERIOR_PEAK", peaks, None, None, None, None, None, PULSE_WIDTH_AMBIGUOUS)
        if len(indices) > 1:
            return PulseDiagnostics("NOT_RESOLVED_MULTI_PEAK", peaks, None, None, None, None, None, PULSE_WIDTH_AMBIGUOUS)
        selected = int(indices[0])
    else:
        if primary_peak_index < 0 or primary_peak_index >= time.size:
            raise ValueError("primary_peak_index is outside the derivative array")
        selected = int(primary_peak_index)
    fwhm = _threshold_width(time, amplitude, selected, 0.5)
    fw1e = _threshold_width(time, amplitude, selected, 1.0 / math.e)
    status = "PASS" if fwhm is not None and fw1e is not None else "NOT_RESOLVED_EVENT_BOUNDS"
    return PulseDiagnostics(
        status=status,
        peaks=peaks,
        primary_peak_index=selected,
        primary_peak_time_s=float(time[selected]),
        primary_peak_abs_A_m_s=float(amplitude[selected]),
        PW_FWHM_s=fwhm,
        PW_FW1E_s=fw1e,
        interpretation=PULSE_WIDTH_AMBIGUOUS,
    )


def compute_esd_spectrum(
    time_s: np.ndarray,
    dMdt_A_m_s: np.ndarray,
    *,
    tukey_alpha: float = 0.5,
    zero_pad_duration_s: float = 12e-9,
) -> tuple[SpectrumResult, np.ndarray]:
    time, derivative = _validate_input(time_s, dMdt_A_m_s)
    intervals = np.diff(time)
    if not np.allclose(intervals, intervals[0], rtol=1e-9, atol=max(1e-30, abs(intervals[0]) * 1e-12)):
        raise ValueError("FFT requires a uniform derivative time grid")
    if not (0.0 <= tukey_alpha <= 1.0):
        raise ValueError("Tukey alpha must be between zero and one")
    dt = float(intervals[0])
    physical_duration = float((time.size - 1) * dt)
    n_fft = max(time.size, int(math.ceil(zero_pad_duration_s / dt)))
    padded_duration = float(n_fft * dt)
    window = tukey(time.size, alpha=tukey_alpha, sym=True)
    windowed = derivative * window
    transform = dt * np.fft.rfft(windowed, n=n_fft)
    frequency = np.fft.rfftfreq(n_fft, dt)
    esd = np.abs(transform) ** 2 / (6.0 * math.pi * epsilon_0 * c**3)
    weights = np.ones_like(frequency)
    if n_fft % 2 == 0:
        weights[1:-1] = 2.0
    else:
        weights[1:] = 2.0
    df = float(1.0 / padded_duration)
    time_integral = float(np.sum(windowed**2) * dt)
    spectral_integral = float(np.sum(weights * np.abs(transform) ** 2) * df)
    relative_error = abs(spectral_integral - time_integral) / max(abs(time_integral), 1e-300)
    result = SpectrumResult(
        frequency_Hz=frequency,
        transform_A_m=transform,
        ESD=esd,
        dt_s=dt,
        n_physical_samples=int(time.size),
        n_fft=int(n_fft),
        physical_record_duration_s=physical_duration,
        zero_padded_duration_s=padded_duration,
        physical_frequency_resolution_Hz=float(1.0 / physical_duration),
        display_fft_bin_spacing_Hz=df,
        nyquist_Hz=float(frequency[-1]),
        window_name="TUKEY",
        tukey_alpha=float(tukey_alpha),
        parseval_time_integral=time_integral,
        parseval_spectral_integral=spectral_integral,
        parseval_relative_error=float(relative_error),
        normalization="X(f)=dt*rFFT[x(t)]",
        one_sided_convention="NONNEGATIVE_FREQUENCIES_STORED_WITHOUT_ESD_DOUBLING",
        dc_handling="DC_RETAINED_AND_EXCLUDED_FROM_RELATIVE_CROSSING_SEARCH",
        nyquist_handling="NYQUIST_STORED_ONCE_FOR_EVEN_NFFT",
    )
    return result, windowed


def _interpolate_log_frequency(f0: float, f1: float, y0: float, y1: float, target: float) -> float:
    if y1 == y0:
        return float(math.sqrt(f0 * f1))
    weight = (target - y0) / (y1 - y0)
    return float(10.0 ** (math.log10(f0) + weight * (math.log10(f1) - math.log10(f0))))


def find_relative_db_crossing(
    frequency_Hz: np.ndarray,
    ESD: np.ndarray,
    *,
    level_dB: float,
    reference_frequency_Hz: float = 1e9,
    stable_bins: int = 3,
) -> CrossingResult:
    frequency = np.asarray(frequency_Hz, dtype=float)
    spectrum = np.asarray(ESD, dtype=float)
    if frequency.ndim != 1 or spectrum.shape != frequency.shape or frequency.size < 3:
        raise ValueError("frequency and ESD must be one-dimensional arrays of equal length")
    if level_dB >= 0.0:
        raise ValueError("relative crossing level must be negative dB")
    if not np.all(np.isfinite(frequency)) or not np.all(np.isfinite(spectrum)) or np.any(spectrum < 0.0):
        raise ValueError("frequency and ESD must be finite and ESD nonnegative")
    if reference_frequency_Hz <= frequency[0] or reference_frequency_Hz >= frequency[-1]:
        return CrossingResult(level_dB, reference_frequency_Hz, math.nan, math.nan, (), None, "REFERENCE_OUT_OF_BAND")
    positive = spectrum[spectrum > 0.0]
    floor = float(np.min(positive) * 1e-12) if positive.size else 1e-300
    reference = float(np.interp(reference_frequency_Hz, frequency, spectrum))
    threshold = reference * 10.0 ** (level_dB / 10.0)
    relative_db = 10.0 * np.log10(np.maximum(spectrum, floor) / max(reference, floor))
    start = int(np.searchsorted(frequency, reference_frequency_Hz, side="left"))
    candidates: list[tuple[int, float]] = []
    for i in range(start, frequency.size - 1):
        if relative_db[i] > level_dB and relative_db[i + 1] <= level_dB:
            candidates.append((i, _interpolate_log_frequency(frequency[i], frequency[i + 1], relative_db[i], relative_db[i + 1], level_dB)))
    raw = tuple(value for _, value in candidates)
    if not candidates:
        return CrossingResult(level_dB, reference_frequency_Hz, reference, threshold, raw, None, "OUT_OF_BAND")
    if len(candidates) > 1:
        return CrossingResult(level_dB, reference_frequency_Hz, reference, threshold, raw, None, "AMBIGUOUS_CROSSING")
    index, value = candidates[0]
    stop = min(relative_db.size, index + 1 + stable_bins)
    if stop - (index + 1) < stable_bins or not np.all(relative_db[index + 1 : stop] <= level_dB):
        return CrossingResult(level_dB, reference_frequency_Hz, reference, threshold, raw, None, "AMBIGUOUS_CROSSING")
    return CrossingResult(level_dB, reference_frequency_Hz, reference, threshold, raw, value, "PASS")


def fit_spectral_slope(frequency_Hz: np.ndarray, ESD: np.ndarray, fit_band_Hz: tuple[float, float]) -> dict[str, object]:
    frequency = np.asarray(frequency_Hz, dtype=float)
    spectrum = np.asarray(ESD, dtype=float)
    low, high = map(float, fit_band_Hz)
    if not (0.0 < low < high):
        raise ValueError("fit band must be finite, positive, and increasing")
    mask = (frequency >= low) & (frequency <= high) & np.isfinite(spectrum) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 3:
        return {"fit_band_Hz": [low, high], "slope_dB_per_decade": None, "R2": None, "number_of_bins": int(np.count_nonzero(mask)), "status": "OUT_OF_BAND"}
    x = np.log10(frequency[mask])
    y = 10.0 * np.log10(spectrum[mask])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return {"fit_band_Hz": [low, high], "slope_dB_per_decade": float(slope), "R2": float(r2), "number_of_bins": int(np.count_nonzero(mask)), "status": "PASS"}


def process_current_moment(
    time_s: np.ndarray,
    M_A_m: np.ndarray,
    *,
    profile: str = "KOILE_COMPATIBLE",
    resampling_strategy: str = "MEDIAN_INPUT_DT",
    user_dt_s: float | None = None,
    reference_dt_s: float | None = None,
    tukey_alpha: float = 0.5,
    zero_pad_duration_s: float = 12e-9,
    reference_frequency_Hz: float = 1e9,
    slope_fit_band_Hz: tuple[float, float] = (1e9, 5e9),
    event_interval_s: tuple[float, float] | None = None,
    primary_peak_index: int | None = None,
) -> ProcessingResult:
    original_time, original_moment = _validate_input(time_s, M_A_m)
    if event_interval_s is not None:
        start, end = event_interval_s
        if not (math.isfinite(start) and math.isfinite(end) and start < end):
            raise ValueError("event interval must be finite and increasing")
        mask = (original_time >= start) & (original_time <= end)
        original_time, original_moment = _validate_input(original_time[mask], original_moment[mask])
    if profile == "RAW":
        intervals = np.diff(original_time)
        if not np.allclose(intervals, intervals[0], rtol=1e-9, atol=max(1e-30, abs(intervals[0]) * 1e-12)):
            raise ValueError("RAW profile requires already-uniform physical input")
        resampled = resample_uniform(original_time, original_moment, strategy="USER_SPECIFIED_DT", user_dt_s=float(intervals[0]))
        processed_time, processed_moment = resampled.time_s, resampled.values
        smoothing_samples = 1
    elif profile in {"NO_SMOOTHING", "KOILE_COMPATIBLE"}:
        resampled = resample_uniform(
            original_time,
            original_moment,
            strategy=resampling_strategy,
            user_dt_s=user_dt_s,
            reference_dt_s=reference_dt_s,
        )
        smoothing_samples = 10 if profile == "KOILE_COMPATIBLE" else 1
        processed_time, processed_moment = _moving_average_valid(resampled.time_s, resampled.values, smoothing_samples)
    else:
        raise ValueError("profile must be RAW, NO_SMOOTHING, or KOILE_COMPATIBLE")
    smoothing_duration = smoothing_samples * resampled.uniform_dt_s
    if profile == "KOILE_COMPATIBLE" and not (1e-12 <= smoothing_duration <= 3e-12):
        smoothing_status = "KOILE_TIMESCALE_MISMATCH"
    else:
        smoothing_status = "PASS"
    derivative_time, derivative = central_difference_second_order(processed_time, processed_moment)
    pulse = analyze_pulses(derivative_time, derivative, primary_peak_index=primary_peak_index)
    original_dt = resampled.original_dt_median_s
    if pulse.PW_FWHM_s is None:
        temporal_resolution = {
            "status": "NOT_RESOLVED_PULSE_WIDTH",
            "original_dt_median_s": original_dt,
            "samples_per_FWHM": None,
            "rule": "USE_ORIGINAL_TIME_RESOLUTION_NOT_INTERPOLATED_DT",
        }
    else:
        samples_per_fwhm = pulse.PW_FWHM_s / original_dt
        if samples_per_fwhm >= 10.0:
            resolution_status = "RESOLVED_AT_10_SAMPLE_TARGET"
        elif samples_per_fwhm >= 5.0:
            resolution_status = "RESOLVED_AT_MINIMUM_5_SAMPLE_TARGET"
        else:
            resolution_status = "NOT_RESOLVED_TEMPORAL_RESOLUTION"
        temporal_resolution = {
            "status": resolution_status,
            "original_dt_median_s": original_dt,
            "samples_per_FWHM": float(samples_per_fwhm),
            "minimum_samples": 5,
            "target_samples": 10,
            "rule": "USE_ORIGINAL_TIME_RESOLUTION_NOT_INTERPOLATED_DT",
        }
    spectrum, windowed = compute_esd_spectrum(
        derivative_time,
        derivative,
        tukey_alpha=tukey_alpha,
        zero_pad_duration_s=zero_pad_duration_s,
    )
    f3 = find_relative_db_crossing(spectrum.frequency_Hz, spectrum.ESD, level_dB=-3.0, reference_frequency_Hz=reference_frequency_Hz)
    f10 = find_relative_db_crossing(spectrum.frequency_Hz, spectrum.ESD, level_dB=-10.0, reference_frequency_Hz=reference_frequency_Hz)
    slope = fit_spectral_slope(spectrum.frequency_Hz, spectrum.ESD, slope_fit_band_Hz)
    return ProcessingResult(
        original_time_s=original_time,
        original_M_A_m=original_moment,
        uniform_time_s=resampled.time_s,
        M_uniform_A_m=resampled.values,
        processed_time_s=processed_time,
        M_processed_A_m=processed_moment,
        derivative_time_s=derivative_time,
        dMdt_A_m_s=derivative,
        windowed_dMdt_A_m_s=windowed,
        resampling=resampled,
        smoothing_metadata={
            "profile": profile,
            "smoothing_samples": smoothing_samples,
            "smoothing_duration_s": float(smoothing_duration),
            "status": smoothing_status,
            "paper_reference_duration_s": [1e-12, 3e-12],
            "window_length_source": PAPER_REPORTED if profile == "KOILE_COMPATIBLE" else PROJECT_OPERATIONAL_CHOICE,
        },
        derivative_metadata={
            "method": "CENTRAL_DIFFERENCE_SECOND_ORDER",
            "source": PROJECT_OPERATIONAL_CHOICE,
            "boundary_policy": "DROP_FIRST_AND_LAST_PROCESSED_SAMPLES",
        },
        temporal_resolution_metadata=temporal_resolution,
        pulse=pulse,
        spectrum=spectrum,
        f3dB=f3,
        f10dB=f10,
        spectral_slope=slope,
        processing_provenance={
            "uniform_resampling_method": PAPER_REPORTED,
            "moving_average_10_samples": PAPER_REPORTED if profile == "KOILE_COMPATIBLE" else "NOT_APPLICABLE",
            "derivative_scheme": PROJECT_OPERATIONAL_CHOICE,
            "tukey_window_type": PAPER_REPORTED,
            "tukey_alpha": PROJECT_OPERATIONAL_CHOICE,
            "zero_pad_12_ns": PAPER_REPORTED if math.isclose(zero_pad_duration_s, 12e-9) else PROJECT_OPERATIONAL_CHOICE,
            "fft_discrete_normalization": PROJECT_OPERATIONAL_CHOICE,
            "ESD_formula": PAPER_REPORTED,
            "pulse_width_interpretation": PULSE_WIDTH_AMBIGUOUS,
        },
    )
