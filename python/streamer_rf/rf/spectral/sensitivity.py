from __future__ import annotations

import numpy as np

from .fft import compute_one_sided_spectrum


def trust_frequency_from_error(
    frequency_Hz: np.ndarray,
    relative_error: np.ndarray,
    *,
    tolerance: float = 0.10,
    rolling_bins: int = 3,
) -> float:
    f = np.asarray(frequency_Hz, dtype=float)
    err = np.asarray(relative_error, dtype=float)
    if rolling_bins > 1:
        kernel = np.ones(rolling_bins) / rolling_bins
        smooth = np.convolve(err, kernel, mode="same")
    else:
        smooth = err
    bad = smooth > tolerance
    if not np.any(bad):
        return float(f[-1])
    first = int(np.argmax(bad))
    return float(f[max(first - 1, 0)])


def spectral_relative_error(reference_time: np.ndarray, reference_signal: np.ndarray, candidate_time: np.ndarray, candidate_signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref = compute_one_sided_spectrum(reference_time, reference_signal, window="hann")
    cand = compute_one_sided_spectrum(candidate_time, candidate_signal, window="hann")
    cand_esd = np.interp(ref.frequency_Hz, cand.frequency_Hz, cand.one_sided_esd[:, 0], left=np.nan, right=np.nan)
    scale = np.maximum(ref.one_sided_esd[:, 0], np.nanmax(ref.one_sided_esd[:, 0]) * 1e-12)
    err = np.abs(cand_esd - ref.one_sided_esd[:, 0]) / scale
    err[~np.isfinite(err)] = np.inf
    return ref.frequency_Hz, err


def compare_rf_mesh_levels(
    coarse_time: np.ndarray,
    coarse_waveform: np.ndarray,
    fine_time: np.ndarray,
    fine_waveform: np.ndarray,
    *,
    tolerance: float = 0.10,
    rolling_bins: int = 3,
) -> dict[str, object]:
    f, err = spectral_relative_error(fine_time, fine_waveform, coarse_time, coarse_waveform)
    return {
        "frequency_Hz": f,
        "relative_error": err,
        "mesh_trust_frequency_Hz": trust_frequency_from_error(f, err, tolerance=tolerance, rolling_bins=rolling_bins),
        "tolerance": tolerance,
    }
