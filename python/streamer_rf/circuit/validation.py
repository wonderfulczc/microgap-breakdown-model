from __future__ import annotations

import numpy as np


def positive_peak_indices(values: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=float)
    idx = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0.0))[0] + 1
    return idx.astype(int)


def estimate_oscillation_frequency(time_s: np.ndarray, values: np.ndarray) -> float:
    t = np.asarray(time_s, dtype=float)
    peaks = positive_peak_indices(values)
    if peaks.size < 3:
        raise ValueError("at least three positive peaks are required")
    periods = np.diff(t[peaks])
    return float(2.0 * np.pi / np.mean(periods))


def estimate_damping_alpha(time_s: np.ndarray, values: np.ndarray) -> float:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    peaks = positive_peak_indices(y)
    if peaks.size < 3:
        raise ValueError("at least three positive peaks are required")
    amp = y[peaks]
    valid = amp > np.max(amp) * 1e-6
    peaks = peaks[valid]
    amp = amp[valid]
    if peaks.size < 3:
        raise ValueError("not enough well-resolved peaks for damping estimate")
    slope, _ = np.polyfit(t[peaks], np.log(amp), 1)
    return float(-slope)


def rlc_theory(L_H: float, C_F: float, R_ohm: float) -> dict[str, float]:
    omega0 = 1.0 / np.sqrt(L_H * C_F)
    alpha = R_ohm / (2.0 * L_H)
    omega_d = np.sqrt(max(omega0 * omega0 - alpha * alpha, 0.0))
    return {"omega0_rad_s": float(omega0), "alpha_1_s": float(alpha), "omega_d_rad_s": float(omega_d)}
