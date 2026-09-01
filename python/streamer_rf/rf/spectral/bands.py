from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Band:
    name: str
    f_low_Hz: float
    f_high_Hz: float


STANDARD_BANDS = (
    Band("VHF", 30e6, 300e6),
    Band("UHF", 300e6, 3e9),
    Band("SHF", 3e9, 30e9),
)

PROJECT_BANDS = (
    Band("60-90 MHz", 60e6, 90e6),
    Band("100-200 MHz", 100e6, 200e6),
    Band("200-300 MHz", 200e6, 300e6),
    Band("300-500 MHz", 300e6, 500e6),
    Band("500 MHz-1 GHz", 500e6, 1e9),
    Band("1-3 GHz", 1e9, 3e9),
    Band("3-10 GHz", 3e9, 10e9),
)


def integrate_band(frequency_Hz: np.ndarray, density: np.ndarray, band: Band) -> float:
    f = np.asarray(frequency_Hz, dtype=float)
    y = np.asarray(density, dtype=float)
    if f.ndim != 1 or y.shape[0] != f.size:
        raise ValueError("frequency and density must have matching first dimension")
    if band.f_high_Hz <= f[0] or band.f_low_Hz >= f[-1]:
        return 0.0
    low = max(band.f_low_Hz, f[0])
    high = min(band.f_high_Hz, f[-1])
    if low >= high:
        return 0.0
    inner = (f > low) & (f < high)
    x = np.r_[low, f[inner], high]
    vals = np.interp(x, f, y)
    return float(np.trapezoid(vals, x))


def integrate_bands(frequency_Hz: np.ndarray, density: np.ndarray, bands: tuple[Band, ...] = PROJECT_BANDS) -> dict[str, float]:
    return {band.name: integrate_band(frequency_Hz, density, band) for band in bands}


def classify_band(band: Band, trusted_low_Hz: float, trusted_high_Hz: float) -> str:
    if band.f_high_Hz <= trusted_low_Hz or band.f_low_Hz >= trusted_high_Hz:
        return "UNTRUSTED"
    if band.f_low_Hz >= trusted_low_Hz and band.f_high_Hz <= trusted_high_Hz:
        return "TRUSTED"
    return "PARTIALLY_TRUSTED"

