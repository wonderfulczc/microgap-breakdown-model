from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CWTResult:
    time_s: np.ndarray
    frequency_Hz: np.ndarray
    coefficients: np.ndarray
    power: np.ndarray
    coi_mask: np.ndarray


def morlet_cwt(
    time_s: np.ndarray,
    signal: np.ndarray,
    frequency_Hz: np.ndarray,
    *,
    n_cycles: float = 6.0,
) -> CWTResult:
    t = np.asarray(time_s, dtype=float)
    x = np.asarray(signal, dtype=float)
    f = np.asarray(frequency_Hz, dtype=float)
    if t.ndim != 1 or x.shape != t.shape:
        raise ValueError("time and signal must be matching 1D arrays")
    dt = float(np.diff(t)[0])
    if not np.allclose(np.diff(t), dt, rtol=1e-9, atol=abs(dt) * 1e-12):
        raise ValueError("CWT requires uniform time")
    coeff = np.zeros((f.size, t.size), dtype=complex)
    half = t.size // 2
    tau = (np.arange(-half, t.size - half) * dt).astype(float)
    X = np.fft.fft(x, n=2 * t.size)
    for i, freq in enumerate(f):
        sigma_t = n_cycles / (2.0 * np.pi * freq)
        wavelet = np.exp(2j * np.pi * freq * tau) * np.exp(-0.5 * (tau / sigma_t) ** 2)
        wavelet /= np.sqrt(np.sum(np.abs(wavelet) ** 2))
        conv = np.fft.ifft(X * np.fft.fft(np.conj(wavelet[::-1]), n=2 * t.size))
        coeff[i] = conv[half : half + t.size]
    distance_to_edge = np.minimum(t - t[0], t[-1] - t)
    coi = np.zeros((f.size, t.size), dtype=bool)
    for i, freq in enumerate(f):
        sigma_t = n_cycles / (2.0 * np.pi * freq)
        coi[i] = distance_to_edge >= np.sqrt(2.0) * sigma_t
    return CWTResult(time_s=t, frequency_Hz=f, coefficients=coeff, power=np.abs(coeff) ** 2, coi_mask=coi)


def cwt_ridge(result: CWTResult) -> tuple[np.ndarray, np.ndarray]:
    masked = np.where(result.coi_mask, result.power, -np.inf)
    idx = np.argmax(masked, axis=0)
    valid = np.isfinite(masked[idx, np.arange(result.time_s.size)])
    ridge = result.frequency_Hz[idx]
    ridge[~valid] = np.nan
    return ridge, valid


def band_energy_vs_time(result: CWTResult, bands: dict[str, tuple[float, float]]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for name, (low, high) in bands.items():
        mask = (result.frequency_Hz >= low) & (result.frequency_Hz <= high)
        power = np.where(result.coi_mask[mask], result.power[mask], 0.0)
        out[name] = np.trapezoid(power, result.frequency_Hz[mask], axis=0) if np.count_nonzero(mask) > 1 else np.zeros(result.time_s.size)
    return out

