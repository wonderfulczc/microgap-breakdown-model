from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WindowInfo:
    name: str
    coherent_gain: float
    energy_correction: float


@dataclass(frozen=True)
class Spectrum:
    frequency_Hz: np.ndarray
    complex_spectrum: np.ndarray
    one_sided_esd: np.ndarray
    dt_s: float
    n_samples: int
    duration_s: float
    df_Hz: float
    nyquist_Hz: float
    window: WindowInfo
    component_names: tuple[str, ...]

    @property
    def amplitude(self) -> np.ndarray:
        weights = one_sided_weights(self.n_samples)
        return weights[:, None] * np.abs(self.complex_spectrum) / max(self.window.coherent_gain, 1e-300)

    @property
    def phase_rad(self) -> np.ndarray:
        return np.angle(self.complex_spectrum)


def validate_uniform_time(time_s: np.ndarray) -> float:
    t = np.asarray(time_s, dtype=float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("time_s must be a 1D array with at least two samples")
    dt = np.diff(t)
    if not np.allclose(dt, dt[0], rtol=1e-9, atol=max(1e-30, abs(dt[0]) * 1e-12)):
        raise ValueError("FFT requires uniformly sampled physical time")
    return float(dt[0])


def make_window(n: int, name: str = "hann") -> tuple[np.ndarray, WindowInfo]:
    if name == "rectangular":
        w = np.ones(n)
    elif name == "hann":
        w = np.hanning(n)
    else:
        raise ValueError("window must be 'rectangular' or 'hann'")
    coherent_gain = float(np.mean(w))
    mean_square = float(np.mean(w * w))
    energy_correction = 1.0 / mean_square if mean_square > 0 else np.inf
    return w, WindowInfo(name=name, coherent_gain=coherent_gain, energy_correction=energy_correction)


def one_sided_weights(n: int) -> np.ndarray:
    weights = np.ones(n // 2 + 1)
    if n % 2 == 0:
        weights[1:-1] = 2.0
    else:
        weights[1:] = 2.0
    return weights


def compute_one_sided_spectrum(
    time_s: np.ndarray,
    values: np.ndarray,
    *,
    window: str = "hann",
    component_names: tuple[str, ...] | None = None,
    apply_window_energy_correction: bool = True,
) -> Spectrum:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.shape[0] != t.size:
        raise ValueError("values first dimension must match time_s")
    dt = validate_uniform_time(t)
    n = t.size
    w, info = make_window(n, window)
    yw = y * w[:, None]
    z = dt * np.fft.rfft(yw, axis=0) * np.exp(-2j * np.pi * np.fft.rfftfreq(n, dt)[:, None] * t[0])
    freq = np.fft.rfftfreq(n, dt)
    weights = one_sided_weights(n)[:, None]
    esd = weights * np.abs(z) ** 2
    if apply_window_energy_correction:
        esd *= info.energy_correction
    duration = float((n - 1) * dt)
    df = float(freq[1] - freq[0]) if freq.size > 1 else np.inf
    names = component_names or tuple(f"c{i}" for i in range(y.shape[1]))
    return Spectrum(
        frequency_Hz=freq,
        complex_spectrum=z,
        one_sided_esd=esd,
        dt_s=dt,
        n_samples=n,
        duration_s=duration,
        df_Hz=df,
        nyquist_Hz=float(freq[-1]),
        window=info,
        component_names=tuple(names),
    )


def parseval_time_energy(time_s: np.ndarray, values: np.ndarray) -> np.ndarray:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    dt = validate_uniform_time(t)
    return np.sum(y * y, axis=0) * dt


def parseval_spectral_energy(spectrum: Spectrum) -> np.ndarray:
    return np.sum(spectrum.one_sided_esd, axis=0) * spectrum.df_Hz
