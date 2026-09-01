from __future__ import annotations

import numpy as np


def multitone(time_s: np.ndarray, tones: dict[float, float]) -> np.ndarray:
    t = np.asarray(time_s, dtype=float)
    y = np.zeros_like(t)
    for freq, amp in tones.items():
        y += amp * np.sin(2.0 * np.pi * freq * t)
    return y


def chirp_signal(time_s: np.ndarray, f0: float, f1: float) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time_s, dtype=float)
    tau = (t - t[0]) / (t[-1] - t[0])
    inst = f0 + (f1 - f0) * tau
    phase = 2.0 * np.pi * (f0 * (t - t[0]) + 0.5 * (f1 - f0) * (t - t[0]) * tau)
    return np.sin(phase), inst


def staged_bursts(time_s: np.ndarray) -> tuple[np.ndarray, dict[str, tuple[float, float]], dict[str, float]]:
    t = np.asarray(time_s, dtype=float)
    centers = {"VHF": 15e-9, "UHF": 35e-9, "GHz": 55e-9}
    freqs = {"VHF": 80e6, "UHF": 400e6, "GHz": 2e9}
    y = np.zeros_like(t)
    sigma = 4e-9
    for name, center in centers.items():
        env = np.exp(-0.5 * ((t - center) / sigma) ** 2)
        y += env * np.sin(2.0 * np.pi * freqs[name] * t)
    bands = {"VHF": (60e6, 100e6), "UHF": (300e6, 500e6), "GHz": (1.5e9, 2.5e9)}
    return y, bands, centers

