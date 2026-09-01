from __future__ import annotations

import numpy as np


def analytic_signal(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = x.size
    X = np.fft.fft(x)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1 : n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1 : (n + 1) // 2] = 2.0
    return np.fft.ifft(X * h)


def polarization_ellipse(e1: np.ndarray, e2: np.ndarray, *, far_field_valid: bool = True) -> dict[str, float | str | bool]:
    if not far_field_valid:
        return {"polarization_valid": False, "reason": "observer is not far-field valid"}
    a1 = analytic_signal(e1)
    a2 = analytic_signal(e2)
    S0 = float(np.mean(np.abs(a1) ** 2 + np.abs(a2) ** 2))
    S1 = float(np.mean(np.abs(a1) ** 2 - np.abs(a2) ** 2))
    cross = np.mean(a1 * np.conj(a2))
    S2 = float(2.0 * np.real(cross))
    S3 = float(-2.0 * np.imag(cross))
    psi = 0.5 * np.arctan2(S2, S1)
    s3n = np.clip(S3 / max(S0, 1e-300), -1.0, 1.0)
    chi = 0.5 * np.arcsin(s3n)
    dolp = float(np.sqrt(S1 * S1 + S2 * S2) / max(S0, 1e-300))
    if abs(S3) < 1e-12 * S0:
        handedness = "linear"
    elif S3 > 0:
        handedness = "left"
    else:
        handedness = "right"
    return {
        "polarization_valid": True,
        "major_axis_angle_rad": float(psi),
        "ellipticity": float(np.tan(chi)),
        "degree_of_linear_polarization": dolp,
        "handedness": handedness,
        "S0": S0,
        "S1": S1,
        "S2": S2,
        "S3": S3,
    }

