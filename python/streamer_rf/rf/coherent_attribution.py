from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from streamer_rf.rf.jefimenko.diagnostics import current_moment_radiation_approx


@dataclass(frozen=True)
class CoherentAttributionResult:
    P_reconstructed: float
    coherent_contributions: dict[str, float]
    eta: dict[str, float]
    eta_sum: float
    eta_sum_error: float
    interference_role: dict[str, str]
    trusted_count: int
    trusted_fraction: float
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def mechanism_radiative_fields(dMdt_terms: dict[str, np.ndarray], observer_vector_m: np.ndarray) -> dict[str, np.ndarray]:
    fields: dict[str, np.ndarray] = {}
    for name, vals in dMdt_terms.items():
        arr = np.asarray(vals, dtype=float)
        if arr.ndim == 1:
            arr = np.column_stack((np.zeros(arr.size), np.zeros(arr.size), arr))
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError("dMdt terms must be vector time series")
        fields[name] = np.vstack([current_moment_radiation_approx(row, observer_vector_m) for row in arr])
    return fields


def uniform_resample(
    time_s: np.ndarray,
    values: np.ndarray,
    valid: np.ndarray,
    *,
    dt_s: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    mask = np.asarray(valid, dtype=bool)
    if t.ndim != 1 or y.shape[0] != t.size or mask.shape != (t.size,):
        raise ValueError("time, values, and valid mask dimensions do not match")
    mask &= np.all(np.isfinite(y), axis=tuple(range(1, y.ndim))) if y.ndim > 1 else np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        raise ValueError("at least three valid samples are required")
    tv = t[mask]
    if np.any(np.diff(tv) <= 0.0):
        raise ValueError("valid times must be strictly increasing")
    dt = float(dt_s) if dt_s is not None else float(np.median(np.diff(tv)))
    if not (dt > 0.0 and math.isfinite(dt)):
        raise ValueError("invalid resampling interval")
    n = max(3, int(math.floor((tv[-1] - tv[0]) / dt)) + 1)
    tu = tv[0] + np.arange(n, dtype=float) * ((tv[-1] - tv[0]) / (n - 1))
    if y.ndim == 1:
        yu = np.interp(tu, tv, y[mask])
    else:
        flat = y.reshape((t.size, -1))
        out = np.column_stack([np.interp(tu, tv, flat[mask, i]) for i in range(flat.shape[1])])
        yu = out.reshape((tu.size, *y.shape[1:]))
    return tu, yu


def coherent_attribution(
    mechanism_coefficients: dict[str, np.ndarray],
    reconstructed_coefficients: np.ndarray,
    mask: np.ndarray,
    *,
    p_floor: float = 0.0,
) -> CoherentAttributionResult:
    rec = np.asarray(reconstructed_coefficients, dtype=complex)
    m = np.asarray(mask, dtype=bool)
    if rec.shape != m.shape:
        raise ValueError("mask and reconstructed coefficients must have identical shape")
    coeffs: dict[str, np.ndarray] = {}
    for name, val in mechanism_coefficients.items():
        arr = np.asarray(val, dtype=complex)
        if arr.shape != rec.shape:
            raise ValueError("all mechanism coefficient arrays must share reconstructed shape")
        coeffs[name] = arr
    count = int(np.count_nonzero(m))
    total = int(m.size)
    if count == 0:
        return CoherentAttributionResult(math.nan, {}, {}, math.nan, math.nan, {}, 0, 0.0, "NO_TRUSTED_COEFFICIENTS")
    P = float(np.sum(np.abs(rec[m]) ** 2, dtype=np.float64))
    if not math.isfinite(P) or P <= p_floor:
        return CoherentAttributionResult(P, {}, {}, math.nan, math.nan, {}, count, count / total, "ZERO_OR_NEAR_ZERO_RECONSTRUCTED_POWER")
    contributions = {
        name: float(np.sum(np.real(arr[m] * np.conj(rec[m])), dtype=np.float64))
        for name, arr in coeffs.items()
    }
    eta = {name: val / P for name, val in contributions.items()}
    eta_sum = float(sum(eta.values()))
    roles = {name: _role(val) for name, val in contributions.items()}
    return CoherentAttributionResult(
        P_reconstructed=P,
        coherent_contributions=contributions,
        eta=eta,
        eta_sum=eta_sum,
        eta_sum_error=abs(eta_sum - 1.0),
        interference_role=roles,
        trusted_count=count,
        trusted_fraction=count / total,
        status="OK",
    )


def cwt_linearity_error(sum_coefficients: np.ndarray, direct_coefficients: np.ndarray, mask: np.ndarray) -> float:
    m = np.asarray(mask, dtype=bool)
    a = np.asarray(sum_coefficients, dtype=complex)
    b = np.asarray(direct_coefficients, dtype=complex)
    if a.shape != b.shape or a.shape != m.shape:
        raise ValueError("CWT arrays and mask must have matching shapes")
    if not np.any(m):
        return math.nan
    return float(np.linalg.norm((a - b)[m]) / max(np.linalg.norm(b[m]), 1e-300))


def _role(contribution: float) -> str:
    if contribution > 0.0:
        return "CONSTRUCTIVE_NET"
    if contribution < 0.0:
        return "DESTRUCTIVE_NET"
    return "NEAR_ZERO_NET"
