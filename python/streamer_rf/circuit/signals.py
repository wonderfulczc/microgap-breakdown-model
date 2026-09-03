from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


def _validate_time_values(time_s: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(values, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.shape != y.shape:
        raise ValueError("time and values must be matching 1D arrays")
    if t.size == 0:
        raise ValueError("at least one sample is required")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)):
        raise ValueError("time and values must be finite")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("time_s must be strictly increasing")
    return t, y


@dataclass(frozen=True)
class PiecewiseLinearSignal:
    time_s: np.ndarray
    values: np.ndarray
    extrapolation: str = "invalid"

    def __post_init__(self) -> None:
        t, y = _validate_time_values(self.time_s, self.values)
        if self.extrapolation not in {"invalid", "hold"}:
            raise ValueError("extrapolation must be 'invalid' or 'hold'")
        object.__setattr__(self, "time_s", t)
        object.__setattr__(self, "values", y)

    @classmethod
    def constant(cls, value: float) -> "PiecewiseLinearSignal":
        return cls(np.array([0.0]), np.array([float(value)]), extrapolation="hold")

    @classmethod
    def from_csv(cls, path: str | Path, value_column: str, *, extrapolation: str = "invalid") -> "PiecewiseLinearSignal":
        data = pd.read_csv(path)
        if "time_s" not in data or value_column not in data:
            raise ValueError(f"{path} must contain time_s and {value_column}")
        return cls(data["time_s"].to_numpy(dtype=float), data[value_column].to_numpy(dtype=float), extrapolation=extrapolation)

    @property
    def support(self) -> tuple[float, float]:
        return float(self.time_s[0]), float(self.time_s[-1])

    def require_supports(self, t0: float, t1: float) -> None:
        if self.extrapolation == "hold":
            return
        lo, hi = self.support
        if t0 < lo or t1 > hi:
            raise ValueError(f"time interval [{t0}, {t1}] exceeds signal support [{lo}, {hi}]")

    def __call__(self, time_s: float | np.ndarray) -> float | np.ndarray:
        t = np.asarray(time_s, dtype=float)
        if self.time_s.size == 1:
            out = np.full_like(t, self.values[0], dtype=float)
        else:
            if self.extrapolation == "invalid" and (np.any(t < self.time_s[0]) or np.any(t > self.time_s[-1])):
                raise ValueError("signal evaluation outside support")
            left = self.values[0] if self.extrapolation == "hold" else np.nan
            right = self.values[-1] if self.extrapolation == "hold" else np.nan
            out = np.interp(t, self.time_s, self.values, left=left, right=right)
        return float(out) if np.ndim(time_s) == 0 else out


@dataclass(frozen=True)
class ConductanceProfile:
    signal: PiecewiseLinearSignal
    source_kind: str = "synthetic"

    @classmethod
    def constant(cls, conductance_S: float, *, source_kind: str = "synthetic") -> "ConductanceProfile":
        if conductance_S < 0.0 or not np.isfinite(conductance_S):
            raise ValueError("conductance must be finite and nonnegative")
        return cls(PiecewiseLinearSignal.constant(conductance_S), source_kind=source_kind)

    @classmethod
    def from_samples(
        cls,
        time_s: np.ndarray,
        Gb_S: np.ndarray,
        *,
        extrapolation: str = "invalid",
        source_kind: str = "synthetic",
    ) -> "ConductanceProfile":
        if np.any(np.asarray(Gb_S, dtype=float) < 0.0):
            raise ValueError("Gb_S must be nonnegative")
        return cls(PiecewiseLinearSignal(time_s, Gb_S, extrapolation=extrapolation), source_kind=source_kind)

    @classmethod
    def from_csv(cls, path: str | Path, *, extrapolation: str = "invalid", source_kind: str = "csv") -> "ConductanceProfile":
        data = pd.read_csv(path)
        if "time_s" not in data:
            raise ValueError(f"{path} must contain time_s")
        if "Gb_S" in data:
            g = data["Gb_S"].to_numpy(dtype=float)
        elif "Rb_ohm" in data:
            r = data["Rb_ohm"].to_numpy(dtype=float)
            g = np.divide(1.0, r, out=np.zeros_like(r), where=np.isfinite(r) & (r > 0.0))
        else:
            raise ValueError(f"{path} must contain Gb_S or Rb_ohm")
        return cls.from_samples(data["time_s"].to_numpy(dtype=float), g, extrapolation=extrapolation, source_kind=source_kind)

    @property
    def support(self) -> tuple[float, float]:
        return self.signal.support

    def require_supports(self, t0: float, t1: float) -> None:
        self.signal.require_supports(t0, t1)

    def __call__(self, time_s: float | np.ndarray) -> float | np.ndarray:
        value = self.signal(time_s)
        if np.any(np.asarray(value) < -1e-30):
            raise ValueError("conductance profile evaluated to a negative value")
        return value
