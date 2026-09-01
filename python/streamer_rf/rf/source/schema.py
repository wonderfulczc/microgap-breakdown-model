from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


REQUIRED_COLUMNS = (
    "level",
    "x_center",
    "y_center",
    "z_center",
    "dx",
    "dy",
    "dz",
    "cell_volume",
    "rho",
    "Jx",
    "Jy",
    "Jz",
)


@dataclass(frozen=True)
class SourceMetadata:
    case_id: str
    solver: str
    solver_version: str
    time_s: float
    coordinate_system: str
    pressure_Pa: float
    temperature_K: float
    geometry_id: str
    voltage_state: str
    photoionization: str
    source_definition: str
    units: dict[str, str] = field(default_factory=dict)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "solver": self.solver,
            "solver_version": self.solver_version,
            "time_s": self.time_s,
            "coordinate_system": self.coordinate_system,
            "pressure_Pa": self.pressure_Pa,
            "temperature_K": self.temperature_K,
            "geometry_id": self.geometry_id,
            "voltage_state": self.voltage_state,
            "photoionization": self.photoionization,
            "source_definition": self.source_definition,
            "units": dict(self.units),
            "extra": dict(self.extra),
        }


@dataclass(frozen=True)
class SourceRecord:
    metadata: SourceMetadata
    columns: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        missing = [name for name in REQUIRED_COLUMNS if name not in self.columns]
        if missing:
            raise ValueError(f"missing SourceRecord columns: {missing}")

        n = len(np.asarray(self.columns[REQUIRED_COLUMNS[0]]))
        if n == 0:
            raise ValueError("SourceRecord must contain at least one cell")

        converted: dict[str, np.ndarray] = {}
        for name, values in self.columns.items():
            arr = np.asarray(values)
            if len(arr) != n:
                raise ValueError(f"column {name} has length {len(arr)}; expected {n}")
            if name != "cell_id" and not np.all(np.isfinite(arr)):
                raise ValueError(f"column {name} contains NaN or Inf")
            converted[name] = arr

        if np.any(converted["cell_volume"] <= 0):
            raise ValueError("cell_volume must be positive")
        for axis in ("dx", "dy", "dz"):
            if np.any(converted[axis] <= 0):
                raise ValueError(f"{axis} must be positive")

        object.__setattr__(self, "columns", converted)

    @property
    def n_cells(self) -> int:
        return len(self.columns["cell_volume"])

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        centers = np.column_stack(
            (self.columns["x_center"], self.columns["y_center"], self.columns["z_center"])
        )
        widths = np.column_stack((self.columns["dx"], self.columns["dy"], self.columns["dz"]))
        return centers - 0.5 * widths, centers + 0.5 * widths

    def geometry_like(self) -> "SourceRecord":
        zeros = np.zeros(self.n_cells)
        cols = {
            "level": self.columns["level"],
            "x_center": self.columns["x_center"],
            "y_center": self.columns["y_center"],
            "z_center": self.columns["z_center"],
            "dx": self.columns["dx"],
            "dy": self.columns["dy"],
            "dz": self.columns["dz"],
            "cell_volume": self.columns["cell_volume"],
            "rho": zeros,
            "Jx": zeros,
            "Jy": zeros,
            "Jz": zeros,
        }
        return SourceRecord(self.metadata, cols)

    def with_source_columns(self, rho: np.ndarray, J: np.ndarray, *, time_s: float | None = None) -> "SourceRecord":
        if J.shape != (self.n_cells, 3):
            raise ValueError("J must have shape (n_cells, 3)")
        meta = self.metadata
        if time_s is not None:
            meta = SourceMetadata(**{**meta.to_dict(), "time_s": float(time_s)})
        cols = dict(self.columns)
        cols["rho"] = np.asarray(rho, dtype=float)
        cols["Jx"] = np.asarray(J[:, 0], dtype=float)
        cols["Jy"] = np.asarray(J[:, 1], dtype=float)
        cols["Jz"] = np.asarray(J[:, 2], dtype=float)
        return SourceRecord(meta, cols)


@dataclass(frozen=True)
class SourceSeries:
    records: tuple[SourceRecord, ...]

    def __post_init__(self) -> None:
        if len(self.records) < 1:
            raise ValueError("SourceSeries requires at least one record")
        times = self.times
        if np.any(np.diff(times) <= 0):
            raise ValueError("SourceSeries times must be strictly increasing")

    @property
    def times(self) -> np.ndarray:
        return np.asarray([record.metadata.time_s for record in self.records], dtype=float)

    @property
    def has_three_snapshots(self) -> bool:
        return len(self.records) >= 3

