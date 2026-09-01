from __future__ import annotations

import numpy as np

from streamer_rf.rf.source.derivatives import central_difference_three_point
from streamer_rf.rf.source.schema import SourceSeries


SOURCE_NAMES = ("rho", "Jx", "Jy", "Jz")


def require_matching_geometry(series: SourceSeries) -> None:
    base = series.records[0]
    keys = ("x_center", "y_center", "z_center", "dx", "dy", "dz", "cell_volume")
    for record in series.records[1:]:
        if record.n_cells != base.n_cells:
            raise ValueError("SourceSeries records must share a remapped common partition")
        for key in keys:
            if not np.allclose(record.columns[key], base.columns[key], rtol=1e-12, atol=1e-18):
                raise ValueError("SourceSeries records must share a remapped common partition")


def stacked_sources(series: SourceSeries) -> tuple[np.ndarray, np.ndarray]:
    rho = np.vstack([record.columns["rho"] for record in series.records])
    J = np.stack(
        [
            np.column_stack((record.columns["Jx"], record.columns["Jy"], record.columns["Jz"]))
            for record in series.records
        ],
        axis=0,
    )
    return rho, J


def central_derivative_series(times: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(times) < 3:
        raise ValueError("at least three source snapshots are required for central derivatives")
    derivatives = []
    for i in range(1, len(times) - 1):
        derivatives.append(central_difference_three_point(times[i - 1 : i + 2], values[i - 1 : i + 2]))
    return times[1:-1].copy(), np.asarray(derivatives)


def linear_interpolate_time(times: np.ndarray, values: np.ndarray, query_times: np.ndarray) -> np.ndarray:
    times = np.asarray(times, dtype=float)
    query_times = np.asarray(query_times, dtype=float)
    idx = np.searchsorted(times, query_times, side="right") - 1
    idx = np.clip(idx, 0, len(times) - 2)
    t0 = times[idx]
    t1 = times[idx + 1]
    frac = (query_times - t0) / (t1 - t0)
    if values.ndim == 2:
        return values[idx, np.arange(query_times.size)] * (1.0 - frac) + values[idx + 1, np.arange(query_times.size)] * frac
    if values.ndim == 3:
        return values[idx, np.arange(query_times.size), :] * (1.0 - frac[:, None]) + values[
            idx + 1, np.arange(query_times.size), :
        ] * frac[:, None]
    raise ValueError("values must have shape (nt, ncell) or (nt, ncell, 3)")


class RetardedSourceInterpolator:
    def __init__(self, series: SourceSeries):
        require_matching_geometry(series)
        self.series = series
        self.times = series.times
        self.rho, self.J = stacked_sources(series)
        self.drho_times, self.drho = central_derivative_series(self.times, self.rho)
        self.dJ_times, self.dJ = central_derivative_series(self.times, self.J)

    @property
    def derivative_support(self) -> tuple[float, float]:
        return float(self.drho_times[0]), float(self.drho_times[-1])

    def evaluate(self, retarded_times: np.ndarray) -> dict[str, np.ndarray | float | bool]:
        tr = np.asarray(retarded_times, dtype=float)
        tmin, tmax = self.derivative_support
        valid = (tr >= tmin) & (tr <= tmax)
        valid_fraction = float(np.count_nonzero(valid) / tr.size)
        safe_tr = np.clip(tr, tmin, tmax)
        return {
            "rho": linear_interpolate_time(self.times, self.rho, safe_tr),
            "J": linear_interpolate_time(self.times, self.J, safe_tr),
            "drho_dt": linear_interpolate_time(self.drho_times, self.drho, safe_tr),
            "dJ_dt": linear_interpolate_time(self.dJ_times, self.dJ, safe_tr),
            "valid_mask": valid,
            "valid_fraction": valid_fraction,
            "all_valid": bool(np.all(valid)),
        }

