from __future__ import annotations

import numpy as np

from streamer_rf.rf.jefimenko.constants import C0, K_E
from streamer_rf.rf.source.integrals import current_moment
from streamer_rf.rf.source.schema import SourceSeries


def source_characteristic_size(series: SourceSeries) -> float:
    record = series.records[0]
    xyz = np.column_stack(
        (record.columns["x_center"], record.columns["y_center"], record.columns["z_center"])
    )
    return float(np.max(np.linalg.norm(xyz - np.mean(xyz, axis=0), axis=1)))


def near_far_audit(series: SourceSeries, observer_position: np.ndarray, source_timescale_s: float) -> dict[str, float]:
    record = series.records[0]
    xyz = np.column_stack(
        (record.columns["x_center"], record.columns["y_center"], record.columns["z_center"])
    )
    center = np.average(xyz, axis=0, weights=record.columns["cell_volume"])
    L = source_characteristic_size(series)
    R = float(np.linalg.norm(np.asarray(observer_position, dtype=float) - center))
    return {
        "source_characteristic_size_m": L,
        "observer_distance_m": R,
        "source_timescale_s": source_timescale_s,
        "epsilon_geom_L_over_R": L / R,
        "epsilon_EM_L_over_ctau": L / (C0 * source_timescale_s),
    }


def current_moment_radiation_approx(dMdt: np.ndarray, observer_vector: np.ndarray) -> np.ndarray:
    R = float(np.linalg.norm(observer_vector))
    rhat = observer_vector / R
    transverse = dMdt - np.dot(dMdt, rhat) * rhat
    return -K_E * transverse / (C0**2 * R)


def current_moment_series(series: SourceSeries) -> np.ndarray:
    return np.vstack([current_moment(record) for record in series.records])

