from __future__ import annotations

import numpy as np

from .schema import SourceRecord, SourceSeries


def total_charge(record: SourceRecord) -> float:
    return float(np.sum(record.columns["rho"] * record.columns["cell_volume"], dtype=np.float64))


def current_moment(record: SourceRecord) -> np.ndarray:
    volume = record.columns["cell_volume"]
    return np.asarray(
        [
            np.sum(record.columns["Jx"] * volume, dtype=np.float64),
            np.sum(record.columns["Jy"] * volume, dtype=np.float64),
            np.sum(record.columns["Jz"] * volume, dtype=np.float64),
        ],
        dtype=float,
    )


def current_profile_z(record: SourceRecord, z_planes: np.ndarray) -> np.ndarray:
    z_planes = np.asarray(z_planes, dtype=float)
    z0 = record.columns["z_center"] - 0.5 * record.columns["dz"]
    z1 = record.columns["z_center"] + 0.5 * record.columns["dz"]
    area = record.columns["cell_volume"] / record.columns["dz"]
    out = np.zeros_like(z_planes, dtype=float)
    for i, zp in enumerate(z_planes):
        mask = (z0 <= zp) & (zp <= z1)
        out[i] = float(np.sum(record.columns["Jz"][mask] * area[mask], dtype=np.float64))
    return out


def frequency_metadata(series: SourceSeries) -> dict[str, float | bool | list[float] | None]:
    times = series.times
    duration = float(times[-1] - times[0]) if len(times) > 1 else 0.0
    dts = np.diff(times)
    uniform = bool(len(dts) > 0 and np.allclose(dts, dts[0], rtol=1e-9, atol=1e-18))
    output_dt = float(dts[0]) if uniform and len(dts) else None
    return {
        "actual_output_times_s": [float(t) for t in times],
        "uniform_output_interval": uniform,
        "output_dt_s": output_dt,
        "record_duration_s": duration,
        "nyquist_frequency_Hz": (0.5 / output_dt) if output_dt and output_dt > 0 else None,
        "raw_frequency_resolution_Hz": (1.0 / duration) if duration > 0 else None,
    }

