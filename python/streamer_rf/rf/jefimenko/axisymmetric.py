from __future__ import annotations

import numpy as np

from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord, SourceSeries


def rotate_axisymmetric_record(record: SourceRecord, n_phi: int) -> SourceRecord:
    if record.metadata.coordinate_system != "axisymmetric_rz":
        raise ValueError("record must use coordinate_system='axisymmetric_rz'")
    phi = (np.arange(n_phi) + 0.5) * (2.0 * np.pi / n_phi)
    r = record.columns["x_center"]
    z = record.columns["z_center"]
    volume = record.columns["cell_volume"] / n_phi
    cols: dict[str, list[np.ndarray]] = {k: [] for k in ("x", "y", "z", "dx", "dy", "dz", "vol", "rho", "Jx", "Jy", "Jz")}
    for p in phi:
        cp = np.cos(p)
        sp = np.sin(p)
        cols["x"].append(r * cp)
        cols["y"].append(r * sp)
        cols["z"].append(z)
        equiv = np.cbrt(volume)
        cols["dx"].append(equiv)
        cols["dy"].append(equiv)
        cols["dz"].append(record.columns["dz"])
        cols["vol"].append(volume)
        cols["rho"].append(record.columns["rho"])
        Jr = record.columns["Jx"]
        cols["Jx"].append(Jr * cp)
        cols["Jy"].append(Jr * sp)
        cols["Jz"].append(record.columns["Jz"])
    n = record.n_cells * n_phi
    meta = SourceMetadata(**{**record.metadata.to_dict(), "coordinate_system": "cartesian_3d_phi_quadrature"})
    return SourceRecord(
        meta,
        {
            "cell_id": np.arange(n),
            "level": np.zeros(n, dtype=int),
            "x_center": np.concatenate(cols["x"]),
            "y_center": np.concatenate(cols["y"]),
            "z_center": np.concatenate(cols["z"]),
            "dx": np.concatenate(cols["dx"]),
            "dy": np.concatenate(cols["dy"]),
            "dz": np.concatenate(cols["dz"]),
            "cell_volume": np.concatenate(cols["vol"]),
            "rho": np.concatenate(cols["rho"]),
            "Jx": np.concatenate(cols["Jx"]),
            "Jy": np.concatenate(cols["Jy"]),
            "Jz": np.concatenate(cols["Jz"]),
        },
    )


def rotate_axisymmetric_series(series: SourceSeries, n_phi: int) -> SourceSeries:
    return SourceSeries(tuple(rotate_axisymmetric_record(record, n_phi) for record in series.records))

