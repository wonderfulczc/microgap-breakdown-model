from __future__ import annotations

import numpy as np

from .schema import SourceMetadata, SourceRecord, SourceSeries


def make_cartesian_record(
    edges_x: np.ndarray,
    edges_y: np.ndarray,
    edges_z: np.ndarray,
    *,
    time_s: float,
    rho: np.ndarray,
    J: np.ndarray,
    case_id: str = "synthetic",
) -> SourceRecord:
    cx = 0.5 * (edges_x[:-1] + edges_x[1:])
    cy = 0.5 * (edges_y[:-1] + edges_y[1:])
    cz = 0.5 * (edges_z[:-1] + edges_z[1:])
    dx = np.diff(edges_x)
    dy = np.diff(edges_y)
    dz = np.diff(edges_z)
    xx, yy, zz = np.meshgrid(cx, cy, cz, indexing="ij")
    dxx, dyy, dzz = np.meshgrid(dx, dy, dz, indexing="ij")
    volume = dxx * dyy * dzz
    meta = SourceMetadata(
        case_id=case_id,
        solver="synthetic",
        solver_version="stage-f1",
        time_s=time_s,
        coordinate_system="cartesian_3d_native_amr",
        pressure_Pa=101325.0,
        temperature_K=300.0,
        geometry_id="analytic_box",
        voltage_state="none",
        photoionization="off",
        source_definition="analytic moving Gaussian with J=rho*v",
    )
    cols = {
        "cell_id": np.arange(volume.size),
        "level": np.zeros(volume.size, dtype=int),
        "x_center": xx.ravel(),
        "y_center": yy.ravel(),
        "z_center": zz.ravel(),
        "dx": dxx.ravel(),
        "dy": dyy.ravel(),
        "dz": dzz.ravel(),
        "cell_volume": volume.ravel(),
        "rho": rho.ravel(),
        "Jx": J[..., 0].ravel(),
        "Jy": J[..., 1].ravel(),
        "Jz": J[..., 2].ravel(),
    }
    return SourceRecord(meta, cols)


def moving_gaussian_on_edges(
    edges_x: np.ndarray,
    edges_y: np.ndarray,
    edges_z: np.ndarray,
    *,
    time_s: float,
    velocity: tuple[float, float, float] = (0.2, 0.0, 0.0),
    sigma: float = 0.45,
    rho0: float = 1.0,
) -> SourceRecord:
    cx = 0.5 * (edges_x[:-1] + edges_x[1:])
    cy = 0.5 * (edges_y[:-1] + edges_y[1:])
    cz = 0.5 * (edges_z[:-1] + edges_z[1:])
    xx, yy, zz = np.meshgrid(cx, cy, cz, indexing="ij")
    v = np.asarray(velocity, dtype=float)
    center = np.asarray([-0.25, 0.0, 0.0], dtype=float) + v * time_s
    r2 = (xx - center[0]) ** 2 + (yy - center[1]) ** 2 + (zz - center[2]) ** 2
    rho = rho0 * np.exp(-0.5 * r2 / sigma**2)
    J = np.zeros(rho.shape + (3,), dtype=float)
    J[..., 0] = rho * v[0]
    J[..., 1] = rho * v[1]
    J[..., 2] = rho * v[2]
    return make_cartesian_record(edges_x, edges_y, edges_z, time_s=time_s, rho=rho, J=J)


def moving_gaussian_series(edges: tuple[np.ndarray, np.ndarray, np.ndarray], times: tuple[float, float, float]) -> SourceSeries:
    return SourceSeries(tuple(moving_gaussian_on_edges(*edges, time_s=t) for t in times))

