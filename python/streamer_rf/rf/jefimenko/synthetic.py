from __future__ import annotations

import math

import numpy as np

from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord, SourceSeries


def _cartesian_geometry(edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    x, y, z = np.meshgrid(centers, centers, centers, indexing="ij")
    dx, dy, dz = np.meshgrid(widths, widths, widths, indexing="ij")
    volume = dx * dy * dz
    return np.column_stack((x.ravel(), y.ravel(), z.ravel())), np.column_stack(
        (dx.ravel(), dy.ravel(), dz.ravel())
    ), volume.ravel(), centers


def make_source_record_from_arrays(
    xyz: np.ndarray,
    widths: np.ndarray,
    volume: np.ndarray,
    rho: np.ndarray,
    J: np.ndarray,
    *,
    time_s: float,
    case_id: str,
    source_definition: str,
) -> SourceRecord:
    meta = SourceMetadata(
        case_id=case_id,
        solver="synthetic",
        solver_version="stage-f2",
        time_s=time_s,
        coordinate_system="cartesian_3d_native_amr",
        pressure_Pa=101325.0,
        temperature_K=300.0,
        geometry_id="synthetic_compact_source",
        voltage_state="none",
        photoionization="off",
        source_definition=source_definition,
    )
    return SourceRecord(
        meta,
        {
            "cell_id": np.arange(volume.size),
            "level": np.zeros(volume.size, dtype=int),
            "x_center": xyz[:, 0],
            "y_center": xyz[:, 1],
            "z_center": xyz[:, 2],
            "dx": widths[:, 0],
            "dy": widths[:, 1],
            "dz": widths[:, 2],
            "cell_volume": volume,
            "rho": rho,
            "Jx": J[:, 0],
            "Jy": J[:, 1],
            "Jz": J[:, 2],
        },
    )


def static_gaussian_charge_series(
    *,
    total_charge_C: float = 1e-12,
    sigma_m: float = 1e-3,
    half_width_m: float = 5e-3,
    n: int = 25,
    times_s: tuple[float, ...] = (0.0, 1e-9, 2e-9),
) -> SourceSeries:
    edges = np.linspace(-half_width_m, half_width_m, n + 1)
    xyz, widths, volume, _ = _cartesian_geometry(edges)
    r2 = np.sum(xyz**2, axis=1)
    g = np.exp(-0.5 * r2 / sigma_m**2)
    g /= np.sum(g * volume)
    rho = total_charge_C * g
    J = np.zeros((volume.size, 3))
    records = tuple(
        make_source_record_from_arrays(
            xyz, widths, volume, rho.copy(), J.copy(), time_s=t, case_id="static_charge", source_definition="static Gaussian charge"
        )
        for t in times_s
    )
    return SourceSeries(records)


def steady_current_element_series(
    *,
    current_moment_Am: float = 1e-6,
    half_width_m: float = 5e-4,
    n: int = 5,
    times_s: tuple[float, ...] = (0.0, 1e-9, 2e-9),
) -> SourceSeries:
    edges = np.linspace(-half_width_m, half_width_m, n + 1)
    xyz, widths, volume, _ = _cartesian_geometry(edges)
    rho = np.zeros(volume.size)
    J = np.zeros((volume.size, 3))
    J[:, 2] = current_moment_Am / np.sum(volume)
    records = tuple(
        make_source_record_from_arrays(
            xyz, widths, volume, rho.copy(), J.copy(), time_s=t, case_id="steady_current", source_definition="steady compact current moment"
        )
        for t in times_s
    )
    return SourceSeries(records)


def dipole_radiation_series(
    *,
    p0_Cm: float = 1e-18,
    omega_rad_s: float = 2.0 * math.pi * 1.0e9,
    sigma_m: float = 5e-4,
    half_width_m: float = 3e-3,
    n: int = 17,
    times_s: np.ndarray | None = None,
) -> SourceSeries:
    if times_s is None:
        times_s = np.linspace(-1.5e-9, 1.5e-9, 61)
    edges = np.linspace(-half_width_m, half_width_m, n + 1)
    xyz, widths, volume, _ = _cartesian_geometry(edges)
    r2 = np.sum(xyz**2, axis=1)
    g = np.exp(-0.5 * r2 / sigma_m**2)
    g /= np.sum(g * volume)
    records = []
    for t in np.asarray(times_s, dtype=float):
        p = p0_Cm * math.sin(omega_rad_s * t)
        pdot = p0_Cm * omega_rad_s * math.cos(omega_rad_s * t)
        rho = p * xyz[:, 2] / sigma_m**2 * g
        J = np.zeros((volume.size, 3))
        J[:, 2] = pdot * g
        records.append(
            make_source_record_from_arrays(
                xyz,
                widths,
                volume,
                rho,
                J,
                time_s=float(t),
                case_id="smooth_electric_dipole",
                source_definition="rho=-p(t) grad_z(g), J=dp/dt g zhat",
            )
        )
    return SourceSeries(tuple(records))


def dipole_far_field_E(theta_rad: float, R_m: float, pddot: float) -> float:
    from streamer_rf.rf.jefimenko.constants import K_E, C0

    return K_E * abs(pddot) * abs(math.sin(theta_rad)) / (C0**2 * R_m)

