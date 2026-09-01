from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

from .schema import SourceMetadata, SourceRecord


AFIVO_STAGE_E_SOURCE_DEFINITION = (
    "AFIVO_STAGE_E_CELL_CENTERED_DRIFT_CURRENT: J=e*mu_e(E/N)*ne*E; "
    "diffusive finite-volume electron flux is not included"
)


def _parse_vector_config(path: Path, key: str) -> tuple[float, float, float]:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s+(.+)$")
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            vals = [float(x) for x in match.group(1).split()[:3]]
            if len(vals) != 3:
                raise ValueError(f"{key} in {path} does not contain three values")
            return vals[0], vals[1], vals[2]
    raise ValueError(f"{key} not found in {path}")


def _parse_scalar_config(path: Path, key: str) -> float:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s+(.+)$")
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1).split()[0])
    raise ValueError(f"{key} not found in {path}")


def infer_afivo_cell_widths_from_config(config_path: str | Path, levels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = Path(config_path)
    domain_len = np.asarray(_parse_vector_config(path, "domain_len"), dtype=float)
    box_size = int(round(_parse_scalar_config(path, "box_size")))
    scale = box_size * np.power(2.0, np.asarray(levels, dtype=float) - 1.0)
    return domain_len[0] / scale, domain_len[1] / scale, domain_len[2] / scale


def _infer_widths_from_coordinates(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    volume: np.ndarray,
    levels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx = np.empty_like(volume)
    dy = np.empty_like(volume)
    dz = np.empty_like(volume)
    for level in np.unique(levels):
        mask = levels == level
        widths = []
        for coord in (x[mask], y[mask], z[mask]):
            unique = np.unique(np.round(coord, decimals=18))
            diffs = np.diff(unique)
            diffs = diffs[diffs > 0]
            widths.append(float(np.min(diffs)) if diffs.size else float(np.cbrt(np.median(volume[mask]))))
        dx[mask], dy[mask], dz[mask] = widths
    inferred = dx * dy * dz
    mismatch = np.abs(inferred - volume) / np.maximum(volume, 1e-300)
    bad = mismatch > 1e-6
    if np.any(bad):
        # Preserve the CSV volume exactly for integrals and distribute any
        # residual anisotropic uncertainty to dz, which is only used for
        # geometric audit helpers.
        dz[bad] = volume[bad] / (dx[bad] * dy[bad])
    return dx, dy, dz


def load_afivo_stage_e_source_csv(
    csv_path: str | Path,
    config_path: str | Path,
    *,
    case_id: str,
    solver_version: str,
    geometry_id: str,
    voltage_state: str,
    photoionization: str,
    max_rows: int | None = None,
    row_indices: np.ndarray | None = None,
) -> SourceRecord:
    csv_path = Path(csv_path)
    data: dict[str, list[float]] = {}
    index_set = None if row_indices is None else {int(i) for i in np.asarray(row_indices, dtype=int)}
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for n, row in enumerate(reader):
            if index_set is not None and n not in index_set:
                continue
            if max_rows is not None and n >= max_rows:
                break
            for key, value in row.items():
                data.setdefault(key, []).append(float(value))
            if index_set is not None and len(data.get("time_s", [])) == len(index_set):
                break
    if not data:
        raise ValueError(f"no source rows selected from {csv_path}")

    levels = np.asarray(data["level"], dtype=int)
    volume = np.asarray(data["cell_volume_m3"], dtype=float)
    x = np.asarray(data["x_m"], dtype=float)
    y = np.asarray(data["y_m"], dtype=float)
    z = np.asarray(data["z_m"], dtype=float)
    try:
        dx, dy, dz = infer_afivo_cell_widths_from_config(config_path, levels)
        rel = np.max(np.abs(dx * dy * dz - volume) / np.maximum(volume, 1e-300))
        if rel > 1e-6:
            dx, dy, dz = _infer_widths_from_coordinates(x, y, z, volume, levels)
    except ValueError:
        dx, dy, dz = _infer_widths_from_coordinates(x, y, z, volume, levels)

    time_s = float(np.asarray(data["time_s"])[0])
    meta = SourceMetadata(
        case_id=case_id,
        solver="afivo-streamer",
        solver_version=solver_version,
        time_s=time_s,
        coordinate_system="cartesian_3d_native_amr",
        pressure_Pa=101325.0,
        temperature_K=300.0,
        geometry_id=geometry_id,
        voltage_state=voltage_state,
        photoionization=photoionization,
        source_definition=AFIVO_STAGE_E_SOURCE_DEFINITION,
        units={
            "rho": "C m^-3",
            "J": "A m^-2",
            "cell_volume": "m^3",
            "E": "V m^-1",
            "ne": "m^-3",
        },
        extra={
            "source_file": str(csv_path),
            "config_file": str(config_path),
            "J_transport_available": False,
            "J_RF_reference_definition": "-e * Gamma_e using the finite-volume electron transport flux",
        },
    )
    cols = {
        "cell_id": np.arange(len(volume)),
        "level": levels,
        "x_center": x,
        "y_center": y,
        "z_center": z,
        "dx": dx,
        "dy": dy,
        "dz": dz,
        "cell_volume": volume,
        "rho": np.asarray(data["rho_Cpm3"], dtype=float),
        "Jx": np.asarray(data["Jx_Apm2"], dtype=float),
        "Jy": np.asarray(data["Jy_Apm2"], dtype=float),
        "Jz": np.asarray(data["Jz_Apm2"], dtype=float),
        "ne": np.asarray(data["ne_m3"], dtype=float),
        "Ex": np.asarray(data["Ex_Vpm"], dtype=float),
        "Ey": np.asarray(data["Ey_Vpm"], dtype=float),
        "Ez": np.asarray(data["Ez_Vpm"], dtype=float),
    }
    return SourceRecord(meta, cols)


def source_record_from_petsc_axisymmetric(
    r_centers: np.ndarray,
    z_centers: np.ndarray,
    dr: float,
    dz: float,
    rho: np.ndarray,
    Jr: np.ndarray,
    Jz: np.ndarray,
    metadata: SourceMetadata,
) -> SourceRecord:
    rr, zz = np.meshgrid(r_centers, z_centers, indexing="ij")
    volume = 2.0 * np.pi * rr.ravel() * dr * dz
    cols = {
        "cell_id": np.arange(volume.size),
        "level": np.zeros(volume.size, dtype=int),
        "x_center": rr.ravel(),
        "y_center": np.zeros(volume.size),
        "z_center": zz.ravel(),
        "dx": np.full(volume.size, dr),
        "dy": np.maximum(2.0 * np.pi * rr.ravel(), dr),
        "dz": np.full(volume.size, dz),
        "cell_volume": volume,
        "rho": rho.ravel(),
        "Jx": Jr.ravel(),
        "Jy": np.zeros(volume.size),
        "Jz": Jz.ravel(),
    }
    return SourceRecord(metadata, cols)
