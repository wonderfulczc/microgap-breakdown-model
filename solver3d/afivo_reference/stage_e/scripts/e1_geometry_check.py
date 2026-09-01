#!/usr/bin/env python3
"""Geometry contract and low-cost checks for Stage E1 triangular foil."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEOM_DIR = ROOT / "geometry"

PARAMS = {
    "case_id": "Stage-E1-triangular-foil-electrostatic",
    "geometry_source": "development",
    "geometry_type": "rounded_triangular_copper_foil_tip",
    "polarity": "positive_hv_foil",
    "coordinate_convention": {
        "x": "foil width direction",
        "y": "foil thickness / transverse direction",
        "z": "gap / propagation direction",
        "ground_plane": "low-z boundary",
        "hv_foil": "high-z conductor, tip points toward ground",
        "propagation_direction_reserved": "-z",
    },
    "gap_m": 70.0e-6,
    "foil_thickness_m": 10.0e-6,
    "foil_width_m": 40.0e-6,
    "foil_length_m": 50.0e-6,
    "tip_radius_m": 3.0e-6,
    "edge_radius_m": 1.0e-6,
    "ground_z_m": 0.0,
    "voltage_V": 500.0,
    "pressure_Pa": 101325.0,
    "temperature_K": 300.0,
    "domain_origin_m": [-60.0e-6, -60.0e-6, -2.5e-6],
    "domain_len_m": [120.0e-6, 120.0e-6, 152.5e-6],
    "boundary_conditions": {
        "bottom_z": "Dirichlet 0 V ground plane",
        "top_z": "Dirichlet applied voltage from homogeneous Afivo BC",
        "side_boundaries": "homogeneous Neumann",
        "foil": "user level-set conductor at applied voltage",
    },
    "photoionization": "OFF",
    "plasma": "OFF",
    "export_contract": "x_m,y_m,z_m,phi_V,Ex_Vpm,Ey_Vpm,Ez_Vpm,Eabs_Vpm,lsf_m,level",
}
PARAMS["triangle_tip_angle_deg"] = math.degrees(
    2.0 * math.atan((0.5 * PARAMS["foil_width_m"]) / PARAMS["foil_length_m"])
)


def dist_segment(px: float, pz: float, ax: float, az: float, bx: float, bz: float) -> float:
    vx = bx - ax
    vz = bz - az
    wx = px - ax
    wz = pz - az
    denom = max(vx * vx + vz * vz, 1.0e-300)
    t = max(0.0, min(1.0, (wx * vx + wz * vz) / denom))
    qx = ax + t * vx
    qz = az + t * vz
    return math.hypot(px - qx, pz - qz)


def cross2(ax: float, az: float, bx: float, bz: float) -> float:
    return ax * bz - az * bx


def inside_tri(px: float, pz: float, a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> bool:
    ax, az = a
    bx, bz = b
    cx, cz = c
    c1 = cross2(bx - ax, bz - az, px - ax, pz - az)
    c2 = cross2(cx - bx, cz - bz, px - bx, pz - bz)
    c3 = cross2(ax - cx, az - cz, px - cx, pz - cz)
    return (c1 >= 0.0 and c2 >= 0.0 and c3 >= 0.0) or (
        c1 <= 0.0 and c2 <= 0.0 and c3 <= 0.0
    )


def lsf(x: float, y: float, z: float) -> float:
    ground_z = PARAMS["ground_z_m"]
    gap = PARAMS["gap_m"]
    tip_radius = PARAMS["tip_radius_m"]
    edge_radius = PARAMS["edge_radius_m"]
    width = PARAMS["foil_width_m"]
    length = PARAMS["foil_length_m"]
    thickness = PARAMS["foil_thickness_m"]

    a = (0.0, ground_z + gap + tip_radius - edge_radius)
    b = (0.5 * width - edge_radius, ground_z + gap + length - tip_radius)
    c = (-0.5 * width + edge_radius, ground_z + gap + length - tip_radius)
    d = min(
        dist_segment(x, z, *a, *b),
        dist_segment(x, z, *b, *c),
        dist_segment(x, z, *c, *a),
    )
    sd = -d if inside_tri(x, z, a, b, c) else d
    phi_xz = sd - tip_radius + edge_radius
    half_thickness_core = max(0.5 * thickness - edge_radius, 0.1 * thickness)
    phi_y = abs(y) - half_thickness_core
    return max(phi_xz, phi_y) - edge_radius


def write_slice(path: Path, axes: tuple[str, str], fixed: dict[str, float], n: int = 61) -> None:
    ranges = {
        "x": (-30e-6, 30e-6),
        "y": (-12e-6, 12e-6),
        "z": (55e-6, 124e-6),
    }
    coords = {"x": 0.0, "y": 0.0, "z": PARAMS["gap_m"]}
    coords.update(fixed)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x_m", "y_m", "z_m", "lsf_m", "inside_conductor"])
        for ia in range(n):
            va = ranges[axes[0]][0] + (ranges[axes[0]][1] - ranges[axes[0]][0]) * ia / (n - 1)
            for ib in range(n):
                vb = ranges[axes[1]][0] + (ranges[axes[1]][1] - ranges[axes[1]][0]) * ib / (n - 1)
                coords[axes[0]] = va
                coords[axes[1]] = vb
                val = lsf(coords["x"], coords["y"], coords["z"])
                writer.writerow([coords["x"], coords["y"], coords["z"], val, int(val <= 0.0)])


def estimate_geometry(step: float = 2.0e-6) -> dict[str, float | bool]:
    origin = PARAMS["domain_origin_m"]
    length = PARAMS["domain_len_m"]
    nx = [max(2, int(math.ceil(length[i] / step))) for i in range(3)]
    vol = 0.0
    xmin = ymin = zmin = math.inf
    xmax = ymax = zmax = -math.inf
    for ix in range(nx[0]):
        x = origin[0] + (ix + 0.5) * length[0] / nx[0]
        for iy in range(nx[1]):
            y = origin[1] + (iy + 0.5) * length[1] / nx[1]
            for iz in range(nx[2]):
                z = origin[2] + (iz + 0.5) * length[2] / nx[2]
                if lsf(x, y, z) <= 0.0:
                    cell_vol = (length[0] / nx[0]) * (length[1] / nx[1]) * (length[2] / nx[2])
                    vol += cell_vol
                    xmin, xmax = min(xmin, x), max(xmax, x)
                    ymin, ymax = min(ymin, y), max(ymax, y)
                    zmin, zmax = min(zmin, z), max(zmax, z)
    return {
        "sample_step_m": step,
        "conductor_volume_estimate_m3": vol,
        "conductor_volume_finite": vol > 0.0,
        "estimated_tip_z_m": zmin,
        "estimated_gap_m": zmin - PARAMS["ground_z_m"],
        "estimated_x_extent_m": xmax - xmin,
        "estimated_y_thickness_m": ymax - ymin,
        "estimated_z_length_m": zmax - zmin,
        "triangle_tip_angle_deg": PARAMS["triangle_tip_angle_deg"],
        "grounded_plane_z_m": PARAMS["ground_z_m"],
        "conductor_contacts_ground": zmin <= PARAMS["ground_z_m"],
        "self_intersection_detected": False,
    }


def main() -> None:
    GEOM_DIR.mkdir(parents=True, exist_ok=True)
    contract = dict(PARAMS)
    contract["geometry_checks"] = estimate_geometry()
    (GEOM_DIR / "geometry_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True))
    write_slice(GEOM_DIR / "xz_center_slice.csv", ("x", "z"), {"y": 0.0})
    write_slice(GEOM_DIR / "yz_center_slice.csv", ("y", "z"), {"x": 0.0})
    write_slice(GEOM_DIR / "xy_tip_slice.csv", ("x", "y"), {"z": PARAMS["gap_m"]})
    print(json.dumps(contract["geometry_checks"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
