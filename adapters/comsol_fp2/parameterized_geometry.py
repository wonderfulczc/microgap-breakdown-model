"""Read-only FP-2B geometry contract; not an injected Stage-C backend."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np

from target_poisson import Grid


class Role(IntEnum):
    GAS = 0
    HV = 1
    RETURN = 2
    OUTER_BOUNDARY = 3
    AXIS = 4
    DIELECTRIC = 5


@dataclass(frozen=True)
class NeedlePlane:
    geometry_id: str
    ground_z_m: float
    ground_thickness_m: float
    tip_center_z_m: float
    tip_radius_m: float
    shank_radius_m: float
    dimension: str = "2D_axi"

    def __post_init__(self) -> None:
        _common(self.geometry_id, self.dimension)
        if min(self.ground_thickness_m, self.tip_radius_m, self.shank_radius_m) <= 0 or self.gap_m <= 0:
            raise ValueError("INVALID_NEEDLE_PLANE_GEOMETRY")

    @property
    def gap_m(self) -> float:
        return self.tip_center_z_m - self.tip_radius_m - (self.ground_z_m + self.ground_thickness_m)

    def classify(self, r_m: float, z_m: float) -> Role:
        _point(r_m, z_m)
        if z_m <= self.ground_z_m + self.ground_thickness_m:
            return Role.RETURN
        if np.hypot(r_m, z_m - self.tip_center_z_m) <= self.tip_radius_m:
            return Role.HV
        if r_m <= self.shank_radius_m and z_m >= self.tip_center_z_m:
            return Role.HV
        return Role.GAS


@dataclass(frozen=True)
class FacingNeedles:
    geometry_id: str
    return_tip_z_m: float
    hv_tip_z_m: float
    return_tip_radius_m: float
    hv_tip_radius_m: float
    return_shank_radius_m: float
    hv_shank_radius_m: float
    return_shank_length_m: float
    hv_shank_length_m: float
    dimension: str = "2D_axi"

    def __post_init__(self) -> None:
        _common(self.geometry_id, self.dimension)
        lengths = (self.return_tip_radius_m, self.hv_tip_radius_m,
                   self.return_shank_radius_m, self.hv_shank_radius_m,
                   self.return_shank_length_m, self.hv_shank_length_m)
        if min(lengths) <= 0 or self.gap_m <= 0:
            raise ValueError("INVALID_FACING_NEEDLES_GEOMETRY")

    @property
    def gap_m(self) -> float:
        return self.hv_tip_z_m - self.return_tip_z_m

    @property
    def return_center_z_m(self) -> float:
        return self.return_tip_z_m - self.return_tip_radius_m

    @property
    def hv_center_z_m(self) -> float:
        return self.hv_tip_z_m + self.hv_tip_radius_m

    def classify(self, r_m: float, z_m: float) -> Role:
        _point(r_m, z_m)
        lower_sphere = np.hypot(r_m, z_m - self.return_center_z_m) <= self.return_tip_radius_m
        lower_shank = (r_m <= self.return_shank_radius_m and
                       self.return_center_z_m-self.return_shank_length_m <= z_m <= self.return_center_z_m)
        upper_sphere = np.hypot(r_m, z_m - self.hv_center_z_m) <= self.hv_tip_radius_m
        upper_shank = (r_m <= self.hv_shank_radius_m and
                       self.hv_center_z_m <= z_m <= self.hv_center_z_m+self.hv_shank_length_m)
        if lower_sphere or lower_shank:
            return Role.RETURN
        if upper_sphere or upper_shank:
            return Role.HV
        return Role.GAS


def _common(geometry_id: str, dimension: str) -> None:
    if not geometry_id:
        raise ValueError("MISSING_GEOMETRY_ID")
    if dimension != "2D_axi":
        raise ValueError("UNSUPPORTED_NONAXISYMMETRIC_GEOMETRY")


def _point(r_m: float, z_m: float) -> None:
    if not np.isfinite(r_m) or not np.isfinite(z_m) or r_m < 0:
        raise ValueError("INVALID_RZ_COORDINATES")


def rasterize(geometry: NeedlePlane | FacingNeedles, grid: Grid) -> np.ndarray:
    """Cell-center roles; exterior boundary cells keep native metal-first priority."""
    if grid.rmax_m <= 0 or grid.zmax_m <= grid.zmin_m or grid.nr < 3 or grid.nz < 4:
        raise ValueError("INVALID_GRID")
    role = np.empty((grid.nz, grid.nr), dtype=np.int8)
    for j in range(grid.nz):
        for i in range(grid.nr):
            physical = geometry.classify(grid.r(i), grid.z(j))
            if physical == Role.GAS and (i == grid.nr-1 or j == 0 or j == grid.nz-1):
                physical = Role.OUTER_BOUNDARY
            role[j, i] = physical
    if not np.any(role == Role.GAS) or not np.any(role == Role.HV) or not np.any(role == Role.RETURN):
        raise ValueError("UNRESOLVED_DISCRETE_GEOMETRY")
    return role


def transport_gas_mask(role: np.ndarray) -> np.ndarray:
    """An outer-boundary row is still gas for Stage-C transport semantics."""
    if np.any((role == Role.DIELECTRIC) | (role == Role.AXIS)):
        raise ValueError("UNSUPPORTED_CELL_MATERIAL")
    return (role == Role.GAS) | (role == Role.OUTER_BOUNDARY)


def gas_boundary_faces(role: np.ndarray, grid: Grid) -> list[dict]:
    """Enumerate all gas/solid, gas/outer and axis faces, without Emax inference."""
    if role.shape != (grid.nz, grid.nr):
        raise ValueError("INVALID_ROLE_MASK")
    gas = transport_gas_mask(role)
    faces = []
    for j in range(grid.nz):
        for i in range(grid.nr):
            if not gas[j, i]:
                continue
            if i == 0:
                faces.append({"gas_cell": (i, j), "terminal": Role.AXIS.name,
                              "normal_gas_outward": (-1, 0),
                              "center_rz_m": (0.0, grid.z(j)), "area_m2": 0.0,
                              "gas_metal_center_distance_m": None})
            exterior = []
            if i == grid.nr-1:
                exterior.append(((1, 0), (grid.rmax_m, grid.z(j)), 2*np.pi*grid.rmax_m*grid.dz))
            area_z = np.pi*((i+1)**2-i**2)*grid.dr**2
            if j == 0:
                exterior.append(((0, -1), (grid.r(i), grid.zmin_m), area_z))
            if j == grid.nz-1:
                exterior.append(((0, 1), (grid.r(i), grid.zmax_m), area_z))
            for normal, center, area in exterior:
                faces.append({"gas_cell": (i, j), "terminal": Role.OUTER_BOUNDARY.name,
                              "normal_gas_outward": normal, "center_rz_m": center,
                              "area_m2": float(area), "gas_metal_center_distance_m": None})
            for di, dj, normal in ((1, 0, (1, 0)), (-1, 0, (-1, 0)),
                                    (0, 1, (0, 1)), (0, -1, (0, -1))):
                ni, nj = i+di, j+dj
                if not (0 <= ni < grid.nr and 0 <= nj < grid.nz):
                    continue
                neighbor = Role(int(role[nj, ni]))
                if gas[nj, ni]:
                    continue
                if di:
                    face_r = (i + (di > 0)) * grid.dr
                    center = (face_r, grid.z(j))
                    area = 2*np.pi*face_r*grid.dz
                else:
                    center = (grid.r(i), grid.z(j)+dj*grid.dz/2)
                    area = np.pi*((i+1)**2-i**2)*grid.dr**2
                faces.append({"gas_cell": (i, j), "terminal": neighbor.name,
                              "normal_gas_outward": normal, "center_rz_m": center,
                              "area_m2": float(area), "gas_metal_center_distance_m": grid.dr if di else grid.dz})
    return faces


def gas_terminal_faces(role: np.ndarray, grid: Grid) -> list[dict]:
    return [face for face in gas_boundary_faces(role, grid) if face["terminal"] in (Role.HV.name, Role.RETURN.name)]


def axis_gap_from_mask(role: np.ndarray, grid: Grid) -> dict:
    axis = role[:, 0]
    return_j = np.flatnonzero(axis == Role.RETURN)
    hv_j = np.flatnonzero(axis == Role.HV)
    if len(return_j) == 0 or len(hv_j) == 0:
        raise ValueError("MISSING_AXIS_TERMINAL")
    lower = int(return_j.max())
    upper = int(hv_j.min())
    if upper <= lower+1 or not np.all(axis[lower+1:upper] == Role.GAS):
        raise ValueError("AMBIGUOUS_AXIS_GAP")
    return {"gas_cell_count": upper-lower-1,
            "discrete_face_gap_m": float((upper-lower-1)*grid.dz),
            "return_face_z_m": float(grid.z(lower)+grid.dz/2),
            "hv_face_z_m": float(grid.z(upper)-grid.dz/2)}
