"""Read-only Stage-C-style axisymmetric static operator for explicit role masks.

This reproduces the stencil and boundary-row semantics in cpp/src/poisson.cpp,
but is not a replacement for the PETSc solver or a physical geometry validator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


EPS0 = 8.8541878128e-12
GAS, HV, RETURN, OUTER = 0, 1, 2, 3


@dataclass(frozen=True)
class Grid:
    nr: int
    nz: int
    rmax_m: float
    zmin_m: float
    zmax_m: float

    @property
    def dr(self) -> float:
        return self.rmax_m / self.nr

    @property
    def dz(self) -> float:
        return (self.zmax_m - self.zmin_m) / self.nz

    def r(self, i: int) -> float:
        return (i + 0.5) * self.dr

    def z(self, j: int) -> float:
        return self.zmin_m + (j + 0.5) * self.dz

    def volume(self, i: int) -> float:
        return np.pi * (((i + 1) * self.dr) ** 2 - (i * self.dr) ** 2) * self.dz


def solve_static(grid: Grid, role: np.ndarray, outer_value: float = 0.0,
                 outer_kind: str = "DIRICHLET",
                 charge_density_C_m3: np.ndarray | None = None) -> dict:
    """Solve zero-charge 1 V weight; role is explicit physical cell classification.

    OUTER applies only to top/bottom/radial-outer boundary cells. Axis i=0 is
    intrinsic symmetry, never an OUTER cell. The matrix is not exported as a
    physical port observable; q_B uses only gas-HV links of this same stencil.
    """
    if grid.nr < 3 or grid.nz < 4 or grid.rmax_m <= 0 or grid.zmax_m <= grid.zmin_m:
        raise ValueError("INVALID_GRID")
    role = np.asarray(role, dtype=np.int8)
    if role.shape != (grid.nz, grid.nr) or not {GAS, HV, RETURN, OUTER} >= set(np.unique(role)):
        raise ValueError("INVALID_ROLE_MASK")
    if not np.any(role == HV) or not np.any(role == RETURN):
        raise ValueError("MISSING_TERMINAL")
    if outer_kind not in {"DIRICHLET", "NEUMANN_ZERO"}:
        raise ValueError("UNSUPPORTED_OUTER_BC")
    rho = np.zeros_like(role, dtype=float) if charge_density_C_m3 is None else np.asarray(charge_density_C_m3, dtype=float)
    if rho.shape != role.shape or not np.isfinite(rho).all() or np.any(rho[role != GAS] != 0):
        raise ValueError("INVALID_SPACE_CHARGE_SUPPORT")
    zero_charge = not np.any(rho)
    n = grid.nr * grid.nz
    matrix = lil_matrix((n, n), dtype=float)
    rhs = np.zeros(n)
    def idx(i: int, j: int) -> int:
        return j * grid.nr + i
    for j in range(grid.nz):
        for i in range(grid.nr):
            row = idx(i, j)
            kind = int(role[j, i])
            boundary = i == grid.nr - 1 or j == 0 or j == grid.nz - 1
            if kind in (HV, RETURN):
                matrix[row, row] = 1.0
                rhs[row] = float(kind == HV)
            elif boundary:
                if kind != OUTER:
                    raise ValueError("UNLABELED_OUTER_BOUNDARY")
                if outer_kind == "DIRICHLET":
                    matrix[row, row] = 1.0
                    rhs[row] = outer_value
                else:
                    ni, nj, h = (i - 1, j, grid.dr) if i == grid.nr - 1 else (
                        (i, j + 1, grid.dz) if j == 0 else (i, j - 1, grid.dz))
                    matrix[row, row] = 1 / h
                    matrix[row, idx(ni, nj)] = -1 / h
            elif kind == OUTER:
                raise ValueError("OUTER_ROLE_IN_INTERIOR")
            else:
                r = grid.r(i)
                ar = (i + 1) * grid.dr / (r * grid.dr ** 2)
                al = i * grid.dr / (r * grid.dr ** 2)
                az = 1 / grid.dz ** 2
                matrix[row, row] = -(ar + al + 2 * az)
                matrix[row, idx(i + 1, j)] = ar
                if i > 0:
                    matrix[row, idx(i - 1, j)] = al
                matrix[row, idx(i, j + 1)] = az
                matrix[row, idx(i, j - 1)] = az
                rhs[row] = -rho[j, i] / EPS0
    potential = np.asarray(spsolve(matrix.tocsr(), rhs)).reshape((grid.nz, grid.nr))
    if not np.isfinite(potential).all():
        raise ValueError("STATIC_SOLVE_FAILED")
    q_hv = 0.0
    energy = 0.0
    outer_energy = 0.0
    links = 0
    for j in range(grid.nz):
        for i in range(grid.nr):
            for di, dj in ((1, 0), (0, 1)):
                ni, nj = i + di, j + dj
                if ni >= grid.nr or nj >= grid.nz:
                    continue
                a, b = int(role[j, i]), int(role[nj, ni])
                if a == OUTER or b == OUTER:
                    if outer_kind == "NEUMANN_ZERO":
                        continue
                if a != GAS and b != GAS:
                    continue
                if di:
                    area = 2 * np.pi * (i + 1) * grid.dr * grid.dz
                    h = grid.dr
                else:
                    area = np.pi * (((i + 1) * grid.dr) ** 2 - (i * grid.dr) ** 2)
                    h = grid.dz
                c = EPS0 * area / h
                dv = potential[j, i] - potential[nj, ni]
                contribution = 0.5 * c * dv * dv
                energy += contribution
                if a == OUTER or b == OUTER:
                    outer_energy += contribution
                if a == HV:
                    q_hv += c * dv
                if b == HV:
                    q_hv -= c * dv
                links += 1
    return {"psi": potential, "q_conjugate_C": q_hv, "C_conjugate_F": q_hv if zero_charge else None,
            "U_operator_J": energy, "outer_link_energy_J": outer_energy,
            "conjugacy_relative_error": abs(q_hv - 2 * energy) / max(abs(q_hv), abs(2 * energy), 1e-300) if zero_charge else None,
            "link_count": links, "operator_kind": "STAGE_C_STENCIL_RECONSTRUCTION_READ_ONLY",
            "geometry_provenance": "EXPLICIT_ROLE_MASK_REQUIRED", "physical_port_status": "PENDING",
            "space_charge_status": "ZERO" if zero_charge else "SYNTHETIC_SOURCE_ONLY"}


def route_geometry(family: str) -> str:
    return "SUPPORTED_NATIVE_NEEDLE_PLANE" if family == "axisymmetric_needle_plane" else "UNSUPPORTED_GEOMETRY"
