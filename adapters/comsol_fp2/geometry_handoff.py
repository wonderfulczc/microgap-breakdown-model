"""Explicit, read-only 2D-axisymmetric static source-to-target mapping checks."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "comsol_fp1"))
from comsol_import import FP1Error  # noqa: E402
from field_mapping import map_field  # noqa: E402
from target_poisson import GAS, HV, RETURN, Grid  # noqa: E402


def audit_handoff(grid: Grid, role: np.ndarray, source_points: np.ndarray,
                  source_psi: np.ndarray, gas_mesh: dict, gap_witness: dict) -> dict:
    """Map only gas-cell samples; roles come from geometry, never Emax."""
    if role.shape != (grid.nz, grid.nr):
        raise FP1Error("INVALID_ROLE_MASK", "target-grid dimensions")
    if gap_witness.get("gas_connected") is not True:
        raise FP1Error("GAP_WITNESS_PENDING", "connected gas path required")
    ends = np.asarray([gap_witness["hv_point_m"], gap_witness["return_point_m"]], dtype=float)
    if ends.shape != (2, 2) or not np.isfinite(ends).all():
        raise FP1Error("INVALID_GAP_WITNESS", "finite r,z witness")
    targets = np.asarray([(grid.r(i), grid.z(j)) for j in range(grid.nz) for i in range(grid.nr) if role[j, i] == GAS])
    mapped, metadata = map_field(source_points, source_psi, targets, gas_mesh)
    interfaces = []
    for j in range(grid.nz):
        for i in range(grid.nr):
            if role[j, i] != GAS:
                continue
            if i + 1 < grid.nr and role[j, i+1] in (HV, RETURN):
                interfaces.append((grid.r(i) + grid.dr/2, grid.z(j)))
            if j + 1 < grid.nz and role[j+1, i] in (HV, RETURN):
                interfaces.append((grid.r(i), grid.z(j) + grid.dz/2))
            if j and role[j-1, i] in (HV, RETURN):
                interfaces.append((grid.r(i), grid.z(j) - grid.dz/2))
    if not interfaces:
        raise FP1Error("NO_TARGET_ELECTRODE_FACE", "no gas-terminal interface")
    faces = np.asarray(interfaces)
    distances = np.linalg.norm(faces[:, None, :] - ends[None, :, :], axis=2)
    return {"source_type": "DERIVED", "target_gas_cells": int(len(targets)),
            "mapped_psi_min_max": [float(np.min(mapped)), float(np.max(mapped))],
            "mapping_method": metadata["method"], "mapping_error_status": metadata["mapping_error_status"],
            "gas_containment_checked": metadata["gas_mesh_containment_checked"],
            "gap_witness_distance_m": float(np.linalg.norm(ends[0]-ends[1])),
            "witness_to_nearest_target_face_m": [float(v) for v in distances.min(axis=0)],
            "physical_face_alignment_status": "REQUIRES_GEOMETRY_SPECIFIC_TOLERANCE_AND_MESH_STUDY"}
