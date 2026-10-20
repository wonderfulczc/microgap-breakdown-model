"""Read-only linear field mapping; never extrapolate or infer conductor boundaries."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

from comsol_import import FP1Error


def map_field(points_rz_m: np.ndarray, values: np.ndarray, targets_rz_m: np.ndarray,
              gas_mesh: dict | None = None) -> tuple[np.ndarray, dict]:
    points = np.asarray(points_rz_m, dtype=float)
    targets = np.asarray(targets_rz_m, dtype=float)
    values = np.asarray(values, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or targets.ndim != 2 or targets.shape[1] != 2:
        raise FP1Error("INVALID_COORDINATES", "r,z arrays required")
    if len(points) != len(values) or not np.isfinite(points).all() or not np.isfinite(targets).all() or not np.isfinite(values).all():
        raise FP1Error("INVALID_MAPPING_INPUT", "length or finite-value mismatch")
    if (targets[:, 0] < 0).any():
        raise FP1Error("INVALID_COORDINATES", "r cannot be negative")
    if gas_mesh is not None:
        from matplotlib.tri import Triangulation

        nodes = np.asarray(gas_mesh["nodes"], dtype=float)
        triangles = np.asarray(gas_mesh["triangles"], dtype=int)
        finder = Triangulation(nodes[:, 0], nodes[:, 1], triangles).get_trifinder()
        if np.any(finder(points[:, 0], points[:, 1]) < 0):
            raise FP1Error("SOURCE_NOT_IN_GAS", "field samples outside actual gas triangles")
        if np.any(finder(targets[:, 0], targets[:, 1]) < 0):
            raise FP1Error("TARGET_NOT_IN_GAS", "target lies outside actual gas triangles")
    else:
        # A convex hull cannot detect an electrode hole; only use this for
        # stand-alone analytic fixtures, never a production Stage-C handoff.
        if len(targets) and len(points) >= 3 and np.any(Delaunay(points).find_simplex(targets) < 0):
            raise FP1Error("MESH_EXTRAPOLATION_REJECTED", "target outside source convex hull")
    triangulation = Delaunay(points)
    if gas_mesh is not None and len(targets):
        simplices = triangulation.find_simplex(targets)
        if np.any(simplices < 0):
            raise FP1Error("MESH_EXTRAPOLATION_REJECTED", "target outside source sample hull")
        corners = points[triangulation.simplices[simplices]]
        probes = np.concatenate((np.mean(corners, axis=1),
                                 0.5 * (corners[:, 0] + corners[:, 1]),
                                 0.5 * (corners[:, 1] + corners[:, 2]),
                                 0.5 * (corners[:, 2] + corners[:, 0])), axis=0)
        if np.any(finder(probes[:, 0], probes[:, 1]) < 0):
            raise FP1Error("INTERPOLATION_CROSSES_NON_GAS", "source simplex spans excluded region")
    mapped = np.asarray(LinearNDInterpolator(triangulation, values)(targets))
    if not np.isfinite(mapped).all():
        raise FP1Error("MESH_EXTRAPOLATION_REJECTED", "source samples do not cover targets")
    return mapped, {
        "method": "LINEAR_NO_EXTRAPOLATION",
        "source_sample_count": len(points),
        "target_count": len(targets),
        "gas_mesh_containment_checked": gas_mesh is not None,
        "mapping_error_status": "NOT_RESOLVED_WITHOUT_INDEPENDENT_REFERENCE",
    }


def mesh_gas_triangles(mesh: dict, gas_regions: set[int]) -> dict:
    nodes_block = mesh["nodes"]
    cells_block = mesh["cells"]
    nodes = np.asarray(nodes_block["rows"], dtype=float)
    nc = {name: i for i, name in enumerate(nodes_block["columns"])}
    cc = {name: i for i, name in enumerate(cells_block["columns"])}
    if not {"node_id", "r_m", "z_m"} <= nc.keys() or not {"domain_id", "v0", "v1", "v2"} <= cc.keys():
        raise FP1Error("INVALID_MESH", "node/cell columns missing")
    ids = nodes[:, nc["node_id"]].astype(int)
    if len(set(ids)) != len(ids) or not np.array_equal(ids, np.arange(len(ids))):
        raise FP1Error("INVALID_MESH", "node IDs must be contiguous zero-based")
    cells = np.asarray(cells_block["rows"])
    selected = np.isin(cells[:, cc["domain_id"]].astype(int), list(gas_regions))
    tri = cells[selected][:, [cc["v0"], cc["v1"], cc["v2"]]].astype(int)
    if tri.size == 0 or tri.min() < 0 or tri.max() >= len(nodes):
        raise FP1Error("INVALID_MESH", "gas connectivity invalid")
    return {"nodes": nodes[:, [nc["r_m"], nc["z_m"]]], "triangles": tri}


def gradient_audit(points_rz_m: np.ndarray, psi: np.ndarray, ew: np.ndarray,
                   gas_mesh: dict, tolerance: float) -> dict:
    """Compare piecewise-linear -grad(psi) with separately exported Ew."""
    from matplotlib.tri import Triangulation

    points = np.asarray(points_rz_m, dtype=float)
    psi = np.asarray(psi, dtype=float)
    ew = np.asarray(ew, dtype=float)
    if len(points) < 3 or len(psi) != len(points) or ew.shape != (len(points), 2):
        raise FP1Error("INVALID_GRADIENT_INPUT", "sample arrays")
    simplices = Delaunay(points).simplices
    if len(simplices) > 4096:
        simplices = simplices[np.linspace(0, len(simplices) - 1, 4096, dtype=int)]
    nodes = np.asarray(gas_mesh["nodes"], dtype=float)
    triangles = np.asarray(gas_mesh["triangles"], dtype=int)
    finder = Triangulation(nodes[:, 0], nodes[:, 1], triangles).get_trifinder()
    centroids = np.mean(points[simplices], axis=1)
    simplices = simplices[finder(centroids[:, 0], centroids[:, 1]) >= 0]
    if len(simplices) == 0:
        raise FP1Error("GRADIENT_UNRESOLVED", "no interior triangles")
    p = points[simplices]
    mat = np.stack((p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    delta = np.stack((psi[simplices[:, 1]] - psi[simplices[:, 0]],
                      psi[simplices[:, 2]] - psi[simplices[:, 0]]), axis=1)
    det = np.linalg.det(mat)
    valid = np.abs(det) > 1e-25
    grad = -np.linalg.solve(mat[valid], delta[valid, :, None])[..., 0]
    reference = np.mean(ew[simplices[valid]], axis=1)
    scale = max(float(np.percentile(np.linalg.norm(reference, axis=1), 95)), 1e-30)
    error = np.linalg.norm(grad - reference, axis=1) / np.maximum(np.linalg.norm(reference, axis=1), 0.01 * scale)
    p95 = float(np.percentile(error, 95))
    return {"sampled_triangles": int(len(error)), "relative_error_p95": p95,
            "tolerance": tolerance, "status": "PASS" if p95 <= tolerance else "FAIL"}
