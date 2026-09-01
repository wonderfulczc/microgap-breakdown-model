from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .derivatives import central_difference_three_point
from .integrals import current_moment, total_charge
from .schema import SourceRecord, SourceSeries


@dataclass(frozen=True)
class RemapAudit:
    charge_before: float
    charge_after: float
    relative_charge_error: float
    M_before: list[float]
    M_after: list[float]
    relative_M_error: float
    rho_min: float
    rho_max: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _relative_error(after: float, before: float) -> float:
    scale = max(abs(before), abs(after), 1e-300)
    return abs(after - before) / scale


def remap_audit(before: SourceRecord, after: SourceRecord) -> RemapAudit:
    q0 = total_charge(before)
    q1 = total_charge(after)
    m0 = current_moment(before)
    m1 = current_moment(after)
    m_scale = max(float(np.linalg.norm(m0)), float(np.linalg.norm(m1)), 1e-300)
    return RemapAudit(
        charge_before=q0,
        charge_after=q1,
        relative_charge_error=_relative_error(q1, q0),
        M_before=[float(x) for x in m0],
        M_after=[float(x) for x in m1],
        relative_M_error=float(np.linalg.norm(m1 - m0) / m_scale),
        rho_min=float(np.min(after.columns["rho"])),
        rho_max=float(np.max(after.columns["rho"])),
    )


def boundary_flux(record: SourceRecord, atol: float = 1e-15) -> float:
    lo, hi = record.bounds()
    domain_lo = np.min(lo, axis=0)
    domain_hi = np.max(hi, axis=0)
    J = (record.columns["Jx"], record.columns["Jy"], record.columns["Jz"])
    widths = (record.columns["dx"], record.columns["dy"], record.columns["dz"])
    flux = 0.0
    for axis in range(3):
        area = record.columns["cell_volume"] / widths[axis]
        low = np.isclose(lo[:, axis], domain_lo[axis], rtol=0.0, atol=atol)
        high = np.isclose(hi[:, axis], domain_hi[axis], rtol=0.0, atol=atol)
        flux += float(np.sum(-J[axis][low] * area[low], dtype=np.float64))
        flux += float(np.sum(J[axis][high] * area[high], dtype=np.float64))
    return flux


def global_continuity_audit(series: SourceSeries, boundary_flux_at_center: float = 0.0) -> dict[str, float]:
    if len(series.records) < 3:
        raise ValueError("global continuity audit requires at least three records")
    q = np.asarray([total_charge(record) for record in series.records[:3]], dtype=float)
    dQdt = float(central_difference_three_point(series.times[:3], q))
    residual = dQdt + boundary_flux_at_center
    scale = max(abs(dQdt), abs(boundary_flux_at_center), 1e-300)
    return {
        "dQdt_Cps": dQdt,
        "surface_flux_A": float(boundary_flux_at_center),
        "continuity_global_abs_A": float(abs(residual)),
        "continuity_global_rel": float(abs(residual) / scale),
    }


def local_continuity_on_matching_grid(series: SourceSeries) -> dict[str, float | str]:
    if len(series.records) < 3:
        raise ValueError("local continuity audit requires at least three records")
    r0, r1, r2 = series.records[:3]
    keys = ("x_center", "y_center", "z_center", "dx", "dy", "dz")
    for key in keys:
        if not np.allclose(r0.columns[key], r1.columns[key]) or not np.allclose(r0.columns[key], r2.columns[key]):
            raise ValueError("local continuity requires matching cell geometry")
    drhodt = central_difference_three_point(
        series.times[:3], np.vstack([r.columns["rho"] for r in (r0, r1, r2)])
    )
    div = _cell_center_divergence(r1)
    residual = drhodt + div
    return {
        "method": "cell-centered finite-difference reconstruction on matching Cartesian partition",
        "continuity_local_L1": float(np.mean(np.abs(residual))),
        "continuity_local_L2": float(np.sqrt(np.mean(residual**2))),
        "continuity_local_Linf": float(np.max(np.abs(residual))),
    }


def _cell_center_divergence(record: SourceRecord) -> np.ndarray:
    x = record.columns["x_center"]
    y = record.columns["y_center"]
    z = record.columns["z_center"]
    ux = np.unique(x)
    uy = np.unique(y)
    uz = np.unique(z)
    if ux.size * uy.size * uz.size != record.n_cells:
        raise ValueError("partition is not a full tensor-product grid")
    shape = (ux.size, uy.size, uz.size)
    order = np.lexsort((z, y, x))
    grid_index = {tuple(np.round([x[i], y[i], z[i]], 15)): n for n, i in enumerate(order)}
    arrs = []
    for name in ("Jx", "Jy", "Jz"):
        arr = np.empty(shape, dtype=float)
        for ix, xv in enumerate(ux):
            for iy, yv in enumerate(uy):
                for iz, zv in enumerate(uz):
                    arr[ix, iy, iz] = record.columns[name][order[grid_index[tuple(np.round([xv, yv, zv], 15))]]]
        arrs.append(arr)
    div = np.gradient(arrs[0], ux, axis=0, edge_order=1)
    div += np.gradient(arrs[1], uy, axis=1, edge_order=1)
    div += np.gradient(arrs[2], uz, axis=2, edge_order=1)
    flat = np.empty(record.n_cells, dtype=float)
    for ix, xv in enumerate(ux):
        for iy, yv in enumerate(uy):
            for iz, zv in enumerate(uz):
                src = order[grid_index[tuple(np.round([xv, yv, zv], 15))]]
                flat[src] = div[ix, iy, iz]
    return flat

