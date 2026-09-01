from __future__ import annotations

import numpy as np

from .schema import SourceRecord


SOURCE_FIELDS = ("rho", "Jx", "Jy", "Jz")


def _overlap_1d(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))


def conservative_remap(source: SourceRecord, target_geometry: SourceRecord) -> SourceRecord:
    """Map source cell averages onto target cells by exact box-overlap volumes."""

    src_lo, src_hi = source.bounds()
    tgt_lo, tgt_hi = target_geometry.bounds()
    target_values = {name: np.zeros(target_geometry.n_cells, dtype=float) for name in SOURCE_FIELDS}

    for ti in range(target_geometry.n_cells):
        tlo = tgt_lo[ti]
        thi = tgt_hi[ti]
        candidates = np.where(
            (src_hi[:, 0] > tlo[0])
            & (src_lo[:, 0] < thi[0])
            & (src_hi[:, 1] > tlo[1])
            & (src_lo[:, 1] < thi[1])
            & (src_hi[:, 2] > tlo[2])
            & (src_lo[:, 2] < thi[2])
        )[0]
        weighted = {name: 0.0 for name in SOURCE_FIELDS}
        covered = 0.0
        for si in candidates:
            ov = (
                _overlap_1d(tlo[0], thi[0], src_lo[si, 0], src_hi[si, 0])
                * _overlap_1d(tlo[1], thi[1], src_lo[si, 1], src_hi[si, 1])
                * _overlap_1d(tlo[2], thi[2], src_lo[si, 2], src_hi[si, 2])
            )
            if ov <= 0:
                continue
            covered += ov
            for name in SOURCE_FIELDS:
                weighted[name] += float(source.columns[name][si]) * ov
        vol = float(target_geometry.columns["cell_volume"][ti])
        if covered > 0:
            for name in SOURCE_FIELDS:
                target_values[name][ti] = weighted[name] / vol

    J = np.column_stack((target_values["Jx"], target_values["Jy"], target_values["Jz"]))
    return target_geometry.with_source_columns(target_values["rho"], J, time_s=source.metadata.time_s)


def nearest_remap(source: SourceRecord, target_geometry: SourceRecord) -> SourceRecord:
    src_xyz = np.column_stack(
        (source.columns["x_center"], source.columns["y_center"], source.columns["z_center"])
    )
    tgt_xyz = np.column_stack(
        (
            target_geometry.columns["x_center"],
            target_geometry.columns["y_center"],
            target_geometry.columns["z_center"],
        )
    )
    values = {name: np.zeros(target_geometry.n_cells, dtype=float) for name in SOURCE_FIELDS}
    for ti, xyz in enumerate(tgt_xyz):
        idx = int(np.argmin(np.sum((src_xyz - xyz) ** 2, axis=1)))
        for name in SOURCE_FIELDS:
            values[name][ti] = source.columns[name][idx]
    J = np.column_stack((values["Jx"], values["Jy"], values["Jz"]))
    return target_geometry.with_source_columns(values["rho"], J, time_s=source.metadata.time_s)

