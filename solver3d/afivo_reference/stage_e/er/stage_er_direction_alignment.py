#!/usr/bin/env python3
"""Stage E-R head-local field direction versus head trajectory analysis."""

from __future__ import annotations

import csv
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results_raw"
SUMMARY = ROOT / "results_summary"
ER_OUT = ROOT / "er"

REFERENCE_AXIS = (0.0, 0.0, -1.0)
POLARITY_SIGN = 1.0
HEAD_DEFINITION = "primary_connected_active_component_leading_edge_ne_ge_e2_head_threshold"
LEGACY_HEAD_DEFINITION = "existing_stage_e_min_z_ne_ge_e2_head_threshold"
FIELD_METHOD = "local_peak_Eabs_within_4_min_dx_of_primary_connected_leading_edge_head"
TRAJECTORY_MIN_DISPLACEMENT_DX = 1.0

STD_USER0 = 33
USER = {
    "head_z_m": STD_USER0 + 10,
    "min_dx_m": STD_USER0 + 6,
}

CASES = {
    "triangular_foil": {
        "prefix": "e2_case_a_pi_off_500V",
        "config": ROOT / "configs" / "e2_triangular_foil_case_a_pi_off_500V.cfg",
        "case_polarity": "positive_hv_foil",
        "geometry": "rounded_triangular_copper_foil_tip",
    },
    "misaligned_needle": {
        "prefix": "e3_misaligned_needle_pair_500V",
        "diagnostic_prefix": "e3_misaligned_needle_pair_500V_er_headlocal",
        "config": ROOT / "configs" / "e3_misaligned_needle_pair_500V.cfg",
        "diagnostic_config": ROOT / "configs" / "e3_misaligned_needle_pair_500V_er_headlocal.cfg",
        "case_polarity": "positive_hv_needle",
        "geometry": "misaligned_needle_pair_offset_10um",
    },
}


def clamp(x: float, lo: float = -1.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def dot(a: Iterable[float], b: Iterable[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm(v: Iterable[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def unit_vector(v: Iterable[float]) -> tuple[tuple[float, float, float] | None, str]:
    vv = tuple(float(x) for x in v)
    if len(vv) != 3 or not all(math.isfinite(x) for x in vv):
        return None, "NONFINITE_INPUT"
    n = norm(vv)
    if n <= 0.0:
        return None, "ZERO_VECTOR"
    return (vv[0] / n, vv[1] / n, vv[2] / n), "VALID"


def vector_angle_deg(a: Iterable[float], b: Iterable[float]) -> float:
    ua, sa = unit_vector(a)
    ub, sb = unit_vector(b)
    if sa != "VALID" or sb != "VALID" or ua is None or ub is None:
        return math.nan
    return math.degrees(math.acos(clamp(dot(ua, ub))))


def direction_alignment(field: Iterable[float], displacement: Iterable[float],
                        reference_axis: Iterable[float] = REFERENCE_AXIS,
                        polarity_sign: float = POLARITY_SIGN) -> dict[str, object]:
    e_raw, field_status = unit_vector(field)
    if field_status == "ZERO_VECTOR":
        field_status = "ZERO_FIELD_MAGNITUDE"
    if e_raw is None:
        return {
            "field_direction_valid": False,
            "trajectory_valid": False,
            "status": field_status,
            "theta_E_deg": math.nan,
            "theta_head_deg": math.nan,
            "delta_theta_deg": math.nan,
        }
    e_prop = tuple(polarity_sign * x for x in e_raw)
    d_head, traj_status = unit_vector(displacement)
    if traj_status == "ZERO_VECTOR":
        traj_status = "ZERO_HEAD_DISPLACEMENT"
    if d_head is None:
        return {
            "field_direction_valid": True,
            "trajectory_valid": False,
            "status": traj_status,
            "theta_E_deg": vector_angle_deg(e_prop, reference_axis),
            "theta_head_deg": math.nan,
            "delta_theta_deg": math.nan,
        }
    return {
        "field_direction_valid": True,
        "trajectory_valid": True,
        "status": "VALID",
        "theta_E_deg": vector_angle_deg(e_prop, reference_axis),
        "theta_head_deg": vector_angle_deg(d_head, reference_axis),
        "delta_theta_deg": vector_angle_deg(e_prop, d_head),
    }


def parse_config_scalar(path: Path, key: str, default: float) -> float:
    needle = f"{key} ="
    for line in path.read_text().splitlines():
        clean = line.split("#", 1)[0].strip()
        if clean.startswith(needle):
            return float(clean.split("=", 1)[1].split()[0])
    return default


def parse_log(prefix: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for line in (RAW / f"{prefix}_log.txt").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("it "):
            continue
        data = line.split()
        if len(data) < STD_USER0:
            continue
        row = {"it": float(data[0]), "time_s": float(data[1])}
        for name, idx in USER.items():
            row[name] = float(data[idx]) if idx < len(data) else math.nan
        rows.append(row)
    return rows


def source_files_by_time(prefix: str, kind: str = "source") -> list[Path]:
    by_time: dict[float, Path] = {}
    for path in sorted(RAW.glob(f"{prefix}_{kind}_*.csv")):
        with path.open() as f:
            reader = csv.DictReader(f)
            first = next(reader, None)
        if first is None:
            continue
        by_time[float(first["time_s"])] = path
    return [by_time[t] for t in sorted(by_time)]


@dataclass
class HeadFieldSample:
    time_s: float
    source_file: str
    min_dx_m: float
    head_valid: bool
    head_status: str
    field_direction_valid: bool
    field_status: str
    head_x_m: float = math.nan
    head_y_m: float = math.nan
    head_z_m: float = math.nan
    Ehead_x_V_m: float = math.nan
    Ehead_y_V_m: float = math.nan
    Ehead_z_V_m: float = math.nan
    Ehead_mag_V_m: float = math.nan
    field_sample_x_m: float = math.nan
    field_sample_y_m: float = math.nan
    field_sample_z_m: float = math.nan
    active_component_count: int = 0
    primary_component_cells: int = 0
    legacy_min_z_m: float = math.nan
    legacy_min_z_in_primary_component: bool = False
    primary_component_min_z_m: float = math.nan


def _row_float(row: dict[str, str], name: str) -> float:
    return float(row[name])


def trajectory_resolution_status(displacement_m: float, dx_m: float,
                                 min_over_dx: float = TRAJECTORY_MIN_DISPLACEMENT_DX) -> tuple[bool, str, float]:
    if not math.isfinite(displacement_m) or not math.isfinite(dx_m):
        return False, "NONFINITE_INPUT", math.nan
    if dx_m <= 0.0:
        return False, "INVALID_SPATIAL_RESOLUTION", math.nan
    if displacement_m <= 0.0:
        return False, "ZERO_HEAD_DISPLACEMENT", 0.0
    over_dx = displacement_m / dx_m
    if over_dx < min_over_dx:
        return False, "INSUFFICIENT_SPATIAL_DISPLACEMENT", over_dx
    return True, "VALID", over_dx


def _component_labels(cells: list[dict[str, float]], min_dx: float) -> tuple[list[int], list[dict[str, float]]]:
    if not cells:
        return [], []
    key_to_index: dict[tuple[int, int, int], int] = {}
    for idx, cell in enumerate(cells):
        key = (
            int(round(cell["x"] / min_dx)),
            int(round(cell["y"] / min_dx)),
            int(round(cell["z"] / min_dx)),
        )
        cell["kx"], cell["ky"], cell["kz"] = key
        key_to_index.setdefault(key, idx)

    labels = [-1] * len(cells)
    components: list[dict[str, float]] = []
    offsets = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    for start in range(len(cells)):
        if labels[start] >= 0:
            continue
        cid = len(components)
        stack = [start]
        labels[start] = cid
        count = 0
        inventory = 0.0
        min_z = math.inf
        while stack:
            idx = stack.pop()
            cell = cells[idx]
            count += 1
            inventory += cell["ne"] * cell["vol"]
            min_z = min(min_z, cell["z"])
            kx, ky, kz = int(cell["kx"]), int(cell["ky"]), int(cell["kz"])
            for ox, oy, oz in offsets:
                nxt = key_to_index.get((kx + ox, ky + oy, kz + oz))
                if nxt is not None and labels[nxt] < 0:
                    labels[nxt] = cid
                    stack.append(nxt)
        components.append({"count": float(count), "inventory": inventory, "min_z": min_z})
    return labels, components


def reconstruct_head_and_field(path: Path, threshold: float, metadata_min_dx_m: float = math.nan) -> HeadFieldSample:
    time_s = math.nan
    source_min_dx = math.inf
    legacy_head_z = math.inf
    finite = True
    active_cells: list[dict[str, float]] = []

    with path.open() as f:
        reader = csv.DictReader(f)
        required = {"time_s", "x_m", "y_m", "z_m", "cell_volume_m3", "ne_m3",
                    "Ex_Vpm", "Ey_Vpm", "Ez_Vpm", "Eabs_Vpm", "lsf_m"}
        if not required <= set(reader.fieldnames or []):
            return HeadFieldSample(math.nan, str(path.relative_to(ROOT)), math.nan, False,
                                   "MISSING_SOURCE_FIELDS", False, "MISSING_SOURCE_FIELDS")
        for row in reader:
            vals = [_row_float(row, k) for k in required]
            if not all(math.isfinite(v) for v in vals):
                finite = False
                continue
            time_s = _row_float(row, "time_s")
            cell_dx = _row_float(row, "cell_volume_m3") ** (1.0 / 3.0)
            source_min_dx = min(source_min_dx, cell_dx)
            if _row_float(row, "lsf_m") > 0.0 and _row_float(row, "ne_m3") >= threshold:
                x = _row_float(row, "x_m")
                y = _row_float(row, "y_m")
                z = _row_float(row, "z_m")
                legacy_head_z = min(legacy_head_z, z)
                active_cells.append({
                    "x": x,
                    "y": y,
                    "z": z,
                    "vol": _row_float(row, "cell_volume_m3"),
                    "ne": _row_float(row, "ne_m3"),
                })

    if not finite:
        dx_out = metadata_min_dx_m if math.isfinite(metadata_min_dx_m) else source_min_dx
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), dx_out, False,
                               "NONFINITE_INPUT", False, "NONFINITE_INPUT")
    min_dx = metadata_min_dx_m if math.isfinite(metadata_min_dx_m) and metadata_min_dx_m > 0.0 else source_min_dx
    if not math.isfinite(legacy_head_z) or not active_cells:
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, False,
                               "NO_VALID_HEAD", False, "NO_VALID_HEAD")

    component_dx = source_min_dx if math.isfinite(source_min_dx) and source_min_dx > 0.0 else min_dx
    labels, components = _component_labels(active_cells, component_dx)
    if not components:
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, False,
                               "NO_VALID_HEAD", False, "NO_VALID_HEAD")
    primary = max(range(len(components)),
                  key=lambda cid: (components[cid]["inventory"], -components[cid]["min_z"]))
    primary_min_z = components[primary]["min_z"]
    legacy_indices = [idx for idx, cell in enumerate(active_cells)
                      if abs(cell["z"] - legacy_head_z) <= 0.5 * min_dx]
    legacy_in_primary = any(labels[idx] == primary for idx in legacy_indices)

    z_band = max(1.25 * min_dx, 1.0e-12)
    sw = sx = sy = sz = 0.0
    primary_count = 0
    for idx, cell in enumerate(active_cells):
        if labels[idx] != primary:
            continue
        primary_count += 1
        if abs(cell["z"] - primary_min_z) > z_band:
            continue
        w = cell["ne"] * cell["vol"]
        sw += w
        sx += cell["x"] * w
        sy += cell["y"] * w
        sz += cell["z"] * w

    if sw <= 0.0:
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, False,
                               "NO_VALID_HEAD", False, "NO_VALID_HEAD")
    hx, hy, hz = sx / sw, sy / sw, sz / sw
    radius = max(4.0 * min_dx, 1.0e-12)
    best: dict[str, float] | None = None
    best_e = -math.inf
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if _row_float(row, "lsf_m") <= 0.0:
                continue
            dx = _row_float(row, "x_m") - hx
            dy = _row_float(row, "y_m") - hy
            dz = _row_float(row, "z_m") - hz
            if math.sqrt(dx * dx + dy * dy + dz * dz) > radius:
                continue
            eabs = _row_float(row, "Eabs_Vpm")
            if eabs > best_e:
                best_e = eabs
                best = {
                    "x": _row_float(row, "x_m"),
                    "y": _row_float(row, "y_m"),
                    "z": _row_float(row, "z_m"),
                    "Ex": _row_float(row, "Ex_Vpm"),
                    "Ey": _row_float(row, "Ey_Vpm"),
                    "Ez": _row_float(row, "Ez_Vpm"),
                    "Eabs": eabs,
                }

    if best is None:
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, True,
                               "VALID", False, "NO_VALID_LOCAL_FIELD", hx, hy, hz,
                               active_component_count=len(components),
                               primary_component_cells=primary_count,
                               legacy_min_z_m=legacy_head_z,
                               legacy_min_z_in_primary_component=legacy_in_primary,
                               primary_component_min_z_m=primary_min_z)
    if best["Eabs"] <= 0.0:
        return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, True,
                               "VALID", False, "ZERO_FIELD_MAGNITUDE", hx, hy, hz,
                               active_component_count=len(components),
                               primary_component_cells=primary_count,
                               legacy_min_z_m=legacy_head_z,
                               legacy_min_z_in_primary_component=legacy_in_primary,
                               primary_component_min_z_m=primary_min_z)
    return HeadFieldSample(time_s, str(path.relative_to(ROOT)), min_dx, True, "VALID", True,
                           "VALID", hx, hy, hz, best["Ex"], best["Ey"], best["Ez"],
                           best["Eabs"], best["x"], best["y"], best["z"],
                           len(components), primary_count, legacy_head_z,
                           legacy_in_primary, primary_min_z)


def analyze_case(case_name: str, info: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object]]:
    prefix = str(info["prefix"])
    diagnostic_prefix = str(info.get("diagnostic_prefix", prefix))
    diagnostic_files = source_files_by_time(diagnostic_prefix, "head_local")
    if diagnostic_files:
        input_prefix = diagnostic_prefix
        input_kind = "head_local"
        config_path = Path(info.get("diagnostic_config", info["config"]))
    else:
        input_prefix = prefix
        input_kind = "source"
        config_path = Path(info["config"])
        diagnostic_files = source_files_by_time(input_prefix, input_kind)
    threshold = parse_config_scalar(config_path, "e2%head_threshold", 1.0e14)
    log_dx_by_time = {float(row["time_s"]): float(row["min_dx_m"]) for row in parse_log(input_prefix)
                      if math.isfinite(float(row["min_dx_m"]))}
    samples = []
    for path in diagnostic_files:
        file_time = math.nan
        with path.open() as f:
            first = next(csv.DictReader(f), None)
            if first is not None:
                file_time = float(first["time_s"])
        metadata_dx = log_dx_by_time.get(file_time, math.nan)
        samples.append(reconstruct_head_and_field(path, threshold, metadata_dx))
    rows: list[dict[str, object]] = []
    for idx, sample in enumerate(samples):
        trajectory_valid = False
        trajectory_status = "INSUFFICIENT_NEXT_POSITION"
        head_dir = (math.nan, math.nan, math.nan)
        theta_head = math.nan
        theta_E = math.nan
        delta = math.nan
        trajectory_displacement_m = math.nan
        trajectory_displacement_over_dx = math.nan
        trajectory_dx_ref_m = math.nan
        trajectory_end_time_s = math.nan
        trajectory_lag_s = math.nan
        if sample.head_valid and sample.field_direction_valid and idx + 1 < len(samples):
            chosen_next: HeadFieldSample | None = None
            chosen_disp: tuple[float, float, float] | None = None
            chosen_status = "INSUFFICIENT_NEXT_POSITION"
            chosen_over_dx = math.nan
            chosen_dx_ref = math.nan
            for nxt in samples[idx + 1:]:
                if not nxt.head_valid:
                    chosen_status = nxt.head_status
                    continue
                disp = (nxt.head_x_m - sample.head_x_m, nxt.head_y_m - sample.head_y_m,
                        nxt.head_z_m - sample.head_z_m)
                trajectory_displacement_m = norm(disp)
                trajectory_dx_ref_m = max(sample.min_dx_m, nxt.min_dx_m)
                resolved, res_status, trajectory_displacement_over_dx = trajectory_resolution_status(
                    trajectory_displacement_m, trajectory_dx_ref_m)
                chosen_status = res_status
                chosen_over_dx = trajectory_displacement_over_dx
                chosen_dx_ref = trajectory_dx_ref_m
                if resolved:
                    chosen_next = nxt
                    chosen_disp = disp
                    break
            if chosen_next is None or chosen_disp is None:
                trajectory_valid = False
                trajectory_status = chosen_status
                trajectory_displacement_over_dx = chosen_over_dx
                trajectory_dx_ref_m = chosen_dx_ref
                theta_E = vector_angle_deg(
                    (POLARITY_SIGN * sample.Ehead_x_V_m,
                     POLARITY_SIGN * sample.Ehead_y_V_m,
                     POLARITY_SIGN * sample.Ehead_z_V_m),
                    REFERENCE_AXIS)
            else:
                trajectory_displacement_m = norm(chosen_disp)
                trajectory_dx_ref_m = chosen_dx_ref
                trajectory_displacement_over_dx = chosen_over_dx
                align = direction_alignment((sample.Ehead_x_V_m, sample.Ehead_y_V_m, sample.Ehead_z_V_m),
                                            chosen_disp)
                trajectory_valid = bool(align["trajectory_valid"])
                trajectory_status = str(align["status"])
                theta_E = float(align["theta_E_deg"])
                theta_head = float(align["theta_head_deg"])
                delta = float(align["delta_theta_deg"])
                if trajectory_valid:
                    head_dir, _ = unit_vector(chosen_disp)
                else:
                    head_dir = (math.nan, math.nan, math.nan)
                trajectory_end_time_s = chosen_next.time_s
                trajectory_lag_s = chosen_next.time_s - sample.time_s
        elif not sample.head_valid:
            trajectory_status = sample.head_status
        elif not sample.field_direction_valid:
            trajectory_status = sample.field_status

        eunit, _ = unit_vector((sample.Ehead_x_V_m, sample.Ehead_y_V_m, sample.Ehead_z_V_m))
        if eunit is None:
            eunit = (math.nan, math.nan, math.nan)
        else:
            eunit = tuple(POLARITY_SIGN * v for v in eunit)
        rows.append({
            "case": case_name,
            "time_s": sample.time_s,
            "source_file": sample.source_file,
            "head_x_m": sample.head_x_m,
            "head_y_m": sample.head_y_m,
            "head_z_m": sample.head_z_m,
            "Ehead_x_V_m": sample.Ehead_x_V_m,
            "Ehead_y_V_m": sample.Ehead_y_V_m,
            "Ehead_z_V_m": sample.Ehead_z_V_m,
            "Ehead_mag_V_m": sample.Ehead_mag_V_m,
            "Eprop_dir_x": eunit[0],
            "Eprop_dir_y": eunit[1],
            "Eprop_dir_z": eunit[2],
            "head_dir_x": head_dir[0],
            "head_dir_y": head_dir[1],
            "head_dir_z": head_dir[2],
            "theta_E_deg": theta_E,
            "theta_head_deg": theta_head,
            "delta_theta_deg": delta,
            "head_valid": sample.head_valid,
            "head_status": sample.head_status,
            "field_direction_valid": sample.field_direction_valid,
            "field_status": sample.field_status,
            "trajectory_valid": trajectory_valid,
            "trajectory_status": trajectory_status,
            "trajectory_displacement_m": trajectory_displacement_m,
            "trajectory_dx_ref_m": trajectory_dx_ref_m,
            "trajectory_displacement_over_dx": trajectory_displacement_over_dx,
            "trajectory_min_displacement_dx": TRAJECTORY_MIN_DISPLACEMENT_DX,
            "trajectory_start_time_s": sample.time_s,
            "trajectory_end_time_s": trajectory_end_time_s,
            "trajectory_lag_s": trajectory_lag_s,
            "head_definition": HEAD_DEFINITION,
            "legacy_head_definition": LEGACY_HEAD_DEFINITION,
            "field_sampling_method": FIELD_METHOD,
            "case_polarity": info["case_polarity"],
            "reference_axis_vector": "0,0,-1",
            "polarity_sign": POLARITY_SIGN,
            "head_threshold_m3": threshold,
            "field_sampling_radius_m": 4.0 * sample.min_dx_m,
            "min_dx_m": sample.min_dx_m,
            "active_component_count": sample.active_component_count,
            "primary_component_cells": sample.primary_component_cells,
            "legacy_min_z_m": sample.legacy_min_z_m,
            "legacy_min_z_in_primary_component": sample.legacy_min_z_in_primary_component,
            "primary_component_min_z_m": sample.primary_component_min_z_m,
            "field_sample_x_m": sample.field_sample_x_m,
            "field_sample_y_m": sample.field_sample_y_m,
            "field_sample_z_m": sample.field_sample_z_m,
            "input_prefix": input_prefix,
            "input_kind": input_kind,
        })

    valid_delta = [float(r["delta_theta_deg"]) for r in rows if bool(r["trajectory_valid"]) and math.isfinite(float(r["delta_theta_deg"]))]
    valid_theta_e = [float(r["theta_E_deg"]) for r in rows if math.isfinite(float(r["theta_E_deg"]))]
    valid_theta_h = [float(r["theta_head_deg"]) for r in rows if math.isfinite(float(r["theta_head_deg"]))]
    head_positions = [(float(r["head_x_m"]), float(r["head_y_m"]), float(r["head_z_m"])) for r in rows if bool(r["head_valid"])]
    displacement = math.nan
    if len(head_positions) >= 2:
        displacement = norm((head_positions[-1][0] - head_positions[0][0],
                             head_positions[-1][1] - head_positions[0][1],
                             head_positions[-1][2] - head_positions[0][2]))

    summary = {
        "case": case_name,
        "prefix": prefix,
        "input_prefix": input_prefix,
        "input_kind": input_kind,
        "geometry": info["geometry"],
        "case_polarity": info["case_polarity"],
        "reference_axis_vector": [0.0, 0.0, -1.0],
        "polarity_sign": POLARITY_SIGN,
        "head_definition": HEAD_DEFINITION,
        "legacy_head_definition": LEGACY_HEAD_DEFINITION,
        "field_sampling_method": FIELD_METHOD,
        "trajectory_min_displacement_dx": TRAJECTORY_MIN_DISPLACEMENT_DX,
        "source_snapshots": len(rows),
        "valid_time_samples": len(valid_delta),
        "time_interval_s": [rows[0]["time_s"], rows[-1]["time_s"]] if rows else [math.nan, math.nan],
        "head_displacement_m": displacement,
        "min_dx_m": min(float(r["min_dx_m"]) for r in rows if math.isfinite(float(r["min_dx_m"]))) if rows else math.nan,
        "invalid_under_resolution_samples": sum(str(r["trajectory_status"]) == "INSUFFICIENT_SPATIAL_DISPLACEMENT" for r in rows),
        "legacy_min_z_outside_primary_component_samples": sum(
            bool(r["head_valid"]) and not bool(r["legacy_min_z_in_primary_component"]) for r in rows),
        "theta_E_range_deg": [min(valid_theta_e), max(valid_theta_e)] if valid_theta_e else [math.nan, math.nan],
        "theta_head_range_deg": [min(valid_theta_h), max(valid_theta_h)] if valid_theta_h else [math.nan, math.nan],
        "delta_theta_median_deg": statistics.median(valid_delta) if valid_delta else math.nan,
        "delta_theta_p25_deg": percentile(valid_delta, 0.25) if valid_delta else math.nan,
        "delta_theta_p75_deg": percentile(valid_delta, 0.75) if valid_delta else math.nan,
        "delta_theta_p90_deg": percentile(valid_delta, 0.90) if valid_delta else math.nan,
        "delta_theta_max_deg": max(valid_delta) if valid_delta else math.nan,
        "invalid_fraction": 1.0 - len(valid_delta) / max(len(rows), 1),
        "field_valid_samples": sum(bool(r["field_direction_valid"]) for r in rows),
        "head_valid_samples": sum(bool(r["head_valid"]) for r in rows),
        "trajectory_valid_samples": sum(bool(r["trajectory_valid"]) for r in rows),
        "temporal_continuity": temporal_continuity(rows),
        "sampling_sensitivity": sampling_sensitivity(rows),
    }
    valid_lags = [float(r["trajectory_lag_s"]) for r in rows
                  if bool(r["trajectory_valid"]) and math.isfinite(float(r["trajectory_lag_s"]))]
    valid_over_dx = [float(r["trajectory_displacement_over_dx"]) for r in rows
                     if bool(r["trajectory_valid"]) and math.isfinite(float(r["trajectory_displacement_over_dx"]))]
    summary["trajectory_lag_median_s"] = statistics.median(valid_lags) if valid_lags else math.nan
    summary["trajectory_lag_p90_s"] = percentile(valid_lags, 0.90) if valid_lags else math.nan
    summary["trajectory_displacement_over_dx_median"] = statistics.median(valid_over_dx) if valid_over_dx else math.nan
    summary["total_head_displacement_over_min_dx"] = (
        displacement / summary["min_dx_m"]
        if math.isfinite(float(displacement)) and math.isfinite(float(summary["min_dx_m"])) and float(summary["min_dx_m"]) > 0.0
        else math.nan
    )
    summary["directional_evidence"] = classify_directional_evidence(summary)
    return rows, summary


def classify_directional_evidence(summary: dict[str, object]) -> str:
    valid = int(summary["valid_time_samples"])
    continuity = summary["temporal_continuity"]
    sampling = summary["sampling_sensitivity"]
    if valid < 3:
        return "NOT_RESOLVED"
    if valid < 6:
        return "RESOLVED_WEAK_ASSOCIATION"
    if isinstance(continuity, dict) and bool(continuity.get("isolated_jump_detected")):
        return "NOT_RESOLVED"
    if isinstance(sampling, dict):
        if sampling.get("status") != "VALID":
            return "NOT_RESOLVED"
        diff = sampling.get("median_difference_deg")
        if isinstance(diff, float) and math.isfinite(diff) and diff > 20.0:
            return "NOT_RESOLVED"
    return "RESOLVED_ASSOCIATION"


def percentile(values: list[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    idx = q * (len(ordered) - 1)
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (idx - lo)


def temporal_continuity(rows: list[dict[str, object]]) -> dict[str, object]:
    deltas = [float(r["delta_theta_deg"]) for r in rows if bool(r["trajectory_valid"]) and math.isfinite(float(r["delta_theta_deg"]))]
    jumps = [abs(b - a) for a, b in zip(deltas, deltas[1:])]
    head_steps = []
    valid = [r for r in rows if bool(r["head_valid"])]
    for a, b in zip(valid, valid[1:]):
        head_steps.append(norm((float(b["head_x_m"]) - float(a["head_x_m"]),
                                float(b["head_y_m"]) - float(a["head_y_m"]),
                                float(b["head_z_m"]) - float(a["head_z_m"]))))
    return {
        "max_delta_theta_jump_deg": max(jumps) if jumps else math.nan,
        "large_delta_theta_jump_count": sum(j > 45.0 for j in jumps),
        "max_head_step_m": max(head_steps) if head_steps else math.nan,
        "isolated_jump_detected": bool(any(j > 45.0 for j in jumps)),
        "diagnosis": "NO_ISOLATED_JUMP" if not any(j > 45.0 for j in jumps) else "CHECK_HEAD_LOCATOR_OR_AMR_OUTPUT",
    }


def sampling_sensitivity(rows: list[dict[str, object]]) -> dict[str, object]:
    full = [float(r["delta_theta_deg"]) for r in rows if bool(r["trajectory_valid"]) and math.isfinite(float(r["delta_theta_deg"]))]
    sparse = [float(r["delta_theta_deg"]) for n, r in enumerate(rows) if n % 2 == 0 and bool(r["trajectory_valid"]) and math.isfinite(float(r["delta_theta_deg"]))]
    if len(full) < 1 or len(sparse) < 1:
        return {"status": "INSUFFICIENT_VALID_ANGLES", "median_difference_deg": math.nan}
    return {
        "status": "VALID",
        "full_median_delta_theta_deg": statistics.median(full),
        "every_second_sample_median_delta_theta_deg": statistics.median(sparse),
        "median_difference_deg": abs(statistics.median(full) - statistics.median(sparse)),
    }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "case", "time_s", "source_file", "head_x_m", "head_y_m", "head_z_m",
        "Ehead_x_V_m", "Ehead_y_V_m", "Ehead_z_V_m", "Ehead_mag_V_m",
        "Eprop_dir_x", "Eprop_dir_y", "Eprop_dir_z",
        "head_dir_x", "head_dir_y", "head_dir_z",
        "theta_E_deg", "theta_head_deg", "delta_theta_deg",
        "head_valid", "head_status", "field_direction_valid", "field_status",
        "trajectory_valid", "trajectory_status", "trajectory_displacement_m",
        "trajectory_dx_ref_m", "trajectory_displacement_over_dx",
        "trajectory_min_displacement_dx", "trajectory_start_time_s",
        "trajectory_end_time_s", "trajectory_lag_s", "head_definition", "legacy_head_definition",
        "field_sampling_method", "case_polarity", "reference_axis_vector",
        "polarity_sign", "head_threshold_m3", "field_sampling_radius_m",
        "min_dx_m", "active_component_count", "primary_component_cells",
        "legacy_min_z_m", "legacy_min_z_in_primary_component", "primary_component_min_z_m",
        "field_sample_x_m", "field_sample_y_m", "field_sample_z_m",
        "input_prefix", "input_kind",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, tuple):
        return [json_safe(v) for v in value]
    return value


def main() -> None:
    ER_OUT.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {
        "case_id": "Stage E-R head-local peak-field direction versus streamer-head trajectory",
        "reference_axis_vector": [0.0, 0.0, -1.0],
        "polarity_sign": POLARITY_SIGN,
        "polarity_convention": "positive HV electrode field points along propagation toward ground; e_E_prop=e_E_raw",
        "coordinate_polarity_verification": {
            "triangular_foil": {
                "hv_electrode_coordinate": "high z triangular foil tip near z=70e-6 m",
                "ground_electrode_coordinate": "low z grounded plane at z=0 m",
                "reference_axis_vector": [0.0, 0.0, -1.0],
                "expected_macroscopic_E_direction": "from positive HV foil toward grounded lower electrode, approximately -z",
                "exported_local_E_direction": "head-local exported Ez is negative in existing snapshots",
                "polarity_sign": POLARITY_SIGN,
            },
            "misaligned_needle": {
                "hv_electrode_coordinate": "upper positive needle, tip surface z=70e-6 m",
                "ground_electrode_coordinate": "lower grounded needle, tip surface z=0 m, x-offset=10e-6 m",
                "reference_axis_vector": [0.0, 0.0, -1.0],
                "expected_macroscopic_E_direction": "from positive upper needle toward grounded lower needle, mostly -z with lateral component",
                "exported_local_E_direction": "head-local exported Ez is negative in existing snapshots",
                "polarity_sign": POLARITY_SIGN,
            },
        },
        "head_definition": HEAD_DEFINITION,
        "legacy_head_definition": LEGACY_HEAD_DEFINITION,
        "field_sampling_method": FIELD_METHOD,
        "trajectory_resolution_rule": "trajectory valid only when the earliest future head displacement reaches >= 1.0 * max(current_min_dx,future_min_dx)",
        "cases": {},
        "afivo_rerun_required": False,
        "afivo_rerun_performed": False,
    }
    combined: list[dict[str, object]] = []
    for case, info in CASES.items():
        rows, case_summary = analyze_case(case, info)
        write_rows(ER_OUT / f"stage_er_direction_alignment_{case}.csv", rows)
        combined.extend(rows)
        summary["cases"][case] = case_summary
        if case_summary.get("input_kind") == "head_local":
            summary["afivo_rerun_performed"] = True
    write_rows(ER_OUT / "stage_er_direction_alignment.csv", combined)
    safe_summary = json_safe(summary)
    (ER_OUT / "stage_er_direction_summary.json").write_text(json.dumps(safe_summary, indent=2, sort_keys=True))
    print(json.dumps(safe_summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
