#!/usr/bin/env python3
"""Summarize Stage E1 electrostatic Afivo outputs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results_raw"
SUMMARY = ROOT / "results_summary"


def parse_log(prefix: str) -> dict[str, float]:
    path = RAW / f"{prefix}_log.txt"
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    data = rows[-1]
    # 3D standard log fixed columns; E1 user variables start after highest(lvl).
    return {
        "Emax_Vpm": float(data[8]),
        "Emax_log_x_m": float(data[9]),
        "Emax_log_y_m": float(data[10]),
        "Emax_log_z_m": float(data[11]),
        "vol_E_gt_0p5_Emax_m3": float(data[33]),
        "vol_E_gt_0p8_Emax_m3": float(data[34]),
        "conductor_volume_afivo_m3": float(data[35]),
        "Emax_user_x_m": float(data[36]),
        "Emax_user_y_m": float(data[37]),
        "Emax_user_z_m": float(data[38]),
        "min_dx_m": float(data[39]),
        "n_cells": float(data[26]),
    }


def parse_line(prefix: str) -> list[dict[str, float]]:
    path = RAW / f"{prefix}_line_000000.txt"
    out: list[dict[str, float]] = []
    with path.open() as f:
        header = f.readline().strip().lstrip("#").split()
        for line in f:
            if not line.strip():
                continue
            values = [float(x) for x in line.split()]
            out.append(dict(zip(header, values)))
    return out


def scan_field_export(prefix: str) -> dict[str, float | bool]:
    path = RAW / f"{prefix}_field_cells_000000.csv"
    max_e = -math.inf
    max_xyz = [math.nan, math.nan, math.nan]
    n_rows = 0
    finite = True
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_rows += 1
            values = [float(row[k]) for k in ("x_m", "y_m", "z_m", "phi_V", "Ex_Vpm", "Ey_Vpm", "Ez_Vpm", "Eabs_Vpm", "lsf_m")]
            if not all(math.isfinite(v) for v in values):
                finite = False
            eabs = values[7]
            if eabs > max_e:
                max_e = eabs
                max_xyz = values[0:3]
    return {
        "field_export_rows": n_rows,
        "field_export_all_finite": finite,
        "field_export_Emax_Vpm": max_e,
        "field_export_Emax_x_m": max_xyz[0],
        "field_export_Emax_y_m": max_xyz[1],
        "field_export_Emax_z_m": max_xyz[2],
    }


def profile_metric(x_prefix: str, y_prefix: str) -> dict[str, float]:
    x_line = parse_line(x_prefix)
    y_line = parse_line(y_prefix)
    n = min(len(x_line), len(y_line))
    ex = [x_line[i]["electric_fld"] for i in range(n)]
    ey = [y_line[i]["electric_fld"] for i in range(n)]
    num = math.sqrt(sum((ex[i] - ey[i]) ** 2 for i in range(n)))
    den = math.sqrt(max(sum(ex[i] ** 2 for i in range(n)), 1.0e-300))
    max_rel = max(abs(ex[i] - ey[i]) / max(abs(ex[i]), 1.0) for i in range(n))
    return {
        "profile_points": n,
        "electrostatic_asymmetry_metric": num / den,
        "profile_max_pointwise_relative_difference": max_rel,
    }


def write_profile_csv(x_prefix: str, y_prefix: str, out_path: Path) -> None:
    x_line = parse_line(x_prefix)
    y_line = parse_line(y_prefix)
    n = min(len(x_line), len(y_line))
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "x_m", "E_x_profile_Vpm", "y_m", "E_y_profile_Vpm", "relative_difference"])
        for i in range(n):
            ex = x_line[i]["electric_fld"]
            ey = y_line[i]["electric_fld"]
            rel = abs(ex - ey) / max(abs(ex), 1.0)
            writer.writerow([i, x_line[i]["x"], ex, y_line[i]["y"], ey, rel])


def main() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    coarse = parse_log("coarse_x")
    medium = parse_log("medium_x")
    coarse.update(scan_field_export("coarse_x"))
    medium.update(scan_field_export("medium_x"))
    coarse_asym = profile_metric("coarse_x", "coarse_y")
    medium_asym = profile_metric("medium_x", "medium_y")
    emax_rel_change = abs(medium["Emax_Vpm"] - coarse["Emax_Vpm"]) / max(abs(coarse["Emax_Vpm"]), 1.0)
    summary = {
        "coarse": coarse,
        "medium": medium,
        "coarse_profile": coarse_asym,
        "medium_profile": medium_asym,
        "coarse_to_medium_Emax_relative_change": emax_rel_change,
        "all_field_exports_finite": bool(coarse["field_export_all_finite"] and medium["field_export_all_finite"]),
    }
    write_profile_csv("coarse_x", "coarse_y", SUMMARY / "e1_profile_comparison_coarse.csv")
    write_profile_csv("medium_x", "medium_y", SUMMARY / "e1_profile_comparison_medium.csv")
    (SUMMARY / "e1_electrostatic_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

