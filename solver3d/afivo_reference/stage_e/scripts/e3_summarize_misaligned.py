#!/usr/bin/env python3
"""Summarize Stage E3 aligned/misaligned needle-pair dynamic benchmarks."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results_raw"
SUMMARY = ROOT / "results_summary"

STD_USER0 = 33
USER = {
    "vol_E_gt_0p5_Emax_m3": STD_USER0 + 0,
    "vol_E_gt_0p8_Emax_m3": STD_USER0 + 1,
    "conductor_volume_afivo_m3": STD_USER0 + 2,
    "Emax_user_x_m": STD_USER0 + 3,
    "Emax_user_y_m": STD_USER0 + 4,
    "Emax_user_z_m": STD_USER0 + 5,
    "min_dx_m": STD_USER0 + 6,
    "x_cm_m": STD_USER0 + 7,
    "y_cm_m": STD_USER0 + 8,
    "r_cm_m": STD_USER0 + 9,
    "head_z_m": STD_USER0 + 10,
    "total_charge_C": STD_USER0 + 11,
    "min_ne_m3": STD_USER0 + 12,
    "ne_moment_asym": STD_USER0 + 13,
    "rho_moment_asym": STD_USER0 + 14,
    "total_electrons_user": STD_USER0 + 15,
    "z_cm_m": STD_USER0 + 16,
}

SOURCE_FIELDS = [
    "time_s",
    "x_m",
    "y_m",
    "z_m",
    "cell_volume_m3",
    "rho_Cpm3",
    "ne_m3",
    "Ex_Vpm",
    "Ey_Vpm",
    "Ez_Vpm",
    "Jx_Apm2",
    "Jy_Apm2",
    "Jz_Apm2",
    "Eabs_Vpm",
    "lsf_m",
    "level",
]

CASES = {
    "aligned": "e3_aligned_needle_pair_500V",
    "misaligned": "e3_misaligned_needle_pair_500V",
}


def parse_log(prefix: str) -> list[dict[str, float]]:
    path = RAW / f"{prefix}_log.txt"
    rows: list[dict[str, float]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("it "):
            continue
        data = line.split()
        if len(data) < STD_USER0:
            continue
        row = {
            "it": float(data[0]),
            "time_s": float(data[1]),
            "dt_s": float(data[2]),
            "velocity_mps": float(data[3]),
            "total_electrons": float(data[4]),
            "total_positive_ions": float(data[5]),
            "sum_charge_number": float(data[6]),
            "Emax_Vpm": float(data[8]),
            "Emax_x_m": float(data[9]),
            "Emax_y_m": float(data[10]),
            "Emax_z_m": float(data[11]),
            "ne_max_m3": float(data[12]),
            "ne_max_x_m": float(data[13]),
            "ne_max_y_m": float(data[14]),
            "ne_max_z_m": float(data[15]),
            "voltage_V": float(data[16]),
            "wc_time_s": float(data[25]),
            "n_cells": float(data[26]),
            "min_dx_m_standard": float(data[27]),
            "highest_level": float(data[32]),
        }
        for name, idx in USER.items():
            row[name] = float(data[idx]) if idx < len(data) else math.nan
        rows.append(row)
    if not rows:
        raise RuntimeError(f"no parsed log rows for {prefix}")
    return rows


def latest_source(prefix: str) -> Path:
    matches = sorted(RAW.glob(f"{prefix}_source_*.csv"))
    if not matches:
        raise FileNotFoundError(f"no source snapshots for {prefix}")
    return max(matches, key=lambda path: path.stat().st_mtime)


def normalized_l2(a: Iterable[float], b: Iterable[float]) -> float:
    av = list(a)
    bv = list(b)
    n = min(len(av), len(bv))
    if n == 0:
        return math.nan
    num = math.sqrt(sum((av[i] - bv[i]) ** 2 for i in range(n)))
    den = math.sqrt(max(sum(av[i] ** 2 for i in range(n)), 1.0e-300))
    return num / den


def rel_diff(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def source_profile_metrics(case: str, prefix: str, min_dx: float) -> dict[str, object]:
    path = latest_source(prefix)
    bin_dx = max(min_dx, 1.0e-12)
    profile_path = SUMMARY / f"e3_profiles_{case}.csv"
    source_rows: list[dict[str, float]] = []
    min_abs_x = math.inf
    min_abs_y = math.inf

    rows = 0
    finite = True
    schema_ok = False
    with path.open() as f:
        reader = csv.DictReader(f)
        schema_ok = reader.fieldnames == SOURCE_FIELDS
        for row in reader:
            rows += 1
            vals = {k: float(row[k]) for k in SOURCE_FIELDS if k != "level"}
            vals["level"] = float(row["level"])
            if not all(math.isfinite(v) for v in vals.values()):
                finite = False
            if vals["lsf_m"] <= 0.0:
                continue
            min_abs_x = min(min_abs_x, abs(vals["x_m"]))
            min_abs_y = min(min_abs_y, abs(vals["y_m"]))
            source_rows.append(vals)

    if not source_rows:
        raise RuntimeError(f"{prefix}: no gas source rows")
    x_slice_tol = min_abs_x + max(0.51 * min_dx, 1.0e-12)
    y_slice_tol = min_abs_y + max(0.51 * min_dx, 1.0e-12)

    plane_data: dict[str, dict[tuple[int, int], dict[str, list[float]]]] = {
        "xz": defaultdict(lambda: defaultdict(list)),
        "yz": defaultdict(lambda: defaultdict(list)),
    }
    for vals in source_rows:
        x = vals["x_m"]
        y = vals["y_m"]
        z = vals["z_m"]
        if abs(y) <= y_slice_tol:
            key = (round(x / bin_dx), round(z / bin_dx))
            bucket = plane_data["xz"][key]
            bucket["s_m"].append(x)
            bucket["z_m"].append(z)
            bucket["Eabs_Vpm"].append(vals["Eabs_Vpm"])
            bucket["ne_m3"].append(vals["ne_m3"])
            bucket["rho_Cpm3"].append(vals["rho_Cpm3"])
        if abs(x) <= x_slice_tol:
            key = (round(y / bin_dx), round(z / bin_dx))
            bucket = plane_data["yz"][key]
            bucket["s_m"].append(y)
            bucket["z_m"].append(z)
            bucket["Eabs_Vpm"].append(vals["Eabs_Vpm"])
            bucket["ne_m3"].append(vals["ne_m3"])
            bucket["rho_Cpm3"].append(vals["rho_Cpm3"])

    def avg(values: list[float]) -> float:
        return sum(values) / len(values)

    compact: dict[str, dict[tuple[int, int], dict[str, float]]] = {"xz": {}, "yz": {}}
    for plane in ("xz", "yz"):
        for key, variables in plane_data[plane].items():
            compact[plane][key] = {name: avg(vals) for name, vals in variables.items()}

    with profile_path.open("w", newline="") as f:
        fields = ["case", "plane", "s_m", "z_m", "Eabs_Vpm", "ne_m3", "rho_Cpm3"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for plane in ("xz", "yz"):
            for key in sorted(compact[plane]):
                row = {"case": case, "plane": plane}
                row.update(compact[plane][key])
                writer.writerow(row)

    common = sorted(set(compact["xz"]) & set(compact["yz"]))
    field_asym = normalized_l2(
        (compact["xz"][k]["Eabs_Vpm"] for k in common),
        (compact["yz"][k]["Eabs_Vpm"] for k in common),
    )
    ne_asym = normalized_l2(
        (compact["xz"][k]["ne_m3"] for k in common),
        (compact["yz"][k]["ne_m3"] for k in common),
    )
    rho_asym = normalized_l2(
        (compact["xz"][k]["rho_Cpm3"] for k in common),
        (compact["yz"][k]["rho_Cpm3"] for k in common),
    )

    return {
        "source_file": str(path.relative_to(ROOT)),
        "profile_file": str(profile_path.relative_to(ROOT)),
        "source_rows": rows,
        "source_schema_ok": schema_ok,
        "source_all_finite": finite,
        "center_slice_min_abs_x_m": min_abs_x,
        "center_slice_min_abs_y_m": min_abs_y,
        "common_profile_bins": len(common),
        "field_plane_asymmetry": field_asym,
        "ne_plane_asymmetry": ne_asym,
        "rho_plane_asymmetry": rho_asym,
    }


def summarize_case(case: str, prefix: str) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows = parse_log(prefix)
    first = rows[0]
    final = rows[-1]
    min_dx = final["min_dx_m"]
    profile = source_profile_metrics(case, prefix, min_dx)
    head_rows = [r for r in rows if r["head_z_m"] > -1.0e90]
    head_delta = 0.0
    head_velocity = 0.0
    if len(head_rows) >= 2:
        head_delta = head_rows[-1]["head_z_m"] - head_rows[0]["head_z_m"]
        dt = head_rows[-1]["time_s"] - head_rows[0]["time_s"]
        if dt > 0.0:
            head_velocity = head_delta / dt

    trajectory = []
    for row in rows:
        trajectory.append(
            {
                "case": case,
                "time_s": row["time_s"],
                "x_cm_m": row["x_cm_m"],
                "y_cm_m": row["y_cm_m"],
                "z_cm_m": row["z_cm_m"],
                "r_cm_m": row["r_cm_m"],
                "Emax_Vpm": row["Emax_Vpm"],
                "ne_max_m3": row["ne_max_m3"],
                "total_electrons": row["total_electrons"],
                "total_charge_C": row["total_charge_C"],
                "head_z_m": row["head_z_m"],
                "ne_moment_asym": row["ne_moment_asym"],
                "rho_moment_asym": row["rho_moment_asym"],
            }
        )

    all_log_finite = True
    for row in rows:
        for name, value in row.items():
            if name == "head_z_m" and value <= -1.0e90:
                continue
            if not math.isfinite(value):
                all_log_finite = False

    ne_growth = final["ne_max_m3"] / max(first["ne_max_m3"], 1.0)
    total_growth = final["total_electrons"] / max(first["total_electrons"], 1.0)
    streamer_propagation = abs(head_delta) > max(2.0 * min_dx, 2.0e-6)
    bridge = final["head_z_m"] > -1.0e90 and final["head_z_m"] <= 2.0e-6
    case_summary: dict[str, object] = {
        "prefix": prefix,
        "n_outputs": len(rows),
        "initial_time_s": first["time_s"],
        "final_time_s": final["time_s"],
        "min_dx_m": min_dx,
        "highest_level": final["highest_level"],
        "n_cells_final": final["n_cells"],
        "runtime_last_wc_s": final["wc_time_s"],
        "Emax_initial_Vpm": first["Emax_Vpm"],
        "Emax_final_Vpm": final["Emax_Vpm"],
        "Emax_max_Vpm": max(r["Emax_Vpm"] for r in rows),
        "ne_max_initial_m3": first["ne_max_m3"],
        "ne_max_final_m3": final["ne_max_m3"],
        "ne_max_growth_factor": ne_growth,
        "total_electrons_initial": first["total_electrons"],
        "total_electrons_final": final["total_electrons"],
        "total_electrons_growth_factor": total_growth,
        "total_charge_final_C": final["total_charge_C"],
        "head_initial_z_m": first["head_z_m"],
        "head_final_z_m": final["head_z_m"],
        "head_delta_m": head_delta,
        "head_velocity_mps": head_velocity,
        "max_abs_x_cm_m": max(abs(r["x_cm_m"]) for r in rows),
        "max_abs_y_cm_m": max(abs(r["y_cm_m"]) for r in rows),
        "max_r_cm_m": max(r["r_cm_m"] for r in rows),
        "final_x_cm_m": final["x_cm_m"],
        "final_y_cm_m": final["y_cm_m"],
        "final_z_cm_m": final["z_cm_m"],
        "final_ne_moment_asym": final["ne_moment_asym"],
        "final_rho_moment_asym": final["rho_moment_asym"],
        "min_ne_m3": min(r["min_ne_m3"] for r in rows),
        "avalanche": bool(ne_growth > 1.05 or total_growth > 1.05),
        "streamer_propagation": bool(streamer_propagation),
        "bridge": bool(bridge),
        "all_log_values_finite": all_log_finite,
        "uncontrolled_negative_density": bool(min(r["min_ne_m3"] for r in rows) < -1.0e8),
        **profile,
    }
    return case_summary, trajectory


def write_case_comparison(cases: dict[str, dict[str, object]]) -> None:
    path = SUMMARY / "e3_case_comparison.csv"
    fields = [
        "quantity",
        "aligned",
        "misaligned",
        "relative_difference_or_ratio",
    ]
    quantities = [
        "Emax_final_Vpm",
        "Emax_max_Vpm",
        "ne_max_final_m3",
        "total_electrons_final",
        "total_charge_final_C",
        "head_final_z_m",
        "max_r_cm_m",
        "field_plane_asymmetry",
        "ne_plane_asymmetry",
        "rho_plane_asymmetry",
        "final_ne_moment_asym",
        "final_rho_moment_asym",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for quantity in quantities:
            a = float(cases["aligned"][quantity])
            b = float(cases["misaligned"][quantity])
            metric = b / max(abs(a), 1.0e-300) if "asymmetry" in quantity else rel_diff(a, b)
            writer.writerow(
                {
                    "quantity": quantity,
                    "aligned": a,
                    "misaligned": b,
                    "relative_difference_or_ratio": metric,
                }
            )


def write_trajectory(rows: list[dict[str, object]]) -> None:
    path = SUMMARY / "e3_centroid_trajectory.csv"
    fields = [
        "case",
        "time_s",
        "x_cm_m",
        "y_cm_m",
        "z_cm_m",
        "r_cm_m",
        "Emax_Vpm",
        "ne_max_m3",
        "total_electrons",
        "total_charge_C",
        "head_z_m",
        "ne_moment_asym",
        "rho_moment_asym",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    cases: dict[str, dict[str, object]] = {}
    trajectory: list[dict[str, object]] = []
    for case, prefix in CASES.items():
        summary, rows = summarize_case(case, prefix)
        cases[case] = summary
        trajectory.extend(rows)

    aligned_field = float(cases["aligned"]["field_plane_asymmetry"])
    misaligned_field = float(cases["misaligned"]["field_plane_asymmetry"])
    aligned_plasma = max(
        float(cases["aligned"]["ne_plane_asymmetry"]),
        float(cases["aligned"]["rho_plane_asymmetry"]),
        float(cases["aligned"]["final_ne_moment_asym"]),
        float(cases["aligned"]["final_rho_moment_asym"]),
    )
    misaligned_plasma = max(
        float(cases["misaligned"]["ne_plane_asymmetry"]),
        float(cases["misaligned"]["rho_plane_asymmetry"]),
        float(cases["misaligned"]["final_ne_moment_asym"]),
        float(cases["misaligned"]["final_rho_moment_asym"]),
    )
    field_increase = misaligned_field / max(aligned_field, 1.0e-300)
    plasma_increase = misaligned_plasma / max(aligned_plasma, 1.0e-300)
    if field_increase > 2.0 and plasma_increase > 1.25:
        conclusion = "STRONG_EVIDENCE"
    elif field_increase > 2.0:
        conclusion = "EARLY_TIME_EVIDENCE"
    else:
        conclusion = "INSUFFICIENT_FIELD_SYMMETRY_BREAKING"

    summary = {
        "case_id": "Stage E3 controlled aligned-vs-misaligned needle pair",
        "geometry_source": "development",
        "voltage_V": 500.0,
        "photoionization": "OFF",
        "offset_cases_m": {"aligned": 0.0, "misaligned": 10.0e-6},
        "seed": {
            "n0_m3": 1.0e16,
            "sigma_m": 3.0e-6,
            "center_m": [0.0, 0.0, 60.0e-6],
            "ne_equals_np": True,
            "nn_zero": True,
        },
        "needle_geometry": {
            "gap_m": 70.0e-6,
            "tip_radius_m": 5.0e-6,
            "rod_radius_m": 5.0e-6,
            "hv_tip_surface_z_m": 70.0e-6,
            "ground_tip_surface_z_m": 0.0,
            "propagation_direction": "-z",
        },
        "resolution": {
            "finest_dx_target_m": 0.7e-6,
            "actual_min_dx_aligned_m": cases["aligned"]["min_dx_m"],
            "actual_min_dx_misaligned_m": cases["misaligned"]["min_dx_m"],
            "omp_threads": 4,
        },
        "matched_output_times_s": [0.0, 1.0e-12, 3.0e-12, 5.0e-12],
        "cases": cases,
        "field_asymmetry_increase_factor": field_increase,
        "plasma_asymmetry_increase_factor": plasma_increase,
        "controlled_symmetry_breaking_conclusion": conclusion,
        "rho_ne_E_J_export_success": all(
            bool(case["source_schema_ok"]) and bool(case["source_all_finite"]) and int(case["source_rows"]) > 0
            for case in cases.values()
        ),
        "source_export_schema": ",".join(SOURCE_FIELDS),
        "blocked": bool(
            conclusion == "INSUFFICIENT_FIELD_SYMMETRY_BREAKING"
            or not all(bool(c["all_log_values_finite"]) for c in cases.values())
            or any(bool(c["uncontrolled_negative_density"]) for c in cases.values())
        ),
    }

    (SUMMARY / "e3_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_case_comparison(cases)
    write_trajectory(trajectory)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
