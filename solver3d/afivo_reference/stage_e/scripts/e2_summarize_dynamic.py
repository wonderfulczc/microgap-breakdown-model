#!/usr/bin/env python3
"""Summarize Stage E2 triangular-foil electrostatic and dynamic runs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


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
    return rows


def parse_line(prefix: str) -> list[dict[str, float]]:
    path = RAW / f"{prefix}_line_000000.txt"
    with path.open() as f:
        header = f.readline().strip().lstrip("#").split()
        return [dict(zip(header, [float(v) for v in line.split()])) for line in f if line.strip()]


def normalized_l2(a: list[float], b: list[float]) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return math.nan
    num = math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(n)))
    den = math.sqrt(max(sum(a[i] ** 2 for i in range(n)), 1.0e-300))
    return num / den


def latest_source(prefix: str) -> Path:
    matches = sorted(RAW.glob(f"{prefix}_source_*.csv"))
    if not matches:
        raise FileNotFoundError(f"no source snapshots for {prefix}")
    return matches[-1]


def source_asymmetry(prefix: str, min_dx: float) -> dict[str, float | bool | int | str]:
    path = latest_source(prefix)
    tol = max(0.51 * min_dx, 1.0e-12)
    bin_dx = max(min_dx, 1.0e-12)
    x_ne: dict[tuple[int, int], list[float]] = {}
    y_ne: dict[tuple[int, int], list[float]] = {}
    x_rho: dict[tuple[int, int], list[float]] = {}
    y_rho: dict[tuple[int, int], list[float]] = {}
    finite = True
    rows = 0
    with path.open() as f:
        reader = csv.DictReader(f)
        expected = {
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
        }
        schema_ok = set(reader.fieldnames or []) == expected
        for row in reader:
            rows += 1
            vals = {k: float(row[k]) for k in expected if k != "level"}
            if not all(math.isfinite(v) for v in vals.values()):
                finite = False
            if vals["lsf_m"] <= 0.0:
                continue
            x = vals["x_m"]
            y = vals["y_m"]
            z = vals["z_m"]
            key_x = (round(abs(x) / bin_dx), round(z / bin_dx))
            key_y = (round(abs(y) / bin_dx), round(z / bin_dx))
            if abs(y) <= tol:
                x_ne.setdefault(key_x, []).append(vals["ne_m3"])
                x_rho.setdefault(key_x, []).append(vals["rho_Cpm3"])
            if abs(x) <= tol:
                y_ne.setdefault(key_y, []).append(vals["ne_m3"])
                y_rho.setdefault(key_y, []).append(vals["rho_Cpm3"])

    def avg_map(data: dict[tuple[int, int], list[float]]) -> dict[tuple[int, int], float]:
        return {k: sum(v) / len(v) for k, v in data.items()}

    x_ne_avg = avg_map(x_ne)
    y_ne_avg = avg_map(y_ne)
    x_rho_avg = avg_map(x_rho)
    y_rho_avg = avg_map(y_rho)
    common_ne = sorted(set(x_ne_avg) & set(y_ne_avg))
    common_rho = sorted(set(x_rho_avg) & set(y_rho_avg))
    ne_metric = normalized_l2([x_ne_avg[k] for k in common_ne], [y_ne_avg[k] for k in common_ne])
    rho_metric = normalized_l2([x_rho_avg[k] for k in common_rho], [y_rho_avg[k] for k in common_rho])
    return {
        "source_file": str(path.relative_to(ROOT)),
        "source_rows": rows,
        "source_schema_ok": schema_ok,
        "source_all_finite": finite,
        "common_ne_bins": len(common_ne),
        "common_rho_bins": len(common_rho),
        "ne_plane_asymmetry": ne_metric,
        "rho_plane_asymmetry": rho_metric,
    }


def electrostatic_summary() -> dict[str, object]:
    e1 = json.loads((SUMMARY / "e1_electrostatic_summary.json").read_text())
    refined_x = parse_log("e2_refined_x")[-1]
    refined_y = parse_log("e2_refined_y")[-1]
    medium_x_line = parse_line("medium_x")
    medium_y_line = parse_line("medium_y")
    refined_x_line = parse_line("e2_refined_x")
    refined_y_line = parse_line("e2_refined_y")
    refined_asym = normalized_l2(
        [r["electric_fld"] for r in refined_x_line],
        [r["electric_fld"] for r in refined_y_line],
    )
    return {
        "medium_Emax_Vpm": e1["medium"]["Emax_Vpm"],
        "refined_Emax_Vpm": refined_x["Emax_Vpm"],
        "medium_to_refined_Emax_relative_change": abs(refined_x["Emax_Vpm"] - e1["medium"]["Emax_Vpm"])
        / max(abs(e1["medium"]["Emax_Vpm"]), 1.0),
        "medium_asymmetry_metric": e1["medium_profile"]["electrostatic_asymmetry_metric"],
        "refined_asymmetry_metric": refined_asym,
        "asymmetry_metric_change": refined_asym - e1["medium_profile"]["electrostatic_asymmetry_metric"],
        "x_profile_medium_to_refined_l2": normalized_l2(
            [r["electric_fld"] for r in medium_x_line],
            [r["electric_fld"] for r in refined_x_line],
        ),
        "y_profile_medium_to_refined_l2": normalized_l2(
            [r["electric_fld"] for r in medium_y_line],
            [r["electric_fld"] for r in refined_y_line],
        ),
        "refined_min_dx_m": refined_x["min_dx_m"],
        "refined_highest_level": refined_x["highest_level"],
        "refined_n_cells": refined_x["n_cells"],
        "refined_high_field_location_m": [
            refined_x["Emax_user_x_m"],
            refined_x["Emax_user_y_m"],
            refined_x["Emax_user_z_m"],
        ],
        "refined_high_field_volumes_m3": {
            "E_gt_0p5_Emax": refined_x["vol_E_gt_0p5_Emax_m3"],
            "E_gt_0p8_Emax": refined_x["vol_E_gt_0p8_Emax_m3"],
        },
        "field_topology_stable": bool(
            abs(refined_x["Emax_user_z_m"] - 70.0e-6) < 5.0e-6
            and math.hypot(refined_x["Emax_user_x_m"], refined_x["Emax_user_y_m"]) < 12.0e-6
        ),
    }


def dynamic_case(prefix: str) -> dict[str, object]:
    rows = parse_log(prefix)
    if len(rows) < 2:
        raise RuntimeError(f"{prefix}: expected at least two log rows")
    first = rows[0]
    final = rows[-1]
    ne_growth = final["ne_max_m3"] / max(first["ne_max_m3"], 1.0)
    total_growth = final["total_electrons"] / max(first["total_electrons"], 1.0)
    head_rows = [r for r in rows if r["head_z_m"] > -1.0e90]
    head_delta = 0.0
    head_velocity = 0.0
    if len(head_rows) >= 2:
        head_delta = head_rows[-1]["head_z_m"] - head_rows[0]["head_z_m"]
        dt = head_rows[-1]["time_s"] - head_rows[0]["time_s"]
        if dt > 0:
            head_velocity = head_delta / dt
    src = source_asymmetry(prefix, final["min_dx_m"])
    all_finite_log = all(
        math.isfinite(v) for row in rows for k, v in row.items() if k != "head_z_m" or v > -1.0e90
    )
    negative_min = min(r["min_ne_m3"] for r in rows)
    return {
        "n_outputs": len(rows),
        "final_time_s": final["time_s"],
        "final_Emax_Vpm": final["Emax_Vpm"],
        "final_ne_max_m3": final["ne_max_m3"],
        "final_total_electrons": final["total_electrons"],
        "final_total_charge_C": final["total_charge_C"],
        "initial_ne_max_m3": first["ne_max_m3"],
        "initial_total_electrons": first["total_electrons"],
        "ne_max_growth_factor": ne_growth,
        "total_electron_growth_factor": total_growth,
        "min_ne_m3": negative_min,
        "max_abs_x_cm_m": max(abs(r["x_cm_m"]) for r in rows),
        "max_abs_y_cm_m": max(abs(r["y_cm_m"]) for r in rows),
        "max_r_cm_m": max(r["r_cm_m"] for r in rows),
        "final_x_cm_m": final["x_cm_m"],
        "final_y_cm_m": final["y_cm_m"],
        "final_head_z_m": final["head_z_m"],
        "head_delta_m": head_delta,
        "head_velocity_mps": head_velocity,
        "final_ne_moment_asym": final["ne_moment_asym"],
        "final_rho_moment_asym": final["rho_moment_asym"],
        "avalanche": bool(ne_growth > 1.05 or total_growth > 1.05),
        "streamer_propagation": bool(abs(head_delta) > max(2.0 * final["min_dx_m"], 2.0e-6)),
        "bridge": bool(final["head_z_m"] > -1.0e90 and final["head_z_m"] <= 2.0e-6),
        "all_log_values_finite": all_finite_log,
        "uncontrolled_negative_density": bool(negative_min < -1.0e8),
        **src,
    }


def write_dynamic_csv(cases: dict[str, dict[str, object]]) -> None:
    path = SUMMARY / "e2_dynamic_diagnostics.csv"
    fields = [
        "case",
        "final_time_s",
        "final_Emax_Vpm",
        "final_ne_max_m3",
        "final_total_electrons",
        "final_total_charge_C",
        "final_head_z_m",
        "head_velocity_mps",
        "final_x_cm_m",
        "final_y_cm_m",
        "ne_plane_asymmetry",
        "rho_plane_asymmetry",
        "avalanche",
        "streamer_propagation",
        "bridge",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name, data in cases.items():
            row = {"case": name}
            row.update({k: data.get(k, "") for k in fields if k != "case"})
            writer.writerow(row)


def main() -> None:
    SUMMARY.mkdir(parents=True, exist_ok=True)
    electro = electrostatic_summary()
    cases = {
        "case_a_pi_off_500V": dynamic_case("e2_case_a_pi_off_500V"),
        "case_b_pi_on_500V": dynamic_case("e2_case_b_pi_on_500V"),
    }
    non_axis = "NO"
    if any(float(case["ne_plane_asymmetry"]) > 0.05 or float(case["rho_plane_asymmetry"]) > 0.05 for case in cases.values()):
        non_axis = "YES"
    elif any(float(case["ne_plane_asymmetry"]) > 0.01 or float(case["rho_plane_asymmetry"]) > 0.01 for case in cases.values()):
        non_axis = "WEAK"
    summary = {
        "electrostatic_refined": electro,
        "dynamic_cases": cases,
        "final_dynamic_voltage_V": 500.0,
        "seed": {
            "n0_m3": 1.0e16,
            "sigma_m": 3.0e-6,
            "center_m": [0.0, 0.0, 60.0e-6],
            "ne_equals_np": True,
            "nn_zero": True,
        },
        "NON_AXISYMMETRIC_STREAMER_RESPONSE": non_axis,
        "source_export_schema": "time_s,x_m,y_m,z_m,cell_volume_m3,rho_Cpm3,ne_m3,Ex_Vpm,Ey_Vpm,Ez_Vpm,Jx_Apm2,Jy_Apm2,Jz_Apm2,Eabs_Vpm,lsf_m,level",
    }
    (SUMMARY / "e2_dynamic_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_dynamic_csv(cases)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

