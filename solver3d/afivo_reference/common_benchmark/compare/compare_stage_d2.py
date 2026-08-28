#!/usr/bin/env python3
"""Small Stage D2 comparison helpers for PETSc/Afivo common benchmarks."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines))


def read_afivo_lineout(path: Path, gap: float) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open() as f:
        header = f.readline().strip().lstrip("#").split()
        for line in f:
            if not line.strip():
                continue
            vals = [float(x) for x in line.split()]
            data = dict(zip(header, vals))
            rows.append(
                {
                    "z_m": data["z"],
                    "z_over_gap": data["z"] / gap,
                    "phi_V": data["phi"],
                    "Eabs_Vpm": data["electric_fld"],
                }
            )
    return rows


def interp(rows: list[dict[str, float]], xkey: str, ykey: str, x: float) -> float:
    data = sorted(rows, key=lambda r: r[xkey])
    if x <= data[0][xkey]:
        return data[0][ykey]
    if x >= data[-1][xkey]:
        return data[-1][ykey]
    for a, b in zip(data, data[1:]):
        if a[xkey] <= x <= b[xkey]:
            t = (x - a[xkey]) / (b[xkey] - a[xkey])
            return (1.0 - t) * a[ykey] + t * b[ykey]
    raise RuntimeError(f"interpolation failed at {x}")


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def max_column(path: Path, names: list[str]) -> float:
    rows = read_csv(path)
    out = 0.0
    for row in rows:
        for name in names:
            out = max(out, abs(float(row[name])))
    return out


def electrostatic(args: argparse.Namespace) -> int:
    petsc = [
        {
            "z_m": float(r["z_m"]),
            "z_over_gap": float(r["z_over_gap"]),
            "phi_V": float(r["phi_V"]),
            "Eabs_Vpm": float(r["Eabs_Vpm"]),
        }
        for r in read_csv(args.petsc)
        if r["cell_class"] == "gas" and 0.02 <= float(r["z_over_gap"]) <= 0.98
    ]
    afivo = [
        r
        for r in read_afivo_lineout(args.afivo, args.gap)
        if 0.02 <= r["z_over_gap"] <= 0.98
    ]
    samples = [0.1, 0.25, 0.5, 0.75, 0.9]
    selected = {}
    for s in samples:
        pp = interp(petsc, "z_over_gap", "phi_V", s)
        ap = interp(afivo, "z_over_gap", "phi_V", s)
        pe = interp(petsc, "z_over_gap", "Eabs_Vpm", s)
        ae = interp(afivo, "z_over_gap", "Eabs_Vpm", s)
        selected[f"{s:.2f}"] = {
            "phi_rel_err": rel(pp / args.voltage, ap / args.voltage),
            "E_rel_err": rel(pe, ae),
            "petsc_E_Vpm": pe,
            "afivo_E_Vpm": ae,
        }

    grid = [0.05 + 0.9 * i / 90 for i in range(91)]
    phi_sq = 0.0
    e_sq = 0.0
    e_ref_sq = 0.0
    for x in grid:
        pp = interp(petsc, "z_over_gap", "phi_V", x) / args.voltage
        ap = interp(afivo, "z_over_gap", "phi_V", x) / args.voltage
        pe = interp(petsc, "z_over_gap", "Eabs_Vpm", x)
        ae = interp(afivo, "z_over_gap", "Eabs_Vpm", x)
        phi_sq += (pp - ap) ** 2
        e_sq += (pe - ae) ** 2
        e_ref_sq += pe**2
    petsc_emax = max(r["Eabs_Vpm"] for r in petsc)
    afivo_emax = max(r["Eabs_Vpm"] for r in afivo)
    result = {
        "phi_normalized_rmse": math.sqrt(phi_sq / len(grid)),
        "E_profile_normalized_rmse": math.sqrt(e_sq / max(e_ref_sq, 1e-300)),
        "selected_point_errors": selected,
        "Emax_petsc_Vpm": petsc_emax,
        "Emax_afivo_Vpm": afivo_emax,
        "Emax_rel_err": rel(petsc_emax, afivo_emax),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


def b0(args: argparse.Namespace) -> int:
    result = {
        "transport_max_rel_err": max_column(
            args.transport,
            [
                "mu_rel_err",
                "D_rel_err",
                "alpha_rel_err",
                "eta2_rel_err",
                "eta3_rel_err",
                "eta_total_rel_err",
            ],
        ),
        "chemistry_max_rel_err": max_column(args.chemistry, ["max_rel_err"]),
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0 if result["transport_max_rel_err"] < 5e-3 and result["chemistry_max_rel_err"] < 1e-12 else 2


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_b0 = sub.add_parser("b0")
    p_b0.add_argument("--transport", type=Path, required=True)
    p_b0.add_argument("--chemistry", type=Path, required=True)
    p_b0.add_argument("--out")
    p_b0.set_defaults(func=b0)
    p_es = sub.add_parser("electrostatic")
    p_es.add_argument("--petsc", type=Path, required=True)
    p_es.add_argument("--afivo", type=Path, required=True)
    p_es.add_argument("--gap", type=float, default=70e-6)
    p_es.add_argument("--voltage", type=float, default=500.0)
    p_es.add_argument("--out")
    p_es.set_defaults(func=electrostatic)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
