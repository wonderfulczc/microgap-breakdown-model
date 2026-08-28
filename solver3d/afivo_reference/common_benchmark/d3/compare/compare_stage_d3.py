#!/usr/bin/env python3
"""Stage D3 PETSc2D/Afivo3D dynamic benchmark comparison."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
from pathlib import Path


HEAD_THRESHOLD = 1.0e14


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines))


def interp(rows: list[dict[str, float]], key: str, t: float) -> float:
    data = sorted(rows, key=lambda r: r["time_s"])
    if t <= data[0]["time_s"]:
        return data[0][key]
    if t >= data[-1]["time_s"]:
        return data[-1][key]
    for a, b in zip(data, data[1:]):
        if a["time_s"] <= t <= b["time_s"]:
            w = (t - a["time_s"]) / (b["time_s"] - a["time_s"])
            return (1.0 - w) * a[key] + w * b[key]
    raise RuntimeError(f"time interpolation failed: {t}")


def read_petsc(path: Path) -> list[dict[str, float]]:
    out = []
    for r in read_csv(path):
        out.append(
            {
                "time_s": float(r["time_s"]),
                "Emax_Vpm": float(r["Emax_Vpm"]),
                "ne_max_m3": float(r["ne_max_m3"]),
                "total_electrons": float(r["total_electrons"]),
                "head_z_m": float(r["head_z_m"]),
                "head_velocity_mps": float(r["head_velocity_mps"]),
                "bridge_flag": float(r["bridge_flag"]),
            }
        )
    return out


def petsc_line_stats(path: Path) -> tuple[float, float]:
    rows = []
    with path.open() as f:
        for r in csv.DictReader(f):
            rows.append(r)
    hits = [float(r["z_m"]) for r in rows if float(r["ne_m3"]) >= HEAD_THRESHOLD]
    head = min(hits) if hits else max(rows, key=lambda r: float(r["ne_m3"]))["z_m"]
    ne_max = max(float(r["ne_m3"]) for r in rows)
    return float(head), ne_max


def read_petsc_with_line(path: Path, line_glob: str) -> list[dict[str, float]]:
    rows = read_petsc(path)
    lineouts = sorted(glob.glob(line_glob))
    prev_head = math.nan
    prev_time = math.nan
    for i, row in enumerate(rows):
        if i < len(lineouts):
            row["head_z_m"], row["ne_max_m3"] = petsc_line_stats(Path(lineouts[i]))
        time = row["time_s"]
        row["head_velocity_mps"] = (
            0.0 if not math.isfinite(prev_head) or time <= prev_time else (row["head_z_m"] - prev_head) / (time - prev_time)
        )
        prev_head = row["head_z_m"]
        prev_time = time
    return rows


def afivo_line_stats(path: Path) -> tuple[float, float]:
    rows = []
    with path.open() as f:
        header = f.readline().strip().lstrip("#").split()
        for line in f:
            vals = [float(x) for x in line.split()]
            rows.append(dict(zip(header, vals)))
    hits = [r["z"] for r in rows if r.get("e", 0.0) >= HEAD_THRESHOLD]
    if hits:
        head = min(hits)
    else:
        head = max(rows, key=lambda r: r.get("e", 0.0))["z"]
    ne_max = max(r.get("e", 0.0) for r in rows)
    return head, ne_max


def read_afivo_log(path: Path, line_glob: str) -> list[dict[str, float]]:
    numeric = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            numeric.append(vals)
    lineouts = sorted(glob.glob(line_glob))
    out = []
    prev_head = math.nan
    prev_time = math.nan
    for i, vals in enumerate(numeric):
        if i < len(lineouts):
            head, ne_max = afivo_line_stats(Path(lineouts[i]))
        else:
            head, ne_max = vals[19], vals[12]
        time = vals[1]
        hv = 0.0 if not math.isfinite(prev_head) or time <= prev_time else (head - prev_head) / (time - prev_time)
        prev_head = head
        prev_time = time
        row = {
            "time_s": time,
            "Emax_Vpm": vals[8],
            "ne_max_m3": ne_max,
            "total_electrons": vals[4],
            "head_z_m": head,
            "head_velocity_mps": hv,
            "bridge_flag": 0.0,
            "x_cm": vals[33] if len(vals) > 33 else 0.0,
            "y_cm": vals[34] if len(vals) > 34 else 0.0,
            "symmetry_metric": vals[35] if len(vals) > 35 else math.nan,
            "min_dx": vals[27] if len(vals) > 27 else math.nan,
        }
        out.append(row)
    return out


def dynamic(args: argparse.Namespace) -> int:
    petsc = read_petsc_with_line(args.petsc, args.petsc_line_glob)
    afivo = read_afivo_log(args.afivo_log, args.afivo_line_glob)
    times = [0.0, 0.5e-12, 1.0e-12, 2.0e-12, 3.0e-12, 5.0e-12]
    rows = []
    metrics = {
        "max_Emax_rel_err": 0.0,
        "max_ne_max_rel_err": 0.0,
        "max_total_electrons_rel_err": 0.0,
        "max_head_position_diff_m": 0.0,
        "max_head_velocity_rel_err": 0.0,
        "afivo_max_cm_offset_m": 0.0,
        "afivo_max_symmetry_metric": 0.0,
        "avalanche": False,
        "streamer_propagation": False,
        "bridge": False,
    }
    for t in times:
        pe = interp(petsc, "Emax_Vpm", t)
        ae = interp(afivo, "Emax_Vpm", t)
        pn = interp(petsc, "ne_max_m3", t)
        an = interp(afivo, "ne_max_m3", t)
        pt = interp(petsc, "total_electrons", t)
        at = interp(afivo, "total_electrons", t)
        ph = interp(petsc, "head_z_m", t)
        ah = interp(afivo, "head_z_m", t)
        pv = interp(petsc, "head_velocity_mps", t)
        av = interp(afivo, "head_velocity_mps", t)
        xcm = interp(afivo, "x_cm", t)
        ycm = interp(afivo, "y_cm", t)
        sym = interp(afivo, "symmetry_metric", t)
        row = {
            "time_s": t,
            "petsc_Emax": pe,
            "afivo_Emax": ae,
            "error_Emax": rel(pe, ae),
            "petsc_ne_max": pn,
            "afivo_ne_max": an,
            "error_ne_max": rel(pn, an),
            "petsc_total_e": pt,
            "afivo_total_e": at,
            "error_total_e": rel(pt, at),
            "petsc_head_z": ph,
            "afivo_head_z": ah,
            "head_difference": abs(ph - ah),
            "petsc_head_v": pv,
            "afivo_head_v": av,
            "afivo_x_cm": xcm,
            "afivo_y_cm": ycm,
            "afivo_symmetry_metric": sym,
        }
        rows.append(row)
        metrics["max_Emax_rel_err"] = max(metrics["max_Emax_rel_err"], row["error_Emax"])
        metrics["max_ne_max_rel_err"] = max(metrics["max_ne_max_rel_err"], row["error_ne_max"])
        metrics["max_total_electrons_rel_err"] = max(metrics["max_total_electrons_rel_err"], row["error_total_e"])
        metrics["max_head_position_diff_m"] = max(metrics["max_head_position_diff_m"], row["head_difference"])
        metrics["max_head_velocity_rel_err"] = max(metrics["max_head_velocity_rel_err"], rel(pv, av) if t > 0 else 0.0)
        metrics["afivo_max_cm_offset_m"] = max(metrics["afivo_max_cm_offset_m"], math.hypot(xcm, ycm))
        metrics["afivo_max_symmetry_metric"] = max(metrics["afivo_max_symmetry_metric"], abs(sym))
    metrics["avalanche"] = max(r["total_electrons"] for r in petsc) > 1.01 * petsc[0]["total_electrons"] and max(
        r["total_electrons"] for r in afivo
    ) > 1.01 * afivo[0]["total_electrons"]
    metrics["streamer_propagation"] = (
        abs(petsc[-1]["head_z_m"] - petsc[0]["head_z_m"]) > 2e-6
        and abs(afivo[-1]["head_z_m"] - afivo[0]["head_z_m"]) > 2e-6
    )
    metrics["bridge"] = bool(max(r["bridge_flag"] for r in petsc) > 0.5)

    with args.out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    args.out_json.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


def chemistry(args: argparse.Namespace) -> int:
    refs = {r["case_id"]: r for r in read_csv(args.reference)}
    max_err = 0.0
    cases = {}
    for case_id, amounts in args.amounts:
        rows = []
        with Path(amounts).open() as f:
            for line in f:
                vals = [float(x) for x in line.split()]
                if vals:
                    rows.append(vals)
        ref = refs[case_id]
        dt = rows[-1][0] - rows[0][0]
        volume = float(ref["volume_m3"])
        dne = (rows[-1][1] - rows[0][1]) / (dt * volume)
        dnp = (rows[-1][2] - rows[0][2]) / (dt * volume)
        dnn = (rows[-1][3] - rows[0][3]) / (dt * volume)
        errs = [
            rel(dne, float(ref["dne_dt_m3s"])),
            rel(dnp, float(ref["dnp_dt_m3s"])),
            rel(dnn, float(ref["dnn_dt_m3s"])),
        ]
        max_err = max(max_err, *errs)
        cases[case_id] = {"dne_dt": dne, "dnp_dt": dnp, "dnn_dt": dnn, "max_rel_err": max(errs)}
    result = {"runtime_chemistry_max_rel_err": max_err, "cases": cases}
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if max_err < 1e-2 else 2


def seed(args: argparse.Namespace) -> int:
    petsc0 = read_petsc_with_line(args.petsc, args.petsc_line_glob)[0]
    afivo0 = read_afivo_log(args.afivo_log, args.afivo_line_glob)[0]
    result = {
        "petsc_ne_max": petsc0["ne_max_m3"],
        "afivo_ne_max": afivo0["ne_max_m3"],
        "ne_max_rel_err": rel(petsc0["ne_max_m3"], afivo0["ne_max_m3"]),
        "petsc_total_electrons": petsc0["total_electrons"],
        "afivo_total_electrons": afivo0["total_electrons"],
        "total_electrons_rel_err": rel(petsc0["total_electrons"], afivo0["total_electrons"]),
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["ne_max_rel_err"] < 0.05 and result["total_electrons_rel_err"] < 0.05 else 2


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    pc = sub.add_parser("chemistry")
    pc.add_argument("--reference", type=Path, required=True)
    pc.add_argument("--amounts", nargs=2, action="append", metavar=("CASE", "FILE"), required=True)
    pc.add_argument("--out", type=Path, required=True)
    pc.set_defaults(func=chemistry)
    ps = sub.add_parser("seed")
    ps.add_argument("--petsc", type=Path, required=True)
    ps.add_argument("--petsc-line-glob", required=True)
    ps.add_argument("--afivo-log", type=Path, required=True)
    ps.add_argument("--afivo-line-glob", required=True)
    ps.add_argument("--out", type=Path, required=True)
    ps.set_defaults(func=seed)
    pd = sub.add_parser("dynamic")
    pd.add_argument("--petsc", type=Path, required=True)
    pd.add_argument("--petsc-line-glob", required=True)
    pd.add_argument("--afivo-log", type=Path, required=True)
    pd.add_argument("--afivo-line-glob", required=True)
    pd.add_argument("--out-csv", type=Path, required=True)
    pd.add_argument("--out-json", type=Path, required=True)
    pd.set_defaults(func=dynamic)
    args = p.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True) if hasattr(args, "out") else None
    if hasattr(args, "out_csv"):
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
