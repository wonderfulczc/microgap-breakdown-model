#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime
import hashlib
import shlex
import sys
from pathlib import Path

import h5py
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.stage4 import detect_collision_event, read_csv_clean


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_to_h5(path: Path) -> Path:
    target = path.with_suffix(".h5")
    if target.exists() and target.stat().st_mtime >= path.stat().st_mtime:
        return target
    d = pd.read_csv(path)
    with h5py.File(target, "w") as h:
        for c in d.columns:
            h.create_dataset(c, data=d[c].to_numpy(), compression="gzip", shuffle=True)
        h.attrs["source_csv_sha256"] = sha(path)
    return target


def find_spec(config: Path, run_name: str) -> dict:
    specs = yaml.safe_load(config.read_text())["runs"]
    for spec in specs:
        if spec["name"] == run_name:
            return spec
    raise ValueError(f"unknown run {run_name}")


def build_command(spec: dict, exe: Path, dest: Path) -> list[str]:
    dr = float(spec["dr_m"])
    nr = round(float(spec["r_max_m"]) / dr)
    nz = round(float(spec["z_max_m"]) / dr)
    nref = float(spec["n0"]) * float(spec["n_ref_ratio"])
    cmd = [
        "mpirun", "-np", str(int(spec["mpi_ranks"])), str(exe), str(dest),
        str(nr), str(nz), str(float(spec["t_end_s"])), str(float(spec["field_value_V_m"])),
        str(nref), str(int(bool(spec["sp3"]))), str(float(spec["eta"])),
        str(float(spec["r_max_m"])), str(float(spec["z_max_m"])), str(float(spec.get("dt_scale", 1.0))),
        str(float(spec["n0"])), str(float(spec["sigma_m"])), str(float(spec["z1_m"])),
        str(float(spec["z2_m"])), str(spec["mode"]), str(int(spec.get("max_steps", 90000))),
        str(float(spec.get("minimum_dt_s", 1e-15))), str(float(spec.get("post_event_s", 2e-10))),
    ]
    if "sigma2_m" in spec:
        cmd += ["--sigma2", str(float(spec["sigma2_m"]))]
    if spec.get("resume_checkpoint"):
        cmd += [str(ROOT / spec["resume_checkpoint"]), str(int(spec.get("initial_step", 0)))]
    return cmd


def append_registry(output: Path, entry: dict, files: list[Path]) -> None:
    rows = []
    for p in files:
        if p.is_file():
            rows.append({**entry, "result_files": str(p.relative_to(ROOT)), "result_sha256": sha(p)})
    reg = output / "provenance/run_registry.csv"
    reg.parent.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(reg).to_dict("records") if reg.exists() else []
    old = [r for r in old if r["run_id"] != entry["run_id"]]
    pd.DataFrame(old + rows).to_csv(reg, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/stage5/runs.yaml")
    ap.add_argument("--output", default="results/stage5")
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--run-name", required=True)
    ap.add_argument("--exit-code", type=int, default=130)
    ap.add_argument("--allow-no-collision", action="store_true")
    args = ap.parse_args()

    output = ROOT / args.output
    config = ROOT / args.config
    spec = find_spec(config, args.run_name)
    run_id = spec["run_id"]
    dest = output / "runs" / args.run_name
    scalar = read_csv_clean(dest / "scalar_history.csv")
    metrics = read_csv_clean(dest / "collision_metrics.csv")
    if len(scalar) == 0 or len(metrics) == 0:
        raise SystemExit("missing clean scalar or collision rows")
    ev = detect_collision_event(metrics) if spec["mode"] == "collision" else None
    actual_end = float(max(scalar.time_s.max(), metrics.time.max()))
    collision_margin = actual_end - (ev.t_collision_s if ev is not None else actual_end)
    no_collision_resource = False
    if ev is not None and (ev.collision_status != "PASS" or collision_margin < float(spec.get("post_event_s", 2e-10))):
        if not args.allow_no_collision:
            raise SystemExit(f"collision window not satisfied: status={ev.collision_status} margin={collision_margin}")
        no_collision_resource = True

    term = pd.DataFrame([{
        "termination_reason": "NO_COLLISION_WITHIN_RESOURCE_WINDOW" if no_collision_resource else ("resource_stop_after_collision_window" if ev is not None else "resource_stop"),
        "requested_end_time": float(spec["t_end_s"]),
        "actual_end_time": actual_end,
        "completed_steps": int(scalar.step.max()),
        "accepted_steps": int(scalar.step.max()),
        "rejected_steps": int(scalar.rejected_steps.iloc[-1]) if "rejected_steps" in scalar else 0,
        "last_dt": float(scalar.dt_s.iloc[-1]),
        "wall_time": "resource_interrupted",
        "checkpoint_path": str((dest / "checkpoint.bin").resolve()),
        "event_time": -1.0 if no_collision_resource else (ev.t_collision_s if ev is not None else -1.0),
        "mpi_ranks": int(spec["mpi_ranks"]),
    }])
    term.to_csv(dest / "termination.csv", index=False)
    (dest / "termination_reason.txt").write_text(str(term.termination_reason.iloc[0]) + "\n")

    generated = list(dest.glob("*.csv")) + list(dest.glob("*.txt")) + list(dest.glob("checkpoint*.bin"))
    for p in dest.glob("fields_*.csv"):
        generated.append(csv_to_h5(p))
    exe = ROOT / args.build_dir / "bin/stage4_run"
    entry = {
        "run_id": run_id,
        "timestamp_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "executable": str(exe.relative_to(ROOT)),
        "executable_sha256": sha(exe),
        "config": str(config.relative_to(ROOT)),
        "config_sha256": sha(config),
        "mpi_ranks": int(spec["mpi_ranks"]),
        "command": shlex.join(build_command(spec, exe, dest)),
        "exit_code": args.exit_code,
        "stdout": str((output / f"provenance/logs/{run_id}.stdout.log").relative_to(ROOT)),
        "stderr": str((output / f"provenance/logs/{run_id}.stderr.log").relative_to(ROOT)),
    }
    append_registry(output, entry, generated)
    print(term.to_string(index=False))


if __name__ == "__main__":
    main()
