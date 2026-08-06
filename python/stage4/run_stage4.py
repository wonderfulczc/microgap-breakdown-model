#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

import h5py
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_command(run_id: str, cmd: list[str], exe: Path, config: Path, ranks: int, output: Path, allow_nonzero: bool = True) -> dict:
    logs = output / "provenance/logs"
    logs.mkdir(parents=True, exist_ok=True)
    so = logs / f"{run_id}.stdout.log"
    se = logs / f"{run_id}.stderr.log"
    with so.open("w") as stdout_file, se.open("w") as stderr_file:
        p = subprocess.run(cmd, cwd=ROOT, text=True, stdout=stdout_file, stderr=stderr_file)
    if p.returncode and not allow_nonzero:
        raise RuntimeError(f"{run_id} failed with {p.returncode}: {se}")
    return dict(
        run_id=run_id,
        timestamp_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        executable=str(exe.relative_to(ROOT)),
        executable_sha256=sha(exe),
        config=str(config.relative_to(ROOT)),
        config_sha256=sha(config),
        mpi_ranks=ranks,
        command=shlex.join(map(str, cmd)),
        exit_code=p.returncode,
        stdout=str(so.relative_to(ROOT)),
        stderr=str(se.relative_to(ROOT)),
    )


def csv_to_h5(path: Path) -> Path:
    d = pd.read_csv(path)
    target = path.with_suffix(".h5")
    with h5py.File(target, "w") as h:
        for c in d.columns:
            h.create_dataset(c, data=d[c].to_numpy(), compression="gzip", shuffle=True)
        h.attrs["source_csv_sha256"] = sha(path)
    return target


def append_registry(output: Path, entries: list[dict], generated: dict[str, list[Path]]) -> None:
    rows = []
    for base in entries:
        for p in generated[base["run_id"]]:
            if p.is_file():
                rows.append({**base, "result_files": str(p.relative_to(ROOT)), "result_sha256": sha(p)})
    reg = output / "provenance/run_registry.csv"
    reg.parent.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(reg).to_dict("records") if reg.exists() else []
    replaced = {e["run_id"] for e in entries}
    old = [r for r in old if r["run_id"] not in replaced]
    pd.DataFrame(old + rows).to_csv(reg, index=False)


def build_cmd(spec: dict, exe: Path, dest: Path) -> list[str]:
    dr = float(spec["dr_m"])
    nr = round(float(spec["r_max_m"]) / dr)
    nz = round(float(spec["z_max_m"]) / dr)
    n0 = float(spec["n0"])
    nref = n0 * float(spec["n_ref_ratio"])
    cmd = [
        "mpirun", "-np", str(int(spec["mpi_ranks"])),
        str(exe), str(dest), str(nr), str(nz), str(float(spec["t_end_s"])),
        str(float(spec["field_value_V_m"])), str(nref), str(int(bool(spec["sp3"]))),
        str(float(spec["eta"])), str(float(spec["r_max_m"])), str(float(spec["z_max_m"])),
        str(float(spec.get("dt_scale", 1.0))), str(n0), str(float(spec["sigma_m"])),
        str(float(spec["z1_m"])), str(float(spec["z2_m"])), str(spec["mode"]),
        str(int(spec.get("max_steps", 90000))), str(float(spec.get("minimum_dt_s", 1e-15))),
        str(float(spec.get("post_event_s", 2e-10))),
    ]
    if "sigma2_m" in spec:
        cmd += ["--sigma2", str(float(spec["sigma2_m"]))]
    if spec.get("resume_checkpoint"):
        cmd += [str(ROOT / spec["resume_checkpoint"]), str(int(spec.get("initial_step", 0)))]
    return cmd


def standard_current_copy(output: Path, name: str, mode: str, dest: Path) -> Path | None:
    current_root = output / "current_moment"
    current_root.mkdir(parents=True, exist_ok=True)
    target = current_root / {"collision": "collision.csv", "left": "left_isolated.csv", "right": "right_isolated.csv"}.get(mode, f"{name}.csv")
    src = dest / "current_moment.csv"
    if src.exists():
        target.write_bytes(src.read_bytes())
        return target
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/stage4/runs.yaml")
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--output", default="results/stage4")
    ap.add_argument("--run-name")
    args = ap.parse_args()
    output = ROOT / args.output
    cfg = ROOT / args.config
    specs = yaml.safe_load(cfg.read_text())["runs"]
    specs = [s for s in specs if args.run_name is None or s["name"] == args.run_name]
    if args.run_name and not specs:
        raise ValueError(f"unknown run {args.run_name}")
    exe = ROOT / args.build_dir / "bin/stage4_run"
    registry = []
    generated: dict[str, list[Path]] = {}
    for spec in specs:
        run_id = spec.get("run_id", "S4-" + spec["name"].upper().replace(".", "_").replace("-", "_"))
        dest = output / "runs" / spec["name"]
        dest.mkdir(parents=True, exist_ok=True)
        cmd = build_cmd(spec, exe, dest)
        entry = run_command(run_id, cmd, exe, cfg, int(spec["mpi_ranks"]), output, allow_nonzero=True)
        registry.append(entry)
        # Convert sparse full-field outputs only; current moment and scalar data remain CSV first-class evidence.
        for p in dest.glob("fields_*.csv"):
            csv_to_h5(p)
        extra = []
        if spec.get("standard_current", True):
            copied = standard_current_copy(output, spec["name"], spec["mode"], dest)
            if copied is not None:
                extra.append(copied)
        generated[run_id] = list(dest.glob("*")) + extra
    append_registry(output, registry, generated)
    print(json.dumps({"runs": [e["run_id"] for e in registry]}, indent=2))


if __name__ == "__main__":
    main()
