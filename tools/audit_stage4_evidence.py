#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 4 evidence audit failed: {message}")


def main() -> None:
    root = ROOT / "results/stage4"
    registry_path = root / "provenance/run_registry.csv"
    require(registry_path.exists(), "missing run registry")
    registry = pd.read_csv(registry_path)
    required_runs = {
        "S4-PAPERLIKE-20UM",
        "S4-LEFT-ISOLATED",
        "S4-RIGHT-ISOLATED",
        "S4-MPI-PROBE-1RANK",
        "S4-MPI-PROBE-2RANK",
    }
    require(required_runs.issubset(set(registry.run_id)), "missing required run ids")
    for row in registry.itertuples():
        p = ROOT / row.result_files
        require(p.exists(), f"registered result missing: {row.result_files}")
        require(sha(p) == row.result_sha256, f"sha mismatch: {row.result_files}")
        require(int(row.exit_code) == 0, f"nonzero registered run: {row.run_id}")

    critical = [
        "runs/paperlike_20um/termination.csv",
        "runs/left_isolated_20um/termination.csv",
        "runs/right_isolated_20um/termination.csv",
        "collision/collision_event.csv",
        "current_moment/collision.csv",
        "current_moment/left_isolated.csv",
        "current_moment/right_isolated.csv",
        "current_moment/delta_current_moment.csv",
        "radiation/delta_current_derivative.csv",
        "radiation/delta_current_spectrum.csv",
        "radiation/delta_current_esd.csv",
        "radiation/band_energy.csv",
        "radiation/output_sampling_sensitivity.csv",
        "mpi/rank_consistency.csv",
        "provenance/derived_artifact_manifest.csv",
    ]
    for rel in critical:
        require((root / rel).exists(), f"missing critical artifact {rel}")

    manifest = pd.read_csv(root / "provenance/derived_artifact_manifest.csv")
    for row in manifest.itertuples():
        p = ROOT / row.file
        require(p.exists(), f"derived artifact missing: {row.file}")
        require(sha(p) == row.sha256, f"derived artifact sha mismatch: {row.file}")

    print("Stage 4 evidence audit passed with traceable measured evidence.")


if __name__ == "__main__":
    main()

