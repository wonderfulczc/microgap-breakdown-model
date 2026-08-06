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
        raise SystemExit(f"Stage 5 evidence audit failed: {message}")


def main() -> None:
    root = ROOT / "results/stage5"
    reg_path = root / "provenance/run_registry.csv"
    require(reg_path.exists(), "missing run registry")
    reg = pd.read_csv(reg_path)
    required = {
        "S5-HIGHFIELD-COLLISION",
        "S5-HIGHFIELD-LEFT-ISOLATED",
        "S5-HIGHFIELD-RIGHT-ISOLATED",
        "S5-LARGERGAP-COLLISION",
        "S5-ASYMMETRIC-COLLISION",
        "S5-HIGHFIELD-MPI-PROBE-1RANK",
        "S5-HIGHFIELD-MPI-PROBE-2RANK",
    }
    require(required.issubset(set(reg.run_id)), f"missing run ids {required - set(reg.run_id)}")
    for row in reg.itertuples():
        p = ROOT / row.result_files
        require(p.exists(), f"registered result missing: {row.result_files}")
        require(sha(p) == row.result_sha256, f"sha mismatch: {row.result_files}")
        code = int(row.exit_code)
        if row.run_id in {"S5-HIGHFIELD-COLLISION", "S5-LARGERGAP-COLLISION", "S5-ASYMMETRIC-COLLISION"}:
            require(code in {0, 130}, f"unexpected resource-aware exit code {code}: {row.run_id}")
        else:
            require(code == 0, f"nonzero isolated/probe run: {row.run_id}")

    critical = [
        "runs/highfield_collision_20um/termination.csv",
        "runs/highfield_left_isolated_20um/termination.csv",
        "runs/highfield_right_isolated_20um/termination.csv",
        "runs/largergap_collision_20um/termination.csv",
        "runs/asymmetric_collision_20um/termination.csv",
        "collision/highfield_collision_event.csv",
        "current_moment/highfield_delta_current_moment.csv",
        "current_moment/highfield_delta_current_pulse_metrics.csv",
        "radiation/highfield_delta_current_derivative.csv",
        "radiation/highfield_delta_current_spectrum.csv",
        "radiation/highfield_delta_current_esd.csv",
        "radiation/highfield_band_energy.csv",
        "radiation/highfield_output_sampling_sensitivity.csv",
        "trends/case_comparison.csv",
        "mpi/rank_consistency.csv",
        "provenance/derived_artifact_manifest.csv",
    ]
    for rel in critical:
        require((root / rel).exists(), f"missing critical artifact {rel}")

    manifest = pd.read_csv(root / "provenance/derived_artifact_manifest.csv")
    for row in manifest.itertuples():
        p = ROOT / row.file
        require(p.exists(), f"derived artifact missing: {row.file}")
        require(sha(p) == row.sha256, f"derived sha mismatch: {row.file}")

    comp = pd.read_csv(root / "trends/case_comparison.csv")
    require(set(comp.metric_type).issubset({"strict_delta", "event_local_proxy", "not_available"}), "invalid metric_type")
    require(comp.loc[comp.metric_type == "event_local_proxy", "Delta_I_peak"].isna().all(), "proxy cases contain strict Delta I")
    print("Stage 5 evidence audit passed with traceable measured evidence.")


if __name__ == "__main__":
    main()
