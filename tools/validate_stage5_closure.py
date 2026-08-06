#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str]) -> str:
    p = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    if p.returncode:
        raise SystemExit(f"command failed: {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
    return p.stdout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Stage 5 closure validation failed: {message}")


def main() -> None:
    py = str(ROOT / ".venv/bin/python")
    require("Stage 1 closure validation passed" in run([py, "tools/validate_stage1_closure.py"]), "Stage 1 regression failed")
    require("Stage 2 closure validation passed" in run([py, "tools/validate_stage2_closure.py"]), "Stage 2 regression failed")
    require("Stage 3 resource-aware preparation validation passed" in run([py, "tools/validate_stage3_resource_preparation.py"]), "Stage 3 regression failed")
    require("Stage 4 resource-aware closure validation passed" in run([py, "tools/validate_stage4_closure.py"]), "Stage 4 regression failed")
    run([py, "-m", "pytest", "tests/stage5", "-q"])
    require("Stage 5 evidence audit passed" in run([py, "tools/audit_stage5_evidence.py"]), "Stage 5 evidence audit failed")

    root = ROOT / "results/stage5"
    reg = pd.read_csv(root / "provenance/run_registry.csv")
    for rid in [
        "S5-HIGHFIELD-COLLISION",
        "S5-HIGHFIELD-LEFT-ISOLATED",
        "S5-HIGHFIELD-RIGHT-ISOLATED",
        "S5-LARGERGAP-COLLISION",
        "S5-ASYMMETRIC-COLLISION",
    ]:
        require(rid in set(reg.run_id), f"missing run_id {rid}")

    high = pd.read_csv(root / "collision/highfield_collision_event.csv").iloc[0]
    require(high.collision_status == "PASS", "highfield collision did not pass event detection")
    require(float(high.t_collision_s) < 1.953740091585568e-9, "highfield collision did not occur earlier than baseline")

    pulse = pd.read_csv(root / "current_moment/highfield_delta_current_pulse_metrics.csv").iloc[0]
    require(pulse.causality_status == "PASS", "highfield Delta I causality failed")
    require(float(pulse.Delta_I_peak_abs) > 0.0 and float(pulse.Delta_I_snr) >= 5.0, "highfield Delta I pulse not significant")

    sampling = pd.read_csv(root / "radiation/highfield_output_sampling_sensitivity.csv")
    require(set(sampling.status) == {"PASS"}, "output sampling sensitivity failed")

    mpi = pd.read_csv(root / "mpi/rank_consistency.csv")
    require(set(mpi.status) == {"PASS"}, "MPI local consistency failed")

    comp = pd.read_csv(root / "trends/case_comparison.csv")
    require(set(comp.case_id) >= {"B0", "F", "S5-LARGERGAP-COLLISION", "S5-ASYMMETRIC-COLLISION"}, "case comparison incomplete")
    require(comp.loc[comp.case_id == "F", "metric_type"].iloc[0] == "strict_delta", "highfield is not strict_delta")
    require(set(comp[comp.metric_type == "event_local_proxy"].collision_status) <= {"NO_COLLISION_WITHIN_RESOURCE_WINDOW", "PASS", "FAIL"}, "bad proxy collision status")

    matrix = ROOT / "docs/stage5_validation_matrix.csv"
    require(matrix.exists(), "missing Stage 5 validation matrix")
    vm = pd.read_csv(matrix)
    require(len(vm) >= 12, "validation matrix too small")
    require(set(vm.status).issubset({"PASS", "RESOURCE_LIMITED", "SKIPPED_RESOURCE_LIMIT"}), "validation matrix contains invalid status")
    for row in vm[vm.status == "PASS"].itertuples():
        require(isinstance(row.run_id, str) and row.run_id, f"PASS without run_id {row.validation_id}")
        require(row.evidence_status in {"measured", "derived_from_measured"}, f"bad evidence status {row.validation_id}")

    transfer = ROOT / "docs/project_simulation_transfer.md"
    report = ROOT / "docs/stage5_closure_report.md"
    require(transfer.exists() and "当前代码已经完成上述全部课题模型" not in transfer.read_text(), "transfer document missing or overclaims")
    require(report.exists(), "missing closure report")
    text = report.read_text().lower()
    require("shi 2019 figures 1–4 have been quantitatively reproduced" not in text, "forbidden figure reproduction claim")
    require("shi 2019的1:1定量复刻" in report.read_text() or "1:1" in report.read_text(), "missing non-1:1 disclaimer")

    print("Stage 5 resource-aware trend and transfer validation passed with traceable measured evidence.")


if __name__ == "__main__":
    main()
