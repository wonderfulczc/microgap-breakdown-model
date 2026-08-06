#!/usr/bin/env python3
from __future__ import annotations

import hashlib
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
        raise SystemExit(f"Stage 4 closure validation failed: {message}")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def registry_has(run_id: str) -> bool:
    reg = pd.read_csv(ROOT / "results/stage4/provenance/run_registry.csv")
    return run_id in set(reg.run_id)


def main() -> None:
    py = str(ROOT / ".venv/bin/python")

    require("Stage 1 closure validation passed" in run([py, "tools/validate_stage1_closure.py"]), "Stage 1 regression failed")
    require("Stage 2 closure validation passed" in run([py, "tools/validate_stage2_closure.py"]), "Stage 2 regression failed")
    require("Stage 3 resource-aware preparation validation passed" in run([py, "tools/validate_stage3_resource_preparation.py"]), "Stage 3 resource-aware regression failed")
    require("Stage 4 evidence audit passed" in run([py, "tools/audit_stage4_evidence.py"]), "Stage 4 evidence audit failed")

    root = ROOT / "results/stage4"
    for rid in ["S4-PAPERLIKE-20UM", "S4-LEFT-ISOLATED", "S4-RIGHT-ISOLATED"]:
        require(registry_has(rid), f"missing run_id {rid}")

    for name, target in [("paperlike_20um", 3e-9), ("left_isolated_20um", 2.25e-9), ("right_isolated_20um", 2.25e-9)]:
        term = pd.read_csv(root / f"runs/{name}/termination.csv").iloc[0]
        require(term.termination_reason == "reached_end_time", f"{name} did not reach requested end time")
        require(abs(float(term.actual_end_time) - target) <= max(1e-18, target * 1e-12), f"{name} wrong final time")
        require(int(term.rejected_steps) == 0, f"{name} rejected steps present")

    ev = pd.read_csv(root / "collision/collision_event.csv").iloc[0]
    require(ev.collision_status == "PASS", "collision status not PASS")
    require(ev.bridge_status == "PASS", "bridge criterion failed")
    require(ev.field_status == "PASS", "field-collapse criterion failed")
    require(ev.distance_status in {"PASS", "NOISY_NONCONTRADICTORY"}, "distance criterion contradicts event")
    require(float(ev.field_drop_fraction) >= 0.20, "field drop below threshold")
    require(float(ev.bridge_growth_ratio) >= 10.0, "bridge growth too weak")

    pulse = pd.read_csv(root / "current_moment/delta_current_pulse_metrics.csv").iloc[0]
    require(pulse.causality_status == "PASS", "Delta I causality check failed")
    require(float(pulse.Delta_I_peak_abs) > 0.0, "Delta I pulse is zero")
    require(float(pulse.Delta_I_snr) >= 5.0, "Delta I SNR below 5")

    deriv = pd.read_csv(root / "radiation/derivative_sensitivity.csv").iloc[0]
    require(float(deriv.peak_relative_difference) <= 0.20, "derivative method disagreement above 20%")

    bands = pd.read_csv(root / "radiation/band_energy.csv")
    require(set(bands.band) == {"VHF", "UHF", "SHF"}, "missing frequency bands")
    require(set(bands.trusted_status).issubset({"trusted", "partial", "untrusted"}), "invalid band trust label")
    spec = pd.read_csv(root / "radiation/delta_current_spectrum.csv")
    require("trusted" in spec.columns and spec.trusted.any(), "spectrum trust labels missing")

    sampling = pd.read_csv(root / "radiation/output_sampling_sensitivity.csv")
    require(float(sampling.Delta_I_peak_rel_diff.max()) <= 0.20, "output sampling Delta I sensitivity failed")
    require(float(sampling.derivative_peak_rel_diff.max()) <= 0.20, "output sampling derivative sensitivity failed")

    mpi = pd.read_csv(root / "mpi/rank_consistency.csv")
    require(set(mpi.status) == {"PASS"}, "MPI rank consistency failed")

    matrix = ROOT / "docs/stage4_validation_matrix.csv"
    require(matrix.exists(), "missing Stage 4 validation matrix")
    vm = pd.read_csv(matrix)
    require(len(vm) >= 16, "validation matrix too small")
    require(set(vm.status) <= {"PASS", "SKIPPED_RESOURCE_LIMIT"}, "validation matrix contains failures")
    for row in vm[vm.status == "PASS"].itertuples():
        require(isinstance(row.run_id, str) and row.run_id, f"PASS without run_id: {row.validation_id}")
        require(isinstance(row.result_sha256, str) and len(row.result_sha256) == 64, f"PASS without result sha: {row.validation_id}")
        require(row.evidence_status in {"measured", "derived_from_measured"}, f"bad evidence status: {row.validation_id}")

    forbidden_text = ""
    for p in [ROOT / "docs/stage4_closure_report.md"]:
        if p.exists():
            forbidden_text += p.read_text(errors="ignore").lower()
    require("shi 2019 case i has been quantitatively reproduced" not in forbidden_text, "forbidden quantitative reproduction claim")
    # Stage 4 originally checked that Stage 5 had not started. Once Stage 5 is
    # legitimately in progress, Stage 4 closure remains a reusable regression
    # gate; Stage 5 validators check that Stage 5 does not alter Stage 4 values.

    print("Stage 4 resource-aware closure validation passed with traceable measured evidence.")


if __name__ == "__main__":
    main()
