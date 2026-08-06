#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
failures: list[str] = []


def check(condition: bool, message: str) -> None:
    if not condition:
        failures.append(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_registry() -> None:
    reg = ROOT / "results/stage3/provenance/run_registry.csv"
    check(reg.exists(), "run registry missing")
    if not reg.exists():
        return
    rows = list(csv.DictReader(reg.open()))
    run_ids = {r["run_id"] for r in rows}
    check("S3-COARSE_ML_RECOVERY_RESUME1" in run_ids, "complete 20um recovery run not registered")
    check("S3-SP3_FORK_OFF_RESUME1" in run_ids, "SP3-off fork run not registered")
    for r in rows:
        p = ROOT / r["result_files"]
        check(p.exists(), f"registered result missing: {p}")
        if p.exists():
            check(sha256(p) == r["result_sha256"], f"registered result hash mismatch: {p}")


def check_complete_20um() -> None:
    term = pd.read_csv(ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/termination.csv").iloc[0]
    check(term.termination_reason == "reached_end_time", "20um recovery did not reach requested end time")
    check(float(term.actual_end_time) >= 2.0e-9, "20um recovery did not reach 2 ns")
    check(int(term.rejected_steps) == 0, "20um recovery had rejected steps")

    summary = pd.read_csv(ROOT / "results/stage3/closure/run_summary.csv")
    row = summary[summary.run == "coarse_ml_recovery_resume1"]
    check(len(row) == 1, "20um recovery missing from independent summary")
    if len(row):
        r = row.iloc[0]
        check(bool(r.formed), "20um connected-front double-headed propagation not formed")
        check(bool(r.opposite_charge), "20um heads do not have opposite charge")
        check(bool(r.head_definitions_agree), "20um head definitions are inconsistent")
        check(float(r.photoionization_ahead) > 0.0, "20um positive-head photoionization missing")
        check(abs(float(r.lower_head_velocity_m_s)) > 1.0e4, "20um lower-head speed too small")
        check(abs(float(r.upper_head_velocity_m_s)) > 1.0e4, "20um upper-head speed too small")
        check(float(r.max_conservation_residual) < 1.0e-6, "20um conservation residual too large")


def check_sp3_fork() -> None:
    term = pd.read_csv(ROOT / "results/stage3/single_seed/sp3_fork_off_resume1/termination.csv").iloc[0]
    check(term.termination_reason == "minimum_dt", "SP3-off fork did not terminate by minimum_dt")
    check(float(term.actual_end_time) < 2.0e-9, "SP3-off fork unexpectedly reached full end time")
    ledger = pd.read_csv(ROOT / "results/stage3/single_seed/sp3_fork_off_resume1/coupling_ledger.csv")
    check((ledger.integral_S_ph == 0.0).all(), "SP3-off fork has nonzero S_ph")

    summary = pd.read_csv(ROOT / "results/stage3/closure/run_summary.csv")
    on = summary[summary.run == "coarse_ml_recovery_resume1"].iloc[0]
    off = summary[summary.run == "sp3_fork_off_resume1"].iloc[0]
    check(bool(on.formed), "SP3-on baseline not formed")
    check(not bool(off.formed), "SP3-off fork formed despite being expected to weaken")
    check(float(on.photoionization_ahead) > 0.0, "SP3-on baseline has no photoionization ahead")
    check(float(off.photoionization_ahead) == 0.0, "SP3-off summary has nonzero photoionization ahead")


def check_partial_10um() -> None:
    d = ROOT / "results/stage3/single_seed/stage3_verification_baseline_10um_resume1"
    check((d / "fields_343.csv").exists(), "10um partial 0.4 ns field missing")
    check((d / "checkpoint.bin").exists(), "10um partial checkpoint missing")
    hist = pd.read_csv(d / "scalar_history.csv")
    check(float(hist.time_s.max()) >= 4.0e-10, "10um partial did not reach 0.4 ns")
    check(int(hist.rejected_steps.max()) == 0, "10um partial had rejected steps")
    check(float(hist.charge_source_residual.abs().max()) < 1.0e-6, "10um partial conservation residual too large")
    strategy = (ROOT / "docs/stage3_resource_strategy.md").read_text()
    check("must not be used to claim strict mesh convergence" in strategy, "resource strategy does not guard partial-grid claims")


def check_forbidden_outputs() -> None:
    for pattern in ("*case_i*", "*collision*", "*current_moment*", "*fft*", "*esd*", "*radiation*"):
        check(not any((ROOT / "results/stage3").rglob(pattern)), f"forbidden output present: {pattern}")


def main() -> None:
    check_registry()
    check_complete_20um()
    check_sp3_fork()
    check_partial_10um()
    check_forbidden_outputs()
    if failures:
        print("\n".join("FAIL: " + f for f in failures))
        raise SystemExit(1)
    print("Stage 3 resource-aware preparation validation passed with measured evidence.")


if __name__ == "__main__":
    main()
