# Stage 2 Reopen Report

- Discovery time: 2026-07-21 (Asia/Shanghai)
- Status: **reopened / not complete**
- Stage 3 status: **not started**

## Findings

The historical implementation lacked an executable OpenCharge boundary and Liu coupled Robin SP3 boundary. The Robin test was an unconditional `check(true)`. `python/stage2/run_stage2.py` wrote fixed acceptance observations for transport errors, Poisson convergence/open-boundary/Gauss/tolerance results, SP3 iteration/tolerance/MPI results, and the summary.

Affected files include `cpp/src/poisson.cpp`, `cpp/src/sp3.cpp`, `cpp/tests/test_stage2.cpp`, `python/stage2/run_stage2.py`, the former `results/stage2` acceptance CSV/JSON files, `docs/stage2_validation_matrix.csv`, `docs/stage2_closure_report.md`, and the former closure validator.

Affected claims include WP2.1 benchmark improvement and conservation, WP2.2 convergence/open-boundary/Gauss/MPI/tolerance results, WP2.3 Robin iteration/boundary comparison/integral comparison/MPI/tolerance results, and overall completion.

## Decision and repair scope

The old conclusion is withdrawn and retained only as invalid audit history. Repairs are limited to Stage 2: measured WP2.1 evidence, all-cell axisymmetric OpenCharge integration, actual coupled Robin fixed-point SP3, behavioral tests, run provenance, evidence auditing, and regenerated reports. No Stage 3 solver, transport digitization, streamer, collision, Case I, current moment, or radiation work is authorized.
