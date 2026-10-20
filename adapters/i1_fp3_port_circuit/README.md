# I1 FP-3 synthetic port engineering

No modules here import native Stage C, Poisson, G0/G1/G2/G3, COMSOL, TENG or RF.
Dependencies are the installed NumPy/SciPy and pytest; no additional installation.
The historical parallel-Cext circuit is read as evidence but never executed or changed.

- `ports.py`: frozen descriptors, enum capacitance/charge/energy scopes, immutable
  accepted state, pure endpoint trials, three synthetic R/C providers and commit guards.
- `circuit.py`: series L-Rs-Cext-port implicit midpoint, actual port trial residual,
  bounded 2-by-2 linear Newton solve and joint accepted-state commit.
- `guards.py`: restricted internal C1 output and explicit FULL/TENG unavailability.
- `test_fp3.py`: isolated elements, signs, ownership, transactions and refusals.
- `validate.py`: the fixed 20 accounting-mode/step-size realizations; independent
  augmented continuous-ODE expm references, full step ledgers and order diagnostics.
- `archive.py`: auditable transaction examples, schema/example, final gate assembly,
  original-file protection and resource/source index; never executes a native solver.
- `execution_plan.json`: frozen input matrix, normalization scales and tolerances.

Actual results: [completed engineering report](results_20261015_final/FP3_PORT_CIRCUIT_ENGINEERING_REPORT.md),
[gate](results_20261015_final/fp3_engineering_gate_summary.json),
[initial clock incident](results_20261015/EXECUTION_INCIDENT.md).
The failed output and failed/successful isolated JUnit files are retained.

Reproduction, **only with explicit new execution authorization**, must use unused
output paths (all result writes are exclusive); do not rerun against these files:

```bash
.venv/bin/python -m pytest adapters/i1_fp3_port_circuit/test_fp3.py -q -o cache_dir=/tmp/fp3-new-cache --basetemp=/tmp/fp3-new-temp --junitxml=/tmp/fp3-new-junit.xml
.venv/bin/python -m adapters.i1_fp3_port_circuit.validate --output /tmp/fp3-new-results
```

All R/C values are synthetic test inputs, not inverse G_R or COMSOL data.
`TrialResult.port_work_J` is total discrete midpoint work;
`component_work_J` excludes the explicit C_p branch in SPLIT mode. Component
and external C_p endpoint energies are separate; add them exactly once.
Trial slopes have units F (Deltaq derivative with respect to endpoint V), while
the outer current residual has units A. Reference-plane signs are differential
H-to-R and never swap conductor identities under polarity reversal.

Port commit requires a validated result from the current token and explicit
convergence acceptance. Coupled commit independently rechecks the current residual
before either accepted state is assigned. Timing/trial-count diagnostics can change
without altering any accepted scientific state. This is an in-process synthetic
transaction, not a durable database transaction or native rollback proof.

See [FULL_PROVIDER_HANDOFF.md](FULL_PROVIDER_HANDOFF.md). Stop after archival;
no automatic physical FULL integration, TENG, cold-hot or RF task follows.
