# Reproducibility and run profiles

## Environment

The frozen reference environment uses Python 3.14, CMake/Ninja, a C++17
compiler, PETSc through `pkg-config`, and Open MPI. Runtime Python versions are
locked in `requirements.txt`; development/test dependencies are in
`requirements-dev.txt`.

Optional external workflows use:

- `AFIVO_STREAMER_ROOT`: pinned external Afivo-streamer checkout.
- `OPENEMS_ROOT`: pinned external openEMS installation prefix.
- `OPENEMS_PROJECT_ROOT`: pinned external openEMS-Project source checkout.
- `OPENEMS_DEPS_ROOT`: optional isolated openEMS dependency prefix.
- `OPENEMS_PYTHON`: Python interpreter for the isolated openEMS environment.
- `PETSC_DIR` and optional `PETSC_ARCH`: non-system PETSc discovery.
- the host MPI launcher; MPI execution is not hidden inside the smoke script.

## Resource profiles

### SMOKE

Runs `scripts/reproduce_smoke.sh`. It configures/builds the C++ project, runs
CTest, imports the Python package, and executes selected Stage-F/G/H/I tests.
It requires neither external backend nor large result data.

### REFERENCE

Runs the full default pytest suite and compact committed generation scripts
where explicitly needed. Frozen small Stage-G/H/I CSV/JSON results are reused;
openEMS and Afivo are not rerun by default.

### FULL_RESEARCH

May regenerate multi-GB Afivo and Stage-F data or execute openEMS. It requires
external backend checkouts, machine-specific environment configuration, and
the configurations referenced by `packaging/large_data_manifest.json`.

## Clean build

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
cmake -S . -B build -G Ninja
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=python .venv/bin/python -m pytest -q
```

## Real experiment reentry

Supply files matching `validation/stage_i/` contracts, then perform:

1. real-data ingestion;
2. provenance gate;
3. quality gate;
4. calibration/reference-plane check;
5. Stage-I comparison;
6. discrepancy-ledger update;
7. validation-status update.

The dry-run pipelines need not be rerun. Scientific status may change only
after real evidence passes these gates.
