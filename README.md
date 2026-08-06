# Microgap Breakdown Model

This repository contains a C++/PETSc/MPI streamer-fluid simulation toolchain prepared for future microgap-breakdown modeling. It was built through a resource-aware reconstruction of the Shi 2019 streamer-collision and broadband-radiation workflow, then reduced for GitHub archival as a reusable solver and workflow foundation.

The repository contains:

- C++17 streamer-fluid solver components.
- PETSc/MPI/CMake build integration.
- Python configuration, run orchestration, postprocessing and plotting tools.
- Unit tests and validation scripts.
- Documentation of the Stage 1–5 reconstruction process.
- Project-transfer notes for adapting the workflow to microgap electrode simulations.

Large historical simulation outputs have been removed. The project should therefore be described as a solver/toolchain/workflow archive, not as a complete raw-result evidence package.

## Purpose

The immediate purpose is to preserve a working numerical framework that can be adapted from simplified free-space streamer cases toward a real microgap discharge model:

实际微间隙电极几何
→ 电极电势边界和实际电压波形
→ 电极附近种子电子/表面发射
→ 微间隙流注与击穿
→ 电极回路时变电流
→ 电流矩与原生宽频辐射
→ 外部 RLC 耦合
→ 接收端信号

Current code covers only part of this chain: streamer-fluid transport, Poisson coupling, SP3 photoionization, ISG-0 electron flux, checkpoint-style run logic, current-moment integration, FFT/ESD postprocessing, and validation infrastructure.

## Recommended platform

Use Linux or Windows + WSL2 Ubuntu. Native Windows builds are not recommended because PETSc, MPI, `pkg-config`, and `mpirun` are used directly by the CMake and Python workflows.

## Dependencies

System packages on Ubuntu/WSL2:

```bash
sudo apt update
sudo apt install -y build-essential cmake ninja-build pkg-config \
  openmpi-bin libopenmpi-dev libpetsc-real-dev python3-venv
```

Python packages:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Build

Configure and build:

```bash
cmake -S . -B build -G Ninja
cmake --build build --parallel
```

Run C++ unit tests:

```bash
ctest --test-dir build --output-on-failure
```

Run a Python import smoke test:

```bash
PYTHONPATH=python PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -c \
  "import streamer_rf.stage4, streamer_rf.stage5; print('python smoke passed')"
```

## How to use the toolchain

The repository is organized around reproducible stages:

- `cpp/`: C++ solver, PETSc/MPI numerical modules, C++ tests.
- `python/`: run drivers, analysis tools, current-moment and radiation postprocessing.
- `config/`: stage and run configuration templates.
- `tests/`: Python regression tests.
- `tools/`: validators, evidence audit tools, cleanup and summary validators.
- `docs/`: method notes, closure reports, project-transfer guidance, and final summaries.

Typical workflow:

1. Edit or add a configuration under `config/`.
2. Build the C++ executables with CMake/Ninja.
3. Run the appropriate Python driver, for example `python/stage4/run_stage4.py` or `python/stage5/finalize_resource_run.py`.
4. Analyze generated outputs with the corresponding Python analysis script.
5. Use the validator tools as consistency checks when result evidence is present.

New simulations will recreate a local `results/` directory. This directory is intentionally ignored by Git.

## Current scientific status

The reconstruction completed the solver, collision, radiation-chain, trend-logic, and research-transfer preparation workflow under resource-aware criteria. It did not complete a 1:1 quantitative reproduction of Shi 2019.

Historical Stage closure validators require the deleted `results/` tree and are retained as reference tools. They will pass again only after the corresponding result evidence is regenerated or restored.

## Limitations

This repository is not yet a complete microgap-breakdown model. The following pieces still need to be implemented or replaced:

- Real electrode geometry.
- Electrode potential boundary conditions.
- Actual voltage waveform input.
- Seed electrons near electrodes or surface-emission models.
- Electrode-circuit time-domain current.
- External RLC coupling.
- Receiver/antenna signal model.
- Calibrated gas transport tables for the target gas and pressure.

It also does not contain the original author transport table from Shi 2019 and does not claim point-by-point reproduction of Shi 2019 figures.

## Documentation

- `docs/github_archive_assessment.md`: archive scope and consequences of result deletion.
- `docs/project_simulation_transfer.md`: how to adapt the solver workflow to a microgap project.
- `docs/simulation_toolchain_summary_for_project_book.md`: Chinese project-book style summary.
- `docs/simulation_reconstruction_summary/index.html`: offline visual summary.

## Minimal validation already checked before upload

```bash
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=python PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -c \
  "import streamer_rf.stage4, streamer_rf.stage5; print('python smoke passed')"
```
