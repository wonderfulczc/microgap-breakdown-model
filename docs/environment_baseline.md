# Step 1 Environment Baseline

Date measured: 2026-08-27

All entries below are measured or detected in the current WSL environment. No
packages, compilers, PETSc, MPI, Python, or Afivo-streamer components were
installed or upgraded during Step 1.

## Platform

- Project root: `/home/helianthusczc/projects/streamer-rf-replica`
- Filesystem note: project is not under `/mnt/c` or `/mnt/d`
- Kernel/WSL: `Linux Chelianthus 6.18.33.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 18 21:54:43 UTC 2026 x86_64 GNU/Linux`
- Ubuntu: `Ubuntu 26.04 LTS`

## C++/MPI/PETSc toolchain

- gcc: `gcc (Ubuntu 15.2.0-16ubuntu1) 15.2.0`, path `/usr/bin/gcc`
- g++: `g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0`, path `/usr/bin/g++`
- CMake: `cmake version 4.2.3`
- Ninja: `1.13.2`
- PETSc: `3.24.4`
- PETSc prefix: `/usr/lib/petscdir/petsc3.24/x86_64-linux-gnu-real`
- mpirun: Open MPI `5.0.10`, path `/usr/bin/mpirun`
- MPI note: `mpirun --version` emitted `pmix_ifinit: socket() failed with errno=1` before printing the Open MPI version.

PETSc remains the detected project PETSc and was not upgraded.

## Python environment

- System Python: `Python 3.14.4`, path `/usr/bin/python3`
- Virtual environment Python: `.venv/bin/python`
- Virtual environment prefix: `/home/helianthusczc/projects/streamer-rf-replica/.venv`
- Shell `VIRTUAL_ENV`: not set during measurement

Detected Python packages in `.venv`:

- numpy `2.5.1`
- scipy `1.18.0`
- pandas `3.0.3`
- matplotlib `3.11.1`
- PyYAML `6.0.3`
- pytest `9.1.1`
- h5py `3.16.0`

HDF5:

- Python h5py HDF5 runtime: `2.0.0`
- `h5cc`: not found
- `pkg-config hdf5`: not found

This is not a Step 1 blocker because the retained project workflow uses Python
`h5py` successfully and the current C++ build does not require system HDF5.

## Future Afivo-streamer prerequisites

- git: `git version 2.53.0`
- gfortran: `GNU Fortran (Ubuntu 15.2.0-16ubuntu1) 15.2.0`, path `/usr/bin/gfortran`
- make: `GNU Make 4.4.1`, path `/usr/bin/make`
- Fortran + OpenMP smoke test: PASS
- OpenMP smoke output: `openmp_threads 16`

Afivo basic compatibility: YES

Afivo-streamer was not cloned, installed, configured, or compiled in Step 1.

## Codex

- Codex CLI: `codex-cli 0.147.0`
- Detection note: command printed a warning that PATH aliases could not be created because the target filesystem was read-only.

## Git baseline note

The original local `.git` metadata was unavailable in the working tree at the
start of Step 1. The user identified the remote repository as
`https://github.com/wonderfulczc/microgap-breakdown-model.git`; read-only
remote probing found `refs/heads/main` at
`f752aa290b63e4c9eecbc7a447777b3bea14298a`.

The remote Git history was restored non-destructively before the Step 1 commit.
Historical local Git metadata before this recovery was unavailable in the
working tree, so no older local-only commits were fabricated.
