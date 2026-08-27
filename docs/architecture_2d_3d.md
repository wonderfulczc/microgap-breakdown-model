# 2D/3D Solver Architecture Boundary

Date measured: 2026-08-27

This Step 1 record freezes the engineering responsibility boundary between the
existing two-dimensional solver and a future Afivo-streamer three-dimensional
backend. It defines interfaces only; no adapter is implemented in Step 1.

## Frozen 2D backend

Backend: C++17 + PETSc + MPI + CMake/Ninja

Role: primary scientific solver

Responsibilities:

- axisymmetric micro-gap discharge
- axisymmetric electrodes, including needle and needle-plane style cases
- parameter sweep
- transport / chemistry sensitivity
- streamer inception
- streamer propagation
- streamer bridging / collision
- discharge current
- current moment
- basic native RF diagnostics

Main project service:

- Subdirection 1: electrode structure -> discharge process -> current/channel
- Subdirection 2: native discharge source and basic broadband mechanism

The existing Stage 1-5 two-dimensional workflow is the frozen reference. Future
project-specific model development must derive from it without overwriting the
reference transport, chemistry, numerical-method, closure, or validation record.

```text
frozen 2D reference
    ↓
future project-specific derived branch
    ↓
micro-gap research model
```

## Future 3D backend

Backend: Afivo-streamer

Role: selected high-value 3D validation backend

Expected technology: Fortran + Afivo AMR + multigrid + OpenMP

Future responsibilities:

- non-axisymmetric electrodes
- triangular copper-foil electrodes
- electrode misalignment
- asymmetric discharge paths
- multi-channel competition when scientifically required
- asymmetric streamer interaction
- full Jx / Jy / Jz distribution
- radiation direction / polarization source data

The 3D backend is not assigned to routine large-scale parameter scans.

## Solver relationship

Afivo-streamer does not replace the 2D solver.

```text
2D screening
    ↓
critical cases selected
    ↓
3D validation
    ↓
evaluate validity range of 2D
```

The solvers remain independent. They are connected through common physical
definitions and postprocessing contracts, not by source-level fusion.

```text
2D raw ─┐
        ├─ common physical/output contract ─ RF postprocessing
3D raw ─┘
```

## Common physical interface

The common physical interface defines comparable simulation meaning across
solvers. Step 1 records the contract only and does not implement adapters.

Required shared definitions:

- geometry definition
- gas composition
- pressure
- temperature
- voltage waveform
- initial condition / seed / inception definition
- transport data
- chemistry
- photoionization
- physical boundary conditions
- output timing
- diagnostic definitions

## Common output logic

Required logical variables:

- time
- geometry_id
- gas
- pressure
- temperature
- voltage
- phi
- E
- rho
- ne
- species
- J
- I
- current_moment
- mesh metadata

Raw solver file formats are not forced to match. The existing 2D workflow may
continue to use HDF5 / CSV / JSON. Future Afivo-streamer workflows may keep
native Silo / VTK or other Afivo-supported outputs. Unification happens at the
physical-variable and RF-postprocessing interface layer.
