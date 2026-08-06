# Stage 4 resource plan

Date: 2026-07-23

Stage 4 validates the collision and radiation workflow under resource-aware acceptance criteria. It does not claim a one-to-one Shi 2019 quantitative reproduction.

## Fixed boundaries

- Grid baseline: 20 um.
- Full 10 um and 5 um long runs are not Stage 4 closure gates.
- 5 um is not executed.
- Full-field output is event-driven and sparse.
- Current moment is integrated directly by C++ at every accepted step.
- Radiation spectra are derived from measured current-moment output.

## Bounded collision search

The collision run budget is:

1. `S4-PAPERLIKE-20UM`: 4.8 MV/m, seed spacing 3.0 mm.
2. `S4-BACKUP-FIELD-1P75-20UM`: 5.6 MV/m, seed spacing 3.0 mm, only if paper-like does not collide.
3. `S4-BACKUP-D0-2P5-20UM`: 4.8 MV/m, seed spacing 2.5 mm, only if the first two do not collide.

No additional parameter search is allowed in Stage 4.

## Controls

After a collision baseline is frozen, exactly two isolated controls are run:

- `S4-LEFT-ISOLATED`
- `S4-RIGHT-ISOLATED`

They use the same global time zero, field, domain, grid, timestep logic, and output logic as the collision baseline.

## Validation scope

Stage 4 can close if the workflow is measured end-to-end:

- collision or collision candidate detected by independent metrics;
- isolated controls are real runs;
- drift and electron current moments are directly output;
- collision increment, derivative, FFT, ESD, and band-energy files are generated;
- time-step and MPI local-window checks are measured;
- evidence audit and closure validator pass.

Strict mesh-independent quantitative agreement remains outside this stage.
