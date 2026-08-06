# Stage 3 resource-aware execution strategy

Date: 2026-07-23

This note records a resource-control decision made during Stage 3 recovery.

## What remains unchanged

- The fluid equations are not changed.
- Morrow--Lowke transport and reaction coefficients are not changed.
- Poisson, SP3, ISG-0, adaptive timestep rules, and density update logic are not coarsened.
- The complete physical streamer evidence is taken from measured solver output, not from hard-coded postprocessing values.

## What is allowed to be coarsened

For preparation-oriented evidence and sensitivity scouting, the following are allowed:

- lower full-field output frequency;
- checkpoint/restart segmentation;
- shorter grid-probe windows;
- partial 10 um / 5 um startup probes clearly marked as partial;
- use of the complete 20 um run as the primary physical-baseline evidence.

## What is not allowed to be claimed from partial probes

Partial 10 um / 5 um runs must not be used to claim strict mesh convergence of the final double-headed streamer propagation stage unless they include a measured propagation window after head formation.

## Current resource observation

The 10 um resume run `S3-STAGE3_VERIFICATION_BASELINE_10UM_RESUME1` advanced from the 0.100 ns checkpoint to 0.400 ns and produced measured field snapshots at approximately 0.2009 ns and 0.4002 ns. The run was advancing normally but was interrupted intentionally because a full startup-to-2 ns 10 um run would be too expensive for the current preparation objective.

## Closure implication

The resource-aware route can support solver-preparation reporting and physical trend checks, but it is not the same as the original strict Stage 3 closure gate requiring 20/10/5 um propagation-stage convergence. If strict closure is required later, complete 10 um and 5 um propagation-window evidence must still be generated.
