# Stage 2: Verified Numerical Building Blocks

Stage 2 starts only after `python tools/validate_stage1_closure.py` passes. It contains exactly three work packages and does not start Shi 2019 Case I.

## WP2.1 Electron transport with ISG-0

Freeze the dimensional reference density, implement ordinary SG and ISG-0 in isolation, verify small-alpha and zero-width limits, verify fallback behavior, and pass the Kulikovsky benchmarks.

## WP2.2 Axisymmetric Poisson solver

Implement an axisymmetric solver with explicit boundary conditions, manufactured-solution convergence, open-boundary validation, and documented tolerances.

## WP2.3 Three-group SP3 photoionization

Implement the three-group SP3 equations with correct SP3—not Helmholtz—coefficients, explicit open-boundary treatment, tolerance studies, and the Bourdon/Liu benchmark.

## Unified Stage 2 exit gate

Stage 2 closes only when:

1. ISG-0 passes Kulikovsky validation.
2. Poisson passes manufactured-solution and open-boundary validation.
3. SP3 passes the Bourdon/Liu benchmark.
4. Interfaces and SI units are consistent across all three modules.
5. All automated tests pass.
6. The modules have not yet been coupled into a complete streamer solver.

Only after this gate may Stage 3 begin. This document is a plan; no Stage 2 implementation is included.
