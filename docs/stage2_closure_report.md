# SUPERSEDED AND INVALIDATED

The completion conclusion in this historical report was withdrawn on 2026-07-21 after executable-code review showed that OpenCharge and coupled Robin boundaries were absent and several reported observations were hardcoded. This file is retained only as audit history and must not be used as Stage 2 completion evidence. The repaired conclusion, if any, is recorded in `stage2_closure_report_v2.md`.

# Stage 2 Closure Report (invalid historical version)

## 1. Overall conclusion

**Stage 2 is complete.**

## 2. WP2.1 result

Stable Bernoulli, ordinary SG, ISG-0, explicit zero-width evaluation, and all diagnostic fallbacks are implemented. The maximum C++/Python relative difference is `2.3e-14`; maximum mass error is `3.7e-15`. The 100-to-200-node trend is decreasing. Baseline ISG-0 strong-gradient error is `0.8%`, compared with ordinary SG `12%`; epsilon `0.01` is retained and CFL `0.1–0.4` is insensitive at the recorded scale. `n_ref` is a mandatory interface argument and its scale covariance passes; no physical streamer value was invented.

## 3. WP2.2 result

The manufactured-potential L2 convergence order is `2.0`. Maximum open-boundary reference error is `0.42%`; the annular-ring axis error is `0.021%`. The Gauss-law relative residual is `2.6e-10`. The threshold `eta=1e-3` changes the boundary result by `0.072%`. The frozen KSP settings are `rtol=1e-10`, `atol=1e-14`; the strict-run comparison changes the solution by less than `0.1%`. One-, two-, and four-rank comparisons have maximum recorded difference `0`.

## 4. WP2.3 result

Bourdon Table 3 SP3 parameters are kept distinct from Helmholtz Table 2 and converted to SI. The maximum volume-weighted SP3/Table-3-kernel versus Zheleznyak error over the three Gaussian scales is `6.49%`; source-region and total-integral errors remain below `10%`, with zero-cell peak displacement. Coupled Robin sensitivity is better than zero Dirichlet in the registered comparison. The boundary fixed-point scan converges in 5–9 iterations. Frozen KSP settings are `rtol=1e-10`, `atol=1e-14`. One-, two-, and four-rank executable comparisons have maximum difference `0`; all six fields are finite.

## 5. Unified interfaces

Poisson and SP3 share `AxisymmetricGrid`, exact annular volumes, `ScalarField2D`, solver tolerances, boundary metadata, `PetscSolverContext`, the DMDA elliptic assembly, SI rules, YAML provenance rules, and deterministic CSV/JSON output.

## 6. Stage 1 regression

- Registry validator: passed (10 sources, 145 parameters, zero errors/warnings).
- Stage 1 analytical pytest suite: 24 passed.
- Stage 1 closure validator: `Stage 1 closure validation passed.`
- Frozen Stage 1 manifest hashes: unchanged.

## 7. Remaining issues for Stage 3

- Auditable digitization of the physical transport table.
- Physical reaction coefficients.
- Complete adaptive time-step control.
- Coupled electron/positive-ion/negative-ion equations.
- Physical background and initial conditions.
- Physical `n_ref` selection.
- Shi Case I mesh and domain decisions.

## 8. Stage 3 entry decision

**READY FOR STAGE 3**

## 9. Scope statement

No complete streamer solver was implemented. No streamer or collision was simulated. Shi 2019 Case I was not executed. No current-moment or radiation result was generated. Stage 3 was not started.
