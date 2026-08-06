# Stage 3 Recovery Plan

## Objective

Recover Stage 3 with measured evidence, without entering Stage 4 or using Liu-Pasko digitization/undisclosed Shi transport tables.

## Execution Protocol

The recovery follows the bounded diagnostic protocol requested for this stage:

1. Confirm Stage 2 anti-fraud closure.
2. Audit the early `coarse_ml` termination.
3. Fix long-time advancement, checkpointing, termination metadata, and coupling ledger output.
4. Replace all-domain peak head diagnostics with front-region diagnostics.
5. Rerun original `coarse_ml` before changing field or seed parameters.
6. If needed, run only the bounded six-case verification matrix.
7. Freeze the lowest-field passing verification baseline.
8. Run SP3 on/off fork, convergence, sensitivity, MPI, checkpoint, evidence audit, and closure validation.

## Bounded Diagnostics

Each blocking gate has at most three diagnostic rounds:

- Round 1: configuration, termination, and data processing.
- Round 2: coupling implementation, units, and signs.
- Round 3: independent control case or module isolation.

Repeated runs of the same configuration are capped at two, with the second allowed only for reproducibility.

## Current Recovery Hypotheses

1. The `coarse_ml` run stopped early because the executable exited when the old formation flag became true.
2. The old head diagnostic used all-domain field maxima rather than a connected streamer-front region.
3. The SP3 response must be evaluated from a shared checkpoint and with formation-aligned metrics, not only final position at one absolute time.
