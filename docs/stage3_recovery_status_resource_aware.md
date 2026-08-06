# Stage 3 recovery status: resource-aware preparation

Date: 2026-07-23

This report records the current Stage 3 recovery state after adopting a resource-aware execution strategy. It does not claim strict Stage 3 closure and does not write `READY FOR STAGE 4`.

## Completed measured evidence

- The previous 0.255 ns early termination was not a numerical failure. The recovered 20 um ML run resumed from checkpoint and reached 2.0 ns.
- `coarse_ml_recovery_resume1` reached `termination_reason=reached_end_time` with 2018 accepted continuation steps and 0 rejected steps.
- Independent connected-front postprocessing identifies double-headed propagation in `coarse_ml_recovery_resume1`.
- The measured connected-front speeds are approximately `-4.09e5 m/s` for the lower head and `+3.93e5 m/s` for the upper head.
- The measured 20 um peak field is `1.2454168416886907e7 V/m`.
- The measured 20 um channel density is `8.25340207015659e19 m^-3`.
- The measured positive-head photoionization-ahead value is `1.925034752890143e28 m^-3 s^-1`.

## SP3 response

The SP3-off fork was started from the same 20 um checkpoint used by the SP3-on recovery continuation.

- SP3-on reached 2.0 ns and formed a sustained double-headed streamer.
- SP3-off terminated by `minimum_dt` at approximately `1.34 ns`.
- SP3-off had `integral_S_ph=0` in the coupling ledger.
- SP3-off did not satisfy the connected-front formation criterion.

This satisfies the preparation-level SP3 response check without relying on absolute final-head-position comparison.

## 10 um resource observation

The 10 um run was restarted from the 0.100 ns checkpoint and advanced to approximately 0.400 ns. Full-field output was reduced to 0.2 ns intervals. The run was advancing normally, but a full startup-to-2 ns 10 um simulation remained too expensive for the current preparation objective.

The 10 um partial run is therefore retained as early grid-trend evidence only. It must not be used as strict propagation-stage mesh convergence evidence.

## Current validation state

The following pass:

- Stage 3 C++ checks: 20 checks passed.
- Stage 3 Python tests: 7 tests passed.
- Stage 3 evidence audit: 0 findings and 0 bad registered hashes.
- Resource-aware preparation validator: `tools/validate_stage3_resource_preparation.py`.

The strict closure validator is not expected to pass yet because it still requires complete `20/10/5 um` propagation-stage convergence evidence.

## Explicit exclusions

No Shi 2019 Case I run was executed. No two-seed collision was simulated. No collision current moment, FFT, ESD, or radiation result was generated.
