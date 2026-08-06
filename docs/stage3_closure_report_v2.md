# Stage 3 Closure Report v2

## Conclusion

Stage 3 is not complete.

The resumed implementation corrected the earlier stage-boundary mistake: Liu-Pasko Figure 1 digitization is no longer treated as a Stage 3 blocker, and the implemented Stage 3 baseline uses the reproducible Morrow-Lowke analytic transport model. The run evidence now shows real coupled Stage 3 execution, but two core physics gates fail and the full convergence/sensitivity matrix was not completed after those failures were measured.

## Previous Attempt

The previous Stage 3 attempt stopped before implementation. The blocker was an overly strict Liu-Pasko Figure 1 digitization requirement. No Stage 3 solver existed at that time. The current implementation uses Morrow-Lowke as the reproducible Stage 3 baseline and does not claim Shi et al. (2019) Case I reproduction.

The retained 600 dpi Liu-Pasko Figure 1 rendering remains a Stage 4 preparatory artifact only.

## Implemented Scope

Implemented source includes Morrow-Lowke transport, reaction coupling, and a coupled axisymmetric three-species streamer executable using the Stage 2 numerical modules:

- `cpp/include/streamer_rf/streamer/MorrowLowke.hpp`
- `cpp/src/streamer/MorrowLowke.cpp`
- `cpp/include/streamer_rf/streamer/ReactionModel.hpp`
- `cpp/src/streamer/ReactionModel.cpp`
- `cpp/include/streamer_rf/streamer/StreamerSolver.hpp`
- `cpp/src/streamer/StreamerSolver.cpp`
- `cpp/apps/stage3_transport.cpp`
- `cpp/apps/stage3_run.cpp`
- `python/streamer_rf/streamer/morrow_lowke_reference.py`
- `python/streamer_rf/streamer/diagnostics.py`
- `python/stage3/run_stage3.py`
- `python/stage3/analyze_stage3.py`
- `tools/audit_stage3_evidence.py`
- `tools/validate_stage3_closure.py`

The code path reuses Stage 2 Poisson, SP3, and ISG-0 modules. No second Poisson or SP3 implementation was added.

## Transport Model

The Morrow-Lowke formulas are audited in `docs/stage3_formula_audit.md` from the locally retained Morrow and Lowke (1997) PDF, Appendix A. Stage 3 does not use a Liu-Pasko lookup table and does not claim the Shi et al. (2019) author transport table.

Measured transport outputs:

- Model breakdown field: `E_k,ML = 2633104.8564786837 V/m`
- Paper reference breakdown field: `3.2e6 V/m`
- Relative difference: about `17.72%`
- Maximum C++/Python coefficient relative difference: `4.048341942137943e-15`

## Coupled Solver Evidence

The measured coarse ML run `S3-COARSE_ML` produced a double-headed candidate:

- final analyzed time: `2.5481538299300147e-10 s`
- peak electric field: `6.785380362802061e6 V/m`
- channel electron density: `9.085903244614232e19 m^-3`
- lower head position: `0.00479 m`
- upper head position: `0.00515 m`
- photoionization ahead of positive head: `2.4428275949957653e26 m^-3 s^-1`
- maximum conservation residual: `3.845920966441989e-14`
- minimum time step: `1.349740976053229e-12 s`

Against the Shi 2017 reference quantities listed for Stage 3 triage, the measured coarse ML peak field is about `69.7%` below `2.24e7 V/m`, and the channel density is about `74.8%` below `3.6e20 m^-3`. Because the basic physics gates below fail, this is not treated as a mere transport-model discrepancy.

The n_ref scan completed for `1e-14`, `1e-12`, `1e-10`, `1e-8`, and `1e-6`. All five points stayed within the 1 percent freeze threshold for the required scalar metrics.

MPI 1/2/4 rank runs completed and produced matching scalar outputs at roundoff scale.

## Blocking Failures

The Stage 3 physical closure gates failed on measured outputs:

- Field-head and charge-head definitions do not agree within two grid cells. For `coarse_ml`, the analyzed final field/charge head offsets are 3-5 cells, so `head_definitions_agree=False`.
- Closing SP3 did not measurably weaken the positive head in the current coarse result. `coarse_no_sp3` has `S_ph=0`, but its upper head position remains `0.00515 m`, the same as `coarse_ml`.

Because these are core physics gates, the baseline/fine grid runs, timestep sensitivity, OpenCharge sensitivity, and domain sensitivity were not completed in this resumed pass. Running the remaining matrix would not close Stage 3 while these two gates remain failed.

## Validation Status

Passed:

- Stage 2 closure regression
- C++ Stage 3 unit tests: 20 checks passed
- Python Stage 3 tests: 7 tests passed
- Morrow-Lowke C++/Python comparison
- model breakdown field calculation
- n_ref freeze scan
- MPI 1/2/4 evidence
- evidence audit with zero static findings and zero bad registered hashes

Failed or incomplete:

- field/charge head agreement gate
- SP3-off physical-response gate
- full 20/10/5 um mesh convergence
- timestep sensitivity
- OpenCharge sensitivity
- domain sensitivity
- final Stage 3 closure validator

## Boundary To Stage 4

Stage 4 is not ready. Before Stage 4, one of the following remains required:

- reliable digitization of the Liu-Pasko solid transport curves, or
- an independent Boltzmann-generated air transport table with graphical consistency checking against Liu-Pasko Figure 1.

The Morrow-Lowke Stage 3 results must not be used as Shi et al. (2019) target reproduction results.

## Scope Guard

Shi et al. (2019) Case I was not executed. No two-seed streamer collision was simulated. No collision current moment, FFT, ESD, VHF/UHF radiation, or Stage 4 result was generated.
