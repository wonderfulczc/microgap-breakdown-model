# Stage 5 Closure Report

Stage 5 target: resource-aware trend, logic and project-transfer validation.

This stage does not attempt Shi 2019 Figure 4b pointwise reproduction and does not claim Shi 2019 1:1 quantitative reproduction.

## 1. Resource boundary

Stage 5 reused the Stage 4 baseline `S4-PAPERLIKE-20UM` without rerunning it. New formal runs were limited to:

- collision runs: high field, larger gap, asymmetric seed;
- isolated runs: high-field left and right only;
- short MPI probes: high-field 1-rank and 2-rank windows;
- no 10 um or 5 um runs.

## 2. Stage 4 baseline

Baseline reused from Stage 4:

- run_id: `S4-PAPERLIKE-20UM`
- grid: 20 um
- E0: `1.5 x 3.2e6 V/m`
- d0: 3.0 mm
- t_collision: `1.953740091585568 ns`
- E_gap drop: `66.183%`
- bridge growth ratio: `188.56`
- Delta I peak: `4.471137688342e-4 A m`

## 3. High background field case

Final high-field collision run:

- run_id: `S5-HIGHFIELD-COLLISION`
- E0: `2.0 x 3.2e6 V/m`
- termination: `resource_stop_after_collision_window`
- actual end time: `1.3745538374968044 ns`
- accepted steps: `2124`
- rejected steps: `0`

The run was interrupted after the collision post-window was already satisfied. The registry preserves exit code 130 and the resource-aware termination file records the measured event window; this is not treated as a numerical failure.

Measured high-field event:

- t_collision: `0.9999612456531687 ns`
- E_gap peak: `1.2759101049521811e7 V/m`
- E_gap drop fraction: `0.6187467148400569`
- bridge growth: `427.44075407626406`

Compared with the Stage 4 baseline, the high-field collision occurred earlier, as expected for stronger background field.

## 4. High-field isolated controls

High-field isolated controls were independently run:

- `S5-HIGHFIELD-LEFT-ISOLATED`: reached `1.25 ns`, `1190` accepted steps, exit code 0.
- `S5-HIGHFIELD-RIGHT-ISOLATED`: reached `1.25 ns`, `1190` accepted steps, exit code 0.

The right isolated run had a paused/corrupted partial attempt preserved only as audit material in `archive/audit_only/stage5_interrupted_runs.tar.zst`; the formal result was rerun and registered.

## 5. High-field complete Delta I and radiation chain

Strict high-field subtraction used the same time coordinate:

`Delta I = I_collision - I_left - I_right`

Measured high-field strict products:

- Delta I peak: `0.0066195567164797535 A m`
- Delta I SNR: `3207.65`
- local polynomial derivative peak: `9.63148956131095e10 A m/s`
- pulse FWHM: `1.3599879857055805e-10 s`
- spectral centroid: `2.0390408236225998e9 Hz`
- trusted frequency limit: `2.125303251497472e11 Hz`

Frequency-band energy:

- VHF: `2.867759143417004e-11 J`
- UHF: `1.0231629487904344e-10 J`
- SHF: `3.065507759666655e-11 J`

The band trusted labels mean sampling/Nyquist trusted only. They do not imply spatial-grid, transport-model, or radiation-model high-frequency convergence.

## 6. Larger-gap proxy case

Run:

- run_id: `S5-LARGERGAP-COLLISION`
- d0: 5.0 mm
- termination: `NO_COLLISION_WITHIN_RESOURCE_WINDOW`
- actual end time: `1.0971808960249658 ns`
- accepted steps: `829`

Observed proxy trend:

- d_head remained about millimeter-scale within the resource window.
- E_gap peak diagnostic: `1.2044451415238287e7 V/m`
- E_gap drop fraction diagnostic: `0.14984390096407574`
- bridge growth diagnostic: `15.047375979165714`
- event-local current proxy peak: `9.722135028260161e-05 A m`

This is a valid long-distance trend proxy, not a strict collision result.

## 7. Asymmetric seed proxy case

Run:

- run_id: `S5-ASYMMETRIC-COLLISION`
- sigma1: 0.1 mm
- sigma2: 0.5 mm
- termination: `NO_COLLISION_WITHIN_RESOURCE_WINDOW`
- actual end time: `0.9615816760178067 ns`
- accepted steps: `760`

Observed proxy trend:

- the wide seed produced a strong early bridge tail;
- head-distance tracking became ambiguous for the broad seed;
- bridge growth diagnostic: `46.51037109409177`;
- E_gap drop diagnostic: `0.691885116016763`;
- event-local current proxy peak: `0.0045526100888156 A m`;
- velocity-ratio diagnostic: `0.19485550642198338`.

Because the online head diagnostic did not provide a stable collision event, Stage 5 does not claim a strict asymmetric collision-position offset.

## 8. Cross-case trend comparison

Cross-case output is stored in:

- `results/stage5/trends/case_comparison.csv`

The table explicitly separates:

- `strict_delta`: Stage 4 baseline and Stage 5 high-field complete radiation chain;
- `event_local_proxy`: larger-gap and asymmetric resource-window trend proxies.

Proxy values are not ranked as strict Delta I values.

## 9. Output-sampling sensitivity

High-field every-2/every-4 output extraction was compared against every-step output using a common sampling/Nyquist-trusted low-frequency centroid and a smoothed derivative peak.

Maximum measured differences:

- Delta I peak: `0%`
- smoothed derivative peak: `2.73%`
- trusted low-band spectral centroid: `0.67%`

All pass the 5% Stage 5 resource-aware target.

## 10. MPI

High-field local window MPI comparison:

- 1 rank vs 2 ranks
- maximum relative difference for non-near-zero primary scalars: `1.88415e-11`
- total charge absolute difference: `3.31614e-25 C`

All MPI checks pass. The total-charge relative value is not used alone because the denominator is near zero.

## 11. Evidence audit

Stage 5 evidence chain:

- source and executable hashes in `results/stage5/provenance/run_registry.csv`;
- raw run outputs under `results/stage5/runs/`;
- derived artifact SHA values in `results/stage5/provenance/derived_artifact_manifest.csv`;
- validation matrix in `docs/stage5_validation_matrix.csv`;
- audit tool: `tools/audit_stage5_evidence.py`;
- closure validator: `tools/validate_stage5_closure.py`.

## 12. Project transfer

The project-transfer assessment is documented in:

- `docs/project_simulation_transfer.md`

Reusable parts: PETSc/MPI solver framework, grid/field data structures, checkpointing, adaptive stepping, current moment integration, radiation postprocessing, and evidence audit workflow.

Required replacements for doctoral micro-gap simulation: electrode geometry, electrode Dirichlet boundaries, realistic seed/emission model, nonuniform applied waveform, calibrated gas transport, and external-circuit coupling.

## 13. Final status

Stage 5 is complete under resource-aware trend and project-transfer acceptance.

The Shi 2019 workflow reconstruction project is complete for solver, collision, radiation-chain, trend-logic, and research-transfer preparation.

This does not mean Shi 2019 Figures 1-4 have been quantitatively reproduced.
