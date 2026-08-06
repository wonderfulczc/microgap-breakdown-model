# Stage 4 Closure Report

Stage 4 target: resource-aware streamer collision and radiation workflow validation.

This stage used the verified Stage 3 Morrow--Lowke-based C++17/PETSc/MPI streamer solver. It did not use the unpublished Shi 2019 author transport table, did not perform a full 10/5 µm long-duration grid matrix, and does not claim one-to-one quantitative reproduction of Shi 2019 Case I.

## Boundary and baseline

Final collision baseline:

- run_id: `S4-PAPERLIKE-20UM`
- configuration: Shi-2019-Case-I-style double seed, not author-table quantitative Case I
- grid: 20 µm
- domain: 0 <= r <= 1.0 mm, 0 <= z <= 10.0 mm
- background field: 1.5 × 3.2e6 V/m
- seeds: z = 3.5 mm and 6.5 mm, n0 = 1e20 m^-3, sigma = 0.1 mm
- termination: reached_end_time at 3.0 ns
- accepted steps: 4148
- rejected steps: 0
- wall time: 5624.84 s

Backup collision cases were not run because the paper-like 20 µm case produced sufficient collision evidence under the resource-aware acceptance criteria.

## Propagation and collision evidence

Collision event source: `results/stage4/collision/collision_event.csv`.

- derived collision time: 1.953740091585568 ns
- bridge criterion: PASS
- bridge mean electron-density growth ratio: 188.5607144893954
- bridge threshold: 8.51987394478503e19 m^-3
- field-collapse criterion: PASS
- E_gap peak: 1.2025854587531555e7 V/m at 0.8818284835165764 ns
- E_gap minimum after peak: 4.0667791620277055e6 V/m
- field-drop fraction: 0.6618303395881607
- distance criterion: NOISY_NONCONTRADICTORY
- d_head diagnostic minimum: 1.1 mm

The C++ head-distance diagnostic is treated as noisy after channel interaction. Stage 4 closure relies on the two strong independent indicators: bridge formation and gap-field collapse. The distance diagnostic shows approach before losing clean head identity and is not used as a forced collision position.

## Isolated controls

Measured isolated controls:

- `S4-LEFT-ISOLATED`: reached 2.25 ns, 2487 accepted steps, 0 rejected steps
- `S4-RIGHT-ISOLATED`: reached 2.25 ns, 2487 accepted steps, 0 rejected steps

The right isolated run had a previous user pause at 0.052 ns with exit code 130. That audit-only fragment is preserved in `archive/audit_only/stage4_interrupted_runs.tar.zst`. The formal evidence run was restarted from t = 0 because the existing resume path would truncate per-step history files.

## Current moment and radiation chain

The C++ executable directly wrote current-moment data every accepted step:

- `results/stage4/current_moment/collision.csv`
- `results/stage4/current_moment/left_isolated.csv`
- `results/stage4/current_moment/right_isolated.csv`

Definitions:

- primary current moment: electron drift current, I_CM,drift
- sensitivity current moment: electron drift plus diffusion
- ion drift current: not included
- displacement current: not included

Measured peak current moments:

- collision drift peak: 5.5846079744866e-3 A m
- left isolated drift peak: 5.133989137883e-4 A m
- right isolated drift peak: 5.134887163777e-4 A m
- collision drift-plus-diffusion relative change: 5.87e-5

Collision increment:

- Delta I peak: 4.471137688342788e-4 A m
- Delta I SNR: 209.783698632248
- causality status: PASS

Derivative and spectrum:

- local-polynomial derivative peak: 2.957134248251927e6 A m s^-1
- uniform five-point derivative peak: 2.953442905802061e6 A m s^-1
- derivative-method relative difference: 0.0012482836895377
- conservative trusted frequency limit: 3.3367227397851086e11 Hz

Frequency-band status:

- VHF 30--300 MHz: trusted
- UHF 0.3--3 GHz: trusted
- SHF 3--30 GHz: trusted

Here `trusted` means sampling/Nyquist trusted for the postprocessed time series. It does not mean that the spatial grid, transport parameters, or radiation model have achieved high-frequency quantitative convergence.

## Resource-aware sensitivity and MPI

Output-sampling sensitivity was used as the resource-aware local sensitivity check:

- every-2-sample Delta I peak relative difference: 0
- every-4-sample Delta I peak relative difference: 0
- every-4-sample derivative peak relative difference: 0.0014088607811951

MPI short-window consistency:

- run_ids: `S4-MPI-PROBE-1RANK`, `S4-MPI-PROBE-2RANK`
- all tracked metrics: PASS
- maximum primary relative observed difference: 9.417414926589372e-12
- near-zero charge and residual quantities were judged with absolute tolerances.

Optional 10 µm event probe:

- status: SKIPPED_RESOURCE_LIMIT
- reason: not a Stage 4 hard gate under the resource-aware plan; 5 µm was explicitly excluded.

## Evidence and validation

Key artifacts:

- run registry: `results/stage4/provenance/run_registry.csv`
- derived artifact manifest: `results/stage4/provenance/derived_artifact_manifest.csv`
- validation matrix: `docs/stage4_validation_matrix.csv`
- figures: `results/stage4/figures/`

Validation:

- Stage 1 closure validation passed.
- Stage 2 closure validation passed with traceable measured evidence.
- Stage 3 resource-aware preparation validation passed with measured evidence.
- Stage 4 evidence audit passed with traceable measured evidence.
- Stage 4 resource-aware closure validation passed with traceable measured evidence.

## Limitations

The mutable `config/stage4/runs.yaml` evolved during Stage 4 as isolated and probe configurations were added. For auditability, the run registry preserves the executable hash, command line, stdout/stderr paths, result files, and result SHA256 values. The command line is the authoritative per-run configuration record for completed runs.

The radiation workflow is validated as a reproducible pipeline from measured C++ current moment to Delta I, derivative, FFT, ESD, and band energy. It is not a calibrated prediction against the Shi 2019 author implementation.

Stage 4 is complete under resource-aware workflow acceptance.

READY FOR STAGE 5 PREPARATION.
