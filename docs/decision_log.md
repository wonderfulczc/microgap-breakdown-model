# Decision Log

## DEC-0001

- decision_id: DEC-0001
- date: 2026-07-21
- decision: 2D axisymmetric quasi-electrostatic main model
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0002

- decision_id: DEC-0002
- date: 2026-07-21
- decision: C++17 main language
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0003

- decision_id: DEC-0003
- date: 2026-07-21
- decision: PETSc and MPI numerical libraries
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0004

- decision_id: DEC-0004
- date: 2026-07-21
- decision: Python post-processing
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0005

- decision_id: DEC-0005
- date: 2026-07-21
- decision: three-group SP3 photoionization
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0006

- decision_id: DEC-0006
- date: 2026-07-21
- decision: Kulikovsky ISG-0 transport
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0007

- decision_id: DEC-0007
- date: 2026-07-21
- decision: independent numerical reconstruction if Figshare missing
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0008

- decision_id: DEC-0008
- date: 2026-07-21
- decision: all inferred values labeled separately
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0009

- decision_id: DEC-0009
- date: 2026-07-21
- decision: full Maxwell FDTD not main route
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0010

- decision_id: DEC-0010
- date: 2026-07-21
- decision: COMSOL not main route
- rationale: frozen A1 architecture/audit decision
- supporting_sources: SRC-0001;SRC-0003;SRC-0005 as applicable
- reversible: yes
- revisit_condition: validation failure or original data/code availability

## DEC-0011

- decision_id: DEC-0011
- date: 2026-07-21
- decision: Primary lifecycle reconstruction uses I_CM0 = 0.44 A m.
- rationale: Shi 2019 main text and Figure 4 caption support 0.44 A m.
- supporting_sources: SRC-0001
- reversible: yes
- revisit_condition: Author clarification or corrected source.

## DEC-0012

- decision_id: DEC-0012
- date: 2026-07-21
- decision: I_CM0 = 0.40 A m is retained only as the Figure 3 caption sensitivity variant.
- rationale: It appears in the Figure 3 caption and conflicts with the main value.
- supporting_sources: SRC-0001
- reversible: yes
- revisit_condition: Author clarification or corrected source.

## DEC-0013

- decision_id: DEC-0013
- date: 2026-07-21
- decision: Continuous analytical Fourier transform is the A2 primary spectrum baseline; discrete FFT is numerical validation only.
- rationale: The closed-form transform avoids finite-window bias and directly implements Eq. (5).
- supporting_sources: SRC-0001
- reversible: yes
- revisit_condition: A different documented transform convention is adopted.

## DEC-0014

- decision_id: DEC-0014
- date: 2026-07-21
- decision: Main FFT validation uses no mean removal and no window function.
- rationale: The lifecycle waveform decays at both generated window edges and this preserves direct continuous-transform comparison.
- supporting_sources: SRC-0001; A2 numerical validation decision
- reversible: yes
- revisit_condition: A later task validates alternative preprocessing.

## DEC-0015

- decision_id: DEC-0015
- date: 2026-07-21
- decision: Figure 3 inset is outside A2 because collision-current data are unavailable.
- rationale: Reconstructing it would require fabricated or inferred collision data.
- supporting_sources: SRC-0001; SRC-0002
- reversible: yes
- revisit_condition: Original data become available.

## DEC-0016

- decision_id: DEC-0016
- date: 2026-07-21
- decision: Figure 4a green collision FFT and purple isolated-streamer curves are outside A2.
- rationale: The required simulation waveforms are unavailable; A2 contains analytical lifecycle curves only.
- supporting_sources: SRC-0001; SRC-0002
- reversible: yes
- revisit_condition: Original data become available.

## DEC-0017

- decision_id: DEC-0017
- date: 2026-07-21
- decision: Stage 1 consists of Stage 1.1 Registry and Model Audit and Stage 1.2 Analytical Lifecycle and Radiation Chain.
- rationale: Establishes one unambiguous overall-stage hierarchy while retaining historical filenames.
- supporting_sources: Stage 1 artifact baseline
- reversible: no
- revisit_condition: Project governance is formally revised.

## DEC-0018

- decision_id: DEC-0018
- date: 2026-07-21
- decision: Stage 1 contains no PDE numerical solver implementation.
- rationale: Stage 1 freezes literature definitions and the analytical radiation chain only.
- supporting_sources: SRC-0001; Stage 1 scope
- reversible: no
- revisit_condition: Never within the closed Stage 1 baseline.

## DEC-0019

- decision_id: DEC-0019
- date: 2026-07-21
- decision: The analytical chain for Shi 2019 Eqs. (4)-(6) is verified.
- rationale: Closed-form derivative, Fourier transform, FFT cross-check, ESD conversion, and band integration passed automated acceptance.
- supporting_sources: SRC-0001; results/stage_A2/summary.json
- reversible: yes
- revisit_condition: A mathematical or unit error is demonstrated.

## DEC-0020

- decision_id: DEC-0020
- date: 2026-07-21
- decision: Primary I_CM0 is 0.44 A m; 0.40 A m remains a literature-conflict sensitivity branch.
- rationale: Main text/Figure 4 support 0.44 while Figure 3 caption states 0.40.
- supporting_sources: SRC-0001
- reversible: yes
- revisit_condition: Author correction or original code becomes available.

## DEC-0021

- decision_id: DEC-0021
- date: 2026-07-21
- decision: T0 is not the actual peak time of published Eq. (4) when alpha differs from beta.
- rationale: Differentiation gives t_peak=T0+ln(alpha/beta)/(alpha+beta).
- supporting_sources: SRC-0001; docs/stage_A2_method.md
- reversible: no
- revisit_condition: The published equation is formally corrected.

## DEC-0022

- decision_id: DEC-0022
- date: 2026-07-21
- decision: The continuous analytical Fourier transform is the lifecycle spectrum baseline; discrete FFT is validation only.
- rationale: Avoids finite-window bias and directly represents Eq. (5).
- supporting_sources: SRC-0001; results/stage_A2/validation/fourier_validation.csv
- reversible: yes
- revisit_condition: A different documented transform convention is adopted.

## DEC-0023

- decision_id: DEC-0023
- date: 2026-07-21
- decision: Figure 2 and collision/isolated curves are deferred to Stage 4.
- rationale: They require numerical streamer and control waveforms that Stage 1 does not generate.
- supporting_sources: SRC-0001; SRC-0002
- reversible: yes
- revisit_condition: Stage 4 begins after Stages 2 and 3 pass.

## DEC-0024

- decision_id: DEC-0024
- date: 2026-07-21
- decision: Figure 4b reconstruction is deferred to Stage 5.
- rationale: Original collision data and discrete integration bins are unavailable.
- supporting_sources: SRC-0001; SRC-0002
- reversible: yes
- revisit_condition: Stage 5 begins or original data become available.

## DEC-0025

- decision_id: DEC-0025
- date: 2026-07-21
- decision: Without Figshare, the project follows an independent numerical reconstruction route.
- rationale: Original data absence must not be disguised as pointwise reproduction.
- supporting_sources: SRC-0001; SRC-0002
- reversible: yes
- revisit_condition: Authentic original data become available.

## DEC-0026

- decision_id: DEC-0026
- date: 2026-07-21
- decision: Stage 1 completion does not imply that a streamer model has been implemented.
- rationale: Transport, Poisson, SP3, coupled propagation, and collision remain unimplemented.
- supporting_sources: Stage 1 scope and artifact manifest
- reversible: no
- revisit_condition: Never; later accomplishments belong to later stages.

## DEC-0027

- decision_id: DEC-0027
- date: 2026-07-21
- decision: Stage 2 will implement ISG-0, axisymmetric Poisson, and three-group SP3 as WP2.1-WP2.3 under one overall stage.
- rationale: Prevents proliferation of ambiguous overall stage labels.
- supporting_sources: docs/stage_mapping.md
- reversible: yes
- revisit_condition: Stage 2 governance is revised before implementation.

## DEC-0028

- decision_id: DEC-0028
- date: 2026-07-21
- decision: Stage 2 completion requires all three numerical building blocks and their unified interfaces to pass before Stage 3.
- rationale: Partial module completion is insufficient for a coupled fluid solver entry.
- supporting_sources: docs/stage2_entry_plan.md
- reversible: yes
- revisit_condition: Verified dependency analysis changes the coupling order.

## DEC-0029

- decision_id: DEC-0029
- date: 2026-07-21
- decision: The Stage 1 reproducible baseline is frozen by stage1_artifact_manifest.csv SHA256 values after closure validation.
- rationale: Enables later detection of accidental changes to closed artifacts.
- supporting_sources: docs/stage1_artifact_manifest.csv
- reversible: yes
- revisit_condition: A documented Stage 1 correction requires a new baseline and release note.

## DEC-0030

- decision_id: DEC-0030
- date: 2026-07-21
- decision: The prior Stage 2 completion conclusion is withdrawn; Stage 2 is reopened and not complete, and Stage 3 is stopped and reset to not started.
- rationale: Executable inspection found no OpenCharge boundary implementation, no coupled SP3 Robin iteration, an unconditional Robin test, and hardcoded acceptance observations inconsistent with the report.
- supporting_sources: docs/stage2_reopen_report.md; results/stage2/forensic/report.md
- reversible: yes
- revisit_condition: A traceable-evidence Stage 2 closure validator passes after real implementations and measurements replace the invalid evidence.

## DEC-0031

- decision_id: DEC-0031
- date: 2026-07-21
- decision: Every Stage 2 PASS metric must trace through a run_id to command, configuration, executable hash, logs, exit status, and result hash.
- rationale: Prevents generated constants, placeholders, and report self-reference from being accepted as measurements.
- supporting_sources: results/stage2/provenance/run_registry.csv; tools/audit_stage2_evidence.py
- reversible: no
- revisit_condition: Never; evidence traceability remains mandatory.

## DEC-0032 — Stage 3 reproducible transport baseline (2026-07-22)

The original Liu–Pasko lookup table is unavailable and Figure 1 solid/dashed curves locally overlap at high field, so automatic digitization is not unique. Stage 3 uses the analytic Morrow–Lowke (1997) Appendix A model for solver development and validation. Liu–Pasko digitization or an independently generated Boltzmann table is deferred to Stage 4 calibration. Stage 3 results must not be described as a quantitative Shi et al. (2019) reproduction. The retained 600 dpi Figure 1 rendering is Stage 4 preparatory material only.

## DEC-0033 — Stage 3 measured gate failure blocks closure (2026-07-22)

Measured Stage 3 coupled runs must not be closed by selecting favorable head definitions. The resumed coarse ML outputs form opposite-sign heads and pass transport/provenance checks, but field-head and charge-head positions differ by more than two cells and the SP3-off run does not show the required positive-head weakening. Stage 3 therefore remains not complete until the coupled physics or acceptance evidence satisfies those gates with traceable measured output.
