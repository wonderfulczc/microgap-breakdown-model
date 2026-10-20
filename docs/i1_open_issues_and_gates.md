# I1 long-term open issues and gates

Established: 2026-10-11. Last reviewed: 2026-10-15. Companion: [geometry, boundary, and port contract](i1_geometry_boundary_port_contract.md). This ledger records **what is established, what remains unresolved, and the evidence needed to change each status**. A status of `VERIFIED` means only the stated numerical/theoretical proposition was checked; it never means experimental validation. `NOT_RESOLVED` is not a synonym for `PASS`.

## Status and maintenance rules

Allowed statuses: `OPEN`, `BLOCKER`, `DOMAIN_RESTRICTED`, `CONDITIONAL`, `NOT_RESOLVED`, `VERIFIED`, `CLOSED`. `BLOCKER` is scoped to the *Impact* field, not a halt on unrelated work. Close an issue only when its stated resolution gate is supported by inspectable evidence, not by a report title, test count, or more simulation runs alone. Retain its ID, history, and evidence after closure. When new evidence addresses the same problem, update the existing ID rather than duplicating it; add a separate issue only for a genuinely different mechanism. Do not upgrade `NOT_RESOLVED` to `PASS` by inference. Do not silently replace `q_A` with `q_B`. Frozen C/G/F solvers and reference outputs remain unchanged. At each I1 Codex review, **check** both documents; edit issue status only with explicit user authorization or when document maintenance is included in that task.

Each history below begins at ledger creation; it does not invent earlier ledger-state transitions. Existing reports retain their original scientific status.

## I1-GEO-01: Parameterized geometry and net-gap validation

- **Issue ID:** `I1-GEO-01`.
- **Status:** `BLOCKER` (geometry-dependent cold-state identification and V5 entry; not all canonical development diagnostics).
- **Problem:** There is no approved geometry-family interface that derives and automatically validates the continuous metal-to-metal net gap from declared facing surfaces and distinguishes it from center-coordinate metadata and mesh-face separation. The historical 65/70/75 um discrepancy is its regression case, not the sole target geometry.
- **Evidence:** [V5-P2 final report](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), section 5; [V5 entry decision](../validation/i1_part_v_preflight2/V5_ENTRY_DECISION.md); [electrode.cpp](../cpp/src/electrode.cpp); [Stage C3](stage_c3_electrode_current_resistance.md); [FP-1 synthetic check and B0 audit](../adapters/comsol_fp1/FP1_INTERFACE_REPORT.md); [B0 recovery](../adapters/comsol_fp1/b0_recovery/B0_RECOVERY_REPORT.md), [FP-2A static routing](../adapters/comsol_fp2/FP2_INTERFACE_REPORT.md), [FP-2B derived audit](../adapters/comsol_fp2/FP2B_ENGINEERING_REPORT.md), [FP-2B native integration](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md), and [FP-2 engineering closeout](../adapters/comsol_fp2/FP2_ENGINEERING_CLOSURE_REPORT.md) (synthetic geometry only; B0 actual geometry not validated).
- **Root cause:** **Confirmed:** coordinate/solid-surface definitions and a 5 um spherical tip radius account for 65/70/75; cell-center mask placement adds mesh dependence. **Unproven:** whether the historical shape matches any actual fabricated electrode.
- **Impact:** Blocks unqualified geometry-dependent K2/V5 comparisons. Lack of actual-electrode dimensions does **not** by itself block verification of a generic parameterized geometry interface; actual-electrode scientific claims separately need source provenance.
- **Owner stage:** Stage B (parameterized geometry and descriptors), Stage C (supported axisymmetric reconstruction), Stage E (3D-required routing); Part V reviews adoption.
- **Resolution gate:** Specify supported geometry families, parameter units, conductor surfaces, coordinate system, gap region, and geometry revision; compute positive `d_net` from the declared facing surfaces and automatically check nonoverlap, metadata consistency, discrete face placement, and 2D/3D routing. Pass independent geometry instances plus the unchanged historical 65 um regression, explicitly preserving its 70 um driver and 75 um metadata definitions. Ambiguous multi-gap cases must be marked unresolved, not silently collapsed. Real dimensions are required only for a later actual-electrode claim, not for this interface gate.
- **Dependencies:** I1-BC-02 for a joint electrostatic problem; I1-NUM-04 before dynamic identification. Actual geometry/metrology only where an actual-electrode claim is made.
- **Next authorized action:** Review the new native view against the unchanged historical 65/70/75 um geometry witness and at least one independently specified actual geometry; quantify cell-mask and net-gap error before geometry-dependent identification. Do not change frozen outputs.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-13.
- **Status history:** 2026-10-11: initialized `BLOCKER` from V5-P2 geometry audit. 2026-10-11: scope and closure gate generalized to parameterized geometry; status unchanged. 2026-10-12: FP-1 synthetic gap check added; B0 CAD provenance pending, status unchanged. 2026-10-12: B0 exported-surface 70 um witness distinguished from historical 65/70/75 um regression; native geometry mismatch confirmed, status unchanged. 2026-10-12: [FP-2B derived geometry audit](../adapters/comsol_fp2/FP2B_ENGINEERING_REPORT.md) passed compiled needle-plane mask regression and synthetic dual-needle classification, but no injectable native geometry interface existed. 2026-10-12: [native integration](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md) added a tested needle-plane/dual-needle view and net-surface-gap calculation; historical 65 um instance and actual geometry/mesh error still need independent acceptance, so `BLOCKER` remains. 2026-10-12: FP-2 closeout retained this gate; no actual-geometry acceptance was added.

- **2026-10-13 evidence review / status history:** [B0 exported-profile static comparison](../adapters/comsol_fp2/FP2_B0_GEOMETRY_BOUNDARY_REPORT.md) established the rounded-cap/cone/cylinder profile and 70 um continuous gap against the native export. [Local cold-model validation](../adapters/comsol_fp2/FP2_LOCAL_NUMERICAL_VALIDATION_REPORT.md) reused it unchanged on A/B domains. This is exported-geometry evidence, not CAD/metrology or dynamic accuracy acceptance; `BLOCKER` remains for unqualified identification. No geometry reconstruction is requested again by this review.

## I1-BC-02: Configurable electrostatic outer-boundary classes

- **Issue ID:** `I1-BC-02`.
- **Status:** `BLOCKER` (global weighting field, capacitance, and physical terminal interpretation).
- **Problem:** Supported geometry families lack a reviewed, configurable set of physically meaningful outer-boundary classes with matching Poisson/weighting conditions and verifiable outer-charge accounting. In the historical case, adjacent HV shank and zero-potential outer gas boundary create a potential jump.
- **Evidence:** [V5-P2 final report](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), sections 3 and 5; [Poisson implementation](../cpp/src/poisson.cpp); [V5 entry decision](../validation/i1_part_v_preflight2/V5_ENTRY_DECISION.md); [FP-1 static interface report](../adapters/comsol_fp1/FP1_INTERFACE_REPORT.md); [B0 recovery](../adapters/comsol_fp1/b0_recovery/B0_RECOVERY_REPORT.md); [FP-2B synthetic face audit](../adapters/comsol_fp2/FP2B_ENGINEERING_REPORT.md); [native integration](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md); [FP-2 closeout](../adapters/comsol_fp2/FP2_ENGINEERING_CLOSURE_REPORT.md) (outer electrostatic correction still unresolved).
- **Root cause:** **Confirmed:** existing Dirichlet boundary arrangement meets the HV shank and produces a mesh-dependent high-field exit; at the finest static grid, the upper region contributes about 49.74% of `q_B`. **Unproven:** which replacement domain/return boundary best represents the experimental fixture.
- **Impact:** Blocks global `Emax`, `C_geom`, and physical terminal interpretation for unreviewed boundary classes; limits K2/Part VI claims. It does not require selecting one universal experimental return boundary before the generic interface can be specified.
- **Owner stage:** Stage B (boundary geometry and class specification), Stage C (Poisson/weighting realization), Part VI (terminal reference).
- **Resolution gate:** Define configurable outer-boundary classes with their physical use conditions, mathematical boundary values/fluxes, return reference, conductor-exit handling, and independently computed `Q_outer`. Verify consistent Poisson/weighting problems and representative domain/mesh sensitivity for each supported class, including a regression that detects the historical shank-exit artifact. A class not yet verified remains unsupported; a changed frozen solver boundary needs separate authorization and baseline comparison.
- **Dependencies:** I1-GEO-01 for geometry identity; I1-PORT-03 for joint terminal interpretation.
- **Next authorized action:** Document candidate boundary classes and read-only verification cases, including the historical failure case; do not silently replace the frozen condition or presume an actual fixture boundary.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-13.
- **Status history:** 2026-10-11: initialized `BLOCKER` from the V5-P2 shank-exit audit. 2026-10-11: closure gate generalized to configurable, testable boundary classes; status unchanged. 2026-10-12: FP-1 synthetic label gate added, no physical boundary-class closure. 2026-10-12: B0 outer zero-charge label retained, but full physical boundary provenance not recovered; status unchanged. 2026-10-12: FP-2B synthetic role/face enumeration distinguished axis, outer and terminals. 2026-10-12: native integration retains outer-row gas and separate exterior faces, but no independent outer electrostatic correction or physical boundary-class validation; `BLOCKER` remains. 2026-10-12: FP-2 closeout confirmed outer **particle** flux is exported but `Q_outer` remains unresolved; status unchanged.
- **2026-10-16 scoped photon-boundary evidence, electrostatic status/Gate unchanged:**
  [B0 photoionization audit](../validation/i1_part_vi/b0_photoionization_20261016/I1_B0_PHOTOIONIZATION_REPLENISHMENT_REPORT.md)
  finds an explicit constructor rejection for geometry_view+SP3 ON. The rectangular SP3 API
  accepts no electrode geometry/internal photon boundary; zeroing metal emission/Sph is not
  such a boundary. This is a distinct photon-boundary dependency recorded here without a
  duplicate Issue or changed electrostatic acceptance definition. DecisionC stops before
  initialization/step/Poisson/SP3; no ON/OFF comparison or weak-photo/failed-bridge conclusion.
  Source defaults (including excitation_ratio0.1) are not device calibration. Archived OFF
  Townsend auxiliary gamma_critical~0.07416 assumes complete ion feedback; actual gamma and
  cathode ion flux remain unknown. Next action is scientific boundary-contract/authorization
  review only, not removing the guard or automatically adding ion emission/transport.

- **2026-10-13 evidence review / status history:** [Local manufactured-charge comparison](../adapters/comsol_fp2/FP2_LOCAL_NUMERICAL_VALIDATION_REPORT.md), sections 1-2, tested zero-Neumann A/B domains at identical 5 um spacing. Maximum common-region space-charge field difference was 1.183%, below the predeclared 5% exploratory target; the source distributions were not tuned. This supports only those compact monopole/dipole sources, not arbitrary later space charge or a global boundary class. `BLOCKER` remains for global weighting/capacitance/physical-terminal claims. The next action is review of the existing evidence, not automatic domain expansion.

## I1-PORT-03: Geometry-family conjugate terminal operator

- **Issue ID:** `I1-PORT-03`.
- **Status:** `BLOCKER` (physical terminal current and K2/Part VI port gate).
- **Problem:** A conjugate terminal-charge operator and signed accepted-step current interface have not been defined and tested across supported geometry/boundary families. Native center-field `q_A` is not conjugate to the historical Poisson operator's Laplace weight; historical same-operator `q_B` has a discrete-network interpretation, but is not yet a family-wide or physical terminal interface.
- **Evidence:** [Part V-A port audit](../validation/i1_part_v/port_mapping_audit.md); [V5 preflight current accounting](../validation/i1_part_v_preflight/p3_current_accounting.md); [V5-P2 final report](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), sections 3-4; [discrete comparison](../validation/i1_part_v_preflight2/discrete_port_comparison.csv); [dynamic budget](../validation/i1_part_v_preflight2/dynamic_current_budget.csv); [Stage C3](stage_c3_electrode_current_resistance.md); [FP-1 synthetic operator check](../adapters/comsol_fp1/FP1_INTERFACE_REPORT.md); [FP-2A read-only operator](../adapters/comsol_fp2/FP2_INTERFACE_REPORT.md); [FP-2B engineering audit](../adapters/comsol_fp2/FP2B_ENGINEERING_REPORT.md); [native port candidate and accepted-step export](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md); [independent three-step charge audit](../adapters/comsol_fp2/FP2_ENGINEERING_CLOSURE_REPORT.md) (physical terminal remains pending).
- **Root cause:** **Confirmed:** center-field surface summation is not the discrete Poisson face flux (the earlier 22.26% mismatch); `q_B` passes same-operator reciprocity/energy and an independently computed source-response budget. The native gas-metal Poisson link uses full `h`, while electron absorption uses `h/2`. **Unproven:** whether `q_B` maps to the physical terminal after geometry/boundary corrections and dynamic convergence.
- **Impact:** Blocks using `q_A` as the unique I1 port charge, equating `i_total_hv` with a G3 external-circuit port, and closing the physical K2/Part VI mapping. The historical `q_B` discrete result remains valid for its tested operator but does not certify other geometries.
- **Owner stage:** Part V/VI (port review) and Stage C (read-only native diagnostics).
- **Resolution gate:** Define an operator construction from each supported geometry's Poisson links, terminal surfaces, and selected outer-boundary class, with an explicit charge/current sign and time-layer contract. Demonstrate static reciprocity and energy conjugacy, independently evaluated outer contributions, and signed accepted-step charge/current budgets (including absorption and numerical corrections) on **at least two independent geometry instances**, not only the historical needle-plane regression. Separately mark physical terminal interpretation `NOT_RESOLVED` until electrode-face/control-surface consistency and per-case dynamic error evidence support it. Any replacement of the native observable needs explicit interface review.
- **Dependencies:** I1-GEO-01, I1-BC-02, I1-NUM-04, I1-MAT-05 for zero-crossing claims.
- **Next authorized action:** Preserve `q_A`; specify a geometry-parameterized `q_B` operator and independent-geometry read-only checks. No automatic interface swap or physical-terminal claim.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-13.
- **Status history:** 2026-10-11: initialized `BLOCKER` for *physical* terminal; same-operator `q_B` reciprocity remains conditionally verified. 2026-10-11: family-wide closure gate added; status unchanged. 2026-10-12: synthetic plate/coax operator bookkeeping passed, but target Poisson links and physical plane still pending; status unchanged. 2026-10-12: independent read-only Stage-C-stencil role-mask checks passed on synthetic masks; B0 mapping and physical port did not pass. 2026-10-12: FP-2B synthetic dual-needle same-operator `q_B`/energy passed. 2026-10-12: [native integration](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md) now exports separate `q_A/q_B`, accepted-step faces and corrections; physical control surface, `Q_outer` and independent actual-geometry validation are still absent, so `BLOCKER` remains. 2026-10-12: [FP-2 closeout](../adapters/comsol_fp2/FP2_ENGINEERING_CLOSURE_REPORT.md) independently attributed the synthetic legacy current residual to activity-floor charge, but correction current exceeds modeled net transport current by about 2.14e4; mathematical accounting passed, physical-port gate unchanged.

- **2026-10-13 evidence review / status history:** [Local validation](../adapters/comsol_fp2/FP2_LOCAL_NUMERICAL_VALIDATION_REPORT.md) preserves separate local q_A/q_B and never equates either to a complete physical terminal. The new species/face/correction audit passes on the recorded short windows, but no Q_outer or global terminal acceptance was added. `BLOCKER` and Part VI dependencies remain unchanged.

## I1-NUM-04: Stage C dynamic field and mechanism convergence

- **Issue ID:** `I1-NUM-04`.
- **Status:** `BLOCKER` (unqualified V5 parameter fitting; not separately authorized bounded state screening, local prediction or the explicitly authorized single-point V6 test).
- **Problem:** Exact weighted finite-volume budgets do not establish numerical accuracy of physical `G_R` and reaction/transport/boundary/mobility trajectories. Limited local pairwise and prediction evidence now exists, but geometry-family-wide applicability and certified continuous error bounds remain unresolved; each new run still needs its own observable/error labels.
- **Evidence:** [Part V-A report](../validation/i1_part_v/PART_V_A_VALIDATION_REPORT.md); [V5 preflight report](../validation/i1_part_v_preflight/V5_PREFLIGHT_REPORT.md); [P2 convergence data](../validation/i1_part_v_preflight/p2_stage_c_convergence.csv); [V5-P2 final report](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), section 5; [V5 entry decision](../validation/i1_part_v_preflight2/V5_ENTRY_DECISION.md); [FP-2B fixture limits](../adapters/comsol_fp2/FP2B_ENGINEERING_REPORT.md); [native three-step synthetic trace](../adapters/comsol_fp2/FP2B_NATIVE_INTEGRATION_REPORT.md); [FP-2 closeout](../adapters/comsol_fp2/FP2_ENGINEERING_CLOSURE_REPORT.md) (no representative dual-needle convergence).
- **Root cause:** **Confirmed:** 30x72 to 40x96 dynamic `G_R` differs about 9.47% at common time; halving `dt` changes it about 0.0323%; A/T/B rates differ about 24/42/38%. Static tip-ROI field and seeded `G_R` improve on finer grids, but are not dynamic convergence. **Unproven:** relative shares of mask, outer boundary, weighting-field, and transport-discretization error in a corrected physical geometry.
- **Impact:** Blocks unqualified Stage C K2 cross-case prediction and formal V5 parameter fitting on unresolved trajectories. Does not require every conceivable geometry to converge before any supported subset can be assessed, nor invalidate V1/V2 identities.
- **Owner stage:** Stage C numerical diagnostics; Part V applies the readiness gate.
- **Resolution gate:** Establish representative verification cases for each supported geometry/boundary class, including spatial levels, time-step controls, `G_R` and A/T/B/M mechanism checks, and separated mask/solver/weighting errors. Define an automated **per-run** provenance and error-status record for mesh, time step, geometry and boundary, with `RESOLVED`/`CONDITIONAL`/`NOT_RESOLVED` on each claimed observable; a new run inherits a representative estimate only after a documented applicability check. Formal V5 uses only runs whose required dynamic quantities meet their declared accuracy gate. No requirement that all possible geometries be pre-converged; no promotion from algebraic or static checks alone.
- **Dependencies:** I1-GEO-01, I1-BC-02; I1-MAT-05 for low-field interpretation.
- **Next authorized action:** Review the completed V7 cold-output-only comparison and Part VI theory checklist. No automatic voltage/grid/shift run, refit, additional state, physical-model execution or Part VI implementation follows.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-15.
- **Status history:** 2026-10-11: initialized `BLOCKER` from V5-P/V5-P2 physical-convergence evidence. 2026-10-11: closure gate revised to representative tests plus per-run error flags; status unchanged. 2026-10-12: FP-2A static synthetic checks do not constitute dynamic mesh/time convergence; status unchanged. 2026-10-12: FP-2B native-function face fixture checked a single synthetic budget. 2026-10-12: native dual-needle three-step synthetic smoke and exact discrete ledger passed with correction accounting, but no representative mesh/time convergence or per-run physical accuracy estimate; `BLOCKER` remains. 2026-10-12: FP-2 closeout reproduced the correction-inclusive identity from exports and found floor-current dominance; no physical-accuracy upgrade.

- **2026-10-13 evidence review / status history:** [Local cold-model report](../adapters/comsol_fp2/FP2_LOCAL_NUMERICAL_VALIDATION_REPORT.md) and [internal V5 decision](../adapters/comsol_fp2/FP2_LOCAL_V5_ENTRY_DECISION.md) add four bounded native runs (436 accepted steps, 222.83 s total). Full/common ROI accounting is independently verified; corrections no longer dominate the unsigned transport scale in these short tests. However, all runs stopped before the fixed 10 ps goal, the larger/finer cases reached only 0.90/0.98 ps, and common-ROI spatial A/B mechanism errors are 12.29%/21.71% despite GR error 3.49%. `LOCAL_DYNAMIC_NUMERICAL_CONVERGENCE` and `V5_INTERNAL_REFERENCE_READINESS` remain `NOT_RESOLVED`; this Issue stays `BLOCKER`. Next authorized action is scientific review of the retained results and resource/discretization limits, not automatic continuation, retuning, or fitting.

- **2026-10-13 accuracy/cost review / status history:** [FP2-EXEC-20261013-01 report](../adapters/comsol_fp2/local_accuracy_run/FP2_ACCURACY_AND_COST_REPORT.md) adds offline attribution and nine bounded native performance/probe runs (45 accepted steps, 46.84 s). The old 0.98 ps common-window GR RMS difference is 3.49%, but growth RMS difference is 12.22% and fine response only 0.847%. An order-averaged initial 2x2 attribution gives field/coefficient/weight contribution +0.12403 S/s and density contribution -0.01298 S/s to the reaction difference. Analytic cell-average greatly reduces electron-count sampling differences, but initial A still differs by 13.22%; geometry/gradient/quadrature effects are not individually identified. Default-OFF weight and fixed-native-matrix caches passed exact short-step equivalence, including source/face/correction arrays and observer ON/OFF. Measured coarse sustained time decreased from 0.464 to 0.156 s/step, while fine native Poisson iteration remains dominant. The complete four-case 10 ps group was estimated at 2981.60 s against 1753.16 s remaining, so it was **not started** under the registered RESOURCE_MISMATCH rule. `LOCAL_DYNAMIC_NUMERICAL_CONVERGENCE`, W_rho reference accuracy and `V5_INTERNAL_REFERENCE_READINESS` remain `NOT_RESOLVED`; I1-NUM-04 remains `BLOCKER`. The only next action is scientific review of native linear-solve cost and scaling/field accuracy before a separately authorized backend/resource decision, not automatic continuation or fitting. Raw true residual and KSP stopping status are separately retained. An [artifact-integrity incident](../adapters/comsol_fp2/local_accuracy_run/ARTIFACT_INTEGRITY_INCIDENT.md) records the historical CTest Gate E write-back: its three CSV hashes are unchanged, but its smoke JSON timing identity was overwritten and cannot be claimed preserved. The prior local-validation artifacts themselves retain all 73 protected hashes; no old Gate is upgraded.

- **2026-10-13 completed-window and reconstruction-support review / status history:** The [same-grid linear/window report](../adapters/comsol_fp2/linear_window_run/FP2_LINEAR_ACCURACY_AND_WINDOW_REPORT.md) completed the separately authorized four 10 ps cases without changing the earlier RESOURCE_MISMATCH record. Tested same-grid GR algebraic sensitivity is at most 5.88126e-11; common-region spatial E_G=2.53537%, E_growth=6.57009%, but sum b_j=7.2223% exceeds the unchanged 5% mechanism gate (b_A=6.84871%). The [controlled reconstruction review](../adapters/comsol_fp2/spatial_reconstruction_run/FP2_SPATIAL_RECONSTRUCTION_REPORT.md) then reused t=0/2/10 ps snapshots: 99.83208%/99.80143%/99.55876% of the common-region local-L1 reaction difference lies outside the fixed union of metal-crossing gradient templates. That union contains only 0.0699878% of the initial fine-grid unsigned reaction scale. The authorized near-only reconstruction cannot alter the dominant frozen-state interior support at fixed A_h and phi/psi, so it was stopped **before implementation**; no new dynamic group ran. This does not uniquely distinguish potential/geometry, interior-gradient and quadrature errors or prove a new method failed. All existing source/binary/results are protected; historical smoke pre-incident JSON bytes remain UNKNOWN. Status remains `BLOCKER`, V5 internal readiness `NOT_RESOLVED`; other physical Gates are unchanged. The next action is review of this support-based stop and unresolved spatial discretization evidence, not automatic implementation of a second method, refinement, or fitting.

- **2026-10-13 single-refinement recovery review / status history:** The authorized h=0.625 um, 800x1600 native case was interrupted with scalar exports through 9.02 ps and an empty 9.0 ps checkpoint. With explicit user approval, an independent driver linked the unchanged native solver library and resumed the valid 8.8 ps checkpoint; the old prefix and all interrupted files remain intact. The [single-refinement report](../adapters/comsol_fp2/single_refinement_run/recovery_8p8/FP2_SINGLE_REFINEMENT_REPORT.md) and [machine-readable decision](../adapters/comsol_fp2/single_refinement_run/recovery_8p8/refinement_gate_summary.json) document completion to 10 ps, exact replay of the available overlap scalar records, and no duplicate increment accumulation. For 1.25 to 0.625 um, COMMON_GAP E_G=1.169067%, E_growth=3.044452%, b_A=3.180728%, sum b_j=3.379323%; the fixed pairwise targets pass, with 9.788880% response. Initial A_R progresses 0.80718118, 0.93015567, 0.99902567, 1.03194840 S/s, with observed orders 0.8364 and 1.0648, not a proven continuum error bound. Common-region correction/unsigned transport is 3.17449e-7, while correction/net boundary-plus-source is 0.0111633; both scales are retained. The separate fine-grid timestep test remains NOT_TESTED, and lost original timing/control records remain missing. **Decision A permits only CONDITIONAL preparation of this local 500 V/SP3-OFF/specified-seed internal candidate reference.** Formal general V5 identification remains `BLOCKER`; no ROM fit, C2-Q-specific admission, COMSOL validation, global terminal, low-field-physics or downstream Gate is upgraded. Next authorized action: review these archived results; no further computation is started.

- **2026-10-13 Part V5 bounded state-sufficiency review / status history:** The explicitly authorized [V5 state screen](../validation/i1_part_v/state_sufficiency/V5_STATE_SUFFICIENCY_REPORT.md) and [Gate JSON](../validation/i1_part_v/state_sufficiency/v5_state_screen_gate_summary.json) evaluate neutral 10.8/12/13.2 um Gaussian cell-average states at the same FULL_GAS G target and W=0, without ROM fitting. Five shape/grid configurations and nine fresh native accepted steps include the selected narrow/wide pair at 0.625 um, 20/10 fs, and tighter same-grid Poisson checks. The fine pair differs by 0.00584174785 S/s in complete A/T/B/M dotG; B contributes +0.00832208047 S/s and A -0.00232122430 S/s. The pre-registered conservative sensitivity envelope is 0.03730060408 S/s, dominated by the sum of per-state grid changes, not linear solve or time-step differences. Execution and 11 read-only tests pass, but both C1 and C2-Q closure screens are `NUMERICALLY_UNRESOLVED`. The raw half-gap/S_pair is 0.28143%; the 5% instantaneous-rate target is `NOT_FALSIFIED_IN_TESTED_PAIRS`, not a sufficiency PASS. State matching has no certified uniform Lipschitz bound. Original FP-2 pairwise evidence, physical port, W-specific dynamic admission and all other Issue statuses are unchanged. No new issue, physical coefficient, solver-core edit, fit or extended trajectory was introduced. **Next authorized action after this task:** scientific review of the archived pair and sensitivity evidence only; no automatic fitting, finer grid, alternative state or new experiment. The requested external roadmap attachment was not found; the [repository-grounded roadmap](i1_research_status_and_roadmap.md) identifies that missing input explicitly.

- **2026-10-14 C1-L0 local-prediction review / status history:** The separately authorized [zero-fit local prediction report](../validation/i1_part_v/local_prediction/V5_LOCAL_PREDICTION_REPORT.md) and [decision](../validation/i1_part_v/local_prediction/v5_local_prediction_summary.json) retain the earlier exact-screen `NUMERICALLY_UNRESOLVED` and instantaneous 5% `NOT_FALSIFIED_IN_TESTED_PAIRS` without refining that small pair gap further. The complete reference pre-layer derivative fixes one effective rate, k_star=9.13463628041107e9 s^-1, with no future-data fit. Reference time extrapolation passes; two 1.25 um known-seed future trajectories and the preselected worst-growth-error wide case at 0.625 um all complete 10 ps. FULL_GAS prediction E_G/E_growth are 0.0259082%/0.495502% (narrow coarse), 0.0334317%/0.640991% (wide coarse), and 0.0858242%/1.61356% (wide fine); fine P95 is 0.166285%. However, same-initial-state native grid sensitivities are 1.17186%/3.09045%, exceeding model errors. The result is `LOCAL_PREDICTION_SUPPORTED_CONDITIONALLY`, not exact closure, continuous error certification, general ROM, identifiability or physical-port acceptance. Full late-time instantaneous dotG-kG is 0.0460849 S/s at 9.98 ps, so successful integrated G prediction does not make the effective rate constant in the PDE. Actual correction budgets are retained on both unsigned and net scales; in particular narrow-case net correction/near-cancelled boundary-source TV is 4.0707%, despite negligible G-growth correction. No W-rho-specific admission follows. This issue remains `BLOCKER` for unqualified fitting, with a supported bounded candidate-output baseline added. **Next authorized action:** review these archived results and limits; stop this task. Further parameter training is not shown necessary in this tested range and is not automatically authorized. No second/third state, new grid, V6/V7, Part VI, thermal or RF computation follows. Only derived timing labels were corrected during archiving, with the draft retained; no historical data or solver files changed. Other Issue statuses remain unchanged.

- **2026-10-15 V5 archive / explicitly authorized V6 evidence:** [V5-G model selection](../validation/i1_part_v/V5_MODEL_SELECTION_FINAL.md) accepts C1-L0 only as a bounded engineering model, with C1 exact closure `NOT_PROVEN`, C2-Q `NOT_SELECTED`, C2-S `OBSERVATIONAL_CANDIDATE_ONLY`. The new prospective [centered600V report](../validation/i1_part_v6/centered_600v/V6_CENTERED_VOLTAGE_TRANSFER_REPORT.md) and [decision](../validation/i1_part_v6/centered_600v/v6_gate_summary.json) record ONE fresh h=1.25um/20fs-cap/10ps native trajectory, 515 accepted steps, same geometry/air/SP3-OFF and unrenormalized cell-average seed, plus one read-only initial RHS replay; scientific process wall1348.822s. The fixed k=9134636280.4 s^-1 is never refit. FULL_GAS E_G=4.281991%, E_growth=45.914630%, P95=6.730002%, true response17.919016%, final signed error-7.084670%. The complete discrete transfer fails the original growth target: `FIXED_K_VOLTAGE_TRANSFER_LIMITED`, not a failed/incomplete solver run. Keep the500V model and do not admit600V fixed-k accuracy. Source/initial density/mask identities are checked; initial600V field scales with600/500 while the independent1V weight field is unchanged. G-correction TV/growth=8.39563e-8, species correction/unsigned transport=5.62384e-7, net-charge correction/net boundary-source TV=0.00732756; all scales and actual corrections are retained. No600V-specific mesh/time error bound exists, so numerical qualification stays conditional;500V sensitivity is not transplanted. New user authorization overrides only old automatic Level2 scheduling, not old position failures or physical gates. The +/-10um cases remain a known applicability limit with no further computation, not a claim about their experimental frequency. No new Issue or status closure; this Issue stays `BLOCKER` for unqualified fitting. Stop for review; no extra voltage, grid, refit, state search or downstream action follows.

- **2026-10-15 document-only position-uncertainty archive (status unchanged):** The user-provided review is now preserved in [the position-uncertainty and paper-wording archive](i1_initial_carrier_position_uncertainty.md), with P1–P8 mapped to the actual unchanged project reports. The +/-10um cases shift the neutral initial carrier cloud, not the electrodes; actual post-initialization ze is about +/-9.958335um and W0 is0. The old V5-C1 `C1_INITIAL_STATE_LIMITED` and growth failures13.072578%/14.077755% remain. Current action is separately **KNOWN_SCOPE_LIMITATION / NO_FURTHER_COMPUTE_AUTHORIZED**, not a new scientific Gate or a claim that off-center inception is rare or absent. No usable device-specific experimental initial-state probability law is supplied (`EXPERIMENTAL_P0=MISSING`); conditional dynamics plus an initial probability law, probability identities and the centroid moment identity remain theoretical memoranda, not implemented stochastic/space-state models. V5-F did not prove ze impossible as a state or a third state necessary; C2-S stays `OBSERVATIONAL_CANDIDATE_ONLY`. V5 model selection and the one600V V6 task are finished; keep `FIXED_K_VOLTAGE_TRANSFER_LIMITED`,500V engineering scope and k_star=9134636280.4 s^-1 as an initial-reference effective rate, not a time average or universal gas constant. The600V10ps/single-grid/GENERIC_DEV restrictions persist. This archive closes no Issue or physical Gate and creates no verification task: no automatic400/800V, geometry/polarity, random-seed, V7, Part VI/FP-3, cold-hot or RF work. External literature and retained paper paragraphs are method background/drafts, not device probability, breakdown or experimental validation. Next authorized action after this document task: stop and await scientific review; no computation is authorized by this ledger entry.

- **2026-10-15 explicitly authorized V7 offline comparison (status unchanged):** The [cold-output report](../validation/i1_part_v/v7_model_comparison/V7_COLD_BASELINE_COMPARISON_REPORT.md) and [Gate](../validation/i1_part_v/v7_model_comparison/v7_gate_summary.json) compare only six archived FULL_GAS scalar curves: centered500V baseline/fixed-peak widths/densities and the completed600V point. No PDE, Poisson, COMSOL, StreamerSolver, compile, fit, new trajectory or shifted-case processing occurred. CONST uses only actual pre G0; TAN/C1 share the frozen9134636280.4 s^-1 and no future truth. At500V, CONST growth error is100%, TAN3.227%--3.469%, C1 0.493%--0.641%; TAN and C1 both meet old exploratory targets. This supports short-window output evolution but does not identify a unique exponential physical law or certify a continuous error bound. At600V, C1 growth error45.914630% and TAN47.783675% retain `FIXED_K_VOLTAGE_TRANSFER_LIMITED`; historical gates/position failures are unchanged. All six C1 historical scalar checks pass. `V7_COLD_OUTPUT_BASELINES=COMPLETED_ON_AVAILABLE_DISCRETE_DATA`, while physical R/R+Vf/RW/Toepler/Mayr/Cassie are `DEFERRED_NO_COMPATIBLE_PORT_REFERENCE`, not FAIL/cancelled, and full physical V7 remains deferred. The [Part VI handoff](../validation/i1_part_v/v7_model_comparison/PART_VI_HANDOFF.md) is ready only for webpage theory review, not implementation; no need to generate fake I=V*G_R or wait for every incompatible classic-model run before reviewing theory. This Issue remains `BLOCKER` for unqualified accuracy/fitting claims; no duplicate Issue or physical gate closure. Stop the bounded task regardless of ranking, with no automatic new model or scientific computation.

## I1-MAT-05: Morrow-Lowke near-zero-field applicability

- **Issue ID:** `I1-MAT-05`.
- **Status:** `DOMAIN_RESTRICTED`.
- **Problem:** The frozen low-positive-field branch has `mu_e=6.87e24/N+338/|E|`, while exactly zero field uses `6.87e24/N`; classical smooth `partial_t mu_e` and zero-crossing port interpretation cannot be presumed.
- **Evidence:** [V5 preflight report](../validation/i1_part_v_preflight/V5_PREFLIGHT_REPORT.md); [Part V-A report](../validation/i1_part_v/PART_V_A_VALIDATION_REPORT.md); [native Morrow-Lowke scan](../validation/i1_part_v/native_zero_field_mobility.csv). The existing source implementation is `cpp/src/streamer/MorrowLowke.cpp`.
- **Root cause:** **Confirmed:** native piecewise coefficient has a nonsmooth/nonbounded positive-field limit versus the `E=0` branch. **Unproven:** physical accuracy of the low-field parameterization for future actual-electrode histories; a numerically smooth interval alone cannot establish it.
- **Impact:** Restricts continuous `dot G_R` mobility-term interpretation, voltage-zero/low-field K1 use, and K2/Part VI port claims. The discrete finite-step product identity remains a separate algebraic result.
- **Owner stage:** Stage C transport/physics audit and Part V domain definition.
- **Resolution gate:** Source-branch and provenance review plus contribution-weighted low-field/zero-crossing statistics on the intended trajectory; specify a supported applicability domain without selecting an arbitrary `E_cut` to improve fit. Physical-coefficient changes need independent authorization.
- **Dependencies:** I1-NUM-04 for representative trajectories; I1-PORT-03 for voltage-zero terminal use.
- **Next authorized action:** Preserve the existing scan and attach low-field occupancy/contribution metadata when a representative, authorized Stage C trajectory exists; do not regularize the frozen coefficient here.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-15.
- **Status history:** 2026-10-11: initialized `DOMAIN_RESTRICTED` from the P1 native-coefficient audit.
- **2026-10-16 archived spatial-path reuse, status unchanged:**
  [Spatial bridge report](../validation/i1_part_vi/spatial_bridge_20261016/I1_SPATIAL_AVALANCHE_AND_CONDUCTION_BRIDGE_REPORT.md)
  evaluates the existing500V/1.25um pre0 and post10ps near-axis70um gas path. No zero/low-field
  path segments occur; frozen single-pass multiplication is14.0564 with transit247.615ps.
  This is neither a low-field repair nor a dynamic avalanche/physical-coefficient validation.
- **2026-10-16 conditional reaction candidate, status unchanged:**
  [K1/K2 bridge report](../validation/i1_part_vi/conductive_channel_bridge_20261016/I1_CONDUCTIVE_CHANNEL_PHYSICAL_BRIDGE_REPORT.md)
  compares a new independent SI implementation to actual unchanged MorrowLowke/ReactionModel
  pure functions and one archived pre cell (194 scalar checks, max relative1.512e-16).
  Exactly-zero-field1ns uses the native D/ionization/attachment=0 branch and e-p recombination;
  it is not a smooth low-field or bipolar recovery validation. Nominal500/800V uniform10ps
  cases stay at prescribed fields, not self-consistent formed-channel fields. No property
  regularization, coefficient change, spatial step or density floor was introduced.
- **2026-10-15 single600V applicability review (status remains DOMAIN_RESTRICTED):** [V6 diagnostics](../validation/i1_part_v6/centered_600v/v6_initial_rate_diagnostic.csv) replay the native complete initial RHS; active zero-field and lowest-branch cells are both0, initial electron-support minimum E is6.652151e6 V/m. Full-step native low-field G/N contributions and low-field A/M projection contributions are0 in this run; sparse pre/post fields are finite. These are contribution-weighted facts for this one600V numerical test, not physical validation of the low-field formula, a new E_cut, or general voltage-zero readiness. The observed initial k0_test=1.620898248e10 s^-1 is diagnostic only; no coefficient or frozen property was changed. Next action after the completed task is review of the recorded domain restriction, not regularization or new runs.
- **2026-10-13 evidence review / status history:** [Local validation](../adapters/comsol_fp2/FP2_LOCAL_NUMERICAL_VALIDATION_REPORT.md) uses the native branch boundary, not a fitted cutoff. Low-field cells account for about 0.85%-3.37% of gas cells, but their electron-weighted and GR contributions are zero in the observed post-floor short windows. This neither repairs the singular limit nor establishes later-time material accuracy; `DOMAIN_RESTRICTED` remains. No coefficient was modified.

## I1-THERM-06: Two-state hot-channel dynamic closure

- **Issue ID:** `I1-THERM-06`.
- **Status:** `NOT_RESOLVED`.
- **Problem:** Fixed-core `(H_c,M_2)` provides a useful local output description in some G1 cases, but a transferable autonomous model for `g_col`, Joule projections, signed loss/redistribution, and state evolution has not passed numerical and independent-condition gates.
- **Evidence:** [Part IV-R final report](../thermal/iv_dual_state_refinement/IV_REFINEMENT_FINAL_REPORT.md); [Part IV first-pass report](../thermal/iv_dual_state_validation/IV_DUAL_STATE_VALIDATION_REPORT.md); [Stage G3](stage_g3_transient_port.md) for the separate circuit reference plane; [FP-4 read-only handoff review](../adapters/i1_fp4_handoff/FP4_COLD_THERMAL_HANDOFF_REVIEW.md) and [separate entry readiness](../adapters/i1_fp4_handoff/fp4_handoff_readiness.json).
- **Root cause:** **Confirmed:** prior report leaves numerical-reference validity, joint state/current coverage, dynamic closure, and independent generalization unresolved; full-rank fitted design alone does not cure coverage. **Possible, unproven:** radial momentum, pressure, higher thermal moments, or core-boundary flux may carry omitted memory.
- **Impact:** Blocks freezing `(H_c,M_2)` as an I1 hot-stage state law or using it as a validated Part VI channel; does not block independent G1 reference studies or the established cold-stage weighted-budget tests.
- **Owner stage:** Stage G1 thermal validation and Part VI only after a reviewed handoff.
- **Resolution gate:** Credible numerical convergence, covered and identifiable train/test domain, independent-condition output and autonomous free-evolution performance above numerical uncertainty, with physical constraints and no truth-state feedback.
- **Dependencies:** G1 numerical-reference and coverage gates; a future cold-to-hot handoff would additionally depend on I1-CIR-07.
- **Next authorized action:** Stop after the separately authorized FP-4 given-state adapter implementation and two-path native G1 replay; submit the thermal-series handoff list for scientific review. No thermal circuit, fit, automatic transition, constructor change or further trajectory is authorized. Preserve Part IV-R `NOT_RESOLVED`; do not add a third state or change G1 physics by default.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-15.
- **Status history:** 2026-10-11: initialized `NOT_RESOLVED` from Part IV-R; no model-freeze claim.
- **2026-10-15 FP-4 review, status unchanged:** AUTO_COLD_TO_THERMAL is `MISSING_PHYSICAL_CLOSURE`: absent electron/heavy/chemical energy partition, neutral mass/momentum, same-EOS energy zero, conservative radial/axial profile and charge/field retention. G0 carries cold connectivity/electrical work, not a thermal profile; its historical result did not reach handoff. `GIVEN_THERMAL_INITIAL_STATE=READY_FOR_BOUNDED_IMPLEMENTATION` only for an independent provenance-bound LTE T/p reference with explicitly given zero radial velocity and column-conduction drive. General nonzero-velocity/conservative input remains `INPUT_CONTRACT_DEFINED_WAITING_DATA` with no current initializer. This does not certify LTE or close the hot ROM. G1 constructor resets time/ledgers and takes no u/U input; public U assignment is not a validated handoff. Existing G2/G3 implementations are acknowledged despite older docs saying they were not started. The 1mm synthetic reference is not converted to a calibrated70um state. No science computation or historical Gate change.
- **2026-10-15 separately authorized FP-4 given-state implementation, status unchanged:** [Adapter report](../adapters/i1_fp4_handoff/given_thermal_run_20261015/FP4_GIVEN_THERMAL_ADAPTER_REPORT.md) and [scoped gates](../adapters/i1_fp4_handoff/given_thermal_run_20261015/fp4_given_thermal_adapter_gate_summary.json) accept `GIVEN_THERMAL_ADAPTER_ENGINEERING_ACCEPTED`. Only the independent native-grid T/p + explicitly zero initial u subset is implemented, with real provenance/EOS/branch/time-support checks and source-bound sidecars. Direct and adapted G1 each completed the existing150-cell1mm10ns primary reference in103 accepted/0 rejected steps,206 total; every accepted T/p/U/time/history comparison is bitwise equal. Final native radial speed169.93798m/s is preserved, not zeroed. Independent initial mass/energy and real boundary/Joule ledgers pass scoped engineering targets; maximum energy residual/Joule6.267e-10. Thirty isolated checks cover27 recorded refusals plus additional input edge cases, including nonzero synthetic clock evaluation without another trajectory. Raw replay source and post-replay rejection-only hardening are separately hashed. No native G1/G2/physics/history edit, old generator execution, G2/RLC/TENG/RF or fitting occurred. This implements the earlier proposal under NEW authorization, not a retroactive change to its historical readiness. AUTO remains missing physical closure, nonzero/conservative initialization unimplemented, true B0 thermal state waiting, LTE experiment unvalidated and Hc/M2 NOT_RESOLVED.

- **2026-10-15 separately authorized thermal-series engineering, status unchanged:**
  [Report](../adapters/i1_fp4_handoff/thermal_series_run_20261015/execution/FP4_THERMAL_SERIES_COUPLING_REPORT.md)
  accepts the independent given-LTE series coupling only. Four10ns trajectories, native G1
  unchanged, actual evolving radial flow retained. Independent mass/energy, column-only RMS
  heating and rejection/common-commit checks pass. This does not close `(H_c,M2)`, supply a
  cold material-energy partition or establish LTE/B0 calibration. AUTO remains
  `MISSING_PHYSICAL_CLOSURE`; no further thermal run or training is authorized after this task.

- **2026-10-15 bounded axial-length sensitivity, physical status unchanged:**
  [Length review](../adapters/i1_fp4_handoff/channel_length_run_20261015/FP4_CHANNEL_LENGTH_DAMPING_REVIEW.md)
  independently reconstructs the original LTE conductivity integral. Given the SAME radial
  T/p, changing ell1mm to a hypothetical70um changes R241.086 to16.876ohm; both material
  inventories and Joule power use the native length. This is not measured B0 thermal data
  or a cold-to-hot mapping. Two authorized100ns native series runs finish in6000 steps;
  thermal resistance rises to about785ohm, and only one signed reversal is resolved.
  LTE, axial uniformity/end losses, energy partition and hot-ROM closure remain open.
  No old four-case reference or transaction/ownership test is rerun; stop for review.

- **2026-10-15 existing thermal-resistance interpretation only, no new Gate:**
  [Offline review](../adapters/i1_fp4_handoff/resistance_interpretation_run_20261015/FP4_THERMAL_RESISTANCE_AND_RINGDOWN_INTERPRETATION.md)
  reads the existing6000 accepted nodes and initial/final profiles only. H finalR/initialR46.5158,
  section conductance2.1498% of initial; native EOS endpoint integral independently matches ell/R.
  QJ4.557763nJ and total material-plus-kinetic increase coexist with kinetic increase6.140176nJ
  and EOS internal-energy decrease1.582413nJ. This is not global energy loss or proof that
  core heating is insufficient: most initial inventory lies outside the descriptive50um region,
  and local cumulative heating/heat/pressure-work budgets are absent. No mechanism percentages
  or experimental root cause inferred. No intermediate profile, no new trajectory, no thermal
  model/closure change; original Hc/M2, LTE, B0 and AUTO statuses remain unchanged.

- **2026-10-16 separately authorized early driven/free comparison, status NOT_RESOLVED unchanged:**
  [Causal review](../adapters/i1_fp4_handoff/early_relaxation_run_20261016/FP4_EARLY_THERMAL_RELAXATION_CAUSAL_REVIEW.md)
  and [decision](../adapters/i1_fp4_handoff/early_relaxation_run_20261016/fp4_early_relaxation_decision.json)
  support `ARTIFICIAL_INITIAL_RELAXATION_DOMINANT_IN_THIS_REFERENCE`, not physical LTE validation.
  Two25ns native paths,1000 steps each, preserve the given70um T/p/EOS/grid. Driven replay
  matches old25ps history within1e-8; all-cell native/mirrored mass/momentum/energy RHS differences0.
  At25ns R_D76.03758ohm vs zero-Joule R_F77.88528ohm; drive effect -2.37233%,
  maximum absolute effect4.25890% over the window. Electrical feedback is resolved and suppresses
  the growth, but is not necessary for its dominant trend. Core50um QJ/heat/hydro outward sums
  are4.55580/50.4550/78.7468nJ in D; core total energy falls124.646nJ vs127.263nJ in F.
  Actual face budgets pass the1e-6 local exchange target; full-domain tiny boundary mass/free-energy
  ratios are roundoff-limited and explicitly not PASS. Hydro flux includes Rusanov numerical terms,
  not an isolated physical expansion-work percentage. Initial hot-core density is2.95855% of
  ambient under this EOS; no measured B0 density or cold-energy mapping is inferred.
  No extra thermal trajectory, temperature/current scan or physics edit is authorized after completion;
  prioritize sourced pre-breakdown/source energy and cold-to-thermal states for scientific review.

- **2026-10-16 RF-priority research decision, physical status unchanged:**
  [Literature/sensitivity review](../validation/i1_part_vi/spark_rf_priority_20261016/I1_SPARK_RESISTANCE_LITERATURE_AND_RF_PRIORITY_REPORT.md)
  stops extending the artificial10V/high-T/low-density LTE cylinder as real70um resistance.
  Preserve all accepted G1/FP3/FP4 and relaxation evidence. A post-collapse10V column voltage
  is not inherently impossible; actual event/source charge and thermal initialization are missing.
  Detailed thermal chemistry is not an unconditional prerequisite for frequency/damping work;
  quantitative energy partition still needs its separate physical inventories and boundaries.
  No new G1, EOS, cold-to-hot closure or Hc/M2 training is performed. Status NOT_RESOLVED unchanged.

## I1-CIR-07: Cold/thermal state to external gap-port partition

- **Issue ID:** `I1-CIR-07`.
- **Status:** `NOT_RESOLVED`.
- **Problem:** The Stage G3 reference port partitions `Cext`, `Cgap`, and streamer/channel branches, but a single physically verified terminal plane and bidirectional Stage C/G/actual-circuit coupling have not been established.
- **Evidence:** [Stage G3 port contract](stage_g3_transient_port.md); [Stage C3 current diagnostic](stage_c3_electrode_current_resistance.md); [V5-P2 final report](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md); [V5 entry decision](../validation/i1_part_v_preflight2/V5_ENTRY_DECISION.md).
- **Root cause:** **Confirmed:** G3's `I_port=I_L-I_Cext=I_sp+I_Cgap` is a lumped-model definition; Stage C's historical `i_total_hv` uses `q_A`, and same-plane equivalence has not been demonstrated. **Unproven:** whether current source, gap capacitance, and full-wave geometry can be partitioned without double counting in the actual fixture.
- **Impact:** Blocks Part VI physical circuit handoff and K2 terminal prediction; does not revoke G3 as a development reference or require unrelated Stage F work to stop.
- **Owner stage:** Part VI, with Stage G and Stage H port/field partition review.
- **Resolution gate:** Reviewed control-surface/sign convention, independently accounted `Cgap` versus `Cext`, consistent Stage C terminal observable, energy/current budget under a real waveform, and proof that downstream propagation representation does not duplicate gap capacitance.
- **Dependencies:** I1-PORT-03, I1-GEO-01, I1-BC-02; I1-THERM-06 if claiming a hot-channel low-order law.
- **Next authorized action:** Compare existing Stage C and G3 port definitions on paper and register the minimal Part VI interface proposal; no circuit coupling or code changes in this task.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-11.
- **Status history:** 2026-10-11: initialized `NOT_RESOLVED`; G3 development-model status is unchanged.
- **2026-10-16 frozen spatial conduction, status unchanged:**
  [Network and port qualification](../validation/i1_part_vi/spatial_bridge_20261016/spatial_bridge_gate_summary.json)
  reports two static real-snapshot solves plus one small cylinder reference. G_R rises9.5931%,
  while conditional G_frozen falls82.6300%; weak continuum-tail connectivity is not physical
  low-resistance bridging. Native center-gradient GR is NOT_COMPARABLE_OPERATOR_IDENTITY
  as a strict bound; only same-link admissible-psi trial energy supplies the discrete bound.
  FP3 stable_V2 acceptance remains intact but its2.5um/60fs state cannot be spliced into this
  1.25um/10ps history. No new I_port, Rphysical, dynamic feedback or RF is generated. Priority
  remains source-grounded spatial channel formation/maintenance and same-plane charge/current,
  not connecting1/G_R or1/G_frozen as an experimental spark resistance.
- **2026-10-16 conditional conductivity demand, status unchanged:**
  [Physical readiness](../validation/i1_part_vi/conductive_channel_bridge_20261016/physical_port_readiness.json)
  separates three passed numerical/algebraic checks from NOT_READY physical channel/RF.
  The actual masked Gaussian pre inventory is96.2174 electrons; conditional uniform70um
  drift columns at20/60/200ohm need3.74e9–4.19e10, not a measured resistance or proof of
  reachability. Three-state reaction ODE lacks spatial bridging, self-consistent E, transport/
  collection and physical port charge. DecisionB prioritizes a source-grounded spatial
  conducting-branch/port map and energy applicability, not longer constant-field extrapolation.
  FP3 local numerical acceptance is retained; no new circuit/RF/thermal run or experimental fit.
- **2026-10-16 experimental CSV evidence, status unchanged:**
  [Waveform review](../validation/stage_i/i1_waveform_review_20261016/I1_EXPERIMENTAL_WAVEFORM_REVIEW.md)
  registers actual root CSVs. CH3 voltage and CH2 current have different acquisition timestamps
  (825.970939s apart); no pointwise V/I, power, phase, KVL identification or physical formation
  history is computed. Probe factors, voltage reference plane, current branch and source state
  remain unknown. `PARAMETER_NOT_IDENTIFIABLE` is a data-qualification result, not an engineering
  Gate reversal. New raw data are available but not calibrated/synchronized port evidence.
- **2026-10-15 new scoped evidence, physical status unchanged:**
  [Thermal-series gates](../adapters/i1_fp4_handoff/thermal_series_run_20261015/execution/thermal_series_gate_summary.json)
  record `GIVEN_THERMAL_SERIES_COUPLING_ENGINEERING=PASS`. Synthetic10V ideal source,
  Rs3ohm/L100nH/**series** Cext101pF, independent1mm given LTE column;
  COLUMN_ONLY50/25ps and syntheticCp0.2pF TOTAL/SPLIT each complete10ns.
  This does not modify the old G2/G3 parallel topology. Column RMS, signed port/column/
  displacement currents and capacitor ownership are separate. All accepted B1/B2 native
  state/history hashes and common ledgers agree exactly. Peak-current/QJ/finalR/axisT cap
  sensitivity is0.14218%/0.08327%/0.02586%/0.01604%; maximum exchange-normalized interval
  energy residual1.44e-7, not normalized to the large initial thermal inventory.
  Full thermal/circuit trial rejection, controlled failure after native step, token guards
  and common commit have local engineering evidence; no native shortening occurred.
  Four trajectories use1000 accepted steps plus5 discarded native trials,18.67s science;
  14 isolated tests pass. Initial API-call failure had0 native steps and is retained.
  Existing given-adapter and cold FULL/stable_V2 acceptances remain separate. No global
  port, actual TENG, cold-to-hot mapping, electrode losses or RF gate is upgraded.
  Current action: stop for scientific review; do not repeat old synthetic or FULL tests.

- **2026-10-15 length/damping evidence, physical status unchanged:**
  [Damping metrics](../adapters/i1_fp4_handoff/channel_length_run_20261015/thermal_damping_metrics.json)
  show the frozen Cp=0 critical totalR62.931678ohm and critical length248.590294um for
  L100nH/seriesCext101pF/Rs3ohm. A70um hypothetical LTE column has initial frozen
  fd47.516077MHz, but actual thermal feedback produces only one reversal around10.91ns,
  with Rtot crossing the frozen critical value near17.85ns. No full repetitive cycle is
  observed; the time-varying system is not certified by instantaneous frozen poles.
  The nominal70um gap is not an effective-channel-length measurement. Source remains
  synthetic10V, not TENG; global port/parasitic/electrode definitions are not validated.
  Prior series coupling/ownership acceptance remains; no additional circuit is authorized.

- **2026-10-15 frozen TENG-source interface engineering, physical status unchanged:**
  [Source report](../adapters/i1_teng_source/results_20261015/TENG_SERIES_SOURCE_INTEGRATION_REPORT.md)
  and [source gates](../adapters/i1_teng_source/results_20261015/teng_source_gate_summary.json)
  accept only the new frozen-position V-Q-x implementation with explicit finite CT/Voc/qT0,
  separate seriesCext and portCp ownership, and joint source/port/circuit trial/commit.
  Four synthetic ratios and one1V capacitor-preserving event pass independent matrix-exponential,
  charge, source-coordinate energy and continuity tests. Max finest state errors0.05214%/0.01438%,
  min temporal order1.99926, max scaled accounting residual3.17e-14, ownership difference0.
  Twenty-two isolated checks pass; initial collection/roundoff-only assertion records remain.
  No thermal/FULL/PETSc/COMSOL/RF backend called and no historical circuit/physical kernel edited.
  SOURCE_PHYSICAL_CALIBRATION remains NOT_VALIDATED: actual CT(x0), Voc(x0), qT0 and
  high-frequency switch/rectifier topology must be sourced, not inferred from70um or101pF.
  Frozen mechanics is an assumption; dynamic mechanics is rejected, not silently frozen.
  A disconnected source would require a different topology, not the same Ceffective formula.
  Stop for scientific review; no automatic real-source/thermal/RF integration.

- **2026-10-15 observable/topology contract, physical status unchanged:**
  [Observation matrix](../adapters/i1_fp4_handoff/resistance_interpretation_run_20261015/signal_observability_and_topology_matrix.csv)
  separates Icol, actual loop current, differential gap voltage, series101pF voltage and Vrx.
  Existing R_used/pre and R_end/post/time-layer audit is offline only. The instantaneous
  ideal-source/Cp0 Rcrit crossing is not a nonlinear pole or a matched experiment result.
  Stage I contract is a no-data template; scope examples explicitly synthetic. Actual source
  CT/Voc/qT0 and high-frequency rectifier/isolation topology remain unknown, not filled with
  synthetic ratios,1V or nominal800V. First request actual waveform location and wiring;
  source-energy and physical-channel changes remain conditional and need review.
  Existing source and FP3/FP4 engineering acceptances are retained, not retested or upgraded.

- **2026-10-16 early driven/free comparison, physical port status unchanged:**
  The above causal review uses only the existing10V ideal-source/series101pF/Cp0 topology,
  with the original provider's actual column RMS and trial/commit. Free relaxation has no
  circuit or invented zero source-work ledger. No existing FP3/FP4/TENG acceptance is retested
  or revoked. Source CT/Voc/qT0, switch topology, same-plane channel current and experimental
  pre-breakdown state remain missing. There is no automatic real-source/thermal integration.

- **2026-10-16 conditional source/circuit sensitivity, physical status unchanged:**
  [Scoped machine status](../validation/i1_part_vi/spark_rf_priority_20261016/rf_resistance_priority_gate.json)
  records seven complete synthetic100ns events and one stopped S500/R20 attempt, not eight PASS.
  The stopped case is not restarted or voltage-scaled into a missing native curve. Remaining cases
  use registered event-based normalization in a new wrapper, without changing old protocol code,
  component parameters or targets. No measured CT/Voc/qT0/Cp or101pF precharge is inferred.
  Constant20/60/200ohm and prescribed20-to200ohm envelope are mathematical comparisons,
  not Toepler/RW calibration or G1 fits. Independent constant-R expm, local dynamic step checks
  and same-interval TOTAL/SPLIT port trial/commit accounting apply only to this synthetic scope.
  Physical port, actual TENG and classical-model benchmark remain unvalidated/deferred.

- **2026-10-16 bounded formation-port implementation, physical status NOT_RESOLVED unchanged:**
  [Formation report](../validation/i1_part_vi/spark_formation_20261016/SPARK_FORMATION_MODELS_ENGINEERING_REPORT.md)
  and [orthogonal judgments](../validation/i1_part_vi/spark_formation_20261016/spark_formation_gate_summary.json)
  accept formula/domain mathematics and the two synthetic1V/5ns coupled software fixtures only.
  Toepler Q[C] and RW S[A^2*s] are immutable trial/accepted histories, separate from x_V,
  TENG source charge and Cp charge; actual trial solves conduction and history together.
  Independent DOP853 reference maximum error2.555e-5 and saved-endpoint accounting1.56e-13;
  73 isolated checks include negative-current math, singular/time/unit guards and rejection/commit.
  Same-endpoint opposite-owner port replay agrees exactly, not an extra complete circuit event.
  Test K/R0/Q0/S0 are explicitly synthetic, not measured70um priors. Physical parameter
  identifiability remains NOT_CALIBRATED. Signed Toepler reversal is not replaced by abs(I);
  physical ports refuse reverse conduction/expired time. RW has no spontaneous recovery at zero
  current. No recovery law, time constant, floor or new physical inventory is implemented.
  Prior seven-case/PARTIAL failure record and all FP3/FP4/TENG physical limitations remain.
  Next action is scientific review with sourced branch/event/formation and recovery evidence;
  no further event, fit, source calibration or thermal model is authorized automatically.

## I1-RF-08: Native source, port, and receiver evidence separation

- **Issue ID:** `I1-RF-08`.
- **Status:** `DOMAIN_RESTRICTED`.
- **Problem:** A numerical Stage F native spectrum, Stage G port spectrum, and Stage H receiver spectrum have different source and transfer meanings. A developed 350 MHz system response cannot establish a trusted 350 MHz native streamer source or calibrated absolute received amplitude.
- **Evidence:** [scientific status](scientific_status.md); [v2.1 ultrafast decision](decisions/DEC-0001-v2.1-ultrafast-rf.md); [Stage H5 evidence chain](stage_h5_signal_evidence_chain.md); [Stage G3](stage_g3_transient_port.md).
- **Root cause:** **Confirmed:** current `NATIVE_RF_350MHZ=NOT_RESOLVED`; H4 absolute amplitude is `NUMERICAL_REFERENCE_ONLY`, Stage5 full-Maxwell reference remains pending, and gap-capacitance partition is required. **Unproven:** actual-electrode native-to-port-to-receiver attribution and experimental calibration.
- **Impact:** Restricts RF claims downstream of K2/Part VI; does not negate validated numerical processing in its documented trust band or a separately labeled system-development reference.
- **Owner stage:** Stage F/H, with Stage G port partition and Part VI source-interface review.
- **Resolution gate:** Source-to-port-to-receiver provenance on one reviewed geometry/reference plane; frequency-specific native trust analysis, nonduplicated transfer/capacitance partition, absolute receiver calibration where claimed, and independent experimental comparison. A 350 MHz system peak alone cannot close the native-source issue.
- **Dependencies:** I1-GEO-01, I1-PORT-03, I1-CIR-07; Stage F/H trust and measurement gates.
- **Next authorized action:** Preserve separate native/port/receiver labels and existing scientific-status flags in any I1 downstream documentation; no RF rerun or reinterpretation here.
- **First established:** 2026-10-11. **Last reviewed:** 2026-10-11.
- **Status history:** 2026-10-11: initialized `DOMAIN_RESTRICTED`; 350 MHz native and absolute receiver validation remain unresolved.
- **2026-10-16 raw waveform review, status unchanged:**
  [Experimental decision](../validation/stage_i/i1_waveform_review_20261016/EXPERIMENT_TO_MODEL_NEXT_DECISION.md)
  is `EXPERIMENTAL_WAVEFORM_FEATURES_READY_ONLY`. CH1 has fast signed oscillations but its
  measurement location is unknown; do not silently name it receiver voltage. Current is a
  separate acquisition with multiple disturbances, not one calibrated damped mode. CH3 slow
  ripple and CH1 fast oscillation remain separate observables. User-reported Tx voltage/current
  oscillation and additional Rx cycles are retained, not disproved or quantitatively established
  across unmatched events. No receiver power, physical spark R or RF validation is inferred.

- **2026-10-16 formation-only engineering, status unchanged:** no50MHz curve, spectrum,
  radiation power or received voltage is produced. Toepler/RW formation histories have no
  extinction/recovery/secondary-breakdown closure and cannot stand in for the full bipolar
  channel. Keep actual transmitter V/I, propagation/antenna response and additional receiver
  cycles distinct. Earlier synthetic sensitivity already establishes the priority of conduction
  survival; do not repeat arbitrary R(t) envelopes or infer experiment constants from ringdown.

- **2026-10-16 resistance/RF priority, status unchanged:** the prescribed synthetic source cases
  demonstrate sensitivity of current and dt*abs(FFT(di_loop/dt)) to resistance level/history.
  At800V, R20/R200 peak-current ratio4.20268 and band-proxy RMS ratio6.33897; these are
  not radiated-power or receiver-voltage ratios. Only R20 has at least three complete visible
  cycles; largest bins in no-cycle cases do not identify oscillation frequencies.100ns gives
 10MHz independent spacing, no zero padding. Stage F spatial J/current moment and Stage H
  geometry/receiver maps remain distinct; no frequency trust or experimental Gate is upgraded.
  Preserve user-reported transmitter V/I oscillation and additional receiver cycles without
  reverse-engineering R(t) from unavailable raw data. Next scientific review should constrain
  event topology/initial charge, first conduction state and survival/recovery, not auto-run scans.

- **2026-10-16 user observation clarification, status unchanged:** the user explicitly reports
  oscillation in transmitter-side voltage **and current**, with more cycles at the receiver.
  Retain this as USER_REPORTED, not receiver-only evidence. Raw waveforms, current branch,
  reference planes, bandwidth/loading and high-frequency topology have not been supplied.
  The artificial driven/free result cannot deny direct transmitter oscillation or identify its
  experimental cause. No RF processing, fitting or new signal simulation was performed.

- **2026-10-15 thermal-length signal identification, status DOMAIN_RESTRICTED unchanged:**
  The two100ns accepted-current signals have one resolved zero crossing and no complete
  repeating cycle. FFT uses4000 points/25ps, periodic Hann, mean removal, no padding;
  independent bins10MHz, Hann ENBW15MHz. `INSUFFICIENT_CYCLES` prevents a precise50MHz
  inference; the dominant10MHz bin is a finite-window transient diagnostic, not a measured
  oscillation. `EXPERIMENTAL_70UM_50MHZ_RF_VALIDATION=NOT_VALIDATED`. A thermal-column
  circuit current remains distinct from a spatial RF source and receiver signal. No RF run.

- **2026-10-15 conditional source-capacitance mode analysis, status unchanged:**
  In the selected connected series topology, synthetic CT/Cext ratios0.1/1/10/infinite
  give undamped f0 about166.095/70.823/52.524/50.0794MHz. These are conditional circuit
  eigenmodes, not measured TENG or RF peaks. Actual source capacitance, pre-event states,
  high-frequency topology and receiver transfer remain unknown;
  `EXPERIMENTAL_50MHZ_VALIDATION=NOT_VALIDATED`. The earlier70um thermal transient remains
  INSUFFICIENT_CYCLES and is not rerun or reclassified by this source-side calculation.

- **2026-10-15 offline experimental-interpretation boundary, status unchanged:**
  A matched measured multi-cycle loop waveform would imply
  `REFERENCE_EXPERIMENT_DISAGREEMENT_CONDITIONAL` for the given-LTE+10V reference, not
  automatically a G1 PDE numerical defect. Receiver-only ringing means
  `RECEIVER_RINGDOWN_NOT_DIRECT_CHANNEL_RESISTANCE_EVIDENCE`; it cannot uniquely invert
  R(t), nor does the numerical channel's lack of repeated cycles prohibit receiver ringing.
  No matched raw measurement is identified in reviewed I1/H/I sources, so actual experimental
  observable remains `OBSERVABLE_NOT_IDENTIFIED`. These are conditional rules, not claims
  that an experiment has been classified. No waveform fit, receiver/RF computation or new Gate.

## 2026-10-15 Part VI theory-only review / existing Issue cross-reference

The separately authorized [I1-PORT-1.0 contract](i1_part_vi_port_circuit_energy_contract.md),
[actual interface mapping](i1_part_vi_interface_mapping.md) and
[machine status](../validation/i1_part_vi/part_vi_theory_status.json) archive the complete
contract supplied in the user message. The separately named attachment file was not found;
this does not mean the supplied full message text is missing. No equation in that new
contract is promoted to an implemented or physically verified project result.
No Issue ID, resolution gate or scientific status is changed by this review.

| Existing Issue / owner | New documentation evidence and confirmed gap | Status and next authorized action |
|---|---|---|
| I1-GEO-01 / B,C,E | Existing geometry/mesh/reference-plane IDs are reusable; they do not certify the physical global terminal. | BLOCKER unchanged; review mapping only, no geometry work authorized. |
| I1-BC-02 / B,C,VI | Q_o, weighted F_o and P_o,in are distinct from unweighted particle loss. Local zero-Neumann does not resolve the global external field; historical C_A forces OpenCharge while C_B uses configured boundaries. | BLOCKER unchanged; no domain expansion or boundary replacement. |
| I1-PORT-03 / C,VI | q_A/q_B remain separate; native i_disp_hv contains space-charge response. Accepted native faces and absorption are reusable, but physical Q_H and a unified interval port ledger are not established. VI-9 must not add electrode absorption again. | BLOCKER unchanged; new theory and mapping documented, implementation not authorized. |
| I1-NUM-04 / C,V | Existing pairwise/local trajectory qualifications remain case-scoped; no new solves or re-evaluation. | BLOCKER unchanged for unqualified fitting; no new verification task. |
| I1-MAT-05 / C,V | Low-field and unsupported-species restrictions persist; the contract does not create ion transport. | DOMAIN_RESTRICTED unchanged; no coefficient work. |
| I1-THERM-06 / G1,VI | Field work, particle/chemical energy and gas heat are separate. Existing G0 diagnostics cannot initialize a calibrated thermal profile; cold port work cannot all become H_c. | NOT_RESOLVED unchanged; no thermal run or fit. |
| I1-CIR-07 / VI,G,H | TOTAL_CHARGE/SPLIT modes defined but not implemented; G2/G3 parallel Cext differs from the new user series-Cext instance. Stage C checkpoint resets electrode history and is not complete trial/rollback. | NOT_RESOLVED unchanged; review minimal adapter/transaction/energy gaps, no FP-3 implementation. |
| I1-RF-08 / F,H,VI | Stage F native electron face-derived center J, physical terminal I, and receiver voltage remain different. Maxwell displacement must not be injected twice. | DOMAIN_RESTRICTED unchanged; no RF/RLC or source reinterpretation. |

This dated entry records review of all eight existing Issues without creating duplicates
or closing any. C1's 500 V internal-output scope, fixed k_star=9134636280.4 s^-1,
exact closure NOT_PROVEN, C2-Q NOT_SELECTED, C2-S OBSERVATIONAL_CANDIDATE_ONLY,
600 V FIXED_K_VOLTAGE_TRANSFER_LIMITED and the position-case
KNOWN_SCOPE_LIMITATION / NO_FURTHER_COMPUTE_AUTHORIZED all remain unchanged.
V7 cold output comparison is completed; physical spark/arc models remain deferred,
not failed. Theory is defined and mapping documented; physical validation is pending.
The current authorized action ends with document review, not a new numerical task.

## Established results not requiring identical reruns

### 2026-10-15 FP-4 given thermal adapter acceptance (bounded engineering only)

[Actual report](../adapters/i1_fp4_handoff/given_thermal_run_20261015/FP4_GIVEN_THERMAL_ADAPTER_REPORT.md),
[input package](../adapters/i1_fp4_handoff/given_thermal_run_20261015/input/handoff_package.json) and
[series-coupling handoff](../adapters/i1_fp4_handoff/given_thermal_run_20261015/THERMAL_SERIES_COUPLING_HANDOFF.md)
are new evidence under I1-FP4-GIVEN-THERMAL-1.0, not a repeated FP3 test or a new research stage.
The earlier FP4 review/schema/readiness remain historical and unchanged.

- **I1-THERM-06:** standalone given-state adapter engineering accepted; no automatic thermalization,
  no velocity clipping, no PDE/EOS copy or constructor-state overwrite. The cold material-energy
  allocation, composition, profile and energy-zero gaps are unchanged; thermal ROM stays NOT_RESOLVED.
- **I1-CIR-07:** I_col remains the prescribed column-conduction branch. The adapter rejects total
  port/inductor-current substitution, checkpoint history and old-G2-as-new-series requests. No
  series thermal PortProvider or complete thermal trial/commit is implemented or tested here.
- **I1-PORT-03 / I1-BC-02 / I1-GEO-01 / I1-NUM-04 / I1-MAT-05 / I1-RF-08:** no status upgrade.
  Two same-algorithm paths establish wrapper equivalence, not independent physical validation,
  new thermal convergence, global port, low-field, experiment or RF readiness.

Stop after the accepted independent-reference engineering task. Next action is scientific review
of the minimum thermal-series interface list; no automatic coupled run or new fit.

### 2026-10-15 FP-4 cold/thermal contract review (documentation only)

[Handoff review](../adapters/i1_fp4_handoff/FP4_COLD_THERMAL_HANDOFF_REVIEW.md),
[field matrix](../adapters/i1_fp4_handoff/cold_to_thermal_field_map.csv) and
[declarative schema](../adapters/i1_fp4_handoff/cold_thermal_handoff_schema.json)
map actual FULL/G0/G1/G2 sources without calling or importing a solver.
Latest `FULL_LOCAL_NUMERICAL_COUPLING_ACCEPTED` is retained and frozen at its local
serial SP3-OFF scope; no synthetic R/C, transaction or weak-quantity rerun is needed.
Old Gate values remain historical, not overwritten.

- **I1-CIR-07, NOT_RESOLVED unchanged:** native local port charge/work and FULL rollback
  exist; they are not thermal material energy or a cross-model switch. G1/G2 fresh
  constructors reset time and integrals; future handoff needs a global clock anchor,
  prior ledgers, physical capacitor charge/current continuity and explicit topology.
  Old G2 parallel Cext and FP3 series Cext remain different instances. I_col is column
  conduction (RMS for G2 intervals), not automatically i_port or electrode collection.
- **I1-THERM-06, NOT_RESOLVED unchanged:** the two entry decisions and precise missing
  relationships are recorded above; the standalone given-state path can be scoped
  independently without fabricating automatic heating or waiting for hot ROM closure.
- **I1-PORT-03 / I1-BC-02 / I1-GEO-01:** local field energy and remaining energy cannot
  be renamed global port/material energy; no compensating C_global-C_local. Same-domain
  electrical precision does not close external-domain/reference-plane mapping.
- **I1-MAT-05 / I1-NUM-04 / I1-RF-08:** no new property, convergence, RF, source or
  experiment evidence. Full material budgets and missing losses stay unavailable,
  not zero and not inferred from residuals. C1/V5/V6/V7 statuses unchanged.

Next action is scientific review of the explicit GivenThermalInitialStateAdapter
proposal; no implementation or computation follows this document task automatically.

### 2026-10-15 FP-3 small-signal numerical interface evidence (no new Issue)

Explicit I1-FP3-SMALL-SIGNAL-1.0 authorization produced the
[weak-port report](../adapters/i1_fp3_port_circuit/weak_port_run_20261015/FP3_WEAK_PORT_ACCURACY_REPORT.md)
and [new scoped gates](../adapters/i1_fp3_port_circuit/weak_port_run_20261015/fp3_weak_port_gate_summary.json).
Decision: **FULL_LOCAL_NUMERICAL_COUPLING_ACCEPTED**, not global physical acceptance.
Historical FULL_PROVIDER_ACCEPTED_COUPLING_PENDING and its original metrics remain intact.

- **I1-CIR-07 / Part VI,G,H:** local increment-based KVL/current commit now resolves the
  old D*rq amplification. Actual FULL feedback, TOTAL/SPLIT each three20fs macros,
  yields max incremental |RV|2.843e-14V (committed rounded endpoints3.168e-12V),
  |Rq|/Qtr6.083e-9, circuit energy/work2.038e-7.
  Native state/history/total-ledger equivalence is exact; repeat/failure/retry restore
  remains exact, including the actual wrapper's cumulative metadata. The limited
  serial local engineering subtask is accepted, not unimplemented. Overall Issue
  remains NOT_RESOLVED for physical port/source/energy integration.
- **I1-NUM-04 / C,V, BLOCKER unchanged for broader accuracy claims:** stable same-network
  source projection and independently solved zero-voltage response agree; finite
  total/weight residual contractions explain the observed error source without an
  invented current or loss. Weak charge/Qtr1.59e-15 and electrical energy/work3.50e-7
  pass the new1e-3 targets on accepted intervals. Finite weight-energy terms and
  negative remaining offsets are retained. This is not spatial/temporal continuous
  PDE certification. Runtime rtol1e-14 changes no assembly or physical formula.
- **I1-PORT-03 / C,VI, BLOCKER unchanged:** explicitly named stable_V2 is a numerical
  evaluation of the same active-network port in the exact-solve limit; q_A,
  native-links q_B and active_V1 remain recorded. No global/local capacitance
  subtraction or Cp=2U/V^2 replacement; no global physical charge claim.
- **I1-GEO-01 / I1-BC-02 / I1-MAT-05 / I1-THERM-06 / I1-RF-08:** unchanged.
  No new geometry, field-domain, low-field, gas heating, TENG or RF evidence.

Resource/evidence:24 native calls,14 explicit static RHS, about45.83s scientific
execution;16 isolated read-only checks. First static-only comparison used differing
Cp energy scopes and stopped with zero steps; corrected at fixed Cp using existing
references, then reused the archived weight/operator. Original incident retained;
no repeated reference/PDE trajectory. No historical CTest or synthetic matrix rerun.
Old artifact-integrity UNKNOWN bytes remain UNKNOWN. This bounded local engineering
task is closed; next action is scientific review, not automatic additional computation.

### 2026-10-15 native FULL transaction/local coupling evidence (existing Issues only)

New explicit I1-FP3-FULL-1.0 authorization produced the
[native FULL report](../adapters/i1_fp3_port_circuit/native_full_run_20261015/FP3_NATIVE_FULL_PROVIDER_REPORT.md),
[gate record](../adapters/i1_fp3_port_circuit/native_full_run_20261015/fp3_full_provider_gate_summary.json)
and [state manifest](../adapters/i1_fp3_port_circuit/native_full_run_20261015/native_state_manifest.json).
It does not rewrite prior theory-only, synthetic-only or historical physics Gate values.

- **I1-CIR-07 / Part VI,G,H, NOT_RESOLVED unchanged:** actual native FULL transactions
  and local series-circuit trials are implemented. All saved fields/histories restore
  exactly; direct one-step replay is exact, two-step independent recomputation differs
  by at most1.187e-15 relative with equal histories. Controlled failure after actual
  step, stale/foreign/duplicate rejection and observer ON/OFF checks pass. There are
  47 native calls including rejected/reference work, and two modes each complete three
  20fs macros with actual charge-root evaluation and unique joint commit. TOTAL/SPLIT
  fields and total ledgers differ by0. This is not an unimplemented smoke; numerical
  qualification is conditional. Root residual <=5.331e-6 of2e-16C meets the new gate,
  but weak-current/energy accuracy is not certified under L/dt² amplification.
- **I1-PORT-03 / C,VI, BLOCKER unchanged for physical port:** q_A and native-links q_B
  remain unchanged. The explicit q_B_active_operator is derived from the same effective
  zero-Neumann network, omitting replaced gas rows. Independent native face-weight and
  actual correction/reaction ledgers are exported. Constant500V charge residual is
  0.0165% of unsigned transport; coupled accepted steps reach31.13% of that scale.
  Endpoint charge remains useful as a qualified local numerical candidate, not a
  resolved physical plasma-current or full-electrode terminal reference.
- **I1-NUM-04 / C,V, BLOCKER unchanged for unqualified claims:** exact transaction
  restore is not a new PDE accuracy proof. Constant500V electrical residual/actual
  exchange reaches0.5525%; coupled residual/conduction work reaches10.43 times despite
  a small fraction of capacitor-charging work. Finite Poisson/weight accuracy and
  small differences remain unresolved; no more accurate same-grid reference was run.
  Negative SPLIT remaining energy around-3.285e-20J is retained and flagged, not clipped
  or removed by redefining capacitance. The old synthetic energy/order evidence is intact.
- **I1-GEO-01 / I1-BC-02:** statuses unchanged; the active local network is not global
  B0 validation. Qo=0 is justified only for that reduced homogeneous-Neumann model;
  excluded boundary-row rho change is0 in the tested windows, not assumed for all cases.
- **I1-MAT-05 / I1-THERM-06 / I1-RF-08:** statuses unchanged; SP3 ON, new ion transport,
  low-field physical validation, gas heat, TENG, cold-hot and RF are not supplied.

Decision: **FULL_PROVIDER_ACCEPTED_COUPLING_PENDING**. The delivered provider and actual
bounded feedback exist; precision admission of weak charge/electrical work remains pending.
Next authorized action after this task: stop for scientific review. The minimal remaining
numerical question is matched same-grid endpoint charge/energy accuracy versus external
root-residual amplification; do not automatically change tolerance, run longer, add states,
repeat the old R/C matrix or create another research phase. Synthetic engineering acceptance,
C1 scope/k, position limitation and all V5/V6/V7 results remain unchanged.

### 2026-10-15 FP-3 synthetic engineering evidence (no physical Issue closure)

New explicit user authorization allowed the independent
[FP-3 synthetic engineering implementation](../adapters/i1_fp3_port_circuit/results_20261015_final/FP3_PORT_CIRCUIT_ENGINEERING_REPORT.md)
and [gate record](../adapters/i1_fp3_port_circuit/results_20261015_final/fp3_engineering_gate_summary.json).
This does not rewrite the previous Part VI theory-only authorization or any historical Gate.

- **I1-CIR-07 (Part VI/G/H), status NOT_RESOLVED unchanged:** the synthetic series-Cext
  adapter now implements TOTAL/SPLIT ownership, actual port-trial residual assembly,
  immutable accepted state, joint commit and rollback/reject behavior. Final 25 isolated
  tests and 20 fixed accounting-mode/time-step realizations pass. Maximum scaled algebra
  residual 2.36535e-14; cumulative energy residual 1.29887e-14; ownership differences 0;
  finest continuous-expm reference RMS error 0.0114583%; order 1.99479–2.00019.
  These are SYNTHETIC tests, not actual gap capacitances or plasma trajectories.
  An initial roundoff-only endpoint-step failure is retained and its driver repair is
  separately documented; thresholds and component parameters are unchanged.
- **I1-PORT-03 (C/VI), status BLOCKER unchanged:** the synthetic charge/current/energy
  protocol is implemented; q_A/q_B/physical Q_H are not redefined or replaced.
  Actual FULL transaction integrity is NOT_TESTED; plasma-circuit coupling is
  NOT_IMPLEMENTED. Synthetic port acceptance does not await or manufacture physical evidence.
- **I1-NUM-04 (C/V), status BLOCKER unchanged for unqualified fitting:** new synthetic
  midpoint order and energy evidence applies only to the specified non-singular linear
  R/C ODEs. No DAE/PDE accuracy, old trajectory, low-field or V5 Gate is upgraded.
- **I1-THERM-06 / I1-RF-08, statuses unchanged:** no gas heating, cold-hot, TENG or RF
  connection. C1's new refusal wrapper explicitly denies physical current/heat/feedback
  and lifting, and restricts the existing 500 V/10 ps supported centered points.
- **I1-GEO-01 / I1-BC-02 / I1-MAT-05:** unchanged; no geometry or material work.

The synthetic skeleton is `FP3_SYNTHETIC_PORT_ENGINEERING_ACCEPTED` and this bounded
engineering task is closed; do not reopen the same R/C matrix because actual data are pending.
Next action is scientific review of [FULL_PROVIDER_HANDOFF](../adapters/i1_fp3_port_circuit/FULL_PROVIDER_HANDOFF.md).
Native FULL trial/commit/rollback and physical-reference interval budgets require a separate
authorization. No new Issue, state, ROM fit or next-stage computation is authorized here.

- Part V-A established the stated analytical and independent finite-volume weighted-budget identities and an exact native-step algebraic budget in its tested scope. Reuse its tests/evidence; this is **not** K2 closure or physical convergence.
- V5-P established the native low-field branch behavior and a time-step-versus-spatial-sensitivity contrast for its development case. Recheck applicability only for a new operating domain.
- V5-P2 attributed the old `q_A`/Laplace mismatch, and independently established same-operator `q_B` reciprocity, energy conjugacy, and a read-only source-response current budget. Do **not** repeat those checks merely to relabel `q_B` as physical.
- The 65/70/75 um definitions are reconciled for the historical geometry. Retain the 65 um case as an unchanged regression fixture for parameterized gap checks; generic geometry/boundary/port gates and actual-electrode claims remain open.
- Part IV-R did **not** authorize freezing the two-state hot model. Its negative/limited evidence is retained, not erased by cold-stage progress.

### I1-COLD-PHYSICS-V3 implementation evidence (2026-10-19 recovery)

[V3 report](../validation/i1_cold_physics_v3/implementation_20261016/I1_COLD_PHYSICS_V3_IMPLEMENTATION_REPORT.md)
and [capabilities](../validation/i1_cold_physics_v3/implementation_20261016/i1_cold_physics_v3_capabilities.json)
record newly authorized native ion transport, actual ion-impact wall emission,
gas-only coupled SP3 metal Robin boundaries, all-enabled-species local budgets and
the serial COLD_FULL_B0 scenario/provider. Historical evidence above is unchanged.

- **I1-BC-02 / I1-PORT-03:** internal metal photon boundary implementation and local
  species/charge accounting pass their scoped engineering checks. Global boundary,
  physical reference plane and physical port remain open. Raw total-potential work
  residual reaches 19.24% of weak face conduction work; no new weak-energy PASS.
- **I1-NUM-04:** 36 native checks, 23 scenario/ownership checks and independent SP3
  direct-solve comparison pass. The one B0 run has complete paired evidence only
  through 1 ps/50 commits, then interruption with no normal exit record. Additional
  wall fragments are not used as complete intervals. B0 execution CONDITIONAL_PASS;
  no second trajectory or fabricated final resource counter. Source revisions and
  post-interruption logging/guard fixes are archived separately.
- **I1-MAT-05:** ion mobility/diffusion, ion-impact yield and photon excitation retain
  synthetic/conditional provenance. Missing emission energy, metal optical parameters
  and negative-ion transport are explicit, not default zeros or calibrated physics.
- **I1-THERM-06 / I1-CIR-07 / I1-RF-08:** no automatic thermal conversion, new RLC/TENG
  event or validated RF source/receiver connection. Existing engineering Gates stand.

Next dependency is a source- and domain-qualified ion/surface parameter package and
corresponding port capability constraints, with the field-work precision limitation
kept separate. No new Issue, scan or follow-up simulation is opened by this archive.

### I1-COLD-PHYSICS-V3.1 evidence (2026-10-19)

[V3.1 report](../validation/i1_cold_physics_v3_1/implementation_20261019/I1_COLD_PHYSICS_V3_1_PARAMETER_AND_RECOVERY_REPORT.md)
adds scoped evidence to existing issues; no historical Gate is rewritten.

- **I1-MAT-05:** versioned ion/surface source/domain contract and semantic refusal implemented.
  Current mobility, diffusion and gamma remain class C synthetic fixtures; actual ion composition,
  surface/incident energy and high-field calibration remain unknown. Real-device mode rejects
  missing evidence or synthetic defaults. Contract acceptance is not parameter calibration.
- **I1-NUM-04 / I1-PORT-03:** same-build serial accepted-state persistence passes a real
  process restart: checkpoint 20 fs, two subsequent steps to 60 fs, 170 comparisons bitwise
  identical; five B0 calls total. Legacy small-probe results match old/new libraries bitwise.
  Atomic accepted publication, source/version/kind/checksum/output-prefix refusals are tested.
  This is not distributed, cross-build or arbitrary external-circuit recovery certification.
- **I1-PORT-03:** offline attribution of the old 0--1 ps evidence is complete on available
  fields. The 19.2373% raw weak-work residual remains unqualified. Its signed correspondence
  with V1/V2 charge evaluation error is documented, not promoted to an independently proved
  solver-residual cause; spatial residual/correction work arrays were not archived.
- **I1-BC-02 / I1-THERM-06 / I1-CIR-07 / I1-RF-08:** global port, emission material energy,
  automatic cold-hot, new circuit physics and RF validation remain unchanged. No long run,
  new physical scan or history reconstruction was started.

Next scientific dependency remains source-qualified high-field ion and material-surface inputs;
the weak electrical-work attribution gap is separate and does not revoke the recovery result.

### I1-COLD-PHYSICS-V3.2 evidence (2026-10-19)

[Implementation](../validation/i1_cold_physics_v3_2/implementation_20261019/I1_COLD_PHYSICS_V3_2_IMPLEMENTATION_REPORT.md)
adds no new Issue and rewrites no historical Gate.

- **I1-MAT-05:** native source-bound EN interpolation/reduced-mobility conversion and explicit
  surface provenance are implemented. No source-backed high-field ion table was obtained;
  actual mixture/surface calibration remains missing. C-class manufactured table is software
  evidence only. Unequal parallel/perpendicular diffusion refuses the unsupported tensor FV path.
- **I1-PORT-03 / I1-NUM-04:** five native calls cover three unique B0 steps through 60 fs,
  including accepted checkpoint/new-process continuation. State/flux/history/ledger replay is
  bitwise identical. New-window independent residual contractions and zero-voltage charge-response
  work achieve the registered weak-work target (maximum signed-work ratio about 3.00e-11).
  Old V3 19.2% evidence remains unchanged; this does not certify material heating or a global port.
- **I1-RF-08 / I1-BC-02:** actual Stage F reader and current-moment routines consume electron
  plus positive-ion accepted FV source, with original faces retained and no duplicate displacement
  or emission current. Global electrode/lead return sources, outer-domain response, propagation
  and receiver validation remain open. Gas-only current moment is not receiver voltage.
- **I1-THERM-06 / I1-CIR-07:** no new thermal, TENG, circuit-feedback or empirical recovery run.
  Source parameter calibration and automatic cold-hot closure remain separate unresolved work.

### I1-COLD-PHYSICS-V3.3 evidence (2026-10-19)

[Implementation](../validation/i1_cold_physics_v3_3/implementation_20261019/I1_V3_3_GLOBAL_SOURCE_AND_RF_INTERFACE_REPORT.md)
adds no Issue and changes no historical Gate.

- **I1-PORT-03 / I1-NUM-04:** versioned conservative graph/source packet and explicit electrode-storage,
  lead-exchange and surface-mapping checks are implemented. Archived native 40--60fs faces cancel
  pairwise; independent full-gas correction/source ledger closes to about 2.30e-15 of unsigned
  cell-face exchange. Cellwise correction/source arrays are absent: that qualification remains
  conditional, not reconstructed from a residual. No new native step or port algorithm change.
- **I1-RF-08:** two finite manufactured sources actually traverse the existing Stage F time-domain
  Jefimenko kernel and Stage H field projection. Independent compact electric/magnetic dipole
  comparisons are below the registered 1% target (worst about 0.803%). This validates a bounded
  source-interface example, not experimental RF, receiver voltage or arbitrary near-wire convergence.
- **I1-GEO-01 / I1-BC-02 / I1-CIR-07:** real B0 remains LOCAL_GAS_ONLY. Missing external leads,
  return geometry, global surface-charge distribution and matching reference-plane inventories
  explicitly refuse full-source assembly. Same current does not specify radiating geometry;
  closed-loop integral J can vanish while its magnetic moment and directional field are nonzero.
- **I1-MAT-05 / I1-THERM-06:** source calibration, surface energy and automatic cold-hot closure
  unchanged. No thermal, TENG, fullwave or resistance-model run was performed.

Next dependency is a physically specified same-reference-plane electrode/lead/return and surface
charge map, with source/switch topology. Stop this task; no automatic follow-up simulation.

### I1-COLD-PHYSICS-V3.4 evidence (2026-10-19)

[Implementation](../validation/i1_cold_physics_v3_4/implementation_20261019/I1_COLD_PHYSICS_V3_4_IMPLEMENTATION_REPORT.md)
adds evidence to existing Issues only; no historical Gate is rewritten.

- **I1-PORT-03 / I1-NUM-04:** actual reaction RHS and raw/pre/post species plus floor/clip arrays now
  accompany accepted net electron and positive-ion faces. Three new B0 steps to 60 fs close every
  gas-cell charge budget: maximum local unsigned-exchange ratio 4.38e-14, global L1 ratio 1.08e-14.
  Original V3.3 missing-cell-data qualification remains CONDITIONAL_PASS. No transport/chemistry change.
- **I1-GEO-01 / I1-BC-02:** native local H/Return surfaces and source-bound profile parameters are
  mapped. Schur interface preserves simultaneous potential, charge functional and energy including
  nonzero outer-domain charge. Its manufactured operator is not a B0 global operator. Actual lead,
  return, terminal selection and surface-charge distribution remain unavailable.
- **I1-RF-08 / I1-CIR-07:** GlobalSourcePacket 1.1 checks explicit finite surface patches, conductor
  currents, contacts and return conservation. Existing source consumers actually read new native
  and manufactured packets; no propagation rerun. Source labels cannot grant global physical or
  receiver qualification, and total port current is not copied into spatial J.
- **I1-MAT-05 / I1-THERM-06:** C-class ion/surface inputs, material energy and cold-hot closure
  retain their previous restrictions. No thermal/source scan or physical calibration.

Next required input is the actual same-reference-plane external geometry and its boundary/surface
response. The software interface is available; a real device global source remains NOT_READY.

### I1-COLD-PHYSICS-V3.5 evidence (2026-10-20)

[Implementation](../validation/i1_cold_physics_v3_5/implementation_20261020/I1_V3_5_MULTIPHASE_ARCHITECTURE_REPORT.md)
adds evidence to the existing Issues, without rewriting historical Gates.

- **I1-GEO-01 / I1-BC-02:** source-bound physical-assembly input contract reuses B0 local surface mapping.
  Complete synthetic geometry and series connectivity are checked; real external leads/return, probes,
  source state, boundary response and surface distributions remain unknown. A Chinese measurement guide
  defines the next input package. Geometry metadata does not certify a global physical source.
- **I1-THERM-06 / I1-CIR-07:** one unified entry actually calls cold FULL (3 accepted steps to 60 fs)
  and the existing independent given-LTE thermal/series provider (4 accepted steps to 0.2 ns).
  No new PDE/EOS/circuit implementation; initial 1 mm reference is not the actual 70 um channel.
  Thermal Joule uses column RMS, not total current; local/anchor time and separate inventories are explicit.
  AUTO_COLD_TO_THERMAL still refuses with MISSING_PHYSICAL_CLOSURE, without constructing G1.
- **I1-PORT-03 / I1-NUM-04 / I1-RF-08:** native conservative species source reaches existing Stage F
  readers and current-moment diagnostics. Thermal Jz is a conditional axial-column profile, not a full
  charge/lead/return source. No propagation or RF receiver test; all physical parameter restrictions remain.

Next dependency is the source-bound actual assembly/topology/measurement package at matching H/Return
planes. No further thermal trajectories, physical source inference or automatic stage transition is authorized here.

### I1-P1 evidence (2026-10-20)

[Source/device implementation](../validation/i1_p1/implementation_20261020/I1_P1_DEVICE_PARAMETERIZATION_AND_SOURCE_REPORT.md)
adds evidence to existing Issues only; all prior Gates remain unchanged.

- **I1-GEO-01 / I1-BC-02:** device package preserves nominal B0 and separate conditional measurement
  metadata. R1/DAG1 units/meaning and D2 remain unverified; no new native geometry or guessed return.
- **I1-CIR-07 / I1-PORT-03:** ideal contact-separation geometry feeds the existing frozen TENG charge
  source and FP3 series transactions. Alternative amplifier topology keeps Rbleed parallel to Cext,
  not equal to Rs; unknown impedance is explicitly an ideal limit. Two 10 ns synthetic references
  pass independent matrix-exponential and energy/charge checks; 33 isolated checks pass.
  Source coordinates and Cp ownership remain separate; no device source calibration.
- **I1-THERM-06 / I1-MAT-05 / I1-RF-08:** dielectric and sigma inputs are declared synthetic because
  the named parameter assessment file was not found. Existing cold/thermal modes accept shared
  nominal identity but reject unsupported TENG requests; no new cold, thermal or RF execution.
  True precharge, material/thermal mapping and global source qualification remain unresolved.

Next dependency is one source-bound pre-breakdown input package with measured topology and consistent
source/Cext initial states, not another artificial resistance or frequency-matching scan.

### I1-E2E-V1 evidence (2026-10-20)

[Implementation](../validation/i1_e2e/v1_20261020/I1_E2E_V1_IMPLEMENTATION_REPORT.md) adds evidence only.

- **I1-CIR-07 / I1-RF-08:** one entry actually connects frozen geometric TENG, original series
  transactions, conservative compact source, existing F/H and a conditional effective-height receiver.
  A 100ns mathematical channel reference completes the chain. Its initial fast RC voltage is
  underresolved (scaled error 0.818), even though late analytic field agreement is 0.051%.
  Startup/antenna/device accuracy is not accepted. No experimental RF Gate upgrade.
- **I1-PORT-03 / I1-NUM-04 / I1-MAT-05:** new axial FV species/Poisson ROM derives port charge from
  same-network electrode storage and actual charged-particle faces. Two 200ps R comparisons and
  one gap comparison genuinely evolve particle states; local charge/work checks pass. Weak
  resistance feedback has no temporal-convergence bound; low-field coefficients and C-class
  ion/surface inputs retain their qualifications. No GR inversion or prescribed-resistance ROM.
- **I1-GEO-01 / I1-BC-02 / I1-THERM-06:** effective area/planar axial geometry and finite equivalent
  radiation conductors are explicit abstractions, not device CAD. Old cold/thermal providers retain
  independent modes and share a result contract; automatic material-energy/state handoff remains missing.

Priority theory is conservative FULL-to-ROM species/flux/charge/work mapping, not another arbitrary
R(t) fit. Existing physical Gates remain unchanged; all four authorized cases are finished.

There are eight initial issues. Current `BLOCKER` issues are I1-GEO-01, I1-BC-02, I1-PORT-03, and I1-NUM-04, each with the scoped impact above. No experimental-science gate is marked verified by this ledger.
