# I1 parameterized geometry, electrostatic boundary, and dynamic port interface contract

Scope: a reusable interface for **supported parameterized electrode geometry families**, their electrostatic boundary problems, and their dynamic terminal observables. It is not a port definition specialized to the historical needle-plane case, nor an experimentally validated port law. Established: 2026-10-11. Last reviewed: 2026-10-11. Maintained with [I1 issue ledger](i1_open_issues_and_gates.md). The historical case below is a regression fixture and evidence of current limitations, not a universal geometry template. This document does not change a frozen solver, observable, or reference result.

## 1. Parameterized geometry and reference surfaces

Each supported geometry instance must identify its geometry family and revision, dimensionality/coordinate system, HV and return electrode solids and named terminal surfaces, gas domain, material/permittivity regions, units, gap-defining facing surfaces, and source/provenance of its parameters. It must declare its static reference voltage and any geometric approximation. Classify axisymmetric 2D-exact, explicitly justified 2D-approximate, or 3D-required routing before any streamer handoff; a nonaxisymmetric structure is not silently promoted to a physical Stage C result. The geometry-to-mesh mapping must record electrode masks/faces and the difference between designed surfaces and discrete surfaces. An unsupported or ambiguous mapping is `GEOMETRY_MAPPING_UNRESOLVED`, not silently assigned a needle-plane gap.

For a general geometry, `d_net` is the **minimum gas-connected separation of the declared facing metal surfaces within the declared gap region**, with the selected surface pair, region, and distance method recorded. If several physically distinct gaps or paths compete, record each or mark the scalar `d_net` unresolved; a center-coordinate distance is not a substitute. Any analytic formula is family-specific. The parameterized interface must compute `d_net` from the solids, validate that it is positive and finite, check electrode nonoverlap and units, and independently compare the discrete gas-face separation against the continuous definition. Validation must be automated for each supported family and must not silently reinterpret legacy metadata.

### Historical needle-plane regression fixture

For the historical axisymmetric needle-plane case, define the continuous gas gap between the metal surfaces, not between input coordinates:

\[
d_{\mathrm{net}}=z_{\mathrm{tip,center}}-r_{\mathrm{tip}}-z_{\mathrm{ground,top}}.
\]

The existing geometry has tip center at 75 um, spherical tip radius 5 um, and ground upper surface at 5 um: **d_net = 65 um**. The Stage C development driver's `gas_gap` is 70 um because it subtracts the ground top from the *tip center*, without subtracting tip radius. Geometry metadata `tip_z - ground_z` is 75 um. These are distinct definitions, not three measurements of the same gap. A cell-center electrode mask adds mesh-dependent face placement; it does not explain the definitional 5 um difference. Preserve all three historical values and the original outputs as a regression test of definitions and automated checks; do not rewrite the old reports. See [electrode.cpp](../cpp/src/electrode.cpp), [Stage C3](stage_c3_electrode_current_resistance.md), and [V5-P2](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), section 5. None is an experimentally verified actual-electrode geometry.

## 2. Configurable electrostatic boundary problem

The geometry instance must select and record a **boundary class** for every nonterminal exterior surface: for example, a modeled grounded enclosure/return conductor, a justified remote truncation with an explicit far-field condition, or a physical symmetry boundary. Insulating or periodic conditions are admissible only where the geometry and physics justify them. The same conductor identities, permittivity map, and compatible outer-boundary class must be used in the physical electrostatic problem and the 1 V weighting problem. For each class, document its physical meaning, mathematical condition, reference potential, domain extent, and how boundary charge/flux enters `Q_outer`; no class is declared universally correct. A metal/exterior-boundary intersection must be explicitly modeled and checked, not hidden by a cell mask.

Boundary verification requires a well-posed potential reference, independent outer-boundary contribution, and representative domain/mesh-sensitivity checks for field and conjugate charge. Record separately a local discharge ROI and global field/capacitance metrics; a localized convergence result cannot upgrade a global or dynamic one. The historical HV shank meets an outer gas Dirichlet boundary with a potential jump. Its global maximum weighting field occurs near the shank/domain exit, not necessarily the needle tip. This is a **regression case exposing a boundary-class failure**, not an outer-boundary prescription for new geometries. See [V5-P2](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), sections 3 and 5, and [Poisson implementation](../cpp/src/poisson.cpp).

## 3. Dynamic terminal sign and reference plane

For each supported geometry and boundary class, choose and record a fixed HV terminal, return terminal, control surface, and outward normal. Define

\[
V_p=V_{\mathrm{HV}}-V_{\mathrm{return}}.
\]

`I_p` is positive when **conventional current enters the HV terminal from the external circuit**. `Q_HV` is positive charge on that HV terminal in the same reference plane. Under Stage C's existing convention, `I_cond,H` is positive for conventional current from the HV electrode into the gas (electron absorption at HV gives a positive contribution). With a consistent control surface, species scope, and time layer, the candidate balance is

\[
I_p=\dot Q_{\mathrm{HV}}+I_{\mathrm{cond},H}.
\]

For the *historical native diagnostic*, `i_cond_hv = e` times the electron absorption number rate at HV, `i_disp_hv = Delta q_hv / Delta t`, and `i_total_hv = i_cond_hv + i_disp_hv`. The name `disp` does **not** imply that this term is only `C_geom dV_p/dt`. A geometry-family port interface must export charge, electrode absorption, outer-boundary contribution, numerical corrections, and geometry-capacitance contribution with explicit units, sign, surface, and time-layer provenance. It must demonstrate its accepted-step current budget without adding a full weighting-current term or a second capacitance current to an already complete same-surface balance. No ion transport current is computed in this Stage C model. See [Stage C3](stage_c3_electrode_current_resistance.md), [V5-A port audit](../validation/i1_part_v/port_mapping_audit.md), and [V5-P2](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), section 4.

Stage G3 defines a lumped gap port with `V_port=V_gap` and `I_port=I_L-I_Cext=I_sp+I_Cgap`: the external `Cext` branch lies outside the gap port, while `Cgap` lies inside it. This is a reference-model convention, **not** a verified direct equality between Stage C's historical `i_total_hv` and G3's port current. Full-wave geometry that already represents the gap capacitance must not double-count it. See [Stage G3](stage_g3_transient_port.md), including its `H_GAP_CAPACITANCE_PARTITION_REQUIRED` constraint.

## 4. Candidate conjugate charge operator and weighting relation

For **each fixed geometry instance** with linear dielectric permittivity, a consistent pair of electrode reference boundaries, and explicitly accounted outer boundaries, define the unit-voltage Laplace weighting potential `psi` (`psi_HV=1`, `psi_return=0`). It must be formed from that instance's *same discrete Poisson operator*, not borrowed from another geometry or boundary class. The *candidate* relation is

\[
Q_{\mathrm{HV}}=C_{\mathrm{geom}}V_p-\int_{\Omega_g}\psi\rho\,dV+Q_{\mathrm{outer}}.
\]

`Q_outer` denotes an **independently evaluated** outer-boundary contribution under the chosen discrete/continuous problem, not a name assigned to an unexplained residual. The volume and boundary partitions, normal, charge sign, and time layer must match the selected terminal observable. A changing `Q_outer` can contribute current even when particle flux through the outer boundary is zero. This equation is neither a universal identity for the native center-field charge nor experimentally validated. Its continuum interpretation and physical terminal plane remain conditional on geometry and boundary convergence. For a supported family, the candidate conjugate operator maps that geometry's discrete potential and electrode links to terminal charge; its static reciprocity, energy conjugacy, and dynamic signed budget must be verified on **independent geometries** in that family. A successful historical needle-plane test alone cannot establish the family-wide interface.

Three electrode-charge objects must remain separate:

| Object | Definition and current evidence | Permitted use |
|---|---|---|
| `q_A` | Historical native `q_hv`: center-field normal-area sum on the stair-step HV/gas interface. It is **not conjugate** to the current Poisson matrix's Laplace `psi`; the reported 22.26% reciprocity discrepancy is not cured by relabeling the outer term. | Preserve as a historical Stage C diagnostic; do not treat as the sole I1 terminal charge. |
| `q_B` | Candidate face charge from the **same Poisson links**, using `epsilon A/h (V_HV-phi_gas)` for the historical operator. Static reciprocity, energy conjugacy, and a read-only same-operator dynamic source-response budget were demonstrated for that operator. | Candidate *discrete-network* terminal observable only. Generalize by operator construction and independent-geometry tests, not by reusing the historical numeric coefficients. Do not silently replace `q_A` or claim physical terminal calibration. |
| `Q_physical` | Charge at the intended physical electrode/control surface and external return reference. | **NOT_RESOLVED** until electrode face location, conductor/exterior-boundary interaction, and dynamic convergence are audited. |

The companion adjoint weight for `q_A` is a response functional, not a fourth physical charge or the Laplace `psi`. In the frozen implementation, metal potential is assigned at a metal cell center and the Poisson gas-metal link uses distance `h`, whereas electron absorption uses `h/2`. This reference-plane mismatch and the outer-boundary artifact limit the physical interpretation of `q_B` despite its exact discrete-network reciprocity. See [V5-P2](../validation/i1_part_v_preflight2/V5_P2_FINAL_REPORT.md), sections 3-5, and [V5 entry decision](../validation/i1_part_v_preflight2/V5_ENTRY_DECISION.md).

## 5. Cold projected state and per-case numerical scope

With fixed geometry, define `E_w=-grad psi` and

\[
G_R=\int_{\Omega_g} e\mu_e n_e |\mathbf E_w|^2\,dV.
\]

Only gas cells and the physical cell volumes belong in the integral: use three-dimensional rotation volumes for an axisymmetric mesh and actual 3D cell volumes for a 3D mesh. `G_R` has units siemens; `E_w` has units inverse metres. The verified analytical and finite-volume structure, including the electron-source, weighted transport, boundary, and mobility-change terms of `dot G_R`, is documented in [Part V-A](../validation/i1_part_v/PART_V_A_VALIDATION_REPORT.md). That evidence is scoped to its tested discretization; it does not verify an unimplemented 3D handoff. Its exact discrete product budget is distinct from the smooth continuous formula. At the native Morrow-Lowke zero-field branch, mobility is nonsmooth/singular in the positive-field limit, so a classical smooth-coefficient `dot G_R` is not automatically valid. See [V5 preflight](../validation/i1_part_v_preflight/V5_PREFLIGHT_REPORT.md).

`G_R` is a cold-stage *projected electron state*. It is **not** the complete terminal current divided by `V_p`. Drift-diffusion, space-charge fields, boundary absorption, geometric/outer-boundary displacement, and the implemented-species scope must be reconciled at one control surface. At `V_p=0`, `G_R` remains definitionally meaningful but a quotient such as `P_E/V_p^2` does not. Exact native discrete-budget closure, static local-mesh improvement, and `q_B` reciprocity do **not** constitute Stage C physical trajectory convergence, a validated terminal map, K2 cross-condition prediction, or experimental validation.

Numerical readiness is **not** a claim that every possible geometry converges in advance. Define representative verification cases for each supported geometry/boundary class and retain their scope. Every actual run must carry its own mesh, time-step, geometry-mask and boundary-class identifiers, monitored observables, available error/sensitivity estimates, and `RESOLVED`, `CONDITIONAL`, or `NOT_RESOLVED` flags for the conclusions it uses. If a necessary estimate is absent, mark that conclusion `NOT_RESOLVED`; do not inherit a representative case's precision without an applicability check. A dynamic `G_R` or A/T/B/M claim requires corresponding dynamic evidence, not merely a static tip-ROI check.

## 6. Dependency and change control

The interface dependency order is: parameterized Stage B geometry and boundary class -> geometry-compatible Stage C (2D) or Stage E (3D) reconstruction, same-operator weighting and case-scoped numerical diagnostics -> Part V internal-state identification -> Part VI reviewed dynamic terminal/circuit coupling -> Stage G/F/H source/receiver interpretation as applicable. The existing numerical evidence concerns Stage C and cannot be transferred to Stage E without its own verification. An actual-electrode scientific claim additionally requires the actual geometry provenance; that requirement does not define the generic interface. These are **gates**, not permission to start the next stage. Stage F native RF, Stage G port spectrum, and Stage H received spectrum are not interchangeable; 350 MHz native RF remains `NOT_RESOLVED` ([scientific status](scientific_status.md)).

Do not alter frozen C/G/F equations or historical reference data to reconcile these definitions. Do not replace `q_A` with `q_B` without an explicitly reviewed interface change. Each later I1 verification should consult this contract and the [issue ledger](i1_open_issues_and_gates.md); update status only when the task authorizes maintenance and the issue's resolution gate is met.
