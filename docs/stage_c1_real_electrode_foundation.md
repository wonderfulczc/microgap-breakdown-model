# Stage C1 Real Axisymmetric Electrode Foundation

Date: 2026-08-27

Scope: real axisymmetric electrode infrastructure with time-varying voltage and
tip-relative local seed support. Stage C1 does not implement electrode current,
conductance, resistance, RLC coupling, or RF diagnostics.

## Added capability

Stage C1 adds an independent real-electrode path beside the frozen Stage 1-5
background-field workflow:

- structured-grid axisymmetric electrode classification
- internal electrode Dirichlet rows in the Poisson linear system
- constant and sampled voltage waveforms
- Gaussian seed placement relative to the electrode tip
- minimal electrostatic diagnostics for Stage C1 smoke validation

The legacy `solve_potential(..., background_field, ...)` API remains available.
Existing Stage 1-5 cases continue to use the uniform background-field path unless
an electrode geometry and voltage waveform are explicitly provided.

## Geometry representation

The first concrete geometry is `AxisymmetricNeedlePlaneGeometry`. It classifies
cell centers as:

- gas
- high-voltage electrode
- grounded electrode

Parameters use SI units:

- `geometry_id`
- `ground_z_m`
- `tip_z_m`
- `tip_radius_m`
- `shank_radius_m`
- `ground_thickness_m`

The current implementation is a structured-grid cell-center mask. It is simple
and robust, but it is not a cut-cell or immersed-boundary representation. Curved
tips are stair-stepped at cell-center resolution, so local tip-field accuracy
near high curvature is grid dependent.

Before claiming quantitative prediction of experimental needle-tip `Emax`, this
representation must be checked by:

- mesh sensitivity
- Stage B electrostatic comparison against COMSOL or an equivalent trusted
  electrostatic reference

Triangular copper-foil and other non-axisymmetric electrodes are intentionally
excluded from Stage C1 and remain future 3D backend scope.

## Poisson boundary treatment

`solve_potential_with_electrodes` assembles conductor cells directly into the
PETSc linear system:

- high-voltage conductor: `phi = V(t)`
- grounded conductor: `phi = 0`
- gas region: axisymmetric Poisson equation,
  `(1/r) d_r(r d_r phi) + d_zz phi = -rho/epsilon0`

The conductor potentials are not introduced by post-solve overwrite. Remaining
outer domain boundaries use the existing project boundary configuration.

## Voltage waveform

Stage C1 provides:

- `ConstantVoltage`
- `SampledVoltage`

Sampled voltage CSV format:

```text
time_s,voltage_V
0.0,0.0
1.0e-9,1000.0
```

Rules:

- `time_s` must be strictly increasing
- units are SI
- interpolation is piecewise linear
- no smoothing or filtering is applied
- no high-order interpolation is used
- default out-of-range behavior is an error, not silent extrapolation
- explicit endpoint-hold behavior is available only through the named policy

## Local initial electrons

Stage C1 reuses the existing `GaussianSeed` initialization. The new helper
`initialize_gaussian_at_tip_offset` computes:

```text
z_seed = z_tip + z_offset_from_tip
```

The existing quasi-neutral initialization remains unchanged:

```text
ne = np
nn = 0
```

Stage C1 does not implement field emission, secondary emission, stochastic
emission, or surface chemistry.

## Diagnostics

Stage C1 adds only electrostatic diagnostics:

- `time`
- `applied_voltage`
- `Emax`
- `phi_HV_residual`
- `phi_ground_residual`
- `geometry_id`
- `gap`
- `tip_radius`
- Poisson/KSP iterations

It does not add:

- `I_electrode`
- `I_displacement`
- `I_conduction`
- `Gb`
- `Rb`
- RLC coupling
- RF diagnostics

## Stage B validation interface

Stage C1 does not depend on COMSOL Stage B output. The reserved future
electrostatic cross-validation interface is:

- `geometry_id`
- `gap`
- `tip_radius`
- `V_ref`
- `Emax`
- `beta_E`
- axis profile: `s, phi, E`

Stage C may use Stage B data for `t = 0`, `rho = 0` electrostatic comparison
only. COMSOL static fields must not be fixed as the streamer field during time
evolution.
