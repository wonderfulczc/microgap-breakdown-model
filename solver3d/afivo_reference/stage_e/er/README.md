# Stage E-R Direction-Alignment Diagnostic

Stage E-R is a post-processing reinforcement node for existing Stage E 3D
smoke data. It does not modify Afivo core source, Stage E physics, Stage D
validation, Stage C-R physics, or Stage F RF source definitions.

## Existing Stage E Head Definition

The completed Stage E Afivo hook reports `e2_head_z` in
`afivo_user/m_user.f90`. It scans gas cells and sets:

```text
head_z = min(z_cell where ne >= e2_head_threshold)
```

If no such cell exists it writes an invalid sentinel. This is a leading-edge
electron-density head locator. It was already used by Stage E2 and Stage E3
summaries.

E-R retains this definition as `legacy_head_definition`. The corrected
post-processor also records whether the legacy minimum-z cell belongs to the
dominant active-electron component.

## Coordinate and Polarity Verification

The Stage E development cases use `z` as the gap direction. The high-voltage
electrode is located at high `z`, the grounded electrode is located at low
`z`, and the documented propagation direction is toward lower `z`:

```text
reference_axis_vector = (0, 0, -1)
```

For the analyzed positive-HV cases, the exported local electric field has
negative `Ez` near the head-local sample, consistent with `E = -grad(phi)` and
with the macroscopic field direction from the positive electrode to the
grounded electrode. Therefore:

```text
polarity_sign = +1
e_E_prop = e_E_raw
```

This sign is fixed from the geometry and boundary-condition metadata. It is
not selected by minimizing `Delta_theta`.

## Corrected Head-Local Field Definition

For each existing source snapshot, E-R reconstructs a minimal 3D head location
that keeps the existing electron-density threshold semantics while avoiding
single-cell threshold-front switching where possible:

1. select gas cells with `lsf_m > 0`;
2. select cells with `ne_m3 >= e2_head_threshold`;
3. identify connected active-electron components in the sparse source export;
4. select the dominant active component;
5. compute an `ne*cell_volume` weighted centroid in its leading-edge slab;
6. sample the local peak `|E|` gas cell within `4*min_dx` of that centroid.

This local peak is intentionally head-local. E-R does not use the global
computational-domain `Emax`, because the global maximum may remain at an
electrode tip.

## Reference Axis

Both analyzed Stage E geometries use the same propagation convention:

```text
reference_axis_vector = (0, 0, -1)
```

The scalar angles `theta_E` and `theta_head` are measured relative to this
axis:

```text
theta = acos(clamp(unit_vector dot reference_axis, -1, 1))
```

## Polarity Convention

The analyzed cases use positive high-voltage electrodes and grounded lower
electrodes. The propagation direction is toward ground. Under this documented
polarity the propagation-oriented field unit vector is:

```text
e_E_prop = e_E_raw
```

where `e_E_raw = E_head / |E_head|`. The sign is fixed from case metadata, not
chosen to improve agreement with the trajectory.

## Trajectory Direction

For each valid head sample at time `t_n`, E-R uses the forward output
displacement:

```text
d_head(t_n) = (r_h(t_(n+1)) - r_h(t_n)) / |r_h(t_(n+1)) - r_h(t_n)|
```

This compares the local field at `t_n` with the immediately subsequent head
motion. Nonuniform output times are allowed because only the displacement
direction is used. Duplicate or stationary head positions are marked invalid.

No smoothing, spline fitting, or filtering is applied.

## Trajectory Resolution Gate

A trajectory direction is numerically meaningful only when a forward head
displacement is resolved by the existing 3D output grid. For each valid sample,
E-R searches for the earliest future valid sample whose displacement reaches
the resolution gate. E-R records:

```text
trajectory_displacement_m
trajectory_dx_ref_m
trajectory_displacement_over_dx
```

where `trajectory_dx_ref_m` is the maximum of the start and selected future
sample minimum grid spacings from the Afivo log metadata. A direction is valid
only when:

```text
trajectory_displacement_m >= 1.0 * trajectory_dx_ref_m
```

Otherwise it is marked:

```text
INSUFFICIENT_SPATIAL_DISPLACEMENT
```

This is a numerical trajectory-reliability gate, not a physical propagation
criterion.

## Full-Vector Angular Mismatch

The primary E-R metric is the 3D vector mismatch:

```text
Delta_theta = acos(clamp(e_E_prop dot d_head, -1, 1))
```

It is not computed as `abs(theta_E - theta_head)`, because two 3D vectors can
have the same polar angle and different azimuth.

## Invalid States

E-R records explicit invalid states:

```text
NO_VALID_HEAD
NO_VALID_LOCAL_FIELD
ZERO_FIELD_MAGNITUDE
INSUFFICIENT_NEXT_POSITION
ZERO_HEAD_DISPLACEMENT
POLARITY_UNRESOLVED
NONFINITE_INPUT
MISSING_SOURCE_FIELDS
INSUFFICIENT_SPATIAL_DISPLACEMENT
INVALID_SPATIAL_RESOLUTION
```

Undefined directions are not encoded as physical zero angles.

## Corrected Current Results

The output files are:

```text
stage_er_direction_alignment.csv
stage_er_direction_alignment_triangular_foil.csv
stage_er_direction_alignment_misaligned_needle.csv
stage_er_direction_summary.json
```

The triangular-foil E2 source shows early avalanche but no spatially resolved
head displacement over the available source snapshots. The corrected forward
head displacement is only about `0.0052 dx`, so the earlier apparent
`Delta_theta ~ 179 deg` is rejected as under-resolved rather than interpreted as
a physical anti-alignment.

The original misaligned-needle E3 source snapshots were full-domain source
exports and were too sparse for a reliable head trajectory. A final E-R
diagnostic rerun therefore used the same frozen E3 physical configuration but
disabled full source export and enabled lightweight head-local output:

```text
e2%write_source_export = f
e2%write_head_local_export = t
output%dt = 1.0e-13
```

The head-local CSV files contain only active gas cells satisfying the existing
`ne >= e2_head_threshold` semantics, with:

```text
time_s,x_m,y_m,z_m,cell_volume_m3,dx_m,ne_m3,
Ex_Vpm,Ey_Vpm,Ez_Vpm,Eabs_Vpm,lsf_m,level
```

This produced 51 head-local samples over `0-5 ps`, with about `80 MB` of new
diagnostic data instead of another multi-GB full-field source export.

The forward-lag trajectory rule searches for the earliest future sample whose
corrected head displacement exceeds `1 dx`. Only one misaligned-needle sample
satisfies this gate. The corrected total head displacement over the available
5 ps window is about `4.09 dx`, below the E-R diagnostic sufficiency target of
`5 dx`, and most future head positions remain under-resolved. The sole
resolved angle is:

```text
Delta_theta = 90.63 deg
trajectory_lag = 0.1 ps
trajectory_displacement_over_dx = 3.70
```

Because this single resolved sample cannot establish a temporal trend, it is
not interpreted as a physical field-head association.

Both analyzed cases are therefore classified as:

```text
directional_evidence = NOT_RESOLVED
```

## Sampling Sensitivity

E-R performs only a lightweight temporal sampling check: full available source
samples versus every second sample. This detects gross sampling artifacts but is
not a convergence study.

For the final head-local misaligned-needle result, the remaining single valid
angle is not sampling-sensitive, but the number of valid samples is
insufficient for a trajectory-level physical conclusion.

## Limitations

E-R establishes quantitative directional association between head-local peak
field and the next resolved leading-edge displacement. It does not by itself
prove strict causality between electric-field direction and streamer
propagation. With the current Stage E smoke cases, the directional association
remains unresolved because the available non-axisymmetric dynamics are
early-avalanche cases without sustained, spatially resolved streamer-head
propagation. This is now frozen as an E-R limitation: later F-R work may
proceed, but it must not use E-R directional association as established
physical evidence.
