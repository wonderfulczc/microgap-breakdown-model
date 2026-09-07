# Stage C-R3 Streamer Head Tracking

C-R3 adds a passive streamer-head diagnostic for the frozen Stage C
axisymmetric solver. It does not modify the PDE update, Morrow-Lowke
transport, reaction chemistry, SP3, timestep control, electrode current, or
`Gb/Rb`.

## Existing Development Head Diagnostic

The previous Stage C scalar `head_position` diagnostic is an axis-line
electron-density tracker. It scans gas cells on `r=0` and reports the lowest
`z` where

```text
ne >= head_ne_threshold
```

If no axis cell crosses the threshold, it reports the axis gas cell with maximum
`ne`. The corresponding velocity is a finite difference of this scalar between
accepted steps. This remains useful as a development cross-check, but it is not
the formal C-R3 space-charge head definition.

## Primary Head Definition

C-R3 identifies the primary head from the existing net space charge

```text
rho = e * (np - ne - nn)
```

The algorithm is:

1. consider gas cells only;
2. find the dominant absolute space-charge extremum `max |rho|`;
3. set `head_polarity = sign(rho_peak)` in automatic mode, or use the
   documented case polarity supplied by the diagnostic caller;
4. select same-polarity cells satisfying
   `head_polarity * rho >= relative_rho_threshold * |rho_peak|`;
5. retain the deterministic 4-neighbor connected component containing the
   dominant extremum.

The default `relative_rho_threshold = 0.2` is a numerical segmentation
parameter, not a breakdown criterion or LFA criterion. For `stage_c2_dynamic`,
the default requested polarity follows the sign of the applied voltage unless
`--head-polarity` is supplied explicitly. This avoids switching between
opposite charge layers within a single positive- or negative-polarity case.

## Charge and Position

The signed head charge is

```text
q_h = integral_{Omega_h} rho dV
```

using the existing axisymmetric cell volume. The centroid uses the
non-cancelling same-polarity weight `|rho| dV`:

```text
z_h = integral z |rho| dV / integral |rho| dV
r_mean = integral r |rho| dV / integral |rho| dV
r_rms = sqrt(integral r^2 |rho| dV / integral |rho| dV)
```

The primary propagation coordinate for Stage C is `z_h`. The radial quantities
are diagnostics of the axisymmetric head extent.

## Kinematics

The tracker stores only accepted diagnostic states containing time, validity,
`z_h`, and `q_h`.

Velocity uses the accepted-state backward difference:

```text
v_n = (z_n - z_{n-1}) / (t_n - t_{n-1})
```

Acceleration uses the nonuniform three-point finite-difference form:

```text
a_n =
2 / (dt_n + dt_{n-1}) *
[
  (z_n - z_{n-1}) / dt_n
  -
  (z_{n-1} - z_{n-2}) / dt_{n-1}
]
```

The raw value is associated with the current accepted diagnostic sample. No
smoothing, spline fitting, or filtering is applied in C-R3.

## Invalid-State Policy

Undefined quantities are represented by explicit status strings and NaN-valued
numbers. Supported statuses include:

- `VALID`
- `NO_HEAD`
- `INSUFFICIENT_CHARGE`
- `INSUFFICIENT_HISTORY_FOR_VELOCITY`
- `INSUFFICIENT_HISTORY_FOR_ACCELERATION`
- `NONFINITE_INPUT`
- `INVALID_SEGMENTATION_PARAMETER`

Velocity and acceleration are not encoded as physical zero when history is
insufficient.

## Output

`stage_c2_dynamic --head-tracking` writes
`streamer_head_tracking.csv`. The diagnostic is off by default. A small
`streamer_head_segmentation_sensitivity.csv` is also written for relative
thresholds `0.1`, `0.2`, and `0.3`.

No full head-mask fields are written by default.

## Scope and Future Use

C-R3 prepares `q_h(t)`, `r_h(t)`, `v_h(t)`, and `a_h(t)` for future RF mechanism
decomposition, for example

```text
dM/dt = (dq_h/dt) v_h + q_h a_h + dM_redis/dt
```

C-R3 does not perform this decomposition and does not prove which term
dominates RF emission. That belongs to later F-R analysis.

Near bridging, collision, or multi-head branching, a single dominant connected
component can become an incomplete description of the discharge morphology.
