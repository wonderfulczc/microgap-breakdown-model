# F-R2 Quasi-Static / Full-Maxwell Applicability Audit

F-R2 audits whether the frozen Stage F discharge sources are
electromagnetically compact enough for the current modelling architecture:

```text
electrostatic / drift-diffusion plasma evolution
+ one-way retarded Jefimenko radiation post-processing
```

The audit concerns the source-evolution approximation only. It does not
remove retardation from the radiation calculation; the Stage F Jefimenko
solver remains a retarded electromagnetic post-processor.

## Source Spatial Scale

The primary source is the frozen RF current:

```text
J_RF = -e * Gamma_e
```

with provenance:

```text
CONTINUITY_CONSISTENT_FINITE_VOLUME_FLUX
```

For each valid source snapshot:

```text
w_i = |J_i| V_i
r_c = sum(w_i r_i) / sum(w_i)
```

For 3D Cartesian records, `r_i` is the cell center. For axisymmetric PETSc
records, the current-weighted centroid lies on the symmetry axis, and the
distance used for `L95` is:

```text
sqrt(r_i^2 + (z_i - z_c)^2)
```

`r95` is the weighted radius containing 95% of the total current weight, and:

```text
L95 = 2 r95
```

The primary stage-level source size is the temporal p95 of `L95(t)`.
`Lbounding` is also reported as a conservative support extent. It can be
inflated by weak current tails and is therefore not used as the sole primary
scale.

## Current-Moment Timescale

The current moment is the existing Stage F quantity:

```text
M(t) = integral J dV
```

F-R2 defines a robust RMS source-evolution timescale:

```text
tau_M = sqrt( integral |M(t)|^2 dt / integral |dM/dt|^2 dt )
```

Actual timestamps are used. No smoothing is applied.

Constant or zero current moment does not produce an artificial finite fast
timescale; it is reported with an explicit status.

## Frequency-Domain Cross-Check

The frequency-domain cross-check uses the same window and computes the
spectral energy ratio of `M(t)` and the numerically differentiated `dM/dt`
through the existing FFT/Parseval utilities. It is a consistency check for
the `tau_M` definition, not a replacement for the Stage F trusted-band logic.

## epsilon_EM

The primary compactness parameter is:

```text
epsilon_EM = L_source / (c tau_M)
```

The equivalent quantities are:

```text
omega_equiv = 1 / tau_M
f_equiv = 1 / (2*pi*tau_M)
lambda_equiv = c / f_equiv
L_over_lambda = L_source / lambda_equiv
epsilon_EM = 2*pi * L_over_lambda
```

F-R2 also reports a conservative variant:

```text
epsilon_EM_conservative =
  max(L95_max, Lbounding_p95) / (c * tau_conservative)
```

where `tau_conservative` is based on p95 `|M|` over p95 `|dM/dt|` inside the
same trusted time window. Isolated derivative spikes are not used.

## Engineering Decision Bands

The following are engineering audit bands, not universal plasma-discharge
constants:

```text
epsilon_EM <= 0.1
  LOW_EM_FEEDBACK_RISK

0.1 < epsilon_EM <= 0.3
  REVIEW_SELECTED_FULL_MAXWELL_REFERENCE

epsilon_EM > 0.3
  SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED
```

The overall decision applies only to the audited frozen Stage4/Stage5 trusted
sources. It is not a statement about all future geometries, all Stage C cases,
or all Stage E cases.

## Frozen Stage4 Result

The audited Stage4 dataset is:

```text
F4-P-stage4-left-isolated
```

The primary result is low compactness risk for the current trusted window.
VHF/UHF/1-3 GHz bands remain governed by Stage F trust reports and are not
reinterpreted here.

## Frozen Stage5 Result

The audited Stage5 dataset is:

```text
F4-C-stage5-highfield-collision
```

The collision source has a larger spatial scale and faster current-moment
evolution. Its primary and conservative `epsilon_EM` fall in the engineering
review band, so a selected full-Maxwell reference case is recommended for this
trusted collision source before a final production manuscript claim.

## Scope Limitations

F-R2 does not implement:

- a full-Maxwell plasma solver,
- two-way EM/plasma feedback,
- a new RF source representation,
- a new streamer simulation,
- new frequency-band attribution.

It only audits whether a selected full-Maxwell reference should be prepared
later for the currently trusted Stage4/Stage5 RF source cases.
