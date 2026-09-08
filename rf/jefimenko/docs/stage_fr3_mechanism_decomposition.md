# Stage F-R3 Current-Moment Mechanism Decomposition

F-R3 adds a passive, reduced-order time-domain decomposition of the
flux-derived current moment used by the frozen Stage F radiation workflow. It
does not modify Stage C, Stage E, the Stage F source contract, or the
production Jefimenko field kernel.

## Current Moment

The source moment is the existing scientific Stage F quantity:

```text
M(t) = integral J_RF dV
J_RF = -e Gamma_e
```

For the current PETSc axisymmetric sources the physically meaningful component
is `Mz`; radial current cancels azimuthally in the compact far-field moment
interpretation.

## Compact-Head Approximation

The compact head contribution is represented as:

```text
M_head(t) = sum_k q_k(t) v_k(t)
M_redis(t) = M(t) - M_head(t)
```

The present production analysis uses one validated dominant head:

```text
head_representation = PRIMARY_HEAD_ONLY
```

The head is extracted from the exact F4 source snapshots using the frozen C-R3
definition: net space charge `rho`, dominant same-polarity extremum, relative
threshold `0.2`, deterministic four-neighbor connected component, and
`|rho| dV` centroid.

`M_redis` is not assumed to be purely channel current. It contains all moment
content not captured by compact-head translation, including distributed
conduction, channel redistribution, head-shape evolution, unresolved secondary
heads, collision interaction, and approximation error.

## Mechanism Terms

The candidate time-domain terms are:

```text
HEAD_CHARGE_EVOLUTION = sum_k dq_k/dt * v_k
HEAD_ACCELERATION = sum_k q_k * a_k
CURRENT_MOMENT_REDISTRIBUTION = dM_redis/dt
```

The first term is deliberately not named impact ionization radiation. The C-R3
head charge is segmented net space charge, so changes in `q_h` can include
chemistry, transport, moving segmentation boundaries, and redistribution.

## Derivative Operator

The primary derivative is a centered local polynomial derivative for
nonuniform timestamps:

```text
degree = 2
radius_samples = 5
```

The local time coordinate is internally normalized before fitting and converted
back to SI derivative units. Samples lacking a complete centered stencil, or
crossing invalid head states, are invalid. A narrower `radius_samples = 2`
operator is retained as the derivative-sensitivity alternative. No spline,
global smoothing, or frequency-domain filtering is used.

## Closure

F-R3 computes:

```text
D_total = D[M_total]
D_reconstructed =
  HEAD_CHARGE_EVOLUTION
  + HEAD_ACCELERATION
  + CURRENT_MOMENT_REDISTRIBUTION

closure_residual = D_total - D_reconstructed
product_rule_residual =
  D[q v] - (D[q] v + q D[v])
```

The closure residual is therefore a direct diagnostic of the discrete product
rule and synchronization quality. It must not be hidden by redefining
`M_redis`.

## Stage4 Result

For `F4-P-stage4-left-isolated`, F-R3 used 50 common head/current-moment
samples from the existing trusted source window. The primary closure was:

```text
normalized_rms_closure = 0.0434428
valid_samples = 30
```

The mechanism-sum current-moment reconstruction correlated with the total
current-moment derivative at `0.998736` with normalized L2 error `0.0434428`.
The pre-existing full-Jefimenko/current-moment cross-check remains
correlation `0.999631` and normalized waveform error `0.0262439`.

Because the mechanism amplitudes change appreciably under the derivative and
segmentation sensitivity checks, all three Stage4 mechanism terms are marked:

```text
INTERPRET_WITH_CAUTION
```

## Stage5 Result

For `F4-C-stage5-highfield-collision`, F-R3 used 95 common samples and 75
valid primary derivative samples. The primary closure was:

```text
normalized_rms_closure = 1.26792
```

The mechanism-sum current-moment reconstruction had correlation `0.661890`
with normalized L2 error `1.26792`. This is not a small closure residual.
Because only a primary head is represented during a collision and F-R2 already
recommended a selected full-Maxwell reference for the collision source, Stage5
mechanism terms are marked:

```text
NOT_RESOLVED
FULL_MAXWELL_REFERENCE_PENDING
```

## Interface to F-R4

F-R3 exports lightweight mechanism time series and a JSON summary under:

```text
rf/jefimenko/validation/fr3_mechanisms/
```

F-R4 may use these time-domain source terms only with the status flags carried
forward. F-R3 does not compute mechanism-stage-frequency percentages and does
not claim that any frequency band is produced by a specific term. Final
coherent mechanism-stage-frequency attribution is deferred to F-R4.

## Limitations

The current decomposition is a candidate reduced-order model. It should be
described as proposed/developed in this work only after separate novelty
verification. It is not a replacement for the full Jefimenko spatial-source
solution, and it is not valid evidence for Stage5 mechanism dominance while
the closure residual remains large.
