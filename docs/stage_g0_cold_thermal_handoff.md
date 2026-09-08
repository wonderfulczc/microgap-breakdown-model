# Stage G0 Cold-Streamer To Thermal-Spark Handoff Framework

G0 establishes a passive interface between the frozen cold Stage C streamer
model and a future thermal spark-channel model. It does not modify Stage C,
Stage D, Stage E, Stage F-R, C-R4 Joule diagnostics, or any RLC circuit model.

## Existing `bridge_flag`

`StreamerSolver::bridge_flag()` is an axial threshold proxy. In electrode mode
it checks the gas cells on the symmetry axis between the high-voltage needle
tip and the grounded electrode and returns true only if every such axial cell
satisfies:

```text
ne >= bridge_ne_threshold
```

This diagnostic does not explicitly build a two-dimensional connected
conductive path through the gas region. G0 therefore records it as
`bridge_legacy` with classification `AXIAL_THRESHOLD_PROXY`.

## Conductive Percolation Diagnostic

When field snapshots are available, G0 adds an independent structured-grid
connectivity diagnostic. The candidate conductive gas mask is:

```text
sigma >= f_sigma * sigma_max
```

where `sigma` is computed from the existing Stage C diagnostic conductivity
definition:

```text
sigma = e * mu_e(E) * ne
```

`f_sigma` is a numerical connectivity parameter, not a physical thermal
transition threshold. The default is `0.1`, matching the C-R4 numerical channel
diagnostic region. Connectivity uses deterministic von-Neumann neighbours,
gas cells only, and explicit high-voltage-adjacent and ground-adjacent seed
sets. It reports:

```text
percolation_valid
percolation_status
percolation_path_length_m
percolation_bottleneck_sigma_S_m
```

An isolated conductive island or two disconnected half-paths cannot be counted
as a bridge.

## Percolation Sensitivity

G0 can evaluate a tiny set of nearby numerical mask fractions, for example:

```text
0.05, 0.1, 0.2
```

This is a data-quality check for catastrophic dependence of the connectivity
diagnostic. It is not threshold calibration.

## C-R4 Joule-Energy Reuse

G0 reuses the frozen C-R4 compact output rather than redefining energy
diagnostics. The C-R4 quantities remain:

```text
PJ_channel_W
QJ_channel_J
channel_volume_m3
channel_length_m
channel_effective_radius_m
sigma_eff_S_m
E_channel_mean_V_m
ne_channel_mean_m3
ne_channel_max_m3
Gb_S
Rb_ohm
dGb_dt_S_s
tau_sigma_s
tau_evolution_s
Xi_sigma
```

The underlying conductive Joule power is still based on
`J_cond = -e Gamma_e` from the finite-volume electron transport flux. G0 does
not recompute C-R4 `PJ`, `QJ`, `Gb`, `Rb`, or `Xi_sigma` with a different
physical definition.

## `Pi_H`

The thermal progress interface is:

```text
Pi_H = QJ_channel / Q_required
```

No production `Q_required` exists in the frozen cold Stage C model. G0 therefore
sets:

```text
Q_required_J = NaN
Pi_H = NaN
thermal_energy_reference_status = NOT_AVAILABLE
```

It does not estimate a target temperature, spark temperature, LTE composition,
or radial thermal profile.

## `Xi_sigma`

G0 carries forward the C-R4 conductivity-evolution audit quantity:

```text
tau_sigma = epsilon0 / sigma_eff
tau_evolution = |Gb / dGb_dt|
Xi_sigma = tau_sigma / tau_evolution
```

`Xi_sigma` is an audit observable only. G0 does not define `Xi_sigma*`.

## Calibration Object

`HandoffCalibration` is an optional future data contract. It may carry:

```text
source_id
source_reference
Q_required_J
Pi_H_threshold
Xi_sigma_threshold
bridge_requirement
```

Values must be finite, positive where applicable, and provenance metadata must
be present. No production default thresholds are supplied. Synthetic calibration
is used only in unit tests.

## Handoff Status Logic

Without a valid calibration object, production status remains:

```text
HANDOFF_CALIBRATION_PENDING
```

G0 never outputs `THERMAL`, `SPARK`, or `HANDOFF_READY` from bridge/percolation
or `Xi_sigma` alone.

## G1 Initial-State Contract

The future G1 thermal solver can consume a compact cold-state contract:

```text
handoff_time_candidate
channel_length_m
channel_effective_radius_m
channel_volume_m3
sigma_eff_S_m
QJ_channel_J
PJ_channel_W
E_channel_mean_V_m
ne_channel_mean_m3
ne_channel_max_m3
Gb_S
Rb_ohm
thermal_profile_status
```

`thermal_profile_status` is currently `NOT_AVAILABLE`. G0 does not fabricate a
radial temperature profile.

## Production Case Result Policy

If a short existing Stage C case does not form a bridge or percolated conductive
path, the correct status is:

```text
PRODUCTION_HANDOFF_NOT_REACHED
HANDOFF_CALIBRATION_PENDING
```

This does not invalidate the G0 framework. It only means the selected cold
development case is not itself a calibrated thermal handoff event.

## Limitations

G0 does not calibrate `Q_required`, `Pi_H*`, `Xi_sigma*`, or a cold-to-thermal
transition threshold. Those require a later G1 overlap thermal simulation and
experimental validation.
