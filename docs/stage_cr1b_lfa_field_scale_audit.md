# C-R1b LFA Field-Scale Audit

Date: 2026-09-06

Scope: diagnostic-only field-scale audit for the frozen Stage C axisymmetric
streamer solver.

This node does not modify the governing equations, Morrow-Lowke transport,
chemistry, SP3 photoionization, timestep control, electrode current definitions,
or `Gb/Rb` definitions.

## Definitions

For gas cells only, the audit computes:

```text
EoverN_Td = |E| / N * 1e21
L_E       = |E| / |grad |E||
tau_E     = |E| / |d|E|/dt|
```

`L_E` is a spatial variation scale of the local electric-field magnitude.
`tau_E` is a temporal variation scale between consecutive accepted solver
states.

## Spatial Derivative Policy

The Stage C field magnitude `|E|` is cell-centered. `grad |E|` is reconstructed
on the existing axisymmetric structured mesh:

- centered differences are used when both neighboring gas cells are available;
- one-sided differences are used next to domain or electrode-mask boundaries;
- stencils never cross conductor cells;
- the radial coordinate uses the existing cell-centered grid, so the axis side
  is handled by the same one-sided policy;
- gas cells with no usable field stencil are marked invalid.

The diagnostic does not divide by an arbitrary physical epsilon. If
`|grad |E||` is zero or below a numerical roundoff-scale threshold, the cell is
excluded from valid `L_E` statistics and counted as `near_zero_gradE`.

## Temporal Derivative Policy

`tau_E` uses only consecutive accepted states. The implementation keeps a
diagnostic-only copy of the previous accepted `|E|` field in `LfaAuditHistory`.
It is not part of the PDE state and is not read by the timestep controller.

The first sample has no previous accepted field and is reported with:

```text
temporal_status = INVALID_INITIAL_SAMPLE
```

Subsequent samples use the actual accepted step size:

```text
d|E|/dt = (|E|_new - |E|_previous) / dt_accepted
```

Cells with zero or roundoff-scale `d|E|/dt` are excluded from valid `tau_E`
statistics and counted as `near_zero_dEdt`.

## Compact Output

The Stage C2 executable writes `lfa_audit.csv` only when explicitly run with
`--lfa-audit`. Default Stage C output is unchanged.

The compact fields are:

- `EoverN_max_Td`
- `LE_min_valid_m`
- `LE_p05_m`
- `LE_median_m`
- `LE_valid_fraction`
- `tauE_min_valid_s`
- `tauE_p05_s`
- `tauE_median_s`
- `tauE_valid_fraction`
- invalid-cell counters and status strings

No full LFA field is dumped at every step by default.

## Region Masks

The audit accepts an optional boolean cell mask so later C-R nodes can compute
statistics for selected gas-cell regions, such as high-field cells or
streamer-head cells. C-R1b does not define those physics regions and does not
introduce a high-field threshold.

## LFA Applicability Limit

The repository currently does not contain traceable electron-energy relaxation
time or relaxation length data. Therefore C-R1b cannot classify the local-field
approximation as physically valid or invalid.

The formal status is:

```text
lfa_relaxation_data_status = NOT_AVAILABLE
LFA_APPLICABILITY = UNRESOLVED_RELAXATION_DATA
```

A future complete audit may compare the field scales with externally validated
electron-energy relaxation scales:

```text
chi_L = lambda_epsilon / L_E
chi_t = tau_epsilon / tau_E
```

Those quantities require traceable electron-energy relaxation data, such as a
validated Boltzmann/BOLSIG-derived table for the same gas, pressure,
temperature, and reduced-field range. They are not invented in this node.
