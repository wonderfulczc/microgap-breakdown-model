# C-R1c Electron Relaxation Interface for LFA Audit

Date: 2026-09-06

Scope: audit-only extension of the C-R1b field-scale diagnostic.

This node does not modify Stage C governing equations, Morrow-Lowke transport,
chemistry, SP3, timestep control, current definitions, or `Gb/Rb`.

## Repository Relaxation-Data Audit

The repository does not currently contain traceable electron-energy relaxation
time or relaxation length data. Existing candidates are not production relaxation
sources:

- `evaluate_morrow_lowke()` provides mobility, diffusion, ionization, attachment,
  Townsend coefficients, and frequencies, but not `tau_epsilon(E/N)` or
  `lambda_epsilon(E/N)`.
- `electron_temperature(mu, D)` is an Einstein-relation helper used by the
  recombination model. It is not an energy relaxation model.
- the D2 `mean_energy_eV` export is derived from `D/mu`; it is not a relaxation
  frequency, relaxation time, or relaxation length dataset.
- existing `TransportTable` is a bounded scalar interpolation helper and carries
  no electron-energy relaxation metadata.

Therefore the current production status remains:

```text
lfa_relaxation_data_status = NOT_AVAILABLE
LFA_APPLICABILITY = UNRESOLVED_RELAXATION_DATA
```

## ElectronRelaxationTable Contract

`ElectronRelaxationTable` is a generic future-data interface. It requires
traceable metadata:

- `source_id`
- `source_reference`
- `gas_composition`
- `pressure_Pa`
- `temperature_K`
- `neutral_density_m3`

and a strictly increasing `reduced_field_Td` grid. Either of the following
quantities may be supplied independently:

- `tau_epsilon_s`
- `lambda_epsilon_m`

The table performs explicit bounded linear interpolation. It does not silently
extrapolate. Outside the supplied range, lookup returns:

```text
OUTSIDE_RELAXATION_TABLE_RANGE
```

## Validation Policy

Malformed tables are rejected:

- empty or incomplete metadata;
- fewer than two reduced-field points;
- nonfinite values;
- duplicate or nonmonotonic `E/N`;
- nonpositive relaxation quantities;
- size mismatch between `E/N` and supplied relaxation arrays;
- absence of both `tau_epsilon` and `lambda_epsilon`.

No production relaxation dataset is added in C-R1c. Synthetic tables are used
only in unit tests.

## Prepared Nondimensional Quantities

If a valid relaxation table is supplied, the diagnostic can compute:

```text
chi_L = lambda_epsilon / L_E
chi_t = tau_epsilon / tau_E
```

Only cells with valid `L_E` or `tau_E` and in-range relaxation lookup contribute
to the corresponding chi statistic.

C-R1c intentionally imposes no universal threshold on `chi_L` or `chi_t`. Even
with a table present, applicability remains unresolved unless a separately
approved physical criterion is supplied.

## Region-Aware Statistics

Whole-domain `LE_min_valid_m` is retained only as a low-level diagnostic. It can
be dominated by low-field or nearly singular numerical locations and should not
be used as a headline physical gradient scale.

The audit now supports these diagnostic regions:

- `ALL_GAS`: all gas cells;
- `ACTIVE_ELECTRON`: gas cells with `ne > 1e-12 * max(ne_max, n_ref)`;
- `HIGH_FIELD`: gas cells with `|E| >= 0.5 * max_gas(|E|)`.

These masks are diagnostic region definitions, not LFA validity thresholds and
not universal physical criteria.

For each region, compact output can report:

- `LE_p05_m`
- `LE_median_m`
- `tauE_p05_s`
- `tauE_median_s`
- `chiL_p95`
- `chiL_median`
- `chiT_p95`
- `chiT_median`
- `relaxation_coverage_fraction`
- `fraction_outside_relaxation_table`

Streamer-head-specific LFA statistics are deferred until C-R3 provides a frozen
head segmentation definition.
