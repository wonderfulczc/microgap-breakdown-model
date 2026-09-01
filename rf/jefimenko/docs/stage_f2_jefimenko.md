# Stage F2 Jefimenko Field Solver

Stage F2 implements the field waveform solver only. It does not perform FFT,
ESD, CWT, band-energy analysis, antenna coupling, RLC coupling, or RF spectrum
interpretation.

For observer position `r`, source position `r'`,

`R_vec = r - r'`, `R = |R_vec|`, `Rhat = R_vec / R`,
and `t_r = t - R/c`.

The implemented standard Jefimenko form is:

`E(r,t) = 1/(4*pi*eps0) integral [ rho(t_r) Rhat/R^2 + drho_dt(t_r) Rhat/(c R) - dJ_dt(t_r)/(c^2 R) ] dV'`

`B(r,t) = mu0/(4*pi) integral [ J(t_r)/R^2 + dJ_dt(t_r)/(c R) ] x Rhat dV'`

The magnetic cross-product order is `(source vector) x Rhat`. For a positive
`+z` current element observed on the `+x` axis this gives `+y` magnetic field,
consistent with the right-hand rule.

## Decomposed Outputs

Each observer sample stores:

- `E_rho_near`
- `E_drho_induction`
- `E_dJ_radiation`
- `E_total`
- `B_J_near`
- `B_dJ_radiation`
- `B_total`

The CSV row schema uses component names such as `Ex_total`, `Ez_dJ`,
`By_J`, and `By_dJ`.

## Native Source Integration

The solver integrates F1 `SourceRecord` leaf cells directly:

`sum integrand_i * cell_volume_i`

No uniform-grid interpolation is performed. If AMR partitions differ in time,
the input `SourceSeries` must first be conservatively remapped by the F1
framework. Observer points inside source cells are rejected; no arbitrary
small-`R` regularization is used.

## Retarded Time

For every source cell the solver evaluates `rho`, `J`, `drho/dt`, and `dJ/dt`
at `t_r = t - R/c`. Field values use linear interpolation in time. Derivatives
use the F1 three-point central derivative, including the nonuniform formula.

Only samples whose retarded times lie inside the derivative support are marked
`retarded_time_valid = true`. No silent extrapolation is used.

## Scientific Gate

F1 currently records Stage E `J` as drift-only:

`CURRENT_SOURCE_PROVENANCE = AFIVO_STAGE_E_CELL_CENTERED_DRIFT_CURRENT`

The formal RF source remains:

`J_RF = -e * Gamma_e`

with `Gamma_e` from the density-update finite-volume electron transport flux.
The Stage E 5 ps data are therefore valid for pipeline smoke only, not for
scientific VHF/UHF/GHz spectrum claims.

