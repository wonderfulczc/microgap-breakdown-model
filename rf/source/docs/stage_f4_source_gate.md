# Stage F4 RF Source Gate

Stage F4 scientific RF attribution requires a continuity-consistent current
source:

`J_RF = -e Gamma_e`

where `Gamma_e` is the finite-volume electron transport flux used by the
density update. A drift-only current, `J = e mu_e n_e E`, is still acceptable
for pipeline smoke tests but is not acceptable for VHF/UHF/GHz scientific
attribution.

## Afivo Source Export

The Stage E benchmark hook exports the historical drift columns
`Jx_Apm2,Jy_Apm2,Jz_Apm2` and, for future F4 reruns, the additional columns
`Jrf_x_Apm2,Jrf_y_Apm2,Jrf_z_Apm2`.

The `Jrf_*` columns are computed from Afivo's public `flux_elec` face field:
opposing faces are averaged to a cell-centered representation and multiplied by
`-e`. This includes drift and diffusion because it uses the same electron
transport flux object as the finite-volume update. Afivo core source files are
not modified.

Existing Stage E raw source files in this workspace do not contain `Jrf_*`
columns. They remain pipeline-only data.

## PETSc Source Export

`StreamerSolver::electron_transport_current_source()` reconstructs
cell-centered `Jr_RF,Jz_RF` from the same gas-gas, absorbing-electrode, and
outer-boundary electron face flux formulas used by `StreamerSolver::step()`.

`stage4_run` now writes these fields as `Jr_RF_A_m2,Jz_RF_A_m2` in sparse field
snapshots, and its `I_CM_electron` column is based on the flux-derived current
moment.

## Current Gate Outcome

The source export path is implemented, but F4 scientific stage attribution is
blocked until the frozen F4 source windows are rerun with these new exports.
The workspace currently has no reusable `results/stage4` or `results/stage5`
directories, and existing Stage E 5 ps source snapshots are drift-only.
