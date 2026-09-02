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

## Recovery Outcome

The source export path is implemented and the F4 source-gate recovery generated
short PETSc 2D source windows under `results/stage_f4/source_gate_recovery/`.
Small tracked summaries are stored under `rf/production/f4_source_recovery/`.

Recovered windows:

- `F4-A-avalanche-inception`: early Stage4 left-isolated window, source-valid
  but too short for scientific RF interpretation.
- `F4-P-streamer-propagation`: Stage4 left-isolated source-gate window with
  flux-derived `J_RF`.
- `F4-C-interaction-collision`: Stage5 high-field collision parameters in a
  short source-gate window with flux-derived `J_RF`.

The short recovery windows restore the source gate for GHz/SHF pipeline smoke.
They do not resolve VHF/UHF attribution. Existing Stage E 5 ps source snapshots
remain drift-only and pipeline-only.
