# Stage F1 RF Source Contract

Stage F1 defines a post-processing data layer only. It does not implement
Jefimenko fields, FFT, ESD, CWT, RLC, antenna coupling, or RF interpretation.

## Current Provenance

`CURRENT_SOURCE_PROVENANCE = AFIVO_STAGE_E_CELL_CENTERED_DRIFT_CURRENT`

The Stage E Afivo source CSV writes

`J = e * mu_e(E/N) * n_e * E`

from the benchmark-specific `m_user.f90` export hook. This is a conventional
electron drift current density and follows Afivo's existing `Je_1..3` output
convention. Afivo's density update uses a finite-volume electron flux that
contains drift and diffusion. That exact face flux is not exported in Stage E.

For formal RF source calculations the reference definition is:

`J_RF = -e * Gamma_e`

where `Gamma_e` is the same finite-volume electron particle flux used in the
electron-density update. Until that transport flux is exported,
`J_transport_available = false` for Stage E source snapshots.

## Canonical SourceRecord

Each leaf cell stores:

- `cell_id` or stable spatial key
- AMR `level`
- `x_center`, `y_center`, `z_center`
- `dx`, `dy`, `dz` or equivalent bounds
- `cell_volume`
- `rho_Cpm3`
- `Jx_Apm2`, `Jy_Apm2`, `Jz_Apm2`
- optional `ne_m3`, `Ex_Vpm`, `Ey_Vpm`, `Ez_Vpm`

Native AMR is the reference representation:

`NATIVE_AMR_SOURCE = reference`

Uniform grids are allowed only as derived products. Production time derivatives
must use conservative overlap remapping when AMR partitions differ.

## Integral Diagnostics

Total charge:

`Q(t) = sum_i rho_i V_i`

Vector current moment:

`M(t) = integral J dV = sum_i J_i V_i`

`M` has units of `A m`. It is not an electrode terminal current.

Cross-section current:

`I_z(z,t) = integral J_z dx dy`

This is a channel/cross-section diagnostic, not a terminal-current diagnostic.

## Continuity Audit

The audited equation is:

`partial_t rho + div J = 0`

Global audit reports:

`dQ/dt + integral_boundary J dot n dA`

Local audit reports L1, L2, and Linf norms of
`partial_t rho + div J`. For Stage E drift-only exported current, the local
audit is a cell-centered reconstruction and is not expected to be exact.

## Frequency Metadata

F1 records:

- actual output times
- whether output spacing is uniform
- output interval
- record duration
- Nyquist frequency
- raw frequency resolution `1/T`

No FFT is performed in F1.

## Stage E 5 ps Data Status

`USABLE_FOR_RF_PIPELINE_TEST = YES`

`USABLE_FOR_VHF_UHF_SCIENTIFIC_SPECTRUM = NO`

The 5 ps Stage E source windows are useful to test source ingestion,
conservation, remapping, and derivative plumbing. They are too short to support
scientific VHF/UHF/GHz spectral conclusions.

