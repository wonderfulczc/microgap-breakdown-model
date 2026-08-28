# Stage D2 PETSc-Afivo Common Benchmark

Stage D2 builds the common benchmark contract for future PETSc2D/Afivo3D
cross-validation. It does not start dynamic streamer comparison.

## Common Physics

The common gas state is dry air at `101325 Pa` and `300 K`. Photoionization is
off for D2-B0 and D2-B1.

Transport is exported directly from the frozen PETSc
`evaluate_morrow_lowke()` implementation by `generate_common_transport`.
The exporter writes:

- `morrow_lowke_transport_full.csv`: PETSc SI transport columns with separate
  `eta2` and `eta3`.
- `morrow_lowke_transport_afivo_old_style.txt`: Afivo old-style SI input with
  one attachment column equal to `eta2 + eta3` at the common gas density.
- `morrow_lowke_common_chemistry.txt`: Afivo reaction-table candidate with
  tabulated first-order ionization/attachment frequencies and analytic
  recombination terms.

The exact PETSc three-species source mapping is:

- ionization: `e -> e + e + P+`, rate `nu_ion(E)`.
- two-body attachment: `e -> N-`, rate `nu_att2(E)`.
- three-body attachment: `e -> N-`, rate `nu_att3(E)` at the common gas density.
- electron-positive recombination: `e + P+ -> neutral`, rate
  `1.138e-11 * Te^-0.7`.
- positive-negative recombination: `P+ + N- -> neutral`, rate
  `2e-13 * sqrt(300/Tg)`.

For Afivo, an exact runtime common chemistry path is feasible through the
generated reaction-table candidate plus a mean-energy table chosen so Afivo's
`Te` reproduces PETSc `e*D/(mu*kB)`. The old-style transport-only path is not
exact because it cannot retain the `eta2`/`eta3` split or recombination.

`COMMON_TRANSPORT_PATH = feasible`

## Seed Mapping

The common seed is axis-centered:

`ne(r,z) = n0 exp(-(r^2 + (z-z0)^2)/(2 sigma^2))`, `np=ne`, `nn=0`.

In Afivo3D this becomes:

`ne(x,y,z) = n0 exp(-((x^2+y^2) + (z-z0)^2)/(2 sigma^2))`.

Afivo's built-in line-segment seed is not guaranteed to be this exact Gaussian
semantics for all falloff choices. D3 should use a minimal `m_user.f90`
initialization hook if exact seed equality is required in a dynamic run.

## B1 Geometry

The canonical D2-B1 geometry is a centered capped rod over a grounded plane:

- ground surface: `z = 0`.
- HV tip apex: `z = 70 um`.
- HV tip center: `z = 75 um`.
- tip/shank radius: `5 um`.
- domain: `x,y in [-80,80] um`, `z in [-2.5,90] um`.
- voltage: `500 V`.

PETSc represents the ground as a thin grounded slab and the needle as the
Stage C axisymmetric structured mask, solved by a D2-only PETSc sampler with
Afivo-matched homogeneous outer boundaries. Afivo represents the HV electrode
with the standard true-3D Cartesian level-set rod and the ground as the lower
Dirichlet plane. This is a deliberate low-cost common benchmark, not a final
geometry-convergence study. The frozen Stage C Poisson/runtime path is not
changed.

## Photoionization

`PHOTOIONIZATION_COMMONIZATION = PARTIAL`

Both solvers have Helmholtz-style photoionization paths, but D2 leaves
photoionization off. Exact coefficient/source/boundary equivalence must be
audited in D3 before enabling it for cross-validation.

## D3 Blockers

The main D3 risk is not backend availability. It is the exact runtime mapping of
three-species chemistry and Gaussian seed initialization in Afivo without
modifying the pinned solver core. The generated reaction table and documented
`m_user.f90` seed hook route are the intended minimal path.

## D2 Run Summary

Generated on 2026-08-28 from PETSc branch `stage-c-real-electrode` and Afivo
commit `a50b5508775086e90dfe423455fb58d812578410`.

B0:

- transport representative-point maximum relative error: `8.65066e-05`.
- pointwise chemistry maximum relative error: `0`.

B1 coarse:

- PETSc resolution: `32 x 74`.
- Afivo resolution: `box_size=8`, `refine_electrode_dx=8e-6 m`.
- `phi` normalized RMSE: `0.0144925`.
- axis `|E|` normalized RMSE: `0.112994`.
- mid-gap `|E|` relative error: `0.00529144`.
- gas-line `Emax`: PETSc `1.70728e7 V/m`, Afivo `8.80637e6 V/m`.

B1 medium:

- PETSc resolution: `64 x 148`.
- Afivo resolution: `box_size=8`, `refine_electrode_dx=4e-6 m`.
- `phi` normalized RMSE: `0.0157364`.
- axis `|E|` normalized RMSE: `0.0300702`.
- selected `|E|` relative errors at `z/gap = 0.1,0.25,0.5,0.75,0.9`:
  `0.0234853`, `0.0231818`, `0.0216527`, `0.00597997`, `0.0519711`.
- gas-line `Emax`: PETSc `1.94854e7 V/m`, Afivo `1.58374e7 V/m`,
  relative error `0.187216`.

The coarse to medium axis-field RMSE improvement is `0.112994 -> 0.0300702`.
The remaining `Emax` difference is treated as the expected near-tip
cell-mask/level-set discretization sensitivity for this low-cost D2 benchmark,
not a D3 blocker.

Runtime/RSS:

- B0 exporter: `0.01 s`, peak RSS `4972 kB`.
- PETSc B1 medium: `0.22 s`, peak RSS `50404 kB`.
- Afivo B1 medium, `OMP_NUM_THREADS=2`: `0.09 s`, peak RSS `31468 kB`.
- B1 comparison: `0.05 s`, peak RSS `16664 kB`.

Validation:

- `cmake --build build --parallel`: PASS.
- `ctest --test-dir build --output-on-failure`: PASS, `10/10`.
- `.venv/bin/python -m pytest`: PASS, `51 passed, 3 skipped, 5 deselected`.

Notes:

- Afivo prints conservative interpolation-error estimates while creating its
  lookup tables, especially for low/zero reaction columns. D2 acceptance uses
  the explicit PETSc direct-evaluator round-trip comparison above.
- OpenMPI/UCX interface warnings appear in the sandbox but did not affect serial
  PETSc or Afivo benchmark completion.
