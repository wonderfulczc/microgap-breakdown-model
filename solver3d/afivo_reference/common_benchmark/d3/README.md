# Stage D3 Dynamic PETSc-Afivo Cross-Validation

Stage D3 uses the Stage D2 common benchmark to compare the frozen PETSc
axisymmetric 2D solver with the pinned Afivo-streamer true 3D Cartesian backend
under an axisymmetric physical setup.

Scope limits:

- Photoionization is off.
- No RLC, RF, current/Rb extraction, bridge tuning, or Stage E geometry.
- Afivo core source is not modified. The benchmark hook in `afivo_user/` is
  copied temporarily over `programs/standard_3d/m_user.f90` only while running
  D3, then the external checkout is restored.

## Physics Contract

- Geometry: centered capped 5 um rod/needle over grounded plane, 70 um gap.
- Gas: 101325 Pa, 300 K. Afivo config uses `gas%pressure = 1.01325` because
  Afivo expects pressure in bar.
- Voltage: constant 500 V.
- Transport: PETSc frozen Morrow-Lowke exported by
  `build/bin/generate_common_transport`.
- Chemistry: PETSc frozen three-species model mapped to Afivo reaction tables.
- Seed: axis-centered Gaussian, `ne=np=n0*exp(-(x^2+y^2+(z-z0)^2)/(2*sigma^2))`,
  `nn=0`, `n0=1e16 m^-3`, `sigma=3 um`, `z0=65 um`.

Afivo built-in Gaussian seed is used with `seed_width=sqrt(2)*sigma`, which
matches the common-case formula. The D3 hook is only used for homogeneous
runtime-chemistry initialization and compact symmetry diagnostics.

## Diagnostics

The shared Python postprocessor computes dynamic diagnostics at matched physical
times, not step numbers:

- `Emax(t)`
- axis-lineout `ne_max(t)`
- total electron count `N_e(t)`
- leading-edge `head_position(t)` and `head_velocity(t)`
- Afivo lateral center of electron density `x_cm(t), y_cm(t)`
- Afivo moment asymmetry `abs(int (x^2-y^2) ne dV) / int (x^2+y^2) ne dV`

Head definition is the minimum axis `z` satisfying `ne >= 1e14 m^-3`; fallback
is the axis position of maximum `ne`.

## Commands

Generate transport and runtime chemistry reference:

```sh
build/bin/generate_common_transport solver3d/afivo_reference/common_benchmark/transport/generated
build/bin/stage_d3_reaction_reference solver3d/afivo_reference/common_benchmark/d3/results/d3_reaction_reference.csv
```

Run PETSc:

```sh
build/bin/stage_d3_petsc_dynamic --voltage 500 \
  --out solver3d/afivo_reference/common_benchmark/d3/petsc/output/d3_petsc_dynamic_500V.csv
```

Run Afivo after temporarily applying `afivo_user/m_user.f90` to the external
`standard_3d` program and rebuilding:

```sh
OMP_NUM_THREADS=4 /home/helianthusczc/projects/afivo-streamer/programs/standard_3d/streamer \
  solver3d/afivo_reference/common_benchmark/d3/afivo/d3_seed_check.cfg
OMP_NUM_THREADS=4 /home/helianthusczc/projects/afivo-streamer/programs/standard_3d/streamer \
  solver3d/afivo_reference/common_benchmark/d3/afivo/d3_dynamic_500V.cfg
```

Compare:

```sh
.venv/bin/python solver3d/afivo_reference/common_benchmark/d3/compare/compare_stage_d3.py chemistry \
  --reference solver3d/afivo_reference/common_benchmark/d3/results/d3_reaction_reference.csv \
  --amounts ionization solver3d/afivo_reference/common_benchmark/d3/afivo/output/d3_chem_ionization_amounts.txt \
  --amounts attachment solver3d/afivo_reference/common_benchmark/d3/afivo/output/d3_chem_attachment_amounts.txt \
  --amounts recombination solver3d/afivo_reference/common_benchmark/d3/afivo/output/d3_chem_recombination_amounts.txt \
  --out solver3d/afivo_reference/common_benchmark/d3/results/d3_chemistry_runtime.json

.venv/bin/python solver3d/afivo_reference/common_benchmark/d3/compare/compare_stage_d3.py dynamic \
  --petsc solver3d/afivo_reference/common_benchmark/d3/petsc/output/d3_petsc_dynamic_500V.csv \
  --petsc-line-glob 'solver3d/afivo_reference/common_benchmark/d3/petsc/output/axis_profile_*.csv' \
  --afivo-log solver3d/afivo_reference/common_benchmark/d3/afivo/output/d3_dynamic_500V_log.txt \
  --afivo-line-glob 'solver3d/afivo_reference/common_benchmark/d3/afivo/output/d3_dynamic_500V_line_*.txt' \
  --out-csv solver3d/afivo_reference/common_benchmark/d3/results/stage_d3_comparison.csv \
  --out-json solver3d/afivo_reference/common_benchmark/d3/results/stage_d3_metrics.json
```

## Qualification Result

The 500 V benchmark shows avalanche-level electron population evolution without
bridge. Afivo remains axisymmetric to far below one finest cell. The robust
cross-solver dynamic diagnostics satisfy the D3 engineering acceptance targets.

Known debts:

- `TIP_FIELD_DISCRETIZATION_DEBT = OPEN`: PETSc structured axisymmetric mask and
  Afivo Cartesian AMR level set retain different tip discretizations.
- `PHOTOIONIZATION_COMMONIZATION = PARTIAL`: main D3 benchmark keeps
  photoionization off; SP3/Helmholtz coefficient commonization is deferred.
- `NUMERICAL_SCHEME_DIFFERENCE`: finite-volume/PETSc 2D and Cartesian AMR/Afivo
  discretizations are intentionally not made identical in Stage D.
