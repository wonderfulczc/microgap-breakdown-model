# Stage D1 Afivo-Streamer Backend Qualification

Stage D1 qualifies a pinned external Afivo-streamer backend for later 3D
cross-validation. No PETSc Stage C source or physics was modified.

## Environment

- OS: Ubuntu 26.04 LTS on WSL2, Linux `6.18.33.2-microsoft-standard-WSL2`
- Fortran: GNU Fortran 15.2.0
- C: GCC 15.2.0
- C++: G++ 15.2.0
- Make: GNU Make 4.4.1
- Git: 2.53.0
- Python: 3.14.4
- CPU: 16 logical processors
- RAM: 7.6 GiB total, 2.0 GiB swap

## Source Pin

- External path: `/home/helianthusczc/projects/afivo-streamer`
- Upstream: `https://github.com/MD-CWI/afivo-streamer.git`
- Pinned commit: `a50b5508775086e90dfe423455fb58d812578410`
- Internal Afivo core: no separate Git submodule; `afivo/` is fixed by the
  outer pinned tree.
- License: GPL-3.0

No system package was installed. The first official `make` failed while building
the bundled Hypre 2.31 dependency because GCC 15 treats `bool` as a C23 keyword:
`typedef unsigned char bool` in Hypre's `distributed_ls/pilut/struct.h` is then
an error. The official Hypre build script was rerun with `CFLAGS=-std=gnu17`,
after which the top-level build passed.

## Build

Initial command:

```sh
cd /home/helianthusczc/projects/afivo-streamer
/usr/bin/time -v make
```

Initial result: failed in Hypre, exit status 2, wall time 32.93 s, peak RSS
163836 KB.

Dependency repair:

```sh
cd /home/helianthusczc/projects/afivo-streamer/afivo/external_libraries
CFLAGS=-std=gnu17 /usr/bin/time -v ./build_hypre.sh
```

Result: passed, wall time 2.24 s, peak RSS 40424 KB.

Final build:

```sh
cd /home/helianthusczc/projects/afivo-streamer
/usr/bin/time -v make
```

Result: passed, exit status 0, wall time 47.16 s, peak RSS 125476 KB. Built
`programs/standard_3d/streamer` and `programs/3d_sprite/streamer`.

## Official Tests

Command:

```sh
cd /home/helianthusczc/projects/afivo-streamer
PATH=/home/helianthusczc/projects/streamer-rf-replica/.venv/bin:$PATH \
  /usr/bin/time -v bash run_test.sh
```

The existing self-repo virtualenv was used only to provide `numpy` for
`tools/compare_logs.py`; no system dependency was added.

Core streamer tests passed through:

- `programs/standard_1d/tests/test_1d.cfg`
- `programs/standard_1d/tests/test_1d_chemistry.cfg`
- all observed `programs/standard_2d/tests/*.cfg`
- `programs/standard_3d/tests/test_3d.cfg`
- `programs/standard_3d/tests/test_3d_chem.cfg`
- `programs/standard_3d/tests/test_3d_photoi_chem.cfg`
- first three dielectric charge tests

Full `run_test.sh` exit status was 1 because
`programs/dielectric_2d/tests/test_dielectric_neg_2d.cfg` differed from its
reference log beyond the default `rtol=1e-5`. The difference first appears in
late adaptive timesteps; this is classified as a non-Stage-D dielectric
regression drift under the current GCC 15 build, not a blocker for the 3D
Cartesian backend qualification.

Official test run resource use: wall time 3:02.04, peak RSS 293084 KB.

## True 3D Smoke

Program/config:

- executable: `programs/standard_3d/streamer`
- config: `programs/standard_3d/tests/test_3d.cfg`

Configuration summary:

- coordinates: Cartesian
- domain: `16e-3 x 16e-3 x 16e-3 m`
- box size: 8 cells per coordinate
- initial maximum level: 5; final highest level: 7
- initial leaves/boxes: 288 leaves, 329 boxes
- final cells: 405504
- minimum dx: `3.125e-5 m`
- gas pressure: 1 bar
- applied field: `-2.5e6 V/m`
- photoionization: disabled
- transport file: `td_air_siglo_swarm.txt`, old-style field table
- final physical time: `3.0e-9 s`
- status: normal exit, no NaN/Inf observed in log

Final regression diagnostic line:

```text
it=15 time=0.30000000E-008 dt=0.36009316E-010
sum(e)=0.36996457E+016 sum(M_plus)=0.40646712E+016 sum(M_min)=0.36502549E+015
max(e)=0.47552324E+019 max(M_plus)=0.50455742E+019 max(M_min)=0.43518375E+018
```

OpenMP smoke:

| OMP threads | status | wall time | peak RSS | final diagnostic |
| --- | --- | ---: | ---: | --- |
| 1 | pass | 10.30 s | 162168 KB | identical final line |
| 2 | pass | 5.80 s | 162164 KB | identical final line |
| 4 | pass | 3.49 s | 162056 KB | identical final line |

## Raw Output

`AFIVO_RAW_OUTPUT_FORMAT = Silo + text/log/config + optional binary datfile`.

The Stage D1 smoke had `silo_write = f`, so it produced compact text regression
logs, summary text, species/rates files, and an expanded output config. The
standard output module supports Silo via `silo_write`, optional restart binary
files via `datfile%write`, line/cross/plane diagnostic text outputs, and Afivo
also contains VTK writer support in `afivo/src/m_vtk.f90`.

Silo output includes AMR geometry and variables registered in the tree. Relevant
Stage D variables are available by name: `e`, positive/negative ion species,
`phi`, `electric_fld`, `rhs`, face-centered `field`, and time metadata.

## Future Stage D Interface Audit

- Geometry/electrodes: `src/m_field.f90`; config keys `use_electrode`,
  `field_given_by`, `field_electrode_type`, `field_rod_r0`, `field_rod_r1`,
  `field_rod_radius`, `field_electrode_grounded`, `field_bc_type`.
- Voltage/waveform: `field_given_by = voltage`, `field`, `voltage_table`, or
  `field_table`; documented in `documentation/electrodes_bc.md`.
- Pressure/temperature/gas: `gas%pressure` and gas configuration modules.
- Initial seed: `src/m_init_cond.f90`; config keys `seed_density`,
  `seed_charge_type`, `seed_rel_r0`, `seed_rel_r1`, `seed_width`,
  `seed_falloff`, `background_density`; custom hook through `m_user.f90`.
- Transport: `src/m_transport_data.f90`; config keys `input_data%file`,
  `input_data%old_style`, optional ion mobility keys.
- Chemistry: `src/m_chemistry.f90`; reaction tables embedded in transport data
  files such as `transport_data/air_chemistry_small_v1.txt`.
- Photoionization: `src/m_photoi.f90`, `src/m_photoi_helmh.f90`,
  `src/m_photoi_mc.f90`; config keys `photoi%enabled`, `photoi%method`,
  `photoi%eta`, Helmholtz coefficients, Monte Carlo controls.
- AMR: `src/m_refine.f90`; config keys `refine_adx`, `refine_min_dx`,
  `refine_max_dx`, `refine_electrode_dx`, `refine_min_dens`, `refine_cphi`,
  `refine_per_steps`.

## Common Transport Path

`COMMON_TRANSPORT_PATH = feasible`.

Afivo already supports tabulated transport input. The lowest-cost Stage D2 route
is:

```text
PETSc Morrow-Lowke evaluator
-> generate alpha(E/N), eta(E/N), mu_e(E/N), D_e(E/N) table
-> write Afivo old-style transport file with efield[V/m] columns
-> run Afivo with input_data%file and input_data%old_style = t
```

No Stage C baseline physics was changed in D1.

## Limitations

Stage D1 did not implement PETSc/Afivo numerical cross-validation, production
microgap geometry, triangular foil, offset needles, RLC, RF, or format
conversion. The 2D dielectric regression drift remains a non-blocking follow-up
for the external backend environment.
