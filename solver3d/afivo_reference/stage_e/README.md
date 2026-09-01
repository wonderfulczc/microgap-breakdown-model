# Stage E1 Afivo Triangular-Foil Electrostatic Benchmark

Stage E1 qualifies a genuinely non-axisymmetric 3D electrode geometry in the
pinned external Afivo-streamer backend. This is a development geometry, not a
claim of final experimental dimensions.

Afivo source is external and pinned at:

- Path: `/home/helianthusczc/projects/afivo-streamer`
- Upstream: `https://github.com/MD-CWI/afivo-streamer.git`
- Commit: `a50b5508775086e90dfe423455fb58d812578410`
- License: GPL-3.0

The benchmark-specific hook in `afivo_user/m_user.f90` is copied temporarily
over `programs/standard_3d/m_user.f90` for E1 runs and restored afterwards. It
defines only the triangular-foil electrode level set and E1 diagnostics; it does
not modify Afivo core `src/`.

## Coordinate Convention

- `x`: foil width direction
- `y`: foil thickness / transverse direction
- `z`: gap / propagation direction
- Ground plane: low-`z` Dirichlet boundary at `z = 0`
- HV triangular foil: high-`z` conductor, tip points toward ground
- Positive HV voltage: `+V` on the foil
- Streamer propagation direction reserved for later stages: `-z`

## Geometry

The canonical development foil is a finite rounded triangular prism:

- `gap_m = 70e-6`
- `foil_thickness_m = 10e-6`
- `foil_width_m = 40e-6`
- `foil_length_m = 50e-6`
- `triangle_tip_angle_deg = 43.6028189727`
- `tip_radius_m = 3e-6`
- `edge_radius_m = 1e-6`
- `voltage_V = 500`

The level-set model is a rounded 2D triangular tip extruded through a finite
foil thickness with finite transverse edge rounding. It intentionally avoids an
ideal mathematical point.

## Commands

Run local geometry checks:

```bash
.venv/bin/python solver3d/afivo_reference/stage_e/scripts/e1_geometry_check.py
```

Temporarily apply the hook to the external Afivo checkout, build
`programs/standard_3d`, and run the E1 configs with `OMP_NUM_THREADS=4`:

```bash
cd /home/helianthusczc/projects/afivo-streamer
make
programs/standard_3d/streamer3d \
  /home/helianthusczc/projects/streamer-rf-replica/solver3d/afivo_reference/stage_e/configs/e1_triangular_foil_medium_x.cfg
```

After the coarse/medium x/y profile runs, summarize:

```bash
.venv/bin/python solver3d/afivo_reference/stage_e/scripts/e1_compare_electrostatic.py
```

## Output Contract

Raw Afivo outputs are kept out of git under `results_raw/`. The E1 hook writes
compact CSV field-cell exports with:

`x_m,y_m,z_m,phi_V,Ex_Vpm,Ey_Vpm,Ez_Vpm,Eabs_Vpm,lsf_m,level`

The standard Afivo lineout files provide selected field profiles. Summary JSON
and CSV files under `results_summary/` are intended for git tracking and future
COMSOL comparison.

## Limits

Stage E1 is electrostatic only:

- `rho = 0`
- `ne = np = nn = 0`
- `photoionization = OFF`
- no avalanche, streamer propagation, bridge, current, Rb, RLC, RF, or Stage E2
