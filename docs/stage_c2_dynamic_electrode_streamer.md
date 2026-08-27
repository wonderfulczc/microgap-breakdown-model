# Stage C2 Dynamic Real-Electrode Streamer Foundation

Date: 2026-08-27

Scope: dynamic avalanche / streamer infrastructure under real axisymmetric
electrode conditions. Stage C2 does not compute electrode current, conductance,
resistance, RLC coupling, current moment, or RF output.

## Plasma domain

The Stage C1 geometry classes remain authoritative for cell classification:

- `Gas`
- `HighVoltageElectrode`
- `GroundElectrode`

Only `Gas` cells are part of the plasma fluid domain. In electrode cells:

- `ne = 0`
- `np = 0`
- `nn = 0`
- `rho = 0`
- `Sph = 0`

Electrode cells do not participate in plasma reaction, transport update,
photoionization source injection, or charge accumulation.

## Electron transport

Gas-gas faces continue to use the existing validated ISG-0 flux. The solver adds
only a degenerate zero-field guard at the call site: if face diffusion and drift
are both zero, the gas-gas flux is zero instead of calling ISG-0 outside its
valid input range.

Gas-electrode faces use an absorbing metal boundary with no emission. For a gas
cell adjacent to a metal face, the coordinate-direction electron flux is:

```text
positive face: F = max(v, 0) ne + 2 D ne / h
negative face: F = min(v, 0) ne - 2 D ne / h
```

where `v = -mu_e E_component`, `D` is the existing electron diffusion
coefficient, and `h` is the local cell spacing. The diffusion term is the
structured-grid one-sided absorbing-wall loss corresponding to a zero electron
density at the metal wall. No secondary emission, field emission, surface
chemistry, reflection, or empirical sticking coefficient is introduced.

Electron loss to HV and ground electrodes is recorded separately and included in
the electron conservation balance.

## Ion treatment

The frozen Stage 1-5 model does not transport positive or negative ions by
drift-diffusion. Stage C2 preserves that model:

- ions evolve by the existing gas-cell reaction model
- ions are zeroed inside electrode cells
- no ion mobility or ion electrode current is introduced

## Voltage coupling

In electrode-aware mode, each Poisson solve evaluates:

```text
V_applied = waveform.value(time)
```

The linear system then imposes:

- HV electrode: `phi = V_applied`
- ground electrode: `phi = 0`
- gas: Poisson equation with current space charge

The legacy uniform `background_field` is not added in electrode-aware mode.
Legacy Stage 1-5 configurations that do not provide electrode geometry and a
voltage waveform still use the old background-field path.

## Polarity and coordinate direction

The C2 development geometry uses `z = 0` at the grounded plane side and places
the high-voltage needle above it. For the positive-streamer development case:

- HV electrode: needle, `z_tip = 75e-6 m`
- ground electrode: plane starting at `z = 0`
- applied voltage: positive needle, ground plane at zero potential
- propagation direction: `-z`, from the needle tip toward the grounded plane
- seed definition: `z_seed = z_tip + z_tip_offset`

The default `z_tip_offset = -10e-6 m` places the Gaussian seed in the gas region
on the ground-facing side of the needle tip.

## Photoionization mask

SP3 may still be solved on the structured computational box, but the emission
source and resulting plasma source are masked so electrode cells cannot generate
`ne`, `np`, `nn`, or `Sph`.

## Diagnostics

Stage C2 adds dynamic diagnostics:

- time
- applied voltage
- Emax
- ne/np/nn maxima
- total electrons
- total charge
- conservation residual
- sigma_max, using `sigma_e = e mu_e ne`
- head position
- head velocity
- bridge flag
- absorbed electrons at HV and ground electrodes
- final physical time
- accepted step count
- `dt_min`, `dt_max`, and `dt_median`
- timestep controller counts and fractions for drift, diffusion, ionization,
  dielectric relaxation, and reaction limits
- initial and maximum `Emax`
- initial and final maximum `E/N`

`sigma_max` is only a local plasma conductivity diagnostic. Stage C2 does not
integrate conductance `Gb` or resistance `Rb`.

## Head diagnostic

The current C2 head diagnostic is an axisymmetric development diagnostic. It
uses the axis cell line and reports the lowest `z` gas cell where:

```text
ne >= head_ne_threshold
```

If no cell crosses the threshold, it reports the axis gas cell with maximum
`ne`. Head velocity is the finite difference of this diagnostic between accepted
steps. This is a simple numerical tracking definition, not a publication-grade
streamer-front definition.

For the current positive needle geometry this lowest-`z` rule corresponds to the
leading edge toward the grounded plane. Each development run writes axis
profiles at initial, selected intermediate, and final states so this scalar head
diagnostic can be checked against `ne(z)` rather than inferred from the scalar
alone.

## Bridge diagnostic

`bridge_flag` is an axial connectivity diagnostic. It is true only if every gas
axis cell between the grounded surface and HV tip has:

```text
ne >= bridge_ne_threshold
```

The threshold is stored in case metadata and is configurable. It is a numerical
development threshold, not an experimental breakdown criterion.

## Development case

The shared low-cost geometry for C2 development is
`stage-c2-dev-70um-needle-plane`:

- geometry: axisymmetric needle-plane
- gap: `70e-6 m`
- tip radius: `5e-6 m`
- shank radius: `2.5e-6 m`
- gas: Morrow-Lowke air baseline
- pressure: `101325 Pa`
- temperature: `300 K`
- seed: Gaussian, `n0 = 1e16 m^-3`, `sigma = 3e-6 m`, `z_tip_offset = -10e-6 m`
- grid: `20 x 48`, `rmax = 80e-6 m`, `zmax = 90e-6 m`
- output: compact diagnostics CSV plus initial, selected, and final snapshots

Two named development cases are retained:

- `stage-c2-avalanche-control`: SP3 off, early-avalanche numerical control.
- `stage-c2-positive-streamer-propagation`: SP3 on, intended positive-streamer
  propagation and bridging development case.

The old `3000 V` control is retained only as an overdriven diagnostic reference.
For propagation closure, candidate voltages are screened first by zero-charge
electrostatic solve and normalized by the project Morrow-Lowke breakdown field
`Ek = morrow_lowke_breakdown_field(...)`. Static COMSOL fields are not imposed
during dynamic streamer evolution.

## C2 Numerical Closure Repair

The final C2 numerical closure repair on 2026-08-27 identified the previous
yoctosecond timestep failure as a low-field reaction-controller interaction in
numerically inactive Gaussian-tail cells. The diagnostic 500 V SP3-on run
recorded the original limiting gas cell at `(i=19, j=47)`, with
`E = 3.534e-7 V/m`, `E/N = 1.445e-11 Td`, `ne = 9.225e4 m^-3`, and
`nu_att3 = 3.011e20 s^-1`. The resulting coefficient-only reaction limit was
`6.643e-22 s`, and the development app's `dt_scale = 1e-3` reduced the
accepted timestep to `6.643e-25 s`.

Only the electrode-aware Stage C timestep path was changed. The legacy
non-electrode Stage 1-5 path keeps its existing coefficient-based timestep
controller. In Stage C mode:

- conductor cells are excluded from drift, diffusion, ionization, dielectric,
  and reaction timestep scans;
- gas-cell species activity is derived from the existing numerical density
  scales, `1e-12 * max(ne_max, n_ref)` and the existing ISG-0 `n_ref` scale;
- inactive electron cells cannot control the global reaction timestep through a
  bare low-field attachment coefficient;
- reaction limiting is species-aware and uses the actual local
  `evaluate_reactions()` sources: `dt = safety * n_k / (-source_k)` for active
  species with negative source;
- the independent ionization growth limiter remains active for positive
  electron growth.

No Morrow-Lowke transport or reaction formula was changed. No `E/N`,
attachment, ionization, mobility, or diffusion clipping was introduced.

Post-repair retained development runs:

| case | voltage | SP3 | steps | final time | dt_min | dt_median | ne_max | bridge |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `stage-c2-avalanche-control-final` | 3000 V | off | 3000 | `1.562e-12 s` | `1.070e-16 s` | `4.279e-16 s` | `5.548e16 m^-3` | no |
| `stage-c2-500v-sp3-final` | 500 V | on | 500 | `3.099e-12 s` | `6.198e-15 s` | `6.198e-15 s` | `6.875e15 m^-3` | no |

The 500 V SP3-on development case now reaches ps-scale physical time without
NaN/Inf, timestep collapse, or negative-density retry exhaustion. Electron
population and the axis head diagnostic evolve during the run. The development
case does not bridge. Bridging validation is deferred to Stage B calibrated
geometry, measured `V(t)`, and experiment-constrained initial condition or
emission models.

`TRANSPORT_LOW_FIELD_VALIDITY = NEEDS_FUTURE_REVIEW`: the Morrow-Lowke low
`E/N` attachment behavior remains a physics-validity question for later review,
but it was not modified in Stage C2.

## Limitations

- Structured-grid electrode masks remain first-order near curved metal surfaces.
- No cut-cell metal boundary is implemented.
- Absorbing electron walls do not include secondary or field emission.
- Ions are not transported.
- The development bridge detector is an axis-line numerical diagnostic only.
- Quantitative tip-field and propagation claims require mesh sensitivity and
  Stage B electrostatic comparison.
- Quantitative propagation and bridging validation is deferred to Stage B
  calibrated geometry, measured voltage, and experiment-constrained seeds or
  emission.
