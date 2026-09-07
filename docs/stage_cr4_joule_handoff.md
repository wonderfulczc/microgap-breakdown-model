# Stage C-R4 Joule Energy And Cold-Thermal Handoff Contract

Stage C-R4 is a passive diagnostic layer. It does not modify the Stage C
continuity equations, Poisson solve, Morrow-Lowke transport, chemistry, SP3,
boundary conditions, timestep selection, current definitions, or `Gb/Rb`.

## Conductive Current Used For `J dot E`

The Joule diagnostic uses the existing flux-derived electron transport current
from `StreamerSolver::electron_transport_current_source()`:

```text
J_cond = -e Gamma_e
```

`Gamma_e` is reconstructed from the same finite-volume gas-gas and absorbing
gas-electrode electron fluxes used by `StreamerSolver::step()`. It includes the
current electron drift and diffusion transport represented by the Stage C
finite-volume update. The cell-centered diagnostic current is the average of
the adjacent face fluxes in each coordinate direction.

The current Stage C model does not evolve ion drift. C-R4 therefore does not
invent ion conductive current.

## Displacement Current Exclusion

Electrode diagnostics still report:

```text
I_total = I_cond + I_disp
```

Joule heating uses only conductive plasma current. Displacement current is not
used in `p_J` and contributes zero conductive Joule heating.

## Local And Integrated Joule Power

The local diagnostic power density is

```text
p_J = J_cond dot E
```

with units `W/m^3`. In the 2D axisymmetric solver it is evaluated as:

```text
p_J = Jr_RF * Er + Jz_RF * Ez
```

The gas-domain power is integrated only over gas cells:

```text
PJ_gas = integral_gas p_J dV
```

using `AxisymmetricGrid::cell_volume(i)`.

## Positive And Negative Power Policy

C-R4 does not apply `abs(J dot E)` and does not silently clamp local values.
Diffusion and face-to-cell reconstruction can produce local negative diagnostic
contributions. The output therefore separates:

```text
PJ_positive_W
PJ_negative_W
PJ_gas_W = PJ_positive_W + PJ_negative_W
```

`PJ_negative_W` is signed and non-positive.

## Accepted-Time Energy Accumulation

`JouleHandoffAccumulator` stores only the previous accepted diagnostic sample.
It accumulates energy with trapezoidal integration:

```text
QJ_n = QJ_(n-1) + 0.5 * (PJ_n + PJ_(n-1)) * dt_n
```

where `dt_n` is the actual accepted time difference between consecutive solver
states. The accumulator is diagnostic-only and does not affect timestep control
or PDE state.

## Numerical Channel Diagnostic Region

The optional channel region is a numerical diagnostic region, not a thermal
transition criterion:

```text
NUMERICAL_CHANNEL_DIAGNOSTIC_REGION:
sigma >= channel_sigma_relative_threshold * sigma_max
```

where

```text
sigma = e * mu_e(E) * ne
```

is the existing Stage C electron conductivity diagnostic used for `Gb/Rb`.
The default relative threshold is `0.1`. This threshold is not a universal
streamer, LFA, breakdown, or thermal-channel threshold.

For a valid channel region C-R4 reports:

```text
channel_volume_m3
channel_length_m
channel_effective_radius_m = sqrt(channel_volume / (pi * channel_length))
sigma_eff_S_m
E_channel_mean_V_m
ne_channel_mean_m3
ne_channel_max_m3
```

Finite channel geometry is not reported when the channel region is undefined.

## Handoff Timescales

C-R4 prepares audit quantities for a future cold-to-thermal handoff:

```text
tau_sigma = epsilon0 / sigma_eff
tau_evolution = |Gb / dGb_dt|
Xi_sigma = tau_sigma / tau_evolution
```

`dGb_dt` uses consecutive accepted-state `Gb` history. `Gb` and `Rb` remain the
frozen Stage C definitions from `conductance_diagnostics()`:

```text
Gb = Pcond / Vgap^2
Rb = 1 / Gb
```

C-R4 does not redefine `Gb/Rb` from the flux-derived Joule diagnostic.

## `Pi_H` Interface

The future thermal progress quantity is:

```text
Pi_H = QJ_channel / Q_required
```

Stage C does not contain a traceable `Q_required`, target temperature,
gas-heating model, or mature-arc model. C-R4 therefore outputs:

```text
thermal_energy_reference_status = NOT_AVAILABLE
Pi_H = NaN
handoff_status = UNRESOLVED_CALIBRATION
```

It never outputs `HANDOFF_READY`, `THERMAL`, or `SPARK_TRANSITION`.

## Stage G0 Data Contract

With `stage_c2_dynamic --joule-handoff`, C-R4 writes
`cold_thermal_handoff.csv` containing conductive Joule power, cumulative Joule
energy, numerical channel summaries, `Gb/Rb`, `tau_sigma`, `tau_evolution`,
`Xi_sigma`, and unresolved thermal-reference status. The file is compact and
does not dump per-cell Joule fields.

## Limitations

The current Stage C model is a cold-fluid discharge model. It does not include
gas heating, thermal ionization, electrode heating, arc hydrodynamics, mature
spark channel conductivity, or ion drift current. C-R4 therefore does not
determine a cold-to-thermal transition threshold. That calibration belongs to
future Stage G0/G1 work.
