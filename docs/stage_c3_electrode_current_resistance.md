# Stage C3 Electrode Current and Plasma Resistance Diagnostics

Date: 2026-08-27

Scope: Stage C3 extends the Stage C2 dynamic real-electrode streamer tool with
terminal current, electrode charge, vacuum capacitance, effective plasma
conductance, and time-dependent breakdown resistance diagnostics.

Stage C3 does not implement RLC coupling, RF, FFT/ESD, Jefimenko fields,
antenna/receiver models, arc heating, gas heating, thermal ionization, or 3D.

## Sign Convention

Terminal currents use conventional current from electrode into gas as positive.

For each gas-electrode face, `n` points from the electrode into the gas. The
electron conventional current is:

```text
J_e = -e Gamma_e
I_cond = integral J_e . n dA
```

The implementation reuses the exact Stage C2 absorbing electron boundary flux
used by the transport update. The diagnostic does not compute a separate
terminal-current approximation. Electron absorption by an electrode gives a
positive conventional current from that electrode into the gas.

Axisymmetric face areas are:

```text
radial face: A = 2 pi r_face dz
axial face:  A = pi (r_outer^2 - r_inner^2)
```

## Displacement Current and Electrode Charge

Electrode charge is computed from the gas-side surface field:

```text
Q_electrode = epsilon0 integral E_normal dA
```

with `E_normal = E . n` and the same electrode-to-gas normal convention. The
displacement current is:

```text
I_disp = dQ_electrode / dt
```

The first dynamic step uses the initialized electrostatic field as the previous
state. If no previous electrode field exists, displacement-current fields remain
invalid instead of being silently treated as physical zero.

Very small `Delta Q` below the Poisson residual scale is treated as numerical
zero before division by `dt`. This prevents PETSc-level charge noise from being
amplified by very small explicit timesteps. The deadband is derived from the
configured elliptic tolerance and machine precision; it is not a physical
current threshold.

## Vacuum Capacitance

The diagnostic `C_gap_vacuum` is computed from a plasma-free electrostatic solve:

```text
C_gap_vacuum = |Q_HV(rho=0, V=1 V)| / 1 V
```

It is a vacuum geometry capacitance for validation and future Stage B/RLC
interfaces. It is not a plasma dynamic capacitance.

Measured C3 value for the current 70 um development geometry:

```text
C_gap_vacuum = 4.796651401354027e-16 F
```

## Effective Conductance and Resistance

Gas-cell electron conductivity is:

```text
sigma_e = e mu_e(E) ne
```

The dissipative effective plasma conductance is:

```text
P_cond = integral_gas sigma_e |E|^2 dV
Gb = P_cond / V_gap^2
Rb = 1 / Gb
```

`Gb` is not defined from `I_total / V`. `Rb` is not defined from
`V / I_total`. If `Gb = 0`, `Rb = +infinity`. If `|V_gap|` is below the
configured numerical voltage tolerance, `Gb_valid = false` and `Rb_valid =
false`; no floor conductance is introduced.

## Current Continuity Audit

The charge/current audit compares the plasma charge derivative with electron
particle current leaving the gas domain:

```text
dQ_plasma/dt ~= I_cond_HV + I_cond_ground + I_outer
```

where `I_outer = e * electron_particle_outflow_outer_boundary`. Reactions are
charge-balanced in the existing model. The residual is normalized by the larger
of the two compared current scales and a small numerical floor.

## Validation Runs

Vacuum linear ramp:

```text
dV/dt = 1.0e11 V/s
C_gap_vacuum = 4.796651401354027e-16 F
C_gap_vacuum * dV/dt = 4.796651401354027e-05 A
simulated last I_disp_HV = 4.7966514013037904e-05 A
relative error = 1.047324549120555e-11
I_cond_HV = 0
```

Vacuum constant voltage:

```text
V = 100 V
I_cond_HV = 0
I_disp_HV = 0
current_continuity_residual = 0
```

Uniform-conductance analytic control:

```text
ne = 1.0e15 m^-3
E = 2.0e6 V/m
V = 200 V
sigma = 7.7872471572728501e-06 S/m
Gb_theoretical = 3.5228676183806e-09 S
Gb_simulated = 3.5228676183806025e-09 S
relative error = 7.0440961922940607e-16
Rb = 2.8385965875711268e8 ohm
```

Case B current extraction demo:

```text
case = 500 V, SP3 on, 70 um gap, 5 um tip radius
seed = n0 1e16 m^-3, sigma 3 um, z_tip_offset -10 um
grid = 20 x 48
steps = 500
final_time = 3.0990278623269544e-12 s
I_cond_HV range = [5.9792169968588989e-09, 1.276711222837963e-07] A
I_disp_HV range = [0, 0] A after numerical Delta Q deadband
I_total_HV range = [5.9792169968588989e-09, 1.276711222837963e-07] A
P_cond range = [5.4649056368108279e-06, 2.2362590051508193e-05] W
Gb range = [2.185962254724331e-11, 8.9450360206032777e-11] S
Rb range = [1.1179384830834446e10, 4.5746444058619339e10] ohm
bridge = no
conservation_residual = 5.4150120655572395e-05
current_continuity_residual_max = 0.015150662958839685
```

The Case B result demonstrates the current/conductance/resistance extraction
pipeline. It is not a validated experimental 70 um breakdown resistance.

## MPI Consistency

The tiny 20-step current case was run with 1 and 2 MPI ranks. The first 20
current diagnostics agree within combined absolute/relative tolerances:

```text
max |I_cond_HV difference| = 3.621876486512436e-19 A
max |I_disp_HV difference| = 0
max |I_total_HV difference| = 3.621876486512436e-19 A
max |Q_HV difference| = 9.139679339018225e-24 C
max |Gb difference| = 1.7752039115469565e-21 S
max |Rb relative difference| = 8.148298967760852e-11
```

## Numerical Debt

`NUMERICAL_DEBT_CONSERVATION = OPEN`

The Stage C3 Case B conservation residual remains `5.415e-5`, matching the C2
Case B residual scale. Current diagnostics did not introduce a multi-order
degradation. Future production validation still needs mesh and timestep
sensitivity.

## Physical Applicability

Stage C3 `Rb(t)` primarily represents electron-dominated avalanche, streamer,
and early-breakdown conductance. The current fluid model does not include gas
heating, thermal ionization, electrode heating, arc hydrodynamics, mature spark
channel thermal conductivity, or ion drift current. Mature spark low-resistance
prediction requires future experiment-constrained extension and calibration.

## Deferred Scientific Validation

Formal scientific `Rb(t)` validation is deferred to:

```text
Stage B calibrated geometry
+ measured V(t)
+ experiment-constrained seed or emission model
```

`BRIDGE_DEMO = NO` for the current Case B interface demo.

