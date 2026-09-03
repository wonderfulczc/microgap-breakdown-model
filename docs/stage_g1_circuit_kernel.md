# Stage G1 Circuit Kernel

Stage G1 adds a reference Python circuit kernel only. It does not run Stage C
production plasma data, does not perform FFT/RF analysis, and does not start
self-consistent plasma-circuit feedback.

## Topology

The canonical benchmark topology is

```text
Vs(t) -> Rs -> L -> gap node
gap node -> ground through Cext || Cgap || Gb(t)
```

The discharge gap branch is frozen as

```text
i_gap = Gb(t) Vgap + Cgap dVgap/dt
```

where `Gb(t)` is the Stage C effective electron conductance interface and
`Cgap` is a constant gap capacitance. Stage G1 uses synthetic `Gb` and `Cgap`
values only. `Rb = V / Itotal` is not used.

## State Equations

The ODE state is

```text
x = [I_L, Vgap]
```

KVL around the driven series branch gives

```text
Vs(t) - Rs I_L - L dI_L/dt - Vgap = 0
```

therefore

```text
dI_L/dt = [Vs(t) - Rs I_L - Vgap] / L
```

KCL at the gap node gives

```text
I_L = Cext dVgap/dt + Cgap dVgap/dt + Gb(t) Vgap
```

therefore

```text
dVgap/dt = [I_L - Gb(t) Vgap] / [Cext + Cgap]
```

The diagnostic currents are stored separately:

```text
I_gap_cond = Gb(t) Vgap
I_gap_disp = Cgap dVgap/dt
I_gap_total = I_gap_cond + I_gap_disp
I_Cext = Cext dVgap/dt
I_source = I_L
```

This avoids double counting displacement current: `Cgap dVgap/dt` appears only
inside the gap branch, while `Cext dVgap/dt` is the external capacitor current.

## Interfaces

`Gb(t)` accepts:

- constant conductance;
- piecewise-linear samples;
- CSV with `time_s,Gb_S`;
- CSV with `time_s,Rb_ohm`, internally converted to `G=1/Rb`.

Time samples must be strictly increasing. The default extrapolation mode is
`invalid`; `hold` is available only when explicitly configured.

The voltage waveform uses the same piecewise-linear `time_s,voltage_V` semantic
as Stage C.

## Energy Audit

The stored energy is

```text
U = 0.5 L I_L^2 + 0.5 Cext Vgap^2 + 0.5 Cgap Vgap^2
```

The cumulative losses and source work are

```text
W_gap = integral Gb(t) Vgap^2 dt
W_Rs = integral Rs I_L^2 dt
W_source = integral Vs(t) I_L dt
```

The reported balance residual is

```text
epsilon_E = [U(t) - U(0)] + W_gap + W_Rs - W_source
```

## Coupling Validity Hook

If a prescribed/reference gap voltage is available, Stage G1 reports

```text
ONE_WAY_COUPLING_ERROR =
max |Vgap_circuit - Vgap_reference| / max |Vgap_reference|
```

No universal threshold is assigned in Stage G1.

## ngspice

Stage G1 only records `which ngspice` and `ngspice --version`. It does not
install ngspice and does not depend on it.
