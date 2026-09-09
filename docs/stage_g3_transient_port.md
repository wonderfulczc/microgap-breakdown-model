# G3: physics-derived transient port

G3 provides a physics-derived transient port representation for Stage H.
It does not yet perform full-wave propagation and it does not establish a
calibrated experimental post-breakdown waveform. G2 was checkpointed as
`8521aa2`; no G0/G1/G2 or canonical circuit source was modified.

## Reference plane and conventions

The actual equations in `circuit/model.py:SeriesRLCGapCircuit.rhs` are

```
L dI_L/dt = Vs - Rs I_L - Vgap
(Cext+Cgap) dVgap/dt = I_L - Gsp Vgap
```

The topology is `Vs -> Rs -> L -> node`, with parallel `Cext`, `Cgap`, `Gsp`
to ground. The G3 boundary leaves Cext on the external lumped-circuit side;
only `Cgap || Gsp` is inside the gap-side reference representation.

```
reference_plane_id = EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP
V_port = Vgap = node potential - ground potential
positive current = external circuit -> physical gap/electrode structure
I_port = I_L - I_Cext
I_port = I_sp + I_Cgap
I_sp = Gsp Vgap
I_Cext = Cext dVgap/dt
I_Cgap = Cgap dVgap/dt
```

Branch currents come from the unchanged canonical `evaluate_solution`, not
from numerical differentiation of saved voltages. The port implementation
computes total current by subtracting the external-capacitor branch, then
checks it independently against conduction plus gap displacement current.
No sign is changed to force positive power. G2 still heats only with I_sp.

## Power, energy and impedance

```
P_port = V_port I_port
P_conduction = V_port I_sp
E_Cgap = Cgap V_port^2/2
integral P_port dt = integral P_conduction dt + Delta E_Cgap
```

CSV energy checks use trapezoidal quadrature on the exported samples and
compare conduction energy independently with G2's accepted ODE quadrature.
No currents or energies are rescaled. Signed negative port power is retained.
The normalized identity residual divides by the largest absolute magnitude
of integrated port energy, integrated conduction energy and capacitor-energy
change, avoiding cancellation in the denominator.

`Rsp_ohm` is the positive thermal conductive resistance held in the corresponding
G2 circuit interval. `Z_dynamic_ohm=V_port/I_port` is an instantaneous terminal
ratio, **not** that resistance, a small-signal impedance, or an S-parameter.
Its numerical gate is

```
abs(I_port) > max(1e-12 A, 1e-8 * max(abs(I_port)))
```

computed once over each complete waveform. Below it the ratio is NaN with
`Z_dynamic_valid=false`; overflow is also invalid, never Inf. This is a
numerical division gate, not a physical transition threshold. Zero voltage
with meaningful nonzero current gives a valid zero ratio. No spectral port
impedance is added in this node.

## Circuit-resolution recovery and uniform waveform

G2 saved 202 thermal-boundary rows, including its roundoff-sized last interval,
but did not retain internal circuit states. The new generator reads only
`thermal/g2/g2_primary_timeseries.csv` and its small reference summary. Each
interval is replayed with the **same held `interval_Rsp_ohm`**, initial I_L/Vgap,
source voltage, circuit parameters, DOP853 electrical tolerances and maximum
internal step (one fourth of the thermal interval). It calls the original
circuit RHS with dense output. This is bounded circuit diagnostic recovery,
not a new coupled run, thermal replay, parameter study or changed feedback.

The replay omits G2's extra energy-quadrature states, so its accepted internal
nodes need not be identical to the original unsaved nodes. Electrical endpoint
agreement is checked: maximum I_L error `9.41e-17 A`, V error `1.42e-14 V`.
The recovered internal-step data contains 1030 rows, including both sides of
each held-R switch. `sample_side` and `interval_index` make duplicate boundary
timestamps explicit. They contribute zero integration width; do not deduplicate
them before native energy integration or blend across the resistance jump.

Median nonzero resolved internal spacing is 12.5 ps. A uniform 0--10 ns grid
uses this spacing (801 points), evaluating the bounded dense solution of its
own interval. At a switch it uses the right interval; at the final endpoint
it uses the final left limit. There is no temporal extrapolation. The
roundoff-sized final interval is excluded only from selection of typical dt,
not from recovery or the endpoint energy check.

Nyquist is 40 GHz. This is a sampling metadata value, **not** a physically
trusted RF bandwidth or an upgrade to any frozen Stage F trust status. No
full-wave response or anti-aliasing validity claim is made.

## Frozen G2 reference result

| Diagnostic | Uniform waveform result |
|---|---:|
| V_port range, V | -5.316852 to 100 |
| I_port range, A | -0.01803782 to 0.33183159 |
| I_sp range, A | -0.01787689 to 0.41478948 |
| I_Cgap range, A | -0.08295790 to 0.000210514 |
| Signed port power range, W | -0.00148243 to 33.18316 |
| Integrated port energy, J | 3.99939477e-9 |
| Integrated conduction energy, J | 4.99971635e-9 |
| Delta E_Cgap, J | -9.99999006e-10 |
| Port identity residual, J | -3.22566491e-13 |
| Port identity relative residual | 6.45170e-5 |
| Maximum KCL residual, A | 4.16334e-17 |

Native exported quadrature gives port/conduction energies 4.00128704e-9 /
5.00208459e-9 J. Uniform versus native differences are 0.04729% / 0.04735%;
peak absolute voltage and current are unchanged. Uniform conduction energy
differs from the frozen G2 4.99797292e-9 J by 0.03488%. All are below the 1%
engineering acceptance target. The native identity residual is 0.01596%.

All 801 sampled ratios pass the current gate; this does not mean that no
continuous-time zero crossing exists between samples. Z_dynamic p10/median/p90
are 301.358 / 458.154 / 655.117 ohm. No extreme ratio near a crossing is promoted
as physics. Separately, thermal Rsp spans 241.086--624.289 ohm.

## Stage H contract and mandatory limitations

`thermal/port.py:PortTransientContract` validates finite branch waveforms,
uniform increasing timestamps, KCL, reference-plane subtraction, resistance
availability and impedance masks. Its JSON metadata binds the uniform CSV to
the source case, SHA256 hashes of frozen G2 inputs, time interval, dt, sample
count, polarity, reference plane, validity and energy statuses. The CSV has
`waveform_valid` and `Z_dynamic_valid` separately.

Method status: `PHYSICS_DERIVED_PORT_MODEL_VALIDATED`.
Production status: `PRODUCTION_PORT_TRANSIENT_NOT_CALIBRATED`.
G0 remains `HANDOFF_CALIBRATION_PENDING`; G1 remains
`THERMAL_REFERENCE_SOLVER_VALIDATED` with `LTE_APPLICABILITY_PENDING_CALIBRATION`;
G2 remains `BIDIRECTIONAL_REFERENCE_COUPLING_VALIDATED` and
`PRODUCTION_THERMAL_RLC_COUPLING_NOT_CALIBRATED`. Thermal resistance still has
`THERMAL_RSP_GRID_SENSITIVITY_PRESENT`. These are inherited, not repaired by
the small port algebraic residual.

**H_GAP_CAPACITANCE_PARTITION_REQUIRED**: the G2 reference already contains
lumped Cgap. Stage H must explicitly choose whether full-wave electrode-gap
geometry resolves this capacitance or Cgap remains lumped. It must not blindly
add an identical lumped capacitor at the same reference plane when geometry
already represents it.

**FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED**: imposing V_port as an ideal source
is a one-way handoff. V and I describe the same original load-dependent
transient, not two independently enforceable boundary conditions on an
arbitrary new load. A changed full-wave load does not feed back into G2 unless
Stage H later embeds the circuit/source impedance. G3 does not solve this.
No Ryu-style step baseline, openEMS run or Stage H development is included.

## Files, resources and checks

Run `.venv/bin/python thermal/g3_port/generate_g3_port.py`.
Outputs: `g3_port_summary.json`, `g3_port_timeseries.csv` (native side-labelled
nodes), `g3_port_uniform.csv` (handoff waveform), all under `thermal/g3_port/`.
The derived data is 565345 bytes. Extraction/validation took 0.59 s excluding
interpreter import and CSV write overhead; peak process RSS was 143752 KiB.

17 new tests cover pure R/C/parallel RC, energy/KCL identities, signed power,
polarity reversal, zero crossings, impedance gates, invalid contracts, temporal
sampling and unresolved production-status propagation. Build passed; CTest
10/10; full pytest 240 passed, 3 skipped, 5 deselected. The skips remain the
unexecuted measurement runner and sandbox MPI cases. All tracked frozen files
remain unchanged after the G2 checkpoint; G3 code/docs/results are new files.
