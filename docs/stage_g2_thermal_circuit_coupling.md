# G2: bidirectional thermal-channel / circuit reference coupling

## Scope and provenance

G2 validates the bidirectional thermal-channel / circuit coupling architecture.
The current coupled reference run is not yet a calibrated representation of the
experimental microgap post-breakdown waveform. No G3 port model is implemented.

G1 checkpoint: `4a45d92` (`Stage G1: add 1D radial thermal spark-channel solver`).
The existing G1 properties, radial solver and historical circuit files are
unchanged; their hashes are recorded in `thermal/g2/g2_reference_summary.json`.
The G0 production candidate remains `G0_THERMAL_INITIALIZATION_UNRESOLVED`.
The reference retains `NUMERICAL_REFERENCE_INITIAL_CONDITION`,
`LTE_APPLICABILITY_PENDING_CALIBRATION`, and `RADIATION_MODEL_NOT_ENABLED`.

## Actual canonical topology and adapter

`python/streamer_rf/circuit/model.py:SeriesRLCGapCircuit.rhs` implements
`Vs -> Rs -> L -> node`, with `Cext || Cgap || Gsp` from node to ground:

```
L dI_L/dt = Vs - Rs I_L - Vgap
(Cext+Cgap) dVgap/dt = I_L - Gsp Vgap
I_cond = Gsp Vgap
I_disp = Cgap dVgap/dt
I_Cext = Cext dVgap/dt
```

`Cgap` is an explicit parallel gap capacitance, not an extra series capacitor.
The conductive branch alone heats the gas: neither `I_L` nor
`I_gap_total=I_cond+I_disp` is the thermal drive. Signed positive branch current
flows from node to ground. Displacement power is capacitor storage, not heat.

`thermal/coupling.py:resistance_conductance` validates finite positive ohms and
returns the existing `ConductanceProfile` in siemens. `circuit_interval` invokes
the **unchanged canonical RHS**, with DOP853 as in the old solver. It augments
the two circuit states only with quadratures of `Gsp Vgap^2`, `Rs I_L^2`, and
`Vs I_L`. This avoids imposing dense output for energy integration; it is not
a second circuit RHS. The existing `evaluate_solution` supplies current,
KCL/KVL and separate inductor/Cext/Cgap energy diagnostics. The old kernel's
output-trapezoid energy diagnostics remain intact for historical validation.

## Thermal resistance and multirate acceptance

The unchanged G1 `joule_heating` uses radial cell areas per unit length:

```
G_length = sum(sigma_i A_i)                 [S m]
Rsp = L_channel/G_length                   [ohm]
Ez = I_cond/G_length                       [V/m]
qJ_i = sigma_i Ez^2                        [W/m^3]
Pthermal = L_channel sum(qJ_i A_i)          [W]
```

The fixed PRIMARY coupling cap is 50 ps, with only a 25 ps alternative. The
thermal interval also obeys the original G1 hydro/conduction/Joule constraints.
Each interval freezes `Rsp_n` from the accepted thermal state, advances the
circuit with its own adaptive internal steps (maximum one fourth of the
thermal interval), then supplies

```
I_rms^2 = integral(I_cond^2 dt)/dt = Qsp_circuit/(Rsp_n dt)
```

to the original G1 explicit thermal step. This is an energy-consistent RMS
drive, not a mean signed current. Current sign remains in the circuit outputs.
G1 independently evaluates `sigma*(I_rms/G_length)^2` on its original state.
No thermal energy is rescaled or overwritten to match the circuit integral.

If G1 accepts a shorter interval, discard the circuit and thermal trial and
repeat both on that shorter interval. Commit states and all integrals only
together. Tests force both a shortened trial and a failed trial to check this.
The adapter owns a shallow copy of G1; current G1 stepping replaces arrays
rather than modifying them in place. That dependency and `_rhs()` access are
explicit integration debt, not a change to the frozen solver.

The mode is `LAGGED`: halving the cap passes the 5% engineering target, so no
predictor-corrector is justified or implemented. At an accepted boundary the
exported `Rsp_ohm` and `I_A=V_channel/Rsp_circuit` refer to the **updated** state;
`interval_Rsp_ohm`, `I_rms_interval_A`, and `P_sp_*` refer to the preceding
accepted interval. Powers are interval averages, not endpoint samples. The
initial row has zero interval duration and no accumulated interval power.
Resistance is piecewise held during circuit integration; endpoint feedback
therefore has small numerical jumps, not an assertion of continuous substep R.

## Energy identities

```
Ecircuit = L I_L^2/2 + Cext Vgap^2/2 + Cgap Vgap^2/2
residual_circuit = Ecircuit-Ecircuit0 + Qfixed + Qsp - Wsource
residual_thermal = Einternal+Ekinetic-Ethermal0-QJ+Qconduction+Qhydro
residual_cross = Qsp_circuit-Qsp_thermal
```

The present circuit is source-free, but source work is retained and tested with
a driven steady circuit. G1 radiation loss is zero because no radiation model
is enabled. Signed outer fluid energy flux and conduction losses are preserved.

## Reference and bounded comparisons

The actual `conductive_gap` case in
`circuit/validation/generate_stage_g1_validation.py` supplies `L=1 uH`,
`Cext=0.8 pF`, `Cgap=0.2 pF`, `Rs=2 ohm`, `Vs=0`, `I_L(0)=0`, `Vgap(0)=100 V`.
Its prescribed conductance is replaced by the G1 thermal resistance as required.
The analysis uses the first 10 ns; no circuit parameter tuning is performed.

The exact G1 synthetic profile is reused: `Rmax=0.5 mm`, length `1 mm`,
`f=exp(-(r/50 um)^2)`, `T=300+9700f K`, `p=101325(1+f) Pa`, `u=0`.
PRIMARY has 150 cells; the only finer-grid comparison has 225 cells.

| Quantity | PRIMARY | Half interval | Finer grid | Fixed-R replay |
|---|---:|---:|---:|---:|
| Peak absolute conductive current, A | 0.414789 | 0.414789 | 0.414319 | 0.414789 |
| Circuit channel energy, nJ | 4.997973 | 4.997970 | 4.997928 | 4.995359 |
| Final axis T, K | 7855.285 | 7856.554 | 7962.067 | 7855.425 |
| Final actual thermal Rsp, ohm | 624.289 | 624.127 | 576.631 | 624.241 |
| Accepted thermal/coupling steps | 201 | 400 | 234 | 201 |
| Circuit internal steps | 867 | 1647 | 1065 | 861 |

The 201st PRIMARY step is a roundoff-sized final remainder, not another physical
event. No rejected interval occurred in these reference runs.
Halving the interval changes peak current by 0%, Qsp by 0.0000666%, final T by
0.0162%, and final Rsp by 0.0260%. The 225-cell comparison changes these by
0.1134%, 0.000908%, 1.359%, and 7.634%, respectively. Preserve
`THERMAL_RSP_GRID_SENSITIVITY_PRESENT`: effects below this uncertainty cannot
be promoted to new physics. The waveform comparisons use bounded interpolation
to PRIMARY accepted timestamps and the discrete sample-vector L2 norm.

The PRIMARY axis cools from 9989.228 to 7855.285 K while the numerical thermal
radius grows from 75.9425 to 80.2006 um. Axis conductivity falls from 2888.438
to 789.135 S/m; resistance rises from 241.086 to 624.289 ohm. Maximum adjacent
resistance and axis-temperature changes are 2.093 ohm and 14.904 K. Current
ranges from -0.017760 to 0.414789 A and voltage from -5.3163 to 100 V on the
accepted grid. Density, pressure and temperature remain positive and finite.
Joule input is only about 2.4e-5 of initial thermal energy, so this case does
not demonstrate heating-dominated expansion or an experimentally important
feedback magnitude. No physical parameter was changed to force such behavior.

In the `FIXED_RSP_REPLAY` baseline only circuit resistance is fixed; thermal
conductivity still evolves under its replayed current. Its actual thermal R
must not be confused with the fixed `Rsp_circuit_ohm=241.086`. Its waveform
differs from PRIMARY by 4.51% in discrete normalized L2. Its thermal input is
5.145004 nJ versus circuit loss 4.995359 nJ, a 2.996% inconsistency: this
deliberately one-way replay is **not** an accepted energy-coupled solution.
Peak current is unchanged because it occurs at the shared initial state.

PRIMARY maximum independent interval-power mismatch is 6.0e-16 relative.
Final cumulative cross-domain mismatch is zero at printed precision. Maximum
circuit energy residual is 6.34e-24 J (1.27e-15 normalized to 5 nJ initial
circuit energy); maximum thermal residual is 7.25e-20 J. These numerical
identities do not remove the lagged resistance-time discretization error.

## Outputs, verification and boundaries

Run `.venv/bin/python thermal/g2/generate_g2_reference.py`. Compact output is
`thermal/g2/g2_reference_summary.json`, four accepted-time CSV files, and only
initial/final PRIMARY radial profiles; no per-circuit-step profiles are saved.
The 945282 bytes of derived data include all allowed comparisons. PRIMARY
runtime was 4.28 s, peak process RSS 140164 KiB; all four configurations took
23.04 s with overall peak 146168 KiB. Later per-case RSS values are process
high-water marks, not isolated measurements.

`tests/thermal/test_stage_g2_coupling.py` adds 21 tests, including old fixed and
prescribed conductance equivalence, nonzero source work, zero/constant-current
energy, independent power closure, resistance feedback, rejection accounting,
invalid inputs, G0 rejection and coupling refinement. Build passed; CTest 10/10;
full pytest 223 passed, 3 skipped, 5 deselected. Skips are the existing missing
measurement and sandbox MPI cases. All 15 old circuit tests are included.

Framework status: `BIDIRECTIONAL_REFERENCE_COUPLING_VALIDATED`.
Production status: `PRODUCTION_THERMAL_RLC_COUPLING_NOT_CALIBRATED`.
G0 remains `HANDOFF_CALIBRATION_PENDING`; G1 remains
`THERMAL_REFERENCE_SOLVER_VALIDATED` with its frozen LTE/radiation limitations.
Stage C/D/E/F and G0/G1, old circuit sources and trusted outputs were not edited.

The exported signed conductive `I_A`, `V_channel_V`, `Rsp_ohm` and
`Rsp_circuit_ohm` are compact inputs for later G3 work. No differential port
impedance, antenna excitation, external full-wave model or G3 solver is defined
here. Positive finite thermal conductance is required; a nonconductive state
fails explicitly instead of inventing a leakage floor or thermal initialization.
