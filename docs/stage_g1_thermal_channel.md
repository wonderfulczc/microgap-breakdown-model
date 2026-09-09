# G1 Radial Thermal Spark-Channel Reference Solver

G1 independently implements a one-dimensional cylindrical LTE reference model
in `python/streamer_rf/thermal`. The historical `streamer_rf.circuit` kernel and
its older `docs/stage_g1_circuit_kernel.md` name are preserved. G0 checkpoint is
`839dace`. No frozen C/D/E/F/G0 equation or output is changed.

## Property Provenance And Transcription

A. D'Angola, G. Colonna, C. Gorse, M. Capitelli, "Thermodynamic and transport
properties in equilibrium air plasmas in a wide pressure and temperature range",
European Physical Journal D 46, 129-150 (2008),
[DOI 10.1140/epjd/e2007-00305-4](https://doi.org/10.1140/epjd/e2007-00305-4).

The [public page-numbered full text](https://www.researchgate.net/profile/Antonio-Dangola/publication/200702842_Thermodynamic_and_transport_properties_in_equilibrium_air_plasmas_in_a_wide_pressure_and_temperature_range/links/0fcfd509786b669986000000/Thermodynamic-and-transport-properties-in-equilibrium-air-plasmas-in-a-wide-pressure-and-temperature-range.pdf)
resolves the earlier provenance stop; a local PDF is not required. Neither the
article nor its figures are committed. `thermal/g1/dangola_property_provenance.json`
records equation/table/page mapping, conventions and coefficient SHA256.

| Property | Equation | Table | Journal pages | Paper units |
|---|---|---|---|---|
| Molar mass |53|21|137,146|kg/mol|
| Density |54|21|137,146|kg/m3|
| Enthalpy |55|22|137,147|cal/g|
| cp |56|23|137,147-148|cal/g/K|
| Electrical conductivity |58|25|137,149|S/m|
| Thermal conductivity |59|26|137,149-150|W/m/K|

On page137, Eq.48 is `exp(-((T-c)/Delta)^2)`, Eq.49 is
`expit(2*(T-c)/Delta)`, Eq.50 is `a-c*exp(-(T/Delta)^w)`, and Eq.51 is
`a*T^w`. Eq.52 uses an ascending pressure polynomial. The text identifies natural
logarithms: all log operations use ln, with pressure input `ln(p_Pa/101325)`.
Density uses Pa and `R=8.31446261815324 J/(mol K)` in `rho=p*M/(R*T)`.
Enthalpy and cp are multiplied by 4184 using the thermochemical calorie.

- Table22 alpha6..alpha9 continue log(a7), log(c7), log(Delta7), before
  Table23 on page147. They are appended to the same rows, not shifted to Table23.
- Table23 sigmoid3 has alpha5/alpha6 continuation; Gaussian6..11 continue on p148.
- Table25 explicitly gives `a7=a3+a4+a5-a1-a2-a6`. Rows marked -a3,-a4,-a5
  yield negative signed amplitudes; the relation is not inferred by tuning.
- Table26 Gaussian7..10 continue on p150; terms3 and8 have negative amplitudes.
  Its first sigmoid parameters are linear, not exponentiated.
- Eq.58 applies xi to ln(T), but sigmoids to T. Eq.59 uses ln(T) for all terms.

Valid inputs are 50-60000 K and 0.01-100 atm inclusive. Outside inputs return
`PROPERTY_OUT_OF_RANGE` and NaN properties; invalid fit results return
`PROPERTY_FIT_INVALID`. No extrapolation is used. Production properties are
`PUBLISHED_EQUILIBRIUM_AIR_FIT`; no competing species chemistry is implemented.

## EOS And Bounded Inversion

`EquilibriumAirProperties.evaluate(T,p)` returns M,rho,h,cp,k,sigma,e and validity.
Conservative thermochemical energy uses `e=h-p/rho`. Independent cp and h fits
are not forced to have identical derivatives. Their largest discrepancy from
`dh/dT|p` at the seven requested atmospheric checkpoints is about 0.622%.

`invert(rho,e,T_guess,p_guess)` uses vectorized damped Newton, with the previous
accepted state as normal guess, then bounded log(T),log(p) least-squares root
search as fallback. All candidates remain within the published domain.
Success requires a forward EOS residual <=2e-8 and valid properties, independently
of optimizer termination. Failure raises an explicit input/inversion/property
error. Only initial guesses may be projected; conserved rho and energy are
never clipped.

Bounded finite differences of rho(T,p),e(T,p) determine sound speed and effective
constant-density heat capacity. Let rt,rp,et,ep denote their derivatives:

```text
det = rt*ep-rp*et
cv_effective = et-ep*rt/rp
sound_speed^2 = (-et+rt*p/rho^2)/det
```

This follows from `de=p/rho^2 drho` on an adiabat and inversion of the same EOS
Jacobian. No constant gamma or invented cv is used. Nonpositive derivative
results stop the step. Independent approximate fits do not guarantee exact
thermodynamic identities over all possible paths; validation covers tested states.

## Cylindrical Conservative Scheme

Uniform radial faces cover [0,Rmax]. Cell locations are midpoints; per-length
volumes are `V_i=pi*(r_right^2-r_left^2)` and face areas `A_f=2*pi*r_f`.
The state and flux are:

```text
U = [rho, rho*u, rho*(e+u^2/2)]
F = [rho*u, rho*u^2+p, (rho*(e+u^2/2)+p)*u]
dU_i/dt = -(A_right*F_right-A_left*F_left)/V_i + S_i
S_momentum = p_i*(A_right-A_left)/V_i
S_energy = qJ-(H_right-H_left)/V_i
```

Rusanov flux uses the maximum adjacent abs(u)+sound_speed. The momentum
geometry term cancels pressure flux for uniform stationary gas. Axis velocity
is reflected with zero area and heat flux. The reference outer ghost is fixed
ambient gas, with conductive ambient temperature over a half cell. This is not
a general nonreflecting boundary; the disturbance must remain well inside the
domain. Closed tests use reflected outer velocity and zero outer heat flux.

Forward Euler evaluates all terms at the same accepted state, without operator
splitting. dt is limited by geometric CFL=0.3, conduction
`0.1*rho*cv_effective*dr^2/k`, and 5% local internal-energy Joule increment.
Invalid trials retry at half dt, at most 16 trials. Rejections do not accumulate
energy or time. No negative-state clipping is used. These are numerical controls.

Radial heat flux is `H_f=-A_f*k_f*(T_right-T_left)/dr`, with harmonic face k.
Radiation is explicitly `RADIATION_MODEL_NOT_ENABLED`; no empirical loss is added.

## Current Drive And Energy Accounting

Current accepts a callable I(t) or bounded `time_s,current_A` CSV, using the
existing piecewise-linear signal convention with explicit support checks.

```text
G_length = sum(sigma_i*V_i)              [S m]
Rsp = L_channel/G_length                [ohm]
Ez = I/G_length                         [V/m]
qJ_i = sigma_i*Ez^2                     [W/m3]
PJ_volume = L_channel*sum(qJ_i*V_i) = I^2*Rsp
```

Displacement current is excluded. Zero current gives zero qJ; nonzero current
with zero conductance raises an error. Tiny finite conductivity may make dt
prohibitive: this model is not intended to create a thermal channel from cold gas.

Budgets track internal/kinetic energy, accepted Euler Joule input, conductive
and hydrodynamic boundary losses, and zero radiation. The residual is
`E(t)-E(0)-QJ+Q_conduction_out+Q_hydro_out`. Open-boundary flow energy cannot be
omitted. Signed negative boundary loss denotes inward flow. Mass has a separate
boundary-flux residual. Algebraic closure does not imply negligible CFD error.

## G0 Contract And Future Energy Reference

The real G0 result retains `G0_THERMAL_INITIALIZATION_UNRESOLVED`. Neither its
invalid candidate nor its latest cold observables supplies a radial temperature.
Future `validate_calibrated_radial_state` requires T,p arrays, grid and provenance.
`required_thermal_energy(mass,cold_e,target_e)` computes a caller-specified,
matched-material internal-energy difference. It does not choose a target or
include unknown work/losses. Q_required, Pi_H*, Xi_sigma* and T* remain uncalibrated.
Profiles, Rsp, G_length and budgets are available for future G2; no G2 exists here.

## Numerical Reference And Results

Run `.venv/bin/python thermal/g1/generate_g1_reference.py`. One fixed case is
evaluated on 100,150,225 cells; 150 is primary. Only its initial/mid/final profiles
are saved along with compact time series and grid/summary tables:

```text
Rmax=500 um; L_channel=1 mm; t_end=10 ns
ambient: 300 K, 101325 Pa
f(r)=exp(-(r/50 um)^2)
T(r,0)=300+9700*f K; p(r,0)=101325*(1+f) Pa; u(r,0)=0
I(t)=0.05*(1+0.5*sin(pi*t/10 ns)^2) A
```

This is `NUMERICAL_REFERENCE_INITIAL_CONDITION`: an overpressured hot core in
atmospheric ambient air, not an isobaric or calibrated handoff. The diagnostic
radius is the first outward crossing of 10% of instantaneous peak temperature
excess, interpolated between cells. It is not a material surface or spark
threshold. Axis observables refer to the innermost cell, not extrapolation to r=0.

The primary run reaches 10 ns in 103 steps without rejection. T_axis drops from
about 9989 to 7870 K, p_axis to 134965 Pa, and maximum outward speed reaches
169.94 m/s. Radius grows from 75.94 to 80.17 um; Rsp rises from about 240 to
620 ohm. Thermal relaxation dominates this weak-current reference; cooling is
not evidence of an incorrect heating sign. Final QJ=1.7291e-8 J.

Final medium/fine differences: T_axis 1.39%, radius 0.316%, Rsp 8.32%, total
energy 0.00139%. Differences shrink relative to coarse/medium. Rsp is least
converged; these are usable trends, not high-precision resistance predictions.
First-order numerical diffusion remains significant. The outer 20% shell has
no pressure change at stored precision, so no second domain run is needed.
Joule closure is below 6e-16; budget residual divided by initial energy is below
4e-16. The JSON also reports residual divided by QJ, to avoid hiding it behind
ambient stored energy. These are numerical conservation identities.

| T (K) | rho (kg/m3) | cp (J/kg/K) | k (W/m/K) | sigma (S/m) |
|---|---|---|---|---|
|300|1.17033|1027.41|0.026111|3.0965e-57|
|1000|0.351083|1121.80|0.071758|1.3805e-19|
|3000|0.114506|2740.87|0.478017|0.023543|
|5000|0.058431|2896.85|0.756226|24.7193|
|10000|0.017135|4888.09|1.361828|2991.46|
|15000|0.007631|22296.89|3.479949|7332.13|
|30000|0.002382|31110.63|5.515685|12060.4|

These are fit evaluations, not experimental points. Magnitudes agree with the
indexed Figure6 axis scale (up to 15000 S/m) and Figure7 scale (up to 7 W/m/K
through 30000 K). Direct figure images remain inaccessible; no digitized or
pointwise overlay agreement is claimed. The Gaussian/transition structure is
retained exactly. Readable original-curve visual comparison remains a limitation.

Tests cover coefficients, ranges, continuity, poor-guess/ordinary EOS roundtrips,
derivatives, geometry, uniform gas, current/Joule identities, constant-power
deposition, closed energy/mass conservation, Neumann Bessel-mode conduction
decay/convergence, axis symmetry, invalid states, current CSV bounds, G0 rejection,
explicit targets and rejected-step accounting.

The single-ionization Saha test is hydrogen-only and
`EQUILIBRIUM_SANITY_CHECK_ONLY`. It uses the equilibrium expression from
[University of Alabama AY521 notes](https://pages.astronomy.ua.edu/townsley/ay521-421/TA2018_09_17.pdf)
and [NIST H I ionization energy](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&level_out=on&shells_out=on&spectra=H+I&unc_out=on&units=1).
It does not supply production air conductivity or composition.

## Physical Scope

[Korytchenko et al., VANT 2018, 4(116),144-149](https://vant.kipt.kharkov.ua/ARTICLE/VANT_2018_4/article_2018_4_144.pdf)
provides radial spark-expansion modelling context. Its nitrogen/circuit
validation does not validate this air implementation. No external code was
copied. Electrode/axial losses, radiation, magnetic forces, viscosity,
nonequilibrium chemistry and two-temperature effects are omitted.

G1 is a thermal spark-channel reference solver. The current G0 production
handoff state is not yet calibrated to initialize G1. The synthetic run does
not represent an experimentally realized microgap thermal channel.
`LTE_APPLICABILITY_PENDING_CALIBRATION` remains: when a real bridged channel
becomes sufficiently thermalized requires overlap simulation and experiment.

## Completed Verification

The completed run records match the current source SHA256 hashes. Build passed;
CTest passed 10/10 in 2.89 s; the normal Python suite passed 202 tests with
3 skipped and 5 deselected in 11.79 s. The new thermal suite passed 29 tests.
Skips retain the existing measurement-runner and sandbox MPI limitations.
Frozen C/D/E/F/G0 and historical circuit paths have no diff from `839dace`.

All three reference grids together took 6.079 s with peak RSS 115384 KiB.
The primary run took 2.253 s. Derived JSON/CSV results total 109450 bytes.
The maximum primary budget residual divided by final Joule input is 3.992e-12.
No simulations or tests were repeated during interruption recovery; only source
hashes, stored results, and the unchanged frozen paths were checked.

G1 engineering verdict: `PASS`, with the documented first-order/grid, LTE,
radiation and original-figure-access limitations. G2 has not been started.
