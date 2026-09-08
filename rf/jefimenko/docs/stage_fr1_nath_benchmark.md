# F-R1 Nath Electron-Avalanche Field Benchmark

This benchmark validates the existing Stage F Jefimenko electromagnetic
field kernel against the published analytical electron-avalanche field
formulation of:

Debasish Nath, "Electromagnetic Fields Due to an Electron Avalanche",
IEEE Transactions on Electromagnetic Compatibility, 64(3), 623-630, 2022,
DOI: 10.1109/TEMC.2021.3139115. The public TechRxiv preprint DOI is
10.36227/techrxiv.14339756.v1.

The benchmark is a solver benchmark only. It does not assert that the Stage C
streamer solution is equivalent to the idealized Nath avalanche.

## Equation Provenance

The TechRxiv/ResearchGate full-text rendering was inspected for Eqs. 8-13
and the magnetic-field derivation through Eq. 28. The DOI download route was
also attempted, but the automated request returned a Cloudflare challenge
HTML page rather than a PDF, so no PDF is committed or retained in the
repository.

The implemented notation follows Nath's Eq. 11:

```text
K = 1 - Rhat . v / c
```

where `R = r_observer - r_source`, `Rhat = R / |R|`, and `v` and `a` are
evaluated at the retarded time. This is the same vector orientation used by
the Stage F Jefimenko kernel.

The analytical electric field uses Nath Eq. 12/13:

```text
E = q/(4*pi*eps0*K^3) [
      (Rhat - v/c) * (1 - |v|^2/c^2) / R^2
    + Rhat x ((Rhat - v/c) x (a/c)) / (c R)
    ]
  + dq/dt/(4*pi*eps0*c^2*K^2*R) [
      Rhat * (Rhat . v) - v
    ]
  + E_positive_ion_trail
```

The second line is the Eq. 13 triple-product form
`Rhat x (Rhat x v/c)`.

The analytical magnetic field uses Nath Eq. 28:

```text
B = 1/(4*pi*eps0*c^2) [
      q*c/(K^3 R^2) * Rhat x (Rhat - v/c) * (1 - |v|^2/c^2)
    + q/(K^3 R) * Rhat x ((Rhat - v/c) x (a/c))
    + dq/dt/(K^2 R) * Rhat x (Rhat - v/c)
    ]
```

The first magnetic term has asymptotic `1/R^2` near-field character. The
last two terms have asymptotic `1/R` radiative character. The cross-product
ordering is not commuted in the implementation.

## Benchmark Avalanche

The primary source is a straight, nonrelativistic z-directed avalanche:

```text
z_h(t) = z0 + v t
q_e(t) = -q0 exp(gamma t)
dq_e/dt = gamma q_e(t)
a = 0
```

The acceleration contribution is exercised only in a separate limiting-case
test. The benchmark does not reproduce the paper's Figure 1 background-field
curve because a machine-readable field/coefficient history is not part of the
current repository inputs.

## Ion Trail and Conservation

The numerical source contains a finite negative electron head and a stationary
positive ion trail. The trail charge is derived from the electron-head charge
growth. A fixed positive seed charge at the origin balances the initial
electron head. The resulting total source charge is conserved to floating-point
roundoff in the generated `SourceSeries`.

## Finite Source Regularization

Nath's formulation treats the electron head as a point charge and the positive
ions as a line charge. The Stage F production Jefimenko solver integrates
finite source cells. F-R1 therefore uses a fixed Cartesian line of small cells
and a narrow Gaussian electron head. The observer is always outside all source
cells.

The reported regularization metadata include:

- source width,
- source cell spacing,
- number of source cells,
- observer distance/source-width ratio.

A three-level source-discretization check is generated. The expected behavior
is decreasing error or a clear discretization plateau.

## Retardation

The analytical moving-head retarded time solves:

```text
t_obs = t_ret + |r_observer - r_h(t_ret)| / c
```

The numerical comparison samples the existing Stage F kernel at observer times
matched to this analytical retarded time. Pre-arrival samples are invalid and
are not silently extrapolated.

## Near and Far Observers

The benchmark uses the paper's observer anchors:

```text
near = (0.25 mm, 0.25 mm, 1 mm)
far  = (1 m, 1 m, 0.001 m)
```

The simplified avalanche is not tuned to reproduce the published waveform
amplitudes. The comparison is analytical-versus-numerical for the same
benchmark source.

## Acceptance Metrics

For well-resolved nontrivial components, the engineering target is:

- normalized waveform L2 error <= 5%,
- relative peak error <= 5%,
- arrival-time error <= one output interval.

Near-observer finite-source effects may be accepted up to 10% only with an
improving discretization trend. Very small cancellation-dominated components
are reported as regularization-sensitive rather than used as primary acceptance
metrics.

## Limitations

F-R1 validates the Stage F electromagnetic field kernel against an independent
published avalanche-field formulation. It does not validate:

- complete Stage C streamer physics,
- Morrow-Lowke transport,
- streamer-head trajectory extraction,
- RF stage attribution,
- antenna/receiver response.
