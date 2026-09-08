# Stage F-R4 Coherent Mechanism-Stage-Frequency Attribution

F-R4 adds a proposed coherent attribution framework for the frozen Stage F
reduced-order mechanism time series. It does not modify Stage C, Stage E, the
Stage F source contract, or the production Jefimenko kernel, and it does not
repair the unresolved Stage5 mechanism decomposition.

## Motivation

The F-R3 mechanisms are time-domain current-moment source terms:

```text
HEAD_CHARGE_EVOLUTION = dq_h/dt * v_h
HEAD_ACCELERATION = q_h * a_h
CURRENT_MOMENT_REDISTRIBUTION = dM_redis/dt
```

The mechanisms are not independent radiation powers. Their fields can interfere
constructively or destructively. Therefore F-R4 uses complex CWT coefficients
and coherent projection onto the reconstructed field, rather than power
fractions such as `|W_m|^2 / sum |W_m|^2`.

## Mechanism-Specific Field Reconstruction

The existing Stage F current-moment far-field helper
`current_moment_radiation_approx` is applied linearly to each mechanism source
term. For mechanism `m`:

```text
E_m(t) = L[D_m(t)]
E_rec(t) = sum_m E_m(t)
```

The operator `L` is the same reduced-order current-moment operator already used
in Stage F/F4 current-moment versus full-Jefimenko cross-checks. This is not an
exact spatial decomposition of the full Jefimenko field.

## Complex CWT

F-R4 reuses the existing Stage F Morlet implementation:

```text
morlet_cwt(time_s, signal, frequency_Hz, n_cycles=5)
```

The CWT returns complex coefficients. F-R4 verifies:

```text
CWT(E_rec) approximately equals sum_m CWT(E_m)
```

before calculating coherent eta values.

## Masks

The common nonnegative mask is binary in the current implementation:

```text
w_s,b(t,f) = stage_mask(t) * band_mask(f) * RF_trust_mask(t,f) * COI_mask(t,f)
```

Stage masks are the frozen F4 stage labels. Frequency masks come from the
stored RFTrustReport trusted interval and existing Stage F band definitions.
No VHF, UHF, or `1-3 GHz` mechanism attribution is produced when those bands are
not trusted or only unresolved in the frozen source context.

## Coherent Contribution

For a stage `s` and trusted band `b`:

```text
P_s,b = sum_{t,f} w_s,b(t,f) |W_rec(t,f)|^2

C_m,s,b =
sum_{t,f} w_s,b(t,f)
Re[ W_m(t,f) conj(W_rec(t,f)) ]

eta_m,s,b = C_m,s,b / P_s,b
```

Because:

```text
W_rec = sum_m W_m
```

then:

```text
sum_m C_m,s,b = P_s,b
sum_m eta_m,s,b = 1
```

up to numerical roundoff, provided `P_s,b` is nonzero and the same mask is used
for all mechanisms.

## Interference Interpretation

`eta_m` is not clamped to `[0, 1]`.

- `eta_m < 0` means destructive coherent contribution.
- `eta_m > 1` is possible when another mechanism contributes negatively.
- near-zero reconstructed power makes eta invalid.

The output records `CONSTRUCTIVE_NET`, `DESTRUCTIVE_NET`, or `NEAR_ZERO_NET` for
each coherent contribution sign.

## Stage4 Result

Stage4 uses the frozen `F4-P-stage4-left-isolated` F-R3 compact mechanism time
series. The trusted interval is:

```text
2.9414085089091916 GHz to 7.966314711629051 GHz
```

The three-level closure diagnostics are:

```text
F-R3 time-domain mechanism closure = 0.0434428
CWT linearity closure = 1.2053e-15
eta sum error <= 1.11e-16 for the primary rows
```

The primary coherent eta values are written to
`fr4_stage4_attribution.csv`. They include destructive contributions; for
example the avalanche trusted-continuous row gives negative
`HEAD_CHARGE_EVOLUTION` and `HEAD_ACCELERATION` eta and eta greater than one for
`CURRENT_MOMENT_REDISTRIBUTION`.

However, derivative robustness propagation shows sign changes in 21 of 27
mechanism-stage-band combinations. Therefore the formal Stage4 mechanism-band
attribution status is:

```text
ATTRIBUTION_NOT_RESOLVED
```

The eta values are retained as diagnostic outputs, not as paper-level mechanism
percentages.

## Stage5 Policy

Stage5 inherits frozen F-R3 statuses:

```text
HEAD_CHARGE_EVOLUTION = NOT_RESOLVED
HEAD_ACCELERATION = NOT_RESOLVED
CURRENT_MOMENT_REDISTRIBUTION = NOT_RESOLVED
FULL_MAXWELL_REFERENCE_PENDING
```

Because F-R3 closure is `1.26792`, F-R4 does not publish formal Stage5 eta
values. `fr4_stage5_status.json` records:

```text
MECHANISM_ATTRIBUTION = NOT_RESOLVED
reason = F_R3_MECHANISM_DECOMPOSITION_NOT_RESOLVED
```

## Relation to Full Jefimenko

F-R4 uses two validation levels:

1. mechanism fields sum to the reduced-order current-moment reconstructed
   field;
2. the total current-moment reconstructed field is compared against the frozen
   full Jefimenko radiative waveform from Stage F/F4.

The frozen full-Jefimenko cross-check remains:

```text
Stage4 correlation = 0.999631, normalized waveform error = 0.0262439
Stage5 correlation = 0.999058, normalized waveform error = 0.0354088
```

Good Level-2 agreement does not repair a failed F-R3 mechanism decomposition.

## Output Interface

F-R4 writes only lightweight files:

```text
rf/jefimenko/validation/fr4_coherent_attribution/fr4_attribution_summary.json
rf/jefimenko/validation/fr4_coherent_attribution/fr4_stage4_attribution.csv
rf/jefimenko/validation/fr4_coherent_attribution/fr4_stage4_robustness.csv
rf/jefimenko/validation/fr4_coherent_attribution/fr4_stage5_status.json
```

No large CWT coefficient fields are dumped by default.

## Limitations

The framework is a candidate method developed in this work. It should not be
described with priority language such as first or first-ever without a separate
novelty review. F-R4 establishes the coherent attribution machinery and its
failure/success gates, but the current Stage4 and Stage5 data do not yet support
robust mechanism-stage-frequency physical percentages.
