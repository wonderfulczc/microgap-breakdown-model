# Stage I WP-I-D: distance, orientation, repeatability, and uncertainty dry run

## Scope and evidence boundary

WP-I-D validates the analysis software for multi-condition 350 MHz-system
datasets. It combines the frozen WP-I-B development events with the supplied
distance/orientation supplement. Every input and every discrepancy entry is
synthetic development evidence. `SYSTEM_350MHZ_VALIDATION` remains
`NOT_MEASURED`, and no result in this work package is experimental evidence.

## Input integrity and condition index

The supplemental manifest protects file size and SHA-256 for all 64 listed
files. The 60 raw events retain `SYNTHETIC_DEVELOPMENT_INPUT`; uncertainty
assumptions separately retain `SYNTHETIC_UNCERTAINTY_ASSUMPTION`. The merged
index contains 90 unique events: five 0-degree distance conditions and four
orientations at 0.6 m. The frozen 0.6 m/0-degree events are referenced once,
not copied or regenerated.

All new waveforms pass finite-value, monotonic-time, sample-interval, record-
length, baseline, clipping, and metadata checks. RAW waveforms are not
normalized or rewritten. The frozen processing policy independently retains
50--500 MHz, 50--100 MHz, and 200--500 MHz metrics. Formal H3 comparison
support remains 200--500 MHz.

## Distance and orientation analysis

Condition statistics preserve event count, means, standard deviations,
coefficients of variation, and 95% confidence intervals. Distance amplitudes
are fitted to the generic development model `V=A*d^(-n)` without fixing the
exponent. Frequency stability is reported relative to 0.6 m.

Orientation amplitudes are compared with `abs(cos(theta))` and fitted to
`sqrt(A_parallel^2*cos(theta)^2 + A_floor^2)` using deterministic bounded,
nonnegative least squares. Peak and RMS 0-to-90-degree contrasts are labelled
`SYNTHETIC_POLARIZATION_CONTRAST`; they are not measured antenna
cross-polarization discrimination.

## Repeatability and uncertainty

The repeatability matrix exposes every condition, including conditions with
larger scatter. A fixed-seed, 2,000-sample Monte Carlo propagates repeatability,
distance, and voltage-scale terms into the distance exponent. First-order
components are retained separately for representative voltage, frequency,
centroid, and orientation-contrast quantities. All numerical uncertainty
values are development placeholders and are not copied into instrument
capability metadata.

## Simulation-comparison boundary

The frozen H3 full-wave reference exists only for its 1.0 m co-polarized
configuration. Other distance conditions carry
`NO_DIRECT_H3_FULLWAVE_REFERENCE`; orientation conditions carry
`NO_DIRECT_H3_ORIENTATION_REFERENCE`. No inverse-distance extrapolation or
fabricated full-wave orientation response is used.

## Outputs and replacement route

The outputs under `validation/stage_i/wp_i_d/` contain the condition index,
event metrics, statistics, fits, repeatability matrix, uncertainty budget,
discrepancy ledger, and status. Future real data can replace raw files and
their geometry, calibration, uncertainty, provenance, and validation metadata
without changing the core parser, statistics, fitting, or uncertainty APIs.
