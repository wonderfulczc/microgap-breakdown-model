# Stage I WP-I-B synthetic pipeline dry run

WP-I-B exercises the Stage-I ingestion and comparison architecture with the
bundle under `validation/stage_i/system_350mhz/synthetic_input`. Every input is
`SYNTHETIC_DEVELOPMENT_INPUT`; the bundle is software-development evidence and
cannot satisfy a VNA, measurement, or system-validation gate.

## Processing boundary

The SHA-256 manifest is verified before parsing. The existing H2 Touchstone
reader ingests the five VNA files without smoothing. All 30 oscilloscope events
remain individually addressable, and the background record is used only as a
noise-reference diagnostic. RAW inputs are never rewritten.

The `SYSTEM_350MHZ` rectangular-window, mean-removed derived spectrum uses the
frozen complex single-sided `2/N` convention. Metrics are kept separate for
`50–500 MHz`, `50–100 MHz`, and `200–500 MHz`. The weak component near 56 MHz
is retained in the first two products but is explicitly outside H3 comparison
support.

Formal H3 comparison is restricted to the overlap with the frozen H3 result,
nominally `200–500 MHz`; the available H3 transform samples span approximately
`205.99–499.38 MHz`. No H3 values are extrapolated below 200 MHz. Normalized
spectral shape is the primary dry-run comparison. An absolute-amplitude ratio
is retained only as `DEVELOPMENT_ABSOLUTE_DIAGNOSTIC_ONLY` because the inputs
and calibration are synthetic.

The synthetic Tx S11 is converted to input impedance and admittance, then the
frozen H3 source voltage is used with the existing loading diagnostic. This is
one-way diagnostics only; `FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED` remains.

## Status boundary

Successful execution sets `WP_I_B_PIPELINE_DRY_RUN = PASS` and
`STAGE_I_END_TO_END_DRY_RUN = PASS`. It also retains all of the following:

- `STAGE_I_EXPERIMENTAL_DATA = NOT_PROVIDED`
- `SYSTEM_350MHZ_VALIDATION = NOT_MEASURED`
- `H2_SCIENTIFIC_VALIDATION = VNA_MEASUREMENT_PENDING`
- `STAGE_H_EXPERIMENTAL_VALIDATION_PENDING = true`

Every discrepancy-ledger entry is `SYNTHETIC_DRY_RUN`. Real files can replace
the bundle by changing data, geometry, calibration, provenance, uncertainty,
and status metadata without redesigning the parser or comparison API.
