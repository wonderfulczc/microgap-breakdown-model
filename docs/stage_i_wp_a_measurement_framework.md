# Stage I WP-I-A: measurement and validation infrastructure

## Scientific purpose and Stage-H relationship

Stage I establishes traceable simulation/experiment comparison. WP-I-A is an
internal work package that defines measurement, calibration, processing,
uncertainty, and discrepancy contracts; it does not perform either physical
comparison campaign. Stage H remains frozen at
`STAGE_H_TOOL_DEVELOPMENT=PASS` and
`STAGE_H_SCIENTIFIC_VALIDATION=PENDING_STAGE_I`, with all H2--H4 limitations
propagated unchanged.

Two later validation profiles remain independent. `SYSTEM_350MHZ` compares
actual 200--500 MHz hardware and received voltage with the G3/H3
post-breakdown reference pathway. `NATIVE_GHZ` compares measurements only
inside the exact Stage-F masks: `2.9414085089091916--7.966314711629051 GHz`
for Stage4 and `3.0472731051868486--10.233758844919172 GHz` for Stage5.
Stage5 retains `FULL_MAXWELL_REFERENCE_PENDING`; native radiation at 350 MHz
remains `NOT_RESOLVED`, not zero.

## Measurement and instrument contracts

`stage_i_measurement_contract.json` records experiment/time identity,
hardware and actual geometry, distance/orientation/polarization, environment,
reference planes, cables/connectors, instruments/calibration, acquisition
settings, repetition, raw path/hash, and uncertainty. Values absent before a
physical campaign are explicitly `NOT_PROVIDED`.

Separate contracts cover the available VNA capability (500 Hz--67 GHz),
oscilloscope capability (8 GHz, 80 GS/s, 12 bit), and optional spectrum
analyzer capability (2 Hz--67 GHz). Capabilities are not measurements. VNA
inputs remain `RX_S11.s1p`, `TX_S11.s1p`, and `TX_RX_S21.s2p`; WP-I-A reuses
`streamer_rf.fullwave.receiver.read_touchstone` and forbids smoothing raw
Touchstone data. Scope rows are `time_s,voltage_V,channel_id`; unnormalized raw
values and all acquisition metadata are mandatory when data arrive.

## Calibration and reference planes

The explicit planes are `SOURCE_PORT`, `TX_FEED`, `FREE_SPACE_REFERENCE`,
`RX_FEED`, and `INSTRUMENT_INPUT`. Direct absolute comparison requires the
same plane. Cross-plane comparison requires a traceable complex transfer.
Cables/connectors may be `CALIBRATED_OUT`, `MEASURED_TRANSFER_AVAILABLE`,
`MODELED`, or `UNRESOLVED`; no loss is removed automatically. De-embedding
must store the raw parent hash, transfer hash, operation parameters, input and
output planes, and a reversible operation record.

Amplitude eligibility is classified as `ABSOLUTE_AMPLITUDE_VALID`,
`RELATIVE_AMPLITUDE_ONLY`, `NORMALIZED_SHAPE_ONLY`, or `NOT_COMPARABLE`.
Absolute validation requires compatible plane, receiver/load, geometry,
instrument impedance, and cable calibration. Normalization never establishes
absolute agreement.

## Data layers and spectral policy

The only layers are `RAW`, `CALIBRATED`, and `DERIVED`. RAW is immutable and
parentless. CALIBRATED and DERIVED products require a parent SHA-256 and a
complete operation record. Individual raw discharge events are never averaged
before storage.

Processing preserves raw intervals and records any derived subwindow or DC
removal. Pairwise simulation/measurement comparisons use identical processing.
The 350-MHz profile can reproduce H3's complex one-sided `2/N` convention;
the native profile reproduces Stage-F's Hann-windowed `dt*rfft`, absolute-time
phase, and one-sided ESD convention. Zero padding creates no new information.
Any filter records type, passband, order, and phase behavior.

Allowed time alignment modes are `ABSOLUTE_TRIGGER_TIME`,
`PROPAGATION_CORRECTED_TIME`, and `FEATURE_ALIGNED_FOR_SHAPE_ONLY`. Feature
alignment is labelled `SHAPE_COMPARISON_ONLY` and cannot validate absolute
timing. Arbitrary shifts are forbidden.

## Metrics, uncertainty, and repetitions

The common API provides feature-frequency and centroid error, normalized
spectral correlation, complex S-parameter error, amplitude ratio/dB error,
waveform normalized RMSE and correlation, arrival-time error, log-log distance
slope, and polarization contrast. WP-I-A assigns no universal acceptance
thresholds.

Uncertainty components cover repeatability, instrument, geometry, distance,
orientation, calibration, and sampling. They remain `NOT_PROVIDED` until
measured. Repeated events have unique repetition indices and raw hashes;
derived ensembles may report mean, median, standard deviation, confidence
interval, and coefficient of variation while retaining event links.

## Mapping, statuses, and discrepancy ledger

`stage_i_simulation_measurement_mapping.csv` maps measured S11/S21 to H2/H3,
200--500 MHz voltage spectra to H3, trusted GHz spectra to H4, and
orientation/distance series to the appropriate model. H3 and H4 are never
compared directly.

The status vocabulary is `NOT_MEASURED`, `DATA_AVAILABLE`,
`CALIBRATION_PENDING`, `COMPARISON_READY`, `PARTIALLY_VALIDATED`, `VALIDATED`,
`NOT_RESOLVED`, and `MODEL_DISCREPANCY`. Data availability alone cannot jump
to validation.

`StageIModelDiscrepancyLedger` stores both values, difference, uncertainty,
possible source, classification, and action. Its source categories include
source model, loading, Tx/Rx geometry, cable/connector, reference plane, mesh,
instrument, environment, and unresolved. The initial ledger is empty and
assigns no blame.

## Current state and next work package

The theoretical H2 Touchstone fixture supplies a software-only ingestion dry
run. It remains `THEORETICAL_SURROGATE_ONLY`, cannot satisfy the VNA gate, and
is not copied into a Stage-I RAW layer. Parser, hashing, reference-plane logic,
and metric APIs pass without experimental claims.

Current status is `STAGE_I_VALIDATION_FRAMEWORK=PASS` and
`STAGE_I_EXPERIMENTAL_DATA=NOT_PROVIDED`. The next work package must acquire,
inspect, and calibrate one real profile's data before any simulation/experiment
validation decision is made.
