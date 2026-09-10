# Stage H4: trusted native-discharge RF through a receiver response

## Role and pathway separation

H4 validates the trusted native-discharge-field-to-receiver development
pathway. It consumes the total Stage-F Jefimenko observer field and applies a
local receiver transfer. It does not consume the G3 port waveform, the H3
200--500 MHz transfer, or unresolved F-R4 mechanism attribution. H4 does not
establish trusted native discharge radiation at 350 MHz and does not predict
the complete experimental received waveform.

## Frozen Stage-F sources and trust

The Stage4 source is
`rf/production/f4_attribution/F4-P-stage4-left-isolated_jefimenko_waveform.csv`;
its frozen continuous trust interval is
`2.9414085089091916--7.966314711629051 GHz`. The Stage5 source is
`rf/production/f4_attribution/F4-C-stage5-highfield-collision_jefimenko_waveform.csv`;
its interval is `3.0472731051868486--10.233758844919172 GHz`. Each source is
paired with its frozen `RFTrustReport` hash in `h4_result_contract.json`.
Only bins inside the exact interval are formal H4 inputs. Untrusted bins are
stored as unavailable (`NaN`), never as zero and never interpolated through a
trust gap. The 200--500 MHz native path, including 350 MHz, remains
`NOT_RESOLVED`.

The field used is `(Ex_total,Ey_total,Ez_total)` at observer
`F4_far_field_reference`, `(0.2,0,0.005) m`, with nominal source-observer
distance `0.2 m`. Stage-F's retarded Jefimenko evaluation already propagates
the plasma source to this observer. Consequently
`NATIVE_PROPAGATION_ALREADY_INCLUDED=true`: H4 applies no additional `1/R`,
phase delay, Friis factor, or transmitter openEMS propagation.

## Projection and receiver fixture

Both trusted spectra are dominated by `Ez_total`. The fixed receiver axis is
therefore `u_rx=(0,0,1)` and `E_parallel=u_rx dot E_native`. The receiver is
`H4_CANONICAL_SHORT_DIPOLE_FIELD_SENSOR`: a center-fed 2.5 mm straight PEC
CurvePort dipole with a 50-ohm load in free space. It is a numerical
development fixture, not a Vivaldi, commercial antenna, or production
receiver. Its CurvePort wire has no unambiguous finite physical radius, so the
result carries `OPENEMS_CURVEPORT_PHYSICAL_WIRE_RADIUS_UNDEFINED`.

## Calibrated field transfer

The local openEMS model uses the official TFSF plane-wave interface
(`exc_type=10`), propagation along `+x`, and parallel field along `+z`. The
incident amplitude is `1 V/m`. Following the installed official
`RCS_Sphere.py` convention, the excitation spectrum is read from `et`, and

`H_rx,E(f) = V_rx,50ohm(f) / E_inc,parallel(f)`.

The result is complex and has units metres. It is not an S parameter and is
not H2's Tx-to-Rx `S21`. A separate port-excited run supplies complex `Zin`
and `S11`. A `+y` incident field supplies the one orthogonal-polarization
check.

## GHz mesh, boundary, and sanity checks

The full-wave support is `2.8--10.5 GHz`, enclosing both exact trusted
intervals with a small upper margin. All boundaries use `PML_8`. The baseline
local mesh has an openEMS-reported 165731 cells and minimum spacing `0.25 mm`;
the fine mesh has 281799 cells and minimum spacing about `0.179 mm`. No
meter-scale 10-GHz domain is used.

At 3.5, 6.0, and 9.5 GHz, fine-grid `|H_rx,E|` is approximately 20.8, 36.1,
and 58.5 micrometres. These loaded values are below the short-dipole
open-circuit effective-height scale `L/2=1.25 mm`, consistent with the highly
capacitive short sensor and 50-ohm loading. Orthogonal response is suppressed
by 113--130 dB. Baseline-to-fine magnitude change is 15.3--15.6%, above the
preferred 5% development target. The formal status is therefore
`RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT`; the likely contributor is the
one-cell CurvePort feed whose discrete geometry changes with mesh. H4 retains
this uncertainty and does not tune the dipole. The machine-readable flag is
`RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT=true`, and absolute receiver
amplitudes carry `H4_ABSOLUTE_AMPLITUDE_STATUS=NUMERICAL_REFERENCE_ONLY`.

## Stage4 and Stage5 application

The native fields use the frozen Stage-F Hann-windowed one-sided transform
(`dt*rfft` with absolute-time phase). Within each exact trusted mask,

`V_rx_native(f) = H_rx,E(f) E_parallel(f)`.

Stage4 has five trusted DFT bins and Stage5 has seven. Fine-grid transfer is
used for the primary result. The corresponding inverse transforms contain
only trusted bins and are labelled
`TRUSTED_BAND_LIMITED_NATIVE_RESPONSE`; they are not complete experimental
waveforms. Stage5 retains `FULL_MAXWELL_REFERENCE_PENDING`.

The normalized incident/received spectral-shape correlations are about
0.954 for Stage4 and 0.978 for Stage5. Their spectral centroids shift upward
by about 0.632 and 1.001 GHz, respectively. These shifts quantify receiver
selection, not changes in plasma-source physics.

## Status and H5 interface

The development method status is `H4_DEVELOPMENT_GATE=PASS` and
`H4_NATIVE_RF_PATHWAY=TRUSTED_BAND_DEVELOPMENT_VERIFIED`. Absolute results
remain conditional on the numerical receiver fixture and its mesh
sensitivity. `PRODUCTION_NATIVE_RF_RECEIVER=NOT_RESOLVED`,
`H2_SCIENTIFIC_VALIDATION=VNA_MEASUREMENT_PENDING`, and
`STAGE_H_EXPERIMENTAL_VALIDATION_PENDING=true` remain unchanged.

`h4_result_contract.json` exports source and trust hashes, observer geometry,
receiver axis and transfer hash, separate Stage4/Stage5 output paths, exact
trust masks, the 350-MHz negative boundary, and all inherited limitations for
H5. It does not merge native RF with the thermal/circuit pathway.
