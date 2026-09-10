# Pre-H4 frequency compatibility audit

This is a bounded compatibility audit, not a formal Stage node. It does not
create H3.5, start H4, or modify frozen Stage C--G/H1/H2/H3 physics or results.

## Classification table

| file | line/function | value | meaning | classification | active_in_H3_350MHz_path | action_required |
|---|---|---:|---|---|---|---|
| `python/streamer_rf/fullwave/foundation.py` | `wavelength_cell(fmax_Hz)` | caller supplied | wavelength mesh limit | `PARAMETERIZED` | yes, through H2 mesh construction | none |
| `python/streamer_rf/fullwave/foundation.py` | `mesh_report(...,fmax_Hz)` | caller supplied | mesh/CFL audit frequency | `PARAMETERIZED` | yes, frozen H2 result | none |
| `python/streamer_rf/fullwave/receiver.py` | `band_mask`, `bounded_complex_interpolate` | caller supplied | generic band and bounded complex interpolation | `PARAMETERIZED` | yes | none |
| `python/streamer_rf/fullwave/receiver.py` | `read_touchstone` | Hz/kHz/MHz/GHz scale map | input-boundary unit conversion to Hz | `PARAMETERIZED` | indirectly through validated H2 inputs | none |
| `python/streamer_rf/fullwave/transient.py` | `H3_BAND_HZ` | 200e6, 500e6 Hz | formal practical H3 band | `ACTIVE_350MHZ_PATH` | yes | none |
| `python/streamer_rf/fullwave/transient.py` | `exact_band_mask` | defaults to `H3_BAND_HZ`, accepts arguments | reusable band clipping | `PARAMETERIZED` | yes | none |
| `fullwave/h3/generate_h3_results.py` | `H2_PATH` | `h2_350mhz_openems_sparameters.csv` | primary H3 passive transfer input | `ACTIVE_350MHZ_PATH` | yes | none |
| `fullwave/h3/generate_h3_results.py` | H2 support gate | 200e6--500e6 Hz, 601 points | prevents wrong transfer input | `ACTIVE_350MHZ_PATH` | yes | none |
| `fullwave/h3/generate_h3_results.py` | `mask`, `interpolate_transfer`, `apply_transfer` | `H3_BAND_HZ` | G3 spectrum to H2 transfer to received spectrum | `ACTIVE_350MHZ_PATH` | yes | none |
| `fullwave/h3/h3_result_contract.json` | `frequency_band_Hz` | 200e6--500e6 Hz | frozen H3 result contract | `ACTIVE_350MHZ_PATH` | yes | none |
| `fullwave/h2/run_h2_350mhz.py` | `f_start/f_stop/f_center` | 200e6/500e6/350e6 Hz | development reference excitation/output grid | `ACTIVE_350MHZ_PATH` | via frozen result | none |
| `fullwave/h2/run_h2_350mhz.py` | `mesh_report(...,f_stop)` | 500e6 Hz | active mesh target, not old 1.3 GHz | `ACTIVE_350MHZ_PATH` | via frozen result | none |
| `fullwave/h2/h2_frequency_correction.json` | system/H3 fields | 350e6 center, 500e6 max, 200e6--500e6 band | practical frequency policy | `ACTIVE_350MHZ_PATH` | yes, corroborating contract | none |
| `fullwave/h2/run_h2_reference.py` | `f_start/f_stop/f_center` | 0.7e9/1.3e9/1.0e9 Hz | old canonical dipole runner | `LEGACY_REFERENCE_ONLY` | no | none |
| `fullwave/h2/generate_h2_results.py` and legacy H2 JSON/CSV | fixture band | 0.7e9--1.3e9 Hz | canonical S11/S21/NF2FF evidence | `LEGACY_REFERENCE_ONLY` | no | none |
| `fullwave/h1/run_h1_reference.py` | canonical excitation/mesh | 1 GHz reference, 2 GHz mesh audit | backend foundation smoke/reference | `LEGACY_REFERENCE_ONLY` | no | none |
| `fullwave/h1/generate_h1_audit.py` | Stage5 trusted high frequency | 10.2337588449 GHz | multiscale H4 cost context | `STAGE_F_NATIVE_RF_ONLY` | no | none |
| Stage-F trust reports and `h2_frequency_correction.json` | Stage4/Stage5 native ranges | 2.9414--7.9663 and 3.0473--10.2338 GHz | frozen native-discharge trust | `STAGE_F_NATIVE_RF_ONLY` | no | none |
| `tests/fullwave/test_h2_receiver.py` | legacy fixture tests | 0.7e9--1.3e9, 1 GHz | backward regression fixture | `LEGACY_REFERENCE_ONLY` | no | none |
| `tests/fullwave/test_h1_foundation.py` | wavelength/mesh tests | 1--2 GHz | generic helper regression | `PARAMETERIZED` | no; proves GHz support | none |
| `tests/fullwave/test_h3_transient.py` | H3 mask/interpolation/result tests | 200e6--500e6, 350e6 | active development path checks | `ACTIVE_350MHZ_PATH` | yes | none |
| H1/H2 documentation | canonical fixture sections | 0.7--1.3 GHz and 1 GHz | historical infrastructure explanation | `LEGACY_REFERENCE_ONLY` | no | none |
| H2/H3 documentation | system-development sections | 200--500 MHz and 350 MHz | active pathway explanation | `ACTIVE_350MHZ_PATH` | documentation only | none |
| H1/H2/H3 documentation | native-RF sections | approximately 3--10 GHz | frozen Stage-F/H4 explanation | `STAGE_F_NATIVE_RF_ONLY` | no | none |

No occurrence was classified `ACTIVE_HARDCODE_DEFECT`.

## Active H3 trace

The frozen G3 file `thermal/g3_port/g3_port_uniform.csv` supplies 801 samples
at 12.5 ps over 0--10 ns. `zero_padded_rfft` in `transient.py` constructs the
source spectrum without a minimum-frequency assumption. `H3_BAND_HZ` sets the
lower and upper limits to 200e6 and 500e6 Hz; `exact_band_mask` performs the
clipping. `generate_h3_results.py` reads only
`fullwave/h2/h2_350mhz_openems_sparameters.csv`, verifies its exact 601-point
200--500 MHz support, and applies `interpolate_transfer` only to in-band FFT
points. `apply_transfer` computes the received complex spectrum, and the same
band-limited transfer feeds the causal linear-convolution output.

The center frequency is a reporting/reference value of 350e6 Hz in the H2
frequency policy and geometry. It does not replace the full 200--500 MHz mask.
No active call reads `h2_reference_sparameters.csv` or falls back to the old
0.7--1.3 GHz fixture.

## File-level frequency verification

The active H2 transfer has SHA-256
`0fe46aefb531e75c2d9b4f6e68df9a244c2bb73f059690654192c03df85d593c`,
601 rows, minimum 200e6 Hz, and maximum 500e6 Hz. The legacy transfer has 301
rows over 700e6--1300e6 Hz and is not referenced by H3.

Actual H3 output axes are:

| output | samples | minimum Hz | maximum Hz |
|---|---:|---:|---:|
| `h3_transfer_function.csv` | 601 | 200000000 | 500000000 |
| `h3_g3_source_spectrum.csv` | 48 | 205992509.3633 | 499375780.2747 |
| `h3_received_spectrum.csv` | 48 | 205992509.3633 | 499375780.2747 |
| `h3_loading_diagnostic.csv` | 48 | 205992509.3633 | 499375780.2747 |

The FFT-derived files begin above 200 MHz because the 10 ns sampled record and
zero-padding grid contain no bin exactly at 200 MHz. Every retained coefficient
is nevertheless inside the inclusive 200--500 MHz mask.

## Native-RF isolation and units

The Stage4 trust report gives the native interval
2.941408508909--7.966314711629 GHz and marks 300--500 MHz `UNTRUSTED`. The
Stage5 report gives 3.047273105187--10.233758844919 GHz and also marks
300--500 MHz `UNTRUSTED`. These reports are copied into the H4 contract only;
H3 does not import them into its mask, transfer support, or interpolation.
Thus 350 MHz native-plasma radiation remains not resolved and is not promoted
by H2/H3 receiver development.

Internal frequencies use Hz, indicated by `_Hz` names and contracts. The
Touchstone parser converts explicit Hz/kHz/MHz/GHz file units once at input.
Wavelength uses `c/f_Hz`; no 350-to-350 Hz or 350 GHz conversion was found.

## Mesh, sampling, and backward compatibility

`wavelength_cell` and `mesh_report` take `fmax_Hz`. At 500 MHz and 20 cells per
wavelength, the limit is 0.0299792458 m. The frozen H2 mesh reports a 500 MHz
target, approximately 0.010 m maximum cells and 1.91201e-11 s CFL estimate; no
1.3 GHz or 10 GHz constraint enters that path.

The same helper gives 0.00149896229 m at 10 GHz and remains usable for H4.
Synthetic masks and bounded complex interpolation pass both 200--500 MHz and
3--10 GHz configurations. GHz support has not been removed.

G3 `dt=12.5 ps` gives approximately 40 GHz Nyquist, comfortably above the H3
band. The 10 ns window's approximately 100 MHz intrinsic Fourier scale is a
separate resolution limitation. Oversampling does not force a GHz target, and
zero padding is not interpreted as added information.

## Verdict

`FREQUENCY_COMPATIBILITY_AUDIT=PASS`. The old H2 GHz fixture does not affect
current 350 MHz H3 results. Stage-F GHz trust metadata remains isolated for H4.
The reusable mesh, unit, masking and interpolation infrastructure accepts both
the <=500 MHz H3 configuration and future GHz H4 configurations. No production
code or frozen result modification was required.
