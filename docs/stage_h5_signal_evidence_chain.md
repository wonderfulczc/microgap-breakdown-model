# Stage H5: Signal-0/1/2 separation and Stage-I handoff

## Signal levels

`SIGNAL_0` is source-level physics before final receiver selection. It has two
different physical branches: `SIGNAL_0_NATIVE` is the Stage-F trusted native
discharge E/B field or spectrum, while `SIGNAL_0_PORT` is the G3
physics-derived `V_port/I_port` transient. These quantities are not
interchangeable and retain their units and provenance.

`SIGNAL_1` is the propagation/structure level. For the post-breakdown branch,
H3 maps G3 total port voltage through its passive reference Tx/space/Rx
transfer. For the native branch, the Stage-F retarded Jefimenko observer field
already contains plasma-source-to-observer propagation, so H4 only applies a
local receiver response. Reapplying `1/R`, phase, Friis, or an openEMS Tx path
to the native field is forbidden.

`SIGNAL_2` is loaded receiver-terminal voltage after receiver selection. H3's
output is `REFERENCE_RECEIVED_VOLTAGE`; H4's output is
`TRUSTED_BAND_LIMITED_NATIVE_RESPONSE`. Neither is currently an experimentally
calibrated terminal waveform.

## Independent pathways and frequency support

The canonical pathway matrix is `fullwave/h5/h5_pathway_matrix.csv`.
`NATIVE_STAGE4` and `NATIVE_STAGE5` use total Stage-F fields, Stage-F
propagation, and the H4 canonical 2.5-mm short-dipole fixture. Their exact
RFTrustReport intervals remain respectively
`2.9414085089091916--7.966314711629051 GHz` and
`3.0472731051868486--10.233758844919172 GHz`.

`POST_BREAKDOWN_G3` uses the G3 thermal/RLC port transient and H3's 350-MHz
reference structure/receiver in `200--500 MHz`. The practical target is
350 MHz and the system maximum is 500 MHz. Native Stage-F radiation in
200--500 MHz remains `NOT_RESOLVED`; this does not mean it is absent. No trust
is interpolated through the gap between 500 MHz and the lower native trusted
range.

## Transfer ownership

The post-breakdown mapping retains three explicit objects:
`V_G3(f)`, `H_structure+receiver(f)`, and `V_rx,G3(f)`. The native mapping
retains `E_native(f)`, local `H_rx,E(f)`, and `V_rx,native(f)`. H5 does not
collapse these into one undocumented transfer function.

Current H3 and H4 products have different frequency support, source quantity,
receiver fixture, reference plane, and validation state. Consequently
`CROSS_PATH_COHERENT_SUMMATION=NOT_PERMITTED_CURRENT_CONFIGURATION` and H5
does not construct `V_H3+V_H4`. A future coherent result requires the same
physical Tx/Rx geometry and receiver, compatible reference plane, time origin,
frequency and phase conventions, and validated transfers.

## Trust hierarchy and inherited limitations

The common vocabulary is `TRUSTED_PHYSICS`, `DEVELOPMENT_VERIFIED`,
`NUMERICAL_REFERENCE_ONLY`, `NOT_RESOLVED`, and
`EXPERIMENTAL_VALIDATION_PENDING`. Stage-F native spectra are trusted only
inside their RFTrustReport masks. H3 is development verified. H4 absolute
receiver voltage remains numerical-reference-only. Native 350-MHz response is
not resolved, and experimental receiver output remains pending. Downstream
processing never upgrades upstream trust.

H3 retains `FULL_WAVE_LOADING_MISMATCH_HIGH`,
`FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED`, a 10-ns source window, an intrinsic
Fourier scale near 100 MHz, and zero padding only for linear convolution. Its
fine frequency grid is not independent high-resolution physical evidence.

H4 retains `RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT` and
`H4_ABSOLUTE_AMPLITUDE_STATUS=NUMERICAL_REFERENCE_ONLY`. Stage4 has only five
trusted FFT bins and Stage5 seven. Their frequency-domain results are
`PRIMARY_NATIVE_RF_OUTPUT`; inverse-transform waveforms are
`SECONDARY_RECONSTRUCTION_ONLY` and time-domain peaks are not primary
quantitative comparison metrics. Stage5 also retains
`FULL_MAXWELL_REFERENCE_PENDING`.

## Spectral metrics and physical interpretation

`h5_common_spectral_metrics.csv` applies one algorithm to each path while
preserving source and receiver units. It records valid support, source and
receiver peaks/centroids, band-integrated metrics, normalized shape, and
receiver-induced shifts. Absolute E-field and voltage metrics are not compared
to each other. `h5_figure_spectral_shapes.csv` contains normalized shapes only
for future plotting and never replaces absolute data.

At the current evidence level, the practical approximately 350-MHz wireless
path is supported by the G3/H3 post-breakdown-port/structure development
route. The Stage-F native-plasma contribution at 350 MHz remains
`NOT_RESOLVED`; H5 assigns no native/RLC percentage. In the trusted GHz bands,
H4 shows that a receiver can alter spectral shape and centroid without changing
the Stage-F trust masks, but its canonical sensor does not establish production
receiver performance.

## Stage-I handoff and completion status

Stage I must supply measured or surveyed actual Tx/Rx geometry, VNA S11/S21
Touchstone data, cable/connector/reference-plane metadata, distance and
orientation, oscilloscope waveforms, measurement bandwidth/sample rate,
calibration/environment records, and repeat/uncertainty data. No such input is
fabricated in H5. Future comparisons cover feature frequency, spectral shape,
calibrated amplitude, phase where reference planes permit, distance dependence,
and orientation/polarization dependence.

Stage H tool development establishes two independently traceable pathways:
native discharge field -> receiver and physics-derived post-breakdown port ->
full-wave structure -> receiver. Current data do not support calibrated
coherent summation of the two pathways.

The final two-axis status is `STAGE_H_TOOL_DEVELOPMENT=PASS` and
`STAGE_H_SCIENTIFIC_VALIDATION=PENDING_STAGE_I`. Production geometry, VNA
validation, experimental receiver waveforms, and full-wave loading feedback
remain pending.
