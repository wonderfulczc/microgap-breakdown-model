# H3: G3 transient through the reference full-wave receiver path

## Role and provenance

H3 validates the development chain from the frozen G3 physics-derived port
transient, through the frozen H2 passive full-wave response, to a reference
received voltage. The input is `thermal/g3_port/g3_port_uniform.csv`: 801
samples on 0--10 ns with 12.5 ps spacing. Its reference plane is
`EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP`, voltage is
`GAP_NODE_MINUS_GROUND`, and positive current is `EXTERNAL_TO_GAP`. H3 imposes
only `V_port(t)` as the one-way excitation; `I_port(t)` is retained solely for
loading diagnostics.

The structure is `H3_350MHZ_REFERENCE_RADIATING_STRUCTURE`, reusing the H2
center-fed ideal Tx/Rx dipoles at 1 m, parallel co-polarized broadside in free
space. It is not the experimental microgap transmitter. Production status is
`PENDING_ACTUAL_GEOMETRY`, and the receiver remains a development fixture.

## Band and voltage transfer

The system-development band is exactly 200--500 MHz. Frozen Stage-F GHz trust
bands are not modified. For the H2 50-ohm two-port wave convention,

```
V_tx,total = a1 + b1 = a1(1 + S11)
V_rx,50ohm = b2 = a1 S21
H_V = V_rx,50ohm / V_tx,total = S21/(1 + S11)
```

Thus H3 uses total Tx port voltage, not incident-wave voltage. The frozen H2
`Zin` gives `Yin=1/Zin`. Complex real/imaginary interpolation is bounded by
200--500 MHz; extrapolation is rejected. No openEMS rerun is needed.

## Transform and causal convolution

The original 12.5 ps time origin and amplitudes are preserved. The complex
single-sided peak-amplitude DFT convention is `2/N` times the unnormalized
RFFT, with the DC and Nyquist exceptions. Sixteen-fold zero padding gives
6.2422 MHz numerical sampling but is labelled
`ZERO_PADDING_FOR_LINEAR_CONVOLUTION_ONLY`; it does not improve the intrinsic
100 MHz information scale of the 10 ns record.

Only coefficients inside 200--500 MHz are transferred. They retain 0.247455
of the non-DC source spectral squared norm. A real band-limited impulse is
formed from the complex H2 transfer, wrapped negative-time content is removed,
and samples before `R/c=3.33564 ns` are zero. The retained 99.9% impulse-energy
window is convolved linearly with the source, so no circular wrap occurs. This
causalization is a deterministic finite-band numerical representation, not a
claim of additional measured bandwidth.

## Received result and delay

The source spectral peak within the band is at 205.993 MHz and its centroid is
341.311 MHz. Receiver selection moves the reference received peak to
330.836 MHz and the centroid to 337.572 MHz. The reference received waveform
has 2951 samples over 36.875 ns, a 0.381380 V absolute peak, 0.123389 V RMS,
an absolute peak at 7.300 ns, and a 0.384785 V principal envelope peak at
7.8125 ns. Its 1% onset is 3.3375 ns, only
1.86 ps after the 1 m free-space delay; antenna delay and band limitation mean
the waveform peak need not occur at `R/c`.

At the nearest transform point, 349.563 MHz,
`|V_G3|=4.62949 V`, `|H_V|=0.0757322`, phase is `-0.460015 rad`, and
`|V_rx|=0.350601 V`. The interpolated source input is
`Zin=80.9318+j31.2531 ohm` and `Yin=0.0107526-j0.00415229 S`.

## Loading diagnostic

H3 computes `I_H,implied=Yin,H V_G3` and compares it with the band-limited G3
current. Engineering bands are fixed as normalized L2 mismatch no more than
0.25 for low, no more than 0.75 for moderate, and above 0.75 for high. They are
diagnostic levels, not universal physical thresholds.

The mismatch is 1.95362, complex-correlation magnitude is 0.575088 with phase
0.266524 rad, and the stable implied/G3 current-magnitude ratio has median
1.54483 and p10/p90 0.721753/3.96535. Status is therefore
`FULL_WAVE_LOADING_MISMATCH_HIGH`. H3 does not feed this load into G2/G3 and
retains `FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED`. Postprocessed energy is not
rescaled or claimed to preserve frozen G3 port energy.

## Source-shape baseline

No traceable Ryu equation and parameter set exists in the repository. H3 does
not invent one. The baseline is `NUMERICAL_SOURCE_SHAPE_BASELINE_ONLY`: a cubic
smooth step with the same 100 V peak and a 0.4625 ns transition duration taken
from the G3 90-to-10% decay interval. Geometry and transfer are unchanged.

The baseline/G3 received peak ratio is 2.50261, RMS ratio is 2.82146, 350 MHz
ratio is 3.74766, centroid difference is 4.85765 MHz, and time-waveform
correlation is 0.162483. These differences demonstrate source-shape sensitivity
of the development chain; they are not experimental evidence.

## Radiation and interfaces

H3 reuses the frozen H2 NF2FF result at 350 MHz: finite fields,
`Dmax=1.54660`, a z-dipole axial null, and a broadside maximum. No FDTD or
angular sweep is rerun.

`h3_result_contract.json` binds G3/H2 hashes, reference structure, receiver,
band, complex transfer/input impedance files, loading status and validation
limits for later H5 use. Results are `REFERENCE_RECEIVED_VOLTAGE`, with
`RECEIVER_STATUS=THEORY_AND_FULLWAVE_DEVELOPMENT_VERIFIED`,
`H2_SCIENTIFIC_VALIDATION=VNA_MEASUREMENT_PENDING`, and
`STAGE_H_EXPERIMENTAL_VALIDATION_PENDING`.

H3 validates the physics-derived port-to-full-wave-to-receiver development
pathway. It does not yet constitute a calibrated prediction of the experimental
received waveform.
