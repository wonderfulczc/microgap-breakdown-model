# H2: receiver transfer and VNA validation gate

H2 characterizes the passive frequency-selective receiver/propagation path
before applying a physics-derived source. It does not inject the G3 transient,
modify Stage G, or start H3/H4/H5. The numerical framework is complete, but the
formal experimental cross-validation remains blocked by absent receiver
hardware geometry and VNA data.

## Frequency correction and four distinct roles

The original 0.7--1.3 GHz, 150 mm dipole result is retained byte-for-byte as
`H2_LEGACY_INFRASTRUCTURE_REFERENCE`. It demonstrates openEMS S11/S21, NF2FF,
mesh convergence and domain convergence only; it is not the project
transmission-frequency reference.

The system tool-development reference is now
`H2_350MHZ_DEVELOPMENT_REFERENCE`: 200--500 MHz, centered at 350 MHz, using a
0.4111439424 m center-fed dipole pair at 1.0 m separation. This is an ideal
reference model, not production receiver hardware. The specified theoretical
files are `THEORETICAL_SURROGATE_ONLY`; they can never satisfy the VNA gate.
The five supplied files under `fullwave/h2/theory_inputs/` pass the existing
Touchstone parser and exact 200--500 MHz, 601-point, 0.5 MHz-spacing, 50 ohm
grid checks. Their provenance remains simulation-only and cannot be promoted
to `VNA_MEASUREMENT` or `VALIDATED_WITH_VNA`.

Stage-F GHz intervals remain a separate `H4_NATIVE_RF_PHYSICS` contract and
are not practical system-transmission bands. Future actual receiver Touchstone
files remain the only route to `VALIDATED_WITH_VNA`. These four roles must not
be combined.

The 350 MHz openEMS geometry uses two z-directed, parallel, co-polarized
broadside dipoles with total length 0.4111439424 m, arm length 0.2055719712 m,
50 ohm center ports, and 1.0 m center separation in ideal free space. No balun,
connector, substrate, matching network or support is invented. The theoretical
wire radius is 1 mm, while CurvePort uses a thin-wire representation; status is
`OPENEMS_THIN_WIRE_RADIUS_NOT_IDENTICAL_TO_THEORY` and no geometry tuning is
performed to conceal this difference.

The new mesh targets 500 MHz rather than retaining the old 1.3 GHz mesh. It
uses approximately 10 mm cells, `1.91201e-11 s` FDTD timestep, 1,621,851 actual
cells and 1950 timesteps. Runtime is 30.94 s, peak RSS is 375692 KiB, and raw
external data are 3.45 MB.

Without tuning, openEMS gives an S11 minimum at 334 MHz and -15.2559 dB. The
-10 dB interval is 320--348.5 MHz. At 350 MHz,
`S11=0.281677+j0.176587`, `S11=-9.56541 dB`, and
`Zin=81.2794+j32.2726 ohm`. At 1.0 m,
`S21=0.0912885-j0.0332096`, or -20.2519 dB, with unwrapped phase
-6.63209 rad. The theoretical phase contains only `-kR`; it cannot validate
the antenna/feed phase or group delay.

Direct evaluation of the supplied theory files gives the S11 minimum at
346.5 MHz and -18.5049 dB. At 350 MHz,
`Zin=64.9370+j5.82388 ohm` and `S11=-17.1205 dB`. Relative to openEMS, the
minimum-frequency mismatch is 12.5 MHz (3.6075%), the complex impedance
difference is `16.3425+j26.4488 ohm`, and the 350 MHz S11 difference is
7.5551 dB. The unscaled S11 magnitude-spectrum correlation is 0.91570.

The theoretical 350 MHz S21 is `0.0539436-j0.0945377`, or -19.2638 dB.
Its openEMS magnitude differs by 0.98809 dB and the unscaled magnitude-spectrum
correlation is 0.94360. These comparisons complete the software development
pipeline but are not hardware validation. The impedance mismatch is interpreted
with `OPENEMS_THIN_WIRE_RADIUS_NOT_IDENTICAL_TO_THEORY`; the theory assumes a
1 mm radius while CurvePort does not guarantee an identical finite-radius feed.

At 200/350/500 MHz the computed `2D^2/lambda` distances are
0.22554/0.39470/0.56385 m, respectively. The 1.0 m configuration is therefore
far-field throughout the band. The 350 MHz NF2FF result is finite, with
`Dmax=1.54660`, an axis null and a broadside maximum. The Friis theoretical
surrogate must not be reused at 0.3 m; such a case requires direct full-wave
openEMS.

The project frequency contract is now:

```
SYSTEM_TRANSMISSION_TARGET_CENTER = 350 MHz
SYSTEM_TRANSMISSION_MAX = 500 MHz
H3_SYSTEM_FULLWAVE_DEVELOPMENT_BAND = 200--500 MHz
status = PROJECT_TARGET_BAND_FOR_TOOL_DEVELOPMENT
```

This system choice does not reinterpret the frozen G3 FFT. The previous G3
spectrum audit remains valid as source-content metadata, while H3 tool
development focuses on the <=500 MHz transmission pathway.

## Hardware and data inventory

A repository-wide search and a bounded search under
`/home/helianthusczc/projects` found no receiver dimensions, PCB/FPCB files,
cable or connector model, matching network, calibration record, `.s1p`, or
`.s2p` measurement. The single-turn coil, commercial 350--450 MHz antenna, and
Vivaldi/UWB receiver named in the H2 specification therefore remain
`PLANNED_NOT_AVAILABLE`; they are not modeled as production hardware.

H2 uses one explicitly non-production configuration:
`H2_CANONICAL_DIPOLE_RX`, classified `GEOMETRY_AVAILABLE_ONLY`. It is paired
with `H2_REFERENCE_TX`, not the unresolved Stage-B transmitter. Both are
project-owned 150 mm thin-wire dipoles based on the pinned openEMS half-wave
dipole test semantics. This fixture validates the two-port and radiating
pipeline only.

## Bands and provenance

The G3 uniform record has 801 samples over 10 ns with `dt=12.5 ps`,
`df=99.875156 MHz`, and a Nyquist frequency of `39.950062 GHz`. Mean-removed,
Hann-windowed RFFT power is used only for source planning. The first resolved
bin is the peak for both V and I. The 99.9% upper cumulative frequencies are
`499.375780 MHz` for V and `599.250936 MHz` for I. H2 therefore exports the
planning interval `99.875156--599.250936 MHz` as
`H3_CIRCUIT_STRUCTURE_BANDS/G3_RESOLVED_99P9_CONTENT`. Content below one
frequency bin is unresolved by the finite 10 ns record; this interval is not a
receiver-validation claim.

The canonical dipole fixture uses `0.7--1.3 GHz`, fixed by its documented
half-wave geometry and resonance, and is `SIMULATION_ONLY`. There is no
`H2_RECEIVER_VALIDATED_BANDS` entry without VNA data. Frozen native RF bands
are copied exactly, not rounded or merged:

- Stage4: `2.941408508909--7.966314711629 GHz`.
- Stage5: `3.047273105187--10.233758844919 GHz`.

These form `H4_NATIVE_RF_TRUSTED_BANDS`. Existing VHF/UHF/1--3 GHz trust
semantics are not upgraded; receiver characterization in a low band would not
make native Stage-F radiation trusted there.

## Geometry, reference planes and loading

Tx and Rx centers are `(-0.15,0,0) m` and `(0.15,0,0) m`, separated by 0.3 m.
Both axes are global z, giving parallel co-polarization in free space without
supports. Each center-gap `CurvePort` has `Z0=50 ohm`. Receiver conventions are:

```
port_reference = CURVEPORT_CENTER_GAP_CONNECTOR_EQUIVALENT
voltage_reference = POSITIVE_Z_ARM_MINUS_NEGATIVE_Z_ARM
positive_current_direction = PORT_INTO_RECEIVER
```

The port is an ideal connector equivalent and contains no cable. The model has
no 70 um gap and does not reopen H1's `PORT_EQUIVALENT_MICROGAP` decision.

`Vrx_50ohm` is the voltage transfer at the loaded receiver terminal.
`Vrx_open_circuit_equivalent` is exported separately using the Thevenin
relation and the receiver input impedance inferred from fixture exchange
symmetry. Neither is an intrinsic field transfer. Since no calibrated incident
plane-wave field has been established, `H_rx,E=V_rx/E_inc` remains
`RX_INTRINSIC_FIELD_TRANSFER_PENDING`; S21 is explicitly a Tx/propagation/Rx
system quantity.

## Near/far-field classification

H2 reports both common engineering boundaries:

```
r_reactive = 0.62 sqrt(D^3/lambda)
r_far = 2 D^2/lambda
```

With `D=0.15 m`, the 1 GHz values are `0.06578 m` and `0.15010 m`. The 0.3 m
separation is beyond `2D^2/lambda` across 0.7--1.3 GHz, so the canonical
configuration is classified `FAR_FIELD`. This classification applies only to
the fixture and selected band.

## openEMS S11, S21 and NF2FF

The baseline uses PML_8 and a 5 mm Cartesian mesh. Complex S parameters follow
the openEMS convention:

```
S11 = tx.uf_ref / tx.uf_inc
S21 = rx.uf_ref / tx.uf_inc
```

The symmetric fixture gives the same receiver and transmitter self-response.
The baseline minimum S11 is `-14.7203 dB` at `894 MHz`; its contiguous
`S11 <= -10 dB` interval is 858--940 MHz. At 1 GHz,
`Zin=103.217+j66.565 ohm`,
`S21=0.0884701+j0.0212325`,
`S21=-20.8209 dB`, unwrapped phase is `-6.04764 rad`, and group delay is
`2.14680 ns`. The full complex data, impedance, loaded/open-circuit voltage
transfers and group delay are retained without magnitude-only reduction.

The first radiating NF2FF run is finite at 1 GHz. It gives
`Prad=1.94461e-25` in the pulse-normalized openEMS result and `Dmax=1.38717`
(`1.42 dBi`). The axial E magnitude is `1.44e-16`, versus `4.02e-12`
broadside, showing the expected z-dipole axis null and broadside maximum.
Absolute Prad depends on pulse normalization and is not used as the mesh
observable; directivity and port ratios are used instead.

## Mesh and domain checks

The baseline, fine and expanded-domain cases are the only H2 FDTD runs.
Baseline-to-fine spacing changes from 5.0 to 3.75 mm. The S11 feature moves
from 894 to 906 MHz (`1.342%`), S21 magnitude changes `1.895%` at 1 GHz,
S21 changes `0.1631 dB`, phase changes `0.04903 rad`, and Dmax changes
`0.0621%`. These key observables meet the 5% engineering target.

The baseline structure has at least 0.15 m clearance to the PML. Expanding
the inner domain by 4/3 at the same mesh changes S21 by `0.0007965 dB`, phase
by `0.005079 rad`, magnitude by `0.00917%`, and Dmax by `0.504%`. PML_8 is
therefore retained without another PML sweep.

## VNA protocol and validation gate

No experimental S parameter exists, so raw-data validation, S11/S21
comparison, phase-delay correction and the acceptance metrics cannot be
evaluated. No synthetic measurement, smoothing, frequency shift, or
peak-by-peak scaling is introduced. Formal status is:

```
receiver_transfer_status = SIMULATION_ONLY
production_receiver_status = NOT_RESOLVED
formal_node_status = MINIMAL_FIX_REQUIRED
blocker = VNA_MEASUREMENT_REQUIRED
```

`h2_vna_measurement_protocol.json` specifies the physical closure action. Use
the existing 500 Hz--67 GHz VNA and a full two-port SOLT calibration suitable
for the actual connectors, with reference planes at the Tx/Rx feed-cable ends.
Measure receiver S11, reference-Tx S11 and Tx-to-Rx S21 over exactly
200--500 MHz with 601 points (0.5 MHz spacing). Keep the theoretical-development
comparison separation at 1.0 m and both
dipoles parallel along z in the anechoic chamber where available. Record
instrument/kit serials, calibration time, cable and adapter IDs/lengths,
height, environment, source power, IF bandwidth and averaging. Source power
must be selected within the passive linear regime and checked by repeating a
trace at reduced power. Save raw `.s1p`, `.s2p` and metadata JSON; use no
primary smoothing.

When data arrive, validation requires monotonic frequencies, finite complex
values, traceable calibration, `Z0=50 ohm` unless metadata proves otherwise,
and bounded interpolation without extrapolation. Engineering targets are:
major-feature frequency mismatch no more than 5%, normalized S21 magnitude
correlation at least 0.90, and approximately 3 dB absolute S21 agreement in
the controlled principal band. Failure attribution must distinguish cable/
connector, geometry, material, reference-plane, chamber and mesh causes.

## Contract and resources

`h2_receiver_transfer_contract.json` binds receiver identity, geometry hash,
port/load convention, complex-data CSV hash, valid frequency interval,
near/far status, mesh status, measurement status and source provenance.
`h2_frequency_bands.json` exports separate H3, H4 and H2 masks; it does not
create a universal Stage-H band.

The three runs used 1.129/2.352/2.352 million actual FDTD cells and
1872/2106/1701 timesteps. Total runtime was 111.67 s, maximum peak RSS was
479244 KiB, and external raw files total 11.21 MB. Only compact JSON/CSV data
are retained in the repository.

H2 has established a reproducible simulation-only receiver-transfer and VNA
ingestion framework. The theoretical-surrogate comparison closes
`H2_DEVELOPMENT_GATE=PASS` and permits H3 tool development, while
`H2_SCIENTIFIC_VALIDATION=VNA_MEASUREMENT_PENDING`,
`STAGE_H_EXPERIMENTAL_VALIDATION_PENDING=true`, and the production receiver
remains `NOT_RESOLVED`. It has not cross-validated production hardware.
