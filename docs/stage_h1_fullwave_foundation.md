# H1: openEMS full-wave numerical foundation

H1 validates the external openEMS toolchain and freezes the electromagnetic
reference-plane and model-partition strategy. It does not inject the G3
transient, construct a production antenna/receiver, or predict an experimental
received waveform.

## External backend

The only backend source is `https://github.com/thliebig/openEMS-Project.git`.
It remains outside this repository at
`/home/helianthusczc/projects/openEMS-Project`; its source tree is clean. Exact
pins are recorded in `fullwave/h1/openems_backend.json`:

| Component | Commit |
|---|---|
| openEMS-Project | `ea3c012c688c9110973b2fdbf77b7b8aa09a099b` |
| openEMS | `8f480d04e0e17a780df28e0ac12e2086be041f9e` |
| CSXCAD | `0458ee11ad711909cb32390b195079c3cff0606e` |
| fparser | `4b9c845b449b520c4b8c5f23c74cd04820084f81` |
| AppCSXCAD | `9249ab7f084c376822cea30ac9a361c47ed18131` |

The user-local no-GUI installation is
`/home/helianthusczc/opt/openems`, built with the Python interface and without
MPI. Native dependencies were unpacked into the isolated
`/home/helianthusczc/opt/openems-deps`; no system libraries were replaced and
the project `.venv` was not modified. AppCSXCAD was intentionally not built.
The installed Python 3.14.4 modules are in the backend-owned `venv`. Runtime
requires the PATH and library prefixes recorded in the provenance JSON.

## Canonical runtime validation

`fullwave/h1/run_h1_reference.py` defines a project-owned passive port fixture:
a 50-ohm lumped source drives PEC terminal strips and an independent 100-ohm
load. A 1 GHz center/1 GHz bandwidth Gaussian pulse, graded Cartesian mesh and
absorbing boundaries exercise CSXCAD, openEMS, lumped-port postprocessing and a
200 MHz frequency-domain field dump. It has no Stage G input and is a numerical
fixture, not a microgap model. Raw FDTD files remain outside the repository
under `/tmp/h1_openems_validation`; only compact derived results are retained.

At 200 MHz the baseline result is
`Zin = 100.007110 + j0.468961 ohm`, a 0.4690% complex error relative to the
100-ohm load. `|S11| = 0.33337796`, consistent with the ideal magnitude 1/3
for a 50-to-100-ohm mismatch. The spectral identity
`P_inc = P_ref + P_acc` closes to `1.99e-16` relative error. Time integration
gives incident/reflected/accepted energies of `1.14453176e-18`,
`1.28329405e-19`, and `1.01620236e-18 J`, with zero represented residual.
These pulse-normalized values validate signs and bookkeeping; they are not
continuous-wave power levels.

The baseline frequency-domain field dataset has shape `3 x 25 x 1 x 17`, is
finite, and has maximum magnitude `1.35517e-10` in openEMS field units. This
exercises the field-output path. NF2FF is not applied to the predominantly
lumped, intentionally non-radiating fixture; formal radiated-field work is
outside H1.

## Port plane and microgap partition

H1 inherits the frozen G3 convention without sign changes:

```
reference_plane_id = EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP
voltage_reference = GAP_NODE_MINUS_GROUND
positive_current_direction = EXTERNAL_TO_GAP
```

The production pathway is `PORT_EQUIVALENT_MICROGAP`. Stage B remains
responsible for local electrostatic geometry, while Stage C/G remain
responsible for discharge, channel and circuit physics. Stage H starts at the
terminal EM plane and will model the radiating electrode/lead/antenna
structure. It does not resolve the 70-um breakdown-gap field.

The Cgap partition is `UPSTREAM_LUMPED_CGAP`. The G2/G3 upstream model retains
the lumped gap capacitance; the H port-equivalent geometry must not add an
identical openEMS Cgap. This resolves
`H_GAP_CAPACITANCE_PARTITION_REQUIRED` for the selected pathway without double
counting. `FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED` remains: later application
of the G3 waveform is a one-way handoff and does not update G2 for a changed
full-wave load.

## Multiscale and frequency cost audit

The engineering audit uses a representative 0.3 m cubical radiating domain,
a 70 um gap, the exact frozen Stage-F trusted upper frequency
`10.2337588449 GHz`, 20 cells per free-space wavelength, PML_8, and a nominal
256 bytes per cell. These are cost-estimation assumptions, not production
geometry or physical thresholds.

Uniform four-cell gap resolution requires 17.5 um cells, approximately
`5.052e12` cells, `3.37e-14 s` CFL scale and about `1.29 PB` nominal memory.
Even an optimistic graded estimate retains `1.525e7` cells and the same tiny
timestep, requiring roughly `2.97e7` steps for a 1 us low-frequency window.
The port-equivalent estimate uses a 0.5 mm minimum cell, about `1.201e7` cells,
`9.63e-13 s` timestep and `3.07 GB` nominal cell storage before solver/PML
overhead. The latter is still a conservative planning model, not a committed
production mesh.

Combining a representative 20 ns GHz window with a 1 us low-frequency decay
has a time-window cost ratio of 50 while retaining high-frequency spatial
resolution. H1 therefore records `MULTIBAND_FULLWAVE_RUNS_REQUIRED`. H2/H3
must freeze exact bands from receiver and source requirements; H1 does not
invent them. Frozen native-source trust intervals remain
2.9414085089--7.9663147116 GHz for Stage4 and
3.0472731052--10.2337588449 GHz for Stage5.

## Geometry and mesh contracts

`h1_reference_geometry.json` labels the fixture `H1_REFERENCE_GEOMETRY` and
keeps electrode, lead, antenna and receiver production fields unset.
`PRODUCTION_ELECTRODE_GEOMETRY_PENDING_STAGE_B` remains mandatory; no final
Stage-B geometry is fabricated.

The reusable mesh helper enforces SI units, wavelength-based maximum cell
size, explicit geometry edges, growth ratio no larger than 1.4, finite positive
coordinates and rejection of electrically irrelevant microscopic refinement
under the port-equivalent policy. It reports min/max spacing, geometric cell
count and a conservative CFL estimate. Mesh policy is frequency-dependent;
there is no single hard-coded mesh for future bands.

The baseline runtime used PML_8, 171349 actual FDTD cells, a
`9.83851e-13 s` solver timestep and 3810 iterations. Coarse/baseline/fine
inner spacings were 0.75/0.50/0.25 mm. Baseline-to-fine changes are 0.6347%
for complex Zin and 0.03334% for `|S11|`, below the 5% engineering target.
PML_8-to-PML_10 changes are `5.89e-7` for complex Zin and `7.11e-7` for
`|S11|`. No late-time growth, negative accepted spectral diagnostic, NaN or
Inf was observed.

## G3 readability and H3 source handoff

The frozen G3 summary and its 801-sample uniform waveform are readable and
match the inherited reference plane. H1 does not use that waveform as an
excitation. The installed Python API exposes `SetCustomExcite` through an
fparser expression, not a sampled CSV waveform API. H1 therefore freezes
`LINEAR_TRANSFER_RESPONSE_THEN_G3_POSTPROCESSING`: excite the later passive
linear structure broadband, obtain its transfer/impulse/frequency response,
then apply the G3 waveform in H3 postprocessing. No upstream source hack or
simultaneous independent prescription of G3 V and I is allowed.

## Reproduction, resources and limitations

Run the audit generator with:

```
.venv/bin/python fullwave/h1/generate_h1_audit.py
```

It consumes existing external runtime results and never starts FDTD. The four
authorized FDTD configurations consumed 42.44 s total; baseline runtime was
8.45 s and maximum peak RSS was 185152 KiB. Their raw external data total
844923 bytes. Compact repository outputs include backend/geometry/mesh/cost/
partition/summary JSON, convergence and PML CSV tables, and the 81-row baseline
port response.

H1 validates the openEMS full-wave infrastructure and freezes the
electromagnetic reference-plane/model-partition strategy. It does not yet
predict the experimental received waveform. Production geometry, receiver
transfer, exact multiband definitions, G3 convolution, radiated-pattern
validation and full-wave loading feedback remain later bounded work.
