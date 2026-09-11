# Stage-I WP-I-B Synthetic Development Dataset

Status: SYNTHETIC_DEVELOPMENT_INPUT
Use: software/pipeline dry run only
Scientific validation: NOT PERMITTED

Frequency coverage: 50-500 MHz
VNA points: 901
VNA step: 0.5 MHz
System target center: 350 MHz
Distances: 0.6 m and 1.0 m
Scope sampling: 5.0 GS/s
Scope record: 2.0 us
Repetitions: 15 per distance

Source-electrode context:
- Curved-needle source, 70 um nominal microgap.
- Exact outsourced-test tip radius/amplitude/SNR values were not used because the raw outsourced test report is not present in this dataset.
- The generator preserves only qualitative project constraints: stable needle/curved-needle readout, readability to 1.0 m, amplitude decay with distance, and a weak low-frequency component near the historical tens-of-MHz readout region.
- The ~350 MHz dominant component is a deliberate Stage-H/Stage-I RF-structure development choice, not a claim that the earlier needle test measured 350 MHz.

Replace this bundle with real VNA/oscilloscope/geometry inputs before any scientific validation claim.
