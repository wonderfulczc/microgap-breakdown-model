# B-RF1 Actual-Electrode Electrostatic Descriptor Handoff

B-RF1 is a read-only interface around independently produced Stage-B
electrostatic results. It does not run COMSOL and does not replace Stage-C's
self-consistent Poisson solve.

## Input

Place one standardized export in a directory containing:

- `metadata.json`
- `field.csv` or `field.h5`

The required fields and SI units are defined in
`contracts/comsol_export_contract.json`. Every actual geometry needs stable
identity, source hashes, explicit `GAP_ANALYSIS_ROI`, mesh metadata, terminal
definitions, and upstream `Cgap` provenance.

## Engineering fixture

Run:

```bash
PYTHONPATH=python .venv/bin/python rf/b_rf1/generate_b_rf1_artifacts.py
```

The generated BASE/REFINED fixture is marked
`SYNTHETIC_DEVELOPMENT_FIXTURE`. It verifies import, ROI, descriptor, routing,
mesh-comparison, and plotting behavior only. It must not be counted as an
actual electrode descriptor or used for a scientific conclusion.

## Current gate

No actual Stage-B COMSOL export was found during this node. Therefore:

- engineering interface: `PASS`
- actual descriptor count: `0`
- F-R6B: not allowed
- next node: wait for Stage-B actual geometry data
