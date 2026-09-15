# Streamer RF Replica

**Authoritative research-development documentation:**
[简体中文 README](README.zh-CN.md). The English README is a concise public entry;
scientific state definitions and maintenance rules are governed by the Chinese documentation.

Resource-aware scientific software for microgap streamer, native-RF,
thermal-channel/circuit, full-wave receiver, and simulation-to-experiment
validation studies. The repository combines a C++17/PETSc/MPI core with a
Python analysis package and explicit external-backend contracts.

Current candidate version: `0.1.0rc1` (built for acceptance, not published).
See [CHANGELOG](CHANGELOG.md) and the
[RC notes](release/0.1.0rc1_release_notes.md).

## Scientific scope

The software implements the frozen v2.0 A--J architecture. It provides working
development and reference pathways, but it is not an experimentally validated
microgap prediction package. In particular:

- Stage-I tool development is complete.
- Stage-I scientific validation is pending real experiments.
- Native discharge RF near 350 MHz is unresolved, not zero.
- H3 full-wave loading feedback is not coupled.
- H4 absolute receiver amplitude is a numerical reference.
- Stage5 retains a pending full-Maxwell reference.

See [scientific status](docs/scientific_status.md) before interpreting results.

## Architecture

| Stage | Role |
|---|---|
| A | configuration, lifecycle, and numerical foundations |
| B | external COMSOL electrostatic/geometry validation interface |
| C | 2D PETSc streamer, electrode current, and cold/thermal handoff observables |
| D--E | external Afivo 3D cross-validation and geometry references |
| F | Jefimenko native field, spectral trust, and mechanism diagnostics |
| G | cold/thermal handoff, LTE thermal channel, RLC coupling, transient port |
| H | openEMS structure/receiver paths and native-field receiver response |
| I | measurement contracts, dry-run pipelines, and real-data reentry |
| J | packaging, reproducibility, and open-source preparation |

The canonical handoff description is in
[software architecture](docs/software_architecture.md).

## Repository layout

- `cpp/`: C++17/PETSc/MPI solver and C++ tests.
- `python/streamer_rf/`: reusable Python analysis package.
- `config/`: core run configurations.
- `tests/`: Python regression tests.
- `rf/`, `thermal/`, `fullwave/`: compact stage scripts, contracts, and results.
- `validation/stage_i/`: measurement contracts and explicitly synthetic fixtures.
- `solver3d/afivo_reference/`: Afivo adapters/configurations, not Afivo source.
- `docs/`: architecture, methods, status, and reproducibility guidance.
- `packaging/`: Stage-J release and reproducibility manifests.
- `results/`, `build/`, `.venv/`: local/regenerable artifacts excluded from release.

## Quick start

Supported reference environment: Linux or WSL2, Python 3.14, CMake, Ninja,
Open MPI, and PETSc discoverable through `pkg-config`.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/pip install -e .
microgap-rf doctor
microgap-rf --version
microgap-rf run examples/configs/smoke.yaml
microgap-rf validate examples/configs/validation_synthetic.yaml
microgap-rf report rp1_smoke
cmake -S . -B build -G Ninja
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=python .venv/bin/python -m pytest -q
```

If PETSc is not installed in a system search path, set `PETSC_DIR` and, where
applicable, `PETSC_ARCH`. Configure the MPI launcher through the host MPI
installation; do not hard-code it into project source.

## Smoke reproduction

The smoke profile builds the C++ targets, runs CTest, and exercises lightweight
Stage-F/G/H/I checks. It does not run Afivo, openEMS, or multi-GB simulations.

```bash
./scripts/reproduce_smoke.sh
```

Alternative interpreter/build locations are explicit:

```bash
PYTHON_BIN=/path/to/python BUILD_DIR=/path/to/build ./scripts/reproduce_smoke.sh
```

See [reproducibility profiles](docs/reproducibility.md) for `SMOKE`,
`REFERENCE`, and `FULL_RESEARCH` scope.

## External backends

Afivo-streamer and openEMS/CSXCAD remain external source trees. Set
`AFIVO_STREAMER_ROOT`, `OPENEMS_ROOT`, and `OPENEMS_PYTHON` when executing their
workflows. Their pinned commits and interfaces are recorded in
`packaging/external_backends.json`. The core smoke suite does not require them.

COMSOL is a proprietary external Stage-B comparison workflow. No COMSOL binary
or proprietary project file is distributed here.

## Data policy

Small contracts, manifests, reference results, and labelled synthetic fixtures
are retained. Approximately 5.12 GB of raw Afivo data and 1.39 GB of Stage-F
raw fields are excluded from normal source releases and must be regenerated or
obtained from an external archive. See `packaging/large_data_manifest.json`.

Stage-I fixtures are marked `SYNTHETIC_DEVELOPMENT_INPUT` or
`SYNTHETIC_DRY_RUN`. They are not experimental measurements.

## Real-data reentry

Real VNA, oscilloscope, geometry, calibration, and uncertainty inputs use the
existing contracts under `validation/stage_i/`. The ordered reentry workflow is
defined in `validation/stage_i/final/stage_i_real_data_reentry_contract.json`
and [reproducibility documentation](docs/reproducibility.md). Core parser and
comparison APIs do not need redesign.

## Citation and license

The project code is licensed under [Apache-2.0](LICENSE). The software author is
`Zach`; citation metadata are provided in [`CITATION.cff`](CITATION.cff). Release
candidate preparation does not imply experimental validation:
`STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT` remains unchanged.

No release, tag, package publication, or repository-visibility change has been
performed by Stage J or Release Preparation RP-3B.

Release Preparation RP-1 adds only a thin configuration/CLI layer over the
existing modules. See [CLI architecture](docs/en/cli_architecture.md); it does
not revise frozen physical models or scientific validation states.
