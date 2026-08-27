# Step 1 Validation Report

Date: 2026-08-27

Scope: engineering audit and minimal repair for the frozen 2D Stage 1-5 solver
workflow plus future 3D Afivo-streamer compatibility boundary.

## Git

The working tree initially exposed an empty `.git` placeholder, so ordinary
`git status`, branch, and commit detection failed. The user-provided remote
repository was reachable for read-only metadata probing and was restored
non-destructively:

- remote: `https://github.com/wonderfulczc/microgap-breakdown-model.git`
- remote main: `f752aa290b63e4c9eecbc7a447777b3bea14298a`
- recovered branch: `main`
- recovered recent history: `f752aa2 Translate README to Chinese`; `9655264 Archive microgap breakdown model toolchain`

Historical local Git metadata before this recovery was unavailable in the
working tree. This record does not invent older commits.

## Pytest baseline split

Default pytest now represents source/core regression health:

```bash
.venv/bin/python -m pytest
```

Historical production artifact checks are retained under the `evidence` marker:

```bash
.venv/bin/python -m pytest -m evidence
```

If `results/` is absent, evidence tests report:

```text
EVIDENCE_UNAVAILABLE: historical production artifact intentionally removed during validated cleanup; see docs/github_archive_assessment.md and docs/github_archive_deletion_manifest.csv
```

This separates current code health from archived large-result availability.

## Commands and results

| command | exit status | result | runtime | peak RSS |
|---|---:|---|---:|---:|
| `cmake --build build --parallel` | 0 | PASS, `ninja: no work to do` | 0.06 s | 19,724 KB |
| `ctest --test-dir build --output-on-failure` | 0 | PASS, 6/6 tests passed | 2.05 s | 49,584 KB |
| `.venv/bin/python -m pytest` | 0 | PASS, 51 passed, 1 skipped, 5 deselected | 2.07 s | 181,180 KB |
| `.venv/bin/python -m pytest -m evidence` | 0 | EVIDENCE_UNAVAILABLE, 5 skipped, 52 deselected | 1.89 s | 148,692 KB |
| `.venv/bin/python tools/validate_cleanup_integrity.py` | 0 | PASS with explicit EVIDENCE_UNAVAILABLE lines for removed `results/` figure sources | 1.15 s | 112,984 KB |

The default pytest skip is the existing Stage 2 measurement-run skip:
`measurement runner has not executed`.

## Low-cost closure validator

`tools/validate_cleanup_integrity.py` is the low-cost retained closure/integrity
validator identified by the repository documentation. In the current archived
state it validates retained HTML assets, source/documentation references,
Python import smoke, and the retained C++ smoke executable. Historical
`results/` figure sources are reported as explicit `EVIDENCE_UNAVAILABLE`
records because the production artifacts were intentionally removed during
validated cleanup.

## Step 1 status basis

The C++ build and CTest baseline pass. Core Python regression tests pass after
formal evidence/core separation. Historical production evidence is unavailable
because `results/` was intentionally removed during validated cleanup; this is a
documentation and artifact-availability state, not a solver failure.
