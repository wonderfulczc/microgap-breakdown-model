# Stage I final audit and Stage-J handoff

## Dual-axis closure

Stage I keeps tool readiness separate from scientific validation.
`STAGE_I_TOOL_DEVELOPMENT=PASS` means that contracts, parsers, quality gates,
spectral metrics, comparison APIs, uncertainty propagation, the discrepancy
ledger, and the synthetic-to-real replacement path operate end to end. It does
not mean that a physical pathway is experimentally validated.

`STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT` remains mandatory. No
real VNA, oscilloscope, production geometry, calibration, or uncertainty
dataset has been supplied. Synthetic dry runs demonstrate software readiness
only.

## Validation boundaries

The practical system context is 50--500 MHz, the formal H3 overlap is
200--500 MHz, and the target center is approximately 350 MHz. Native Stage-F
radiation near 350 MHz remains `NOT_RESOLVED`; this does not mean zero native
radiation.

Native Stage4 and Stage5 comparisons retain their exact RFTrustReport bands and
sparse trusted bins. Stage5 retains `FULL_MAXWELL_REFERENCE_PENDING`. The
native and post-breakdown paths cannot currently be coherently summed, and no
absolute cross-mechanism contribution ratio is available.

## Inherited debts

The final debt ledger preserves the H2 VNA gate, H3 loading mismatch and absent
full-wave feedback, H4 receiver-transfer mesh sensitivity and numerical-only
absolute amplitude, sparse Stage-F trusted bins, Stage5 full-Maxwell reference,
missing real Stage-I data, and missing production Tx/Rx geometry. These debts
block scientific validation or limit production interpretation, but none blocks
packaging the reproducible development toolchain.

## Reproducibility and open-source boundary

The repository records build/test commands, Python dependencies, environment
documentation, stage contracts, input hashes, result manifests, and explicit
synthetic provenance. Afivo and openEMS remain external backends; their source
trees must not be copied into the project package. Machine-local paths in old
backend configurations require environment configuration or templates during
Stage J. Build trees, virtual environments, temporary openEMS output, and large
raw Afivo/Stage-F fields are excluded or externally archived; nothing is
deleted by this audit.

## Real-data reentry

When physical data arrive, Stage I resumes at real-data ingestion, provenance
and quality gates, calibration/reference-plane checks, comparison, discrepancy
ledger update, and validation-status update. Existing tool-development dry runs
do not need to be repeated. Parser and metric APIs remain unchanged; only raw
data, actual geometry, instrument/calibration/uncertainty metadata, provenance,
and status are replaced.

Stage-J tool packaging is allowed because the software and reproducibility
checks pass. Public scientific validation remains incomplete, and Stage J is
not started by WP-I-E.
