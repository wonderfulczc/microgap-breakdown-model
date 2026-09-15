# Frozen v2.0 software architecture

## Layers

The repository is divided into five explicit layers.

1. **Core software:** C++17/PETSc/MPI streamer components, the
   `streamer_rf` Python package, shared contracts, and tests.
2. **External backends:** Afivo-streamer and openEMS/CSXCAD source trees remain
   external. COMSOL is an external proprietary comparison workflow.
3. **Reference data:** compact frozen outputs, contracts, hashes, and manifests.
4. **Large scientific data:** raw Afivo, Stage-F fields, and temporary full-wave
   products are regenerated or archived externally.
5. **Stage-I synthetic fixtures:** deterministic development inputs with
   provenance guards; they are not measurements.

## Physical and data handoffs

Stage B provides a future electrostatic/geometry comparison contract; Stage C
owns the evolving 2D discharge and electrode observables. Stages D--E compare
against an external 3D Afivo backend. Stage F maps source charge/current data to
Jefimenko observer fields and trust-masked spectra.

Stage G maps cold-streamer observables to an unresolved calibrated handoff,
then provides an LTE thermal reference, bidirectional thermal/RLC reference,
and a transient port contract. Stage H has two separate paths:

- Stage-F native observer field to a local receiver response. Propagation is
  already contained in the Jefimenko field and is not applied twice.
- G3 transient port through a passive full-wave structure to a reference
  receiver.

Stage I maps both paths to measurement contracts, calibration/reference-plane
checks, metrics, uncertainty, and discrepancy ledgers. Current paths are not
eligible for coherent summation. Stage J packages these interfaces without
changing their scientific status.

## Canonical references

Detailed decisions remain in `docs/stage_c*.md`, `docs/stage_d*.md`,
`docs/stage_g*.md`, `docs/stage_h*.md`, and `docs/stage_i*.md`. Machine-readable
handoffs are under `thermal/`, `fullwave/`, and `validation/stage_i/`.
