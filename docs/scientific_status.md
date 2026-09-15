# Scientific status boundary

The following status classes must not be conflated:

- `FROZEN_COMPLETE_TOOL_NODE`: implementation and its tests are frozen.
- `DEVELOPMENT_VERIFIED`: a development/reference workflow operates.
- `TRUSTED_PHYSICS`: supported only inside an explicit trust mask.
- `NUMERICAL_REFERENCE_ONLY`: absolute physical calibration is unavailable.
- `NOT_RESOLVED`: evidence is insufficient; this does not mean a zero effect.
- `PENDING_REAL_EXPERIMENT`: physical validation has not occurred.

## Current truth state

- Stage H tool development: `PASS`.
- Stage I tool development: `PASS`.
- Stage I scientific validation: `PENDING_REAL_EXPERIMENT`.
- Public scientific validation complete: `false`.
- Practical system context: 50--500 MHz; H3 comparison: 200--500 MHz;
  target center: approximately 350 MHz.
- Native Stage-F radiation near 350 MHz: `NOT_RESOLVED`.
- H3: `FULL_WAVE_LOADING_MISMATCH_HIGH` and
  `FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED`.
- H4 absolute amplitude: `NUMERICAL_REFERENCE_ONLY`, with receiver-transfer
  mesh sensitivity present.
- Stage4/Stage5 native GHz physics is trusted only at sparse RFTrustReport
  bins; Stage5 retains `FULL_MAXWELL_REFERENCE_PENDING`.
- Production Tx/Rx geometry and VNA/oscilloscope data: `NOT_PROVIDED`.
- Cross-path coherent summation: `NOT_PERMITTED_CURRENT_CONFIGURATION`.
- Cross-mechanism absolute contribution: `NOT_RESOLVED`.

Stage-J packaging success cannot upgrade any scientific status above.
