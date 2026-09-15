# Stage I WP-I-C native-GHz validation pipeline

WP-I-C prepares the future experimental comparison of the frozen Stage-F/H4
native-discharge path. It is independent of the practical `50–500 MHz` system
profile and does not change the 350 MHz target.

The source contracts, RFTrustReport hashes, observer metadata, receiver
transfer hash, and H4 spectra are verified without recomputation. Stage4 uses
five trusted FFT bins within `2.941408508909–7.966314711629 GHz`; Stage5 uses
seven bins within `3.047273105187–10.233758844919 GHz`. Stage5 retains
`FULL_MAXWELL_REFERENCE_PENDING`.

Future comparisons use only the intersection of measurement support, the exact
Stage-F trusted-bin mask, and H4 receiver-transfer support. Trusted bins are
not synthesized across gaps. Frequency-domain trusted bins are the primary
comparison; H4 inverse-transform waveforms remain
`SECONDARY_RECONSTRUCTION_ONLY`.

An 8 GHz, 80 GS/s oscilloscope can cover the complete Stage4 trust band but
only part of Stage5. The 2 Hz–67 GHz spectrum-analyzer capability can cover
both bands. These are capability statements, not measurement settings or
evidence.

The ingestion hook accepts only explicitly supplied
`SYNTHETIC_NATIVE_GHZ_DRY_RUN` data with Hz frequency units, declared
voltage/field units, exact trusted-bin semantics, repetition indices, and a
false scientific-validation flag. No such fixture is created in WP-I-C.

Completion means `WP_I_C_VALIDATION_FRAMEWORK = PASS` and
`NATIVE_GHZ_PIPELINE_READY = true`. Both native profiles remain
`NOT_MEASURED`; production receiver, absolute amplitude, VNA, mesh, and
full-Maxwell debts remain open.
