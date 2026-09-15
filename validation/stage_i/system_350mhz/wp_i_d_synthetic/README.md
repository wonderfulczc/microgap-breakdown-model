# Stage-I WP-I-D Synthetic Supplemental Dataset

Status: SYNTHETIC_DEVELOPMENT_INPUT
Scientific validation: NOT PERMITTED

This supplement is intended to be combined with the prior WP-I-B synthetic bundle.

It adds:
- Distance sweep at 0.3, 0.5, 0.8 m, 0° polarization, 10 repetitions each.
- Orientation sweep at 0.6 m and 30°, 60°, 90°, 10 repetitions each.
- Synthetic uncertainty assumptions for software uncertainty propagation.

Reuse from prior WP-I-B bundle:
- 0.6 m / 0° / 15 repetitions
- 1.0 m / 0° / 15 repetitions
- VNA synthetic data and background noise record

Important:
- Distance amplitude follows a development-only empirical power law with exponent 1.05.
- Orientation follows a cosine-like linear-polarization response with a synthetic -25 dB amplitude floor.
- These assumptions exist only to test pipeline behavior. They are not experimental evidence and must not be used to infer physical laws.
