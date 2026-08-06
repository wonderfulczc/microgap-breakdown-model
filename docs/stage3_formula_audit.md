# Stage 3 Formula Audit

Stage 3 uses the reproducible Morrow–Lowke model and is not a Shi et al. (2019) quantitative reconstruction. The authoritative transport source is the locally retained Morrow and Lowke (1997) PDF, SHA256 `8095418484042060a50ec1bd7df4037c825833853fa86b9e17219d1141e588fd`. Validator source tag: Morrow–Lowke 1997 PDF p.14, A1–A11.

| ID | Source/page/equation | Program quantity | SI conversion | Implementation | Test/ambiguity |
|---|---|---|---|---|---|
| S3-ML-01 | Morrow–Lowke 1997 PDF p.14, A1–A2 | `ionization_townsend`, `ionization_frequency` | cm² and cm⁻¹ to m⁻¹; frequency = α|W| | `MorrowLowke.cpp` | piecewise junction and Python comparison |
| S3-ML-02 | p.14, A3–A4 | two-body attachment | same Townsend conversion | `MorrowLowke.cpp` | negative fitted tail is clamped to zero and recorded |
| S3-ML-03 | p.14, A5 | three-body attachment | N in cm⁻³, η₃/N² in cm⁵ | `MorrowLowke.cpp` | pressure-density scaling test |
| S3-ML-04 | p.14, A7 | diffusion | original E is V/cm; final D is m²/s | `MorrowLowke.cpp` | independent transcription |
| S3-ML-05 | p.14, A8–A11 | drift speed and mobility | cm/s to m/s; μ=|W|/|E| | `MorrowLowke.cpp` | A10's printed endpoint mismatch is regularized by a cubic Hermite connection to A9/A11 (≤1.3% correction); zero-field uses the A11 slope limit |
| S3-RXN-01 | Liu–Pasko 2004 PDF p.4, text after Figure 1 | βep and electron temperature | SI as printed | `ReactionModel.cpp` | analytic unit test |
| S3-RXN-02 | Liu–Pasko 2004 PDF p.4 | βnp | SI as printed, T=300 K baseline | `ReactionModel.cpp` | analytic unit test |
| S3-EQ-01 | Liu–Pasko 2004 PDF pp.2–3, (1)–(4); Shi 2017 PDF p.3, (1)–(4) | three continuity equations and Poisson | SI; electron flux follows Stage 2 convention | `StreamerSolver.cpp` | charge-source balance and one-step tests |
| S3-EQ-02 | Bourdon 2007 and Liu 2007; Stage 2 audited implementation | SP3 emission and coupled Robin boundary | SI before matrix assembly | reused `sp3.cpp` | Stage 2 closure plus coupling test |
| S3-EQ-03 | Kulikovsky 1995; Stage 2 audited implementation | electron face flux | SI and explicit n_ref | reused `transport.cpp` | Stage 2 closure plus one-step test |
| S3-DT-01 | user-specified Stage 3 numerical policy | drift/diffusion/ionization/dielectric/reaction limits | seconds | `AdaptiveTimeStep.cpp` | each controller branch |

The Liu–Pasko 600 dpi rendering is a Stage 4 preparatory artifact only. No digitized lookup table is used in Stage 3.
