# Stage 1.1 Registry and Model Audit Report

- Sources: 10
- Parameters: 145
- direct: 66
- derived: 67
- inferred_from_figure: 6
- assumption: 1
- unresolved: 5
- Equations: 35
- Unknowns: 18

## Critical unknowns

- Figshare original data unavailable
- original grid spacing unknown
- output sampling interval unknown
- Poisson tolerance unknown
- SP3 tolerance unknown
- electron transport table unavailable
- diffusion current inclusion ambiguous
- displacement current inclusion ambiguous
- isolated-streamer control construction unknown
- SP3 open boundary may differ from original

## Readiness

Registry review can proceed. A2, pointwise Case I reproduction, and original-data FFT reproduction remain blocked.

## Conflicts and cross-paper differences

- Shi 2019 I_CM0 is 0.44 A m in prose/Figure 4 caption but 0.4 A m in Figure 3 caption.
- SP3 Table 3 A_j (cm^-1 Torr^-1) and Helmholtz Table 2 A_j (cm^-2 Torr^-2) are separate.
- Shi 2019 is quasi-electrostatic; Luque 2017 full-Maxwell FDTD/PML is validation only.
- Ihaddadene ~4 ps and Shi ~10 ps collapse times refer to different cases.

## Researcher confirmation

Confirm I_CM0; current composition; transport-table reconstruction; solver tolerances; grid convergence; FFT preprocessing; isolated control; Figure 1 digitization.

## Validator

PASS: `python tools/validate_registry.py` returned exit code 0; 10 sources and 145 parameters; 0 errors, 0 warnings.

Stage 1.2 has not started. No solver, Poisson/SP3/ISG-0/FFT implementation or simulation result was created.

