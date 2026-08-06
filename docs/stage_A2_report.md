# Stage 1.2 Analytical Lifecycle and Radiation Chain Report

## Execution environment

- Date: 2026-07-21
- Python: 3.13.1
- Platform: Windows 11
- Numerical stack: NumPy 2.3.5, SciPy 1.16.3, Matplotlib 3.10.7, PyYAML 6.0.3, pytest 9.1.1

## Files and parameter provenance

The stage adds the `python/streamer_rf` analytical package, reproduction entry point, four test modules, sourced YAML configuration, analytical datasets, validation tables, six PNG/PDF figure pairs, and these method/report documents. The configuration is cross-checked against `docs/parameter_registry.csv` before execution.

Primary parameters are alpha=1.7e9 s^-1, beta=7.5e7 s^-1, I0=0.44 A m, and T0=6 ns from SRC-0001, Shi 2019 PDF pp. 5 and 7. The 0.40 A m Figure 3 caption value is retained only as a sensitivity variant.

## Test and validation results

- Automated tests: 24/24 passed.
- Stage 1.1 registry validator: 10 sources and 145 parameters, 0 errors and 0 warnings.
- Maximum analytical-transform versus independent quad relative error: 1.2375722305025981e-11.
- Maximum FFT relative error in the automatically trusted band: 2.226720457802179e-4 (0.0223%).
- rad/s versus Hz energy-integral relative difference: 0.0 at binary64 precision.
- FFT grid: 524,288 points; true frequency resolution 2.051237509 MHz. No mean removal or window was used.

## Primary peak check

- T0: 6.000000 ns
- Actual t_peak: 7.758250939 ns
- t_peak-T0: 1.758250939 ns
- I0: 0.44 A m
- Actual I_peak: 0.3693459485 A m
- I_peak/I0: 0.8394226102

Shi 2019 states that T0 controls when the current moment reaches its maximum. Under the published Eq. (4), T0 is not the exact maximum time when alpha differs from beta; the equation was not altered.

## I0 conflict sensitivity

Changing I0 from 0.44 to 0.40 A m leaves all waveform timescales and shapes unchanged. Current moment peak, positive derivative peak, negative derivative magnitude, and spectral amplitude are each 10% larger for 0.44 than for 0.40. ESD and every integrated analytical band energy are 21% larger because energy scales as I0 squared. UNK-0018 remains open and does not block Stage 1.2 Analytical Lifecycle and Radiation Chain.

## Analytical band energies for the primary lifecycle model

- VHF, 30-300 MHz: 6.0720198616e-8 J
- UHF, 0.3-3 GHz: 1.8445438289e-9 J
- SHF, 3-30 GHz: 1.1741410333e-33 J

These standard-band analytical lifecycle energies are not a reconstruction of Shi 2019 Figure 4b, whose collision data and original discrete integration bins are unavailable.

## Figure assessment

Figure 3a reproduces the asymmetric lifecycle current-moment shape and separately marks T0 and the actual peak. Figure 3b reproduces the asymmetric bipolar analytical derivative. The partial Figure 4a reconstruction contains only the primary, slow-decay, and later-transition analytical lifecycle curves and reproduces their expected low-frequency sensitivity and high-frequency convergence/trend.

Not reconstructed: Figure 2; Figure 3 inset collision signal; Figure 4a green collision FFT and purple isolated-streamer curves; Figure 4b. No collision data, graphical digitization, or synthetic substitutes were used.

## Remaining issues and stage recommendation

The original I0 caption conflict remains open. High-frequency FFT points below the finite-window numerical noise floor are explicitly untrusted; the continuous analytical transform remains the baseline. Original collision and isolated-streamer waveforms remain unavailable.

Stage 2 Numerical Methods may begin as a separate stage because Stage 1.2 Analytical Lifecycle and Radiation Chain meets its analytical, FFT, ESD, test, and documentation acceptance criteria. This recommendation does not authorize or start Poisson, SP3, ISG-0, PETSc/MPI, or a streamer fluid solver.
