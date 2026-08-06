# Stage 1 Release Notes

Stage 1 freezes the literature/model audit and analytical lifecycle-to-radiation reconstruction. It includes traceable registries, Eq. (4) and its derivative, analytical Fourier spectrum, FFT cross-validation, ESD/band energy, figures, tests, closure gate, and SHA256 artifact manifest.

Validated results: 24/24 tests pass; registries report zero errors and warnings; peak time is 7.758250939 ns; peak current moment is 0.3693459485 A·m; analytical quadrature error is 1.2375722305025981e-11; trusted FFT error is 2.226720457802179e-4; Hz/rad·s⁻¹ energy consistency error is zero. The 0.40 A·m variant changes amplitude-linear results by 10% and ESD by 21% relative to 0.44 A·m.

Known limitations include unavailable Figshare data and absent transport, Poisson, SP3, streamer, and collision implementations. Stage 2 may start only the three verified numerical building blocks described in `docs/stage2_entry_plan.md`.

The current directory is not a Git repository. Suggested commands after placing it under Git are:

```text
git add <stage1 files>
git commit -m "Close Stage 1 model audit and analytical reconstruction"
git tag -a stage1-closed -m "Stage 1 closed"
```
