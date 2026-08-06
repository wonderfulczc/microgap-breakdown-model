# Stage 1 Closure Report

## 1. Overall conclusion

**Stage 1 is complete.** The automated closure gate passes on the frozen manifest baseline.

## 2. Completed scope

Stage 1.1 completed source, parameter, equation, unknown, decision, and automated consistency registries. Stage 1.2 completed the Eq. (4) lifecycle model, analytical derivative and Fourier transform, FFT cross-check, ESD/band energy, Figure 3 main plots, and the three analytical Figure 4a curves.

## 3. Counts

- Sources: 10; parameters: 145; equations: 35; unknowns: 19; decisions: 29.
- Tests: 24 passed; figures: 6 logical figures (12 PNG/PDF files); data files: 5 primary CSV files plus 4 validation CSV files and 2 summaries.
- Manifest: 51 artifacts, including `references/metadata/local_manifest.json`.

## 4. Key validation results

- Registry: 10 sources and 145 parameters, 0 errors, 0 warnings.
- Pytest: 24/24 passed, no failures or skips.
- `t_peak = 7.758250939 ns`; `I_peak = 0.3693459485 A·m`.
- Analytical-versus-quadrature maximum relative error: 1.2375722305025981e-11.
- Trusted FFT-band maximum relative error: 2.226720457802179e-4.
- Hz versus rad/s energy-integral relative error: 0.
- Changing `I_CM0` from 0.44 to 0.40 A·m changes amplitude-linear results by 10% and ESD by 21%.
- Closure validator: `Stage 1 closure validation passed.`

## 5. Physical conclusions

The lifecycle model contributes mainly in VHF, while its propagation-model SHF energy is extremely small. GHz-and-higher collision radiation requires later collision data. `T0 = 6 ns` differs from the actual Eq. (4) peak time. The `I_CM0` conflict changes amplitude, not the lifecycle time scale.

## 6. Unresolved items

Primary transfers are Stage 2: 5 (transport reconstruction, Poisson/SP3 criteria, SP3 boundary treatment, ISG reference density); Stage 3: 0; Stage 4: 9 (grid/domain, sampling/differentiation, current components, isolated control); Stage 5: 3 (Figshare/Cases II–IV/Figure 4b); pointwise-only: 2 (FFT preprocessing). The full dispositions remain in `docs/unknowns.md`.

## 7. Frozen baseline

`docs/stage1_artifact_manifest.csv` records every artifact's SHA256 in its `sha256` column. The manifest intentionally cannot hash itself.

## 8. Stage 2 entry

**READY FOR STAGE 2.** Stage 2 may begin only verified numerical building blocks; it does not begin Case I.

## 9. Scope statement

ISG-0, Poisson, SP3, PETSc/MPI, the streamer solver, and collision results have not been implemented or generated. Figure 2, the Figure 3 inset, Figure 4a's green collision FFT curve and purple isolated-streamer curve, and Figure 4b remain unreconstructed.

The repository is not currently a Git repository. If initialized or moved into one, the suggested archival commands are:

```text
git add <stage1 files>
git commit -m "Close Stage 1 model audit and analytical reconstruction"
git tag -a stage1-closed -m "Stage 1 closed"
```
