# GitHub archive assessment

## 1. Archive positioning

The repository can be reduced and uploaded to a personal GitHub archive if it is positioned as a reusable simulation toolchain and workflow repository.

The correct archive scope is:

- C++17/PETSc/MPI solver source;
- Python configuration, run orchestration and postprocessing tools;
- CMake/Ninja build configuration;
- tests and validators as implementation references;
- stage reports and methodology documents;
- project-transfer notes for future microgap-discharge modeling.

The archive should not be described as a complete, result-bearing evidence package after raw simulation outputs are deleted. Historical Stage 1–5 closure validators depend on `results/`; once `results/` is removed, those validators become historical audit tools rather than runnable closure checks until the simulations or reduced fixtures are regenerated.

## 2. Effect of deleting the 116 uncertain files

The 116 uncertain files recorded in `results/cleanup/retained_uncertain.csv` were not referenced by the current source, configuration, tests or build system in the artifact inventory. They consisted of:

- local literature PDFs;
- one Liu–Pasko Figure 1 preparatory image;
- unreferenced or invalidated result files.

Deleting them does not affect compilation, C++ solver source, Python modules, configuration templates or the future use of the toolchain. It does remove local literature copies and historical raw evidence that would be needed for re-auditing the old reconstruction without re-obtaining those files.

## 3. Effect of deleting current simulation results

Deleting `results/` removes the large historical run outputs, run registries, scalar histories, field snapshots, spectra and result manifests. This is acceptable for a GitHub repository focused on environment, tools and workflow.

Expected consequences:

- the solver and workflow code remain usable;
- new simulations can still write new outputs under `results/`;
- historical closure reports remain as documentation;
- old result-dependent validators and result-dependent Python tests will fail until results are regenerated or replaced by small test fixtures;
- static summary figures in `docs/simulation_reconstruction_summary/assets/` are retained as lightweight documentation, but their original raw source files are no longer present after result deletion.

Post-deletion verification:

- `results/` removed;
- `archive/audit_only/` removed;
- 116 uncertain files removed according to `docs/github_archive_deletion_manifest.csv`;
- C++ build passed;
- CTest passed: 6/6;
- Python import smoke passed for the retained package modules.

## 4. Microgap workflow verification

The proposed research workflow is technically coherent and should be the target workflow for the next project:

实际微间隙电极几何
→ 电极电势边界和实际电压波形
→ 电极附近种子电子/表面发射
→ 微间隙流注与击穿
→ 电极回路时变电流
→ 电流矩与原生宽频辐射
→ 外部RLC耦合
→ 接收端信号

Current repository coverage is partial:

| Workflow item | Current status |
| --- | --- |
| 实际微间隙电极几何 | Not implemented; current model uses simplified axisymmetric domains. |
| 电极电势边界和实际电压波形 | Not implemented as the main workflow; current streamer cases use simplified background fields. |
| 电极附近种子电子/表面发射 | Not implemented; current cases use Gaussian seeds. |
| 微间隙流注与击穿 | Partially covered at the streamer-fluid-solver level, not as a full electrode microgap breakdown model. |
| 电极回路时变电流 | Not implemented; current current moment uses plasma electron drift/current-density integration, not electrode-circuit current. |
| 电流矩与原生宽频辐射 | Partially covered by current moment, derivative, FFT and ESD workflow. |
| 外部RLC耦合 | Not implemented. |
| 接收端信号 | Not implemented. |

Therefore, the repository is ready as a foundation for the microgap project, but the next development step is not more Shi-style replication. The next step is to replace the simplified free-space streamer setup with real electrode geometry, voltage boundary conditions, seed/surface-emission models and circuit coupling.
