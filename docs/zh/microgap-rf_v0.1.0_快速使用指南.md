# microgap-rf v0.1.0 快速使用指南

> 适用版本：`v0.1.0`。v0.1.0 未发布到 PyPI。

> **使用前先区分：** 如果要执行本文的 `examples/`、Stage C–I、350 MHz 或 synthetic validation 示例，请使用 Git 仓库工作副本；仅安装 wheel 时，这些仓库级数据和脚本并不全部随包提供。

# 1. 环境检查

```bash
microgap-rf --version
microgap-rf doctor
```

| 输出 | 含义 |
|---|---|
| `AVAILABLE` | 可使用 |
| `OPTIONAL_MISSING` | 只影响 Afivo/openEMS/COMSOL 等可选工作流 |
| `REQUIRED_MISSING` | **完整科研工具链不完整**：核心 Python 依赖或 CMake/C++/MPI/PETSc 中有缺失；`doctor` 返回退出码 2 |

`REQUIRED_MISSING` 不意味着 `report` 或所有 Python 轻量功能都不可用。但执行 PETSc/C++ 仿真前，必须修复对应依赖。

| 使用模式 | 可以直接做什么 | 需要仓库吗 |
|---|---|---|
| **仅安装 wheel** | `--help`、`doctor`、自定义 smoke config、`report`、Python API、包内 `literature impact`/`release audit` | 否 |
| **Git 仓库 + 环境** | examples、PETSc Stage、350 MHz frozen contract、Stage-I synthetic、Stage-F/G/H | 是 |
| **仓库 + 外部 backend** | Afivo/openEMS/COMSOL 实际计算 | 是，并需外部软件 |

源码环境：

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/pip install -e .
cmake -S . -B build -G Ninja
cmake --build build --parallel
```

外部 backend 按需设置：

```bash
export AFIVO_STREAMER_ROOT=/path/to/afivo-streamer
export OPENEMS_ROOT=/path/to/openEMS
export OPENEMS_PYTHON=/path/to/openems-python
```

# 2. 常用命令

```bash
microgap-rf --help
microgap-rf doctor
microgap-rf run CONFIG.yaml
microgap-rf validate CONFIG.yaml
microgap-rf report CASE_ID
microgap-rf literature impact PAPER_ID
microgap-rf release audit
```

v0.1.0 **没有** `microgap-rf simulate ...`。`run` 只执行轻量 smoke 或冻结 H3 contract 复用，不启动大型物理仿真。
`release audit` 只检查包内发行门且无发布副作用，不查询 GitHub 实时发行状态。

# 3. 新建 case

1. 从 `examples/configs/` 复制最接近的 YAML。
2. 改为唯一 `case_id` 和未存在的 `output.directory`。
3. 不改动 `schema_version/contract_version/model_version`，除非正在执行受控模型变更。
4. 先运行 `doctor`，再运行 case。

```yaml
schema_version: "1.0"
contract_version: rp1.1
model_version: frozen-v2.0
case_id: my_smoke_001
workflow:
  type: smoke
backend:
  core: core
output:
  directory: results/my_smoke_001
scientific_claims:
  native_rf_trusted_at_target: false
```

```bash
microgap-rf run /path/to/my_smoke_001.yaml
microgap-rf report my_smoke_001 --results-root /path/to/results
```

case 目录已存在时 CLI 会拒绝覆盖。
在 wheel 模式下，上述相对 `output.directory` 基于 YAML 所在目录，因此该 case 位于 `/path/to/results/my_smoke_001`。在源码仓库模式下，相对输出基于仓库根目录。

# 4. 三种典型工作流

## 4.1 轻量架构 smoke

```bash
microgap-rf run examples/configs/smoke.yaml
microgap-rf report rp1_smoke
```

看：`results/rp1_smoke/data/result.json`、`status.json`、`report/report.md`。

## 4.2 350 MHz 系统路径

```bash
microgap-rf run examples/configs/system_350mhz_dry_run.yaml
microgap-rf report system_350mhz_development
```

该命令复用 `fullwave/h3/h3_result_contract.json`；不重跑 openEMS。
正式 H3 比较带为 200–500 MHz，target 约 350 MHz。`NATIVE_RF_350MHZ=NOT_RESOLVED`。
v0.1.0 runner 固定输出 `target_Hz=350e6`；不要通过改 `target_MHz` 将这个入口当成通用 H3 重算器。

## 4.3 Stage-I synthetic validation dry-run

```bash
microgap-rf validate examples/configs/validation_synthetic.yaml
microgap-rf report validation_synthetic
```

看：`data/validation_result.json`。结果必须保持 `SYNTHETIC_DRY_RUN` 和 `scientific_validation_allowed=false`，不能作为实验验证。

## 大型/Stage 工作流速查

| 任务 | 实际入口 |
|---|---|
| PETSc 2D 静电 | `build/bin/stage_c1_smoke` |
| PETSc 2D streamer | `build/bin/stage_c2_dynamic OUTPUT [options]` |
| PETSc 电流 | `build/bin/stage_c3_current OUTPUT [options]` |
| Afivo 3D | 外部 `streamer3d CONFIG.cfg`；项目配置在 `solver3d/afivo_reference/` |
| native RF | `rf/**/generate_*.py` |
| thermal/RLC/port | `thermal/g1/generate_g1_reference.py`、`thermal/g2/generate_g2_reference.py`、`thermal/g3_port/generate_g3_port.py` |
| openEMS | `fullwave/h1-h4/run_*.py`（backend）+ `generate_*.py`（后处理） |

# 5. 输入数据

| 任务 | 必备输入 |
|---|---|
| 统一 CLI | `case_id`、`workflow.type`、`backend`、`output.directory` |
| 350 MHz 路径 | `[200,500] MHz`、`target_MHz: 350`、H3 contract |
| PETSc C2/C3 | 电压、种子、steps/dt、SP3/诊断开关；几何/网格在 app 中固定 |
| Afivo | `.cfg`、外部 backend、hook/transport/化学/AMR |
| Stage F | 冻结总 `J_RF`、observer/time grid、RFTrustReport |
| thermal/RLC | thermal 初值、径向网格、`L/C/Rs/Vs`、电导/电流输入 |
| openEMS | Tx/Rx 几何、mesh/PML、port、频带、距离/方向 |
| VNA | `RX_S11.s1p`、`TX_S11.s1p`、`TX_RX_S21.s2p`、校准/reference plane |
| scope | `time_s,voltage_V,channel_id` 波形合同+采样/触发/阻抗/重复元数据 |

v0.1.0 统一 YAML 没有通用 `geometry`/`source` 字段；几何和物理源仍是 backend-specific。

# 6. 输出目录

`results/<case_id>/` 是逻辑结构。相对 `output.directory` 在源码仓库模式下基于 repository root，在 wheel+外部 YAML 模式下基于 YAML 所在目录。

```text
results/<case_id>/
|-- config_input.yaml
|-- config_resolved.yaml
|-- run_manifest.json
|-- provenance.json
|-- status.json
|-- logs/
|-- data/
|-- figures/
`-- report/
```

必保留：两份 config、manifest、provenance、status 和原始 data。
`run` 写 `data/result.json`，`validate` 写 `data/validation_result.json`。
冻结 Stage 结果仍在 `rf/`、`thermal/`、`fullwave/`、`validation/`，不会自动迁移。

# 7. report

```bash
microgap-rf report CASE_ID
microgap-rf report CASE_ID --results-root /other/results
```

`report CASE_ID` 默认在当前执行根的 `results/` 下查找。如果运行时使用了其他目录中的 YAML，请用 `--results-root` 指向实际 `results` 目录。

需要：`run_manifest.json` + `status.json`。

生成：`report/report.md`。

自动包含：backend、退出状态、runtime、科学状态、Git 和配置 hash。

不包含：论文式物理解释、机制定论、自动实验验证。

# 8. Stage-I 真实实验导入

```text
RAW 真实数据
 -> provenance/hash
 -> quality gate
 -> calibration/reference-plane check
 -> simulation comparison
 -> discrepancy ledger
 -> validation status update
```

需要同时提供：

- VNA Touchstone 或 scope waveform；
- 实际 Tx/Rx 几何、距离、方向和极化；
- 仪器、线缆/连接器、校准和 reference plane；
- repetition index 和 uncertainty metadata；
- 原始文件 SHA-256。

合同位置：`validation/stage_i/stage_i_*_contract.json`。

重入顺序：`validation/stage_i/final/stage_i_real_data_reentry_contract.json`。

当前统一 `validate` 只直接运行 synthetic bundle；真实数据使用 Stage-I 脚本/API，不需重跑所有 dry-run。
`validation/stage_i/wp_i_d/generate_wp_i_d.py` 也是冻结 synthetic 专用脚本；真实距离/方向数据应按重入合同复用 `streamer_rf.validation` API，不直接替换脚本内的 synthetic 路径。

# 9. 常见状态

| 状态 | 现在应怎么理解 |
|---|---|
| `STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT` | 工具可用，真实实验尚未闭环 |
| `SYSTEM_350MHZ_VALIDATION=NOT_MEASURED` | 350 MHz 系统路径未有真实测量验证 |
| `NATIVE_RF_350MHZ=NOT_RESOLVED` | 不能认为可信，也不能认为零 |
| `NUMERICAL_REFERENCE_ONLY` | 可用于工具链/数值对照，不是生产幅值预测 |
| `FULL_MAXWELL_REFERENCE_PENDING` | Stage5 仍缺全 Maxwell 参考 |
| `SYNTHETIC_DRY_RUN` | 只证明软件流程通过 |
| `FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED` | H3 加载没有回馈到 G2/G3 |

# 10. 常用文件位置

| 内容 | 位置 |
|---|---|
| 完整工具说明 | `docs/zh/microgap-rf_v0.1.0_工具说明书.md` |
| 统一配置 schema | `config/unified_case_schema.yaml` |
| 配置示例 | `examples/configs/` |
| C++ 应用 | `cpp/apps/` |
| Python API | `python/streamer_rf/` |
| Afivo 接口 | `solver3d/afivo_reference/` |
| native RF | `rf/` |
| thermal/RLC/port | `thermal/` 和 `circuit/` |
| full-wave | `fullwave/h1` 至 `fullwave/h5` |
| Stage-I 合同/结果 | `validation/stage_i/` |
| 真实数据重入 | `validation/stage_i/final/stage_i_real_data_reentry_contract.json` |
| 文献 registry | `literature/literature_registry.yaml` |
| 决策记录 | `docs/decisions/` |
| 科学状态 | `docs/zh/科学状态.md` |
| 软件版本/依赖 | `pyproject.toml` |
| 发行说明 | `release/0.1.0_release_notes.md` |
