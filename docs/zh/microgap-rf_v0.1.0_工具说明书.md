# microgap-rf v0.1.0 工具说明书

> 适用版本：`v0.1.0`
>
> 软件名：Microgap Breakdown RF Simulation & Validation Toolkit
>
> 仓库：<https://github.com/wonderfulczc/microgap-breakdown-model>

# 1. 软件定位

`microgap-rf` 是面向微间隙击穿放电与 RF 接收链路的科研工具包。
它将二维流注、三维外部流注、原生 RF、热通道/RLC、全波传播和实验对比通过明确合同连接。
软件提供四类能力：仿真、后处理、实验验证基础设施、科研模型维护。
统一 CLI 用于轻量 case、冻结合同复用、合成验证 dry-run、报告和环境检查。
大型 PETSc、Afivo、openEMS 及 Stage-F 工作流仍使用已有 Stage 程序或外部 backend。
软件已发布不等于物理结论已经真实实验验证。

# 2. 软件总体架构

```text
用户配置 / 实验数据
          |
          v
统一 CLI: microgap-rf
          |
          v
任务调度 / 数据合同 / 科学状态门
          |
          v
PETSc | Afivo | COMSOL | Jefimenko | thermal/RLC | openEMS | Stage-I
          |
          v
频谱 / 传递 / 统计 / 不确定度 / 差异台账
          |
          v
results/<case_id>/  +  Stage 冻结结果  +  报告 / provenance
```

| Backend/子系统 | 类型 | v0.1.0 边界 |
|---|---|---|
| C++17/PETSc/MPI | 项目内置源码 | 二维轴对称静电、流注、电流/电阻诊断；需本地编译 |
| `streamer_rf` Python | 项目内置 | RF、thermal/RLC、full-wave 后处理、Stage-I、CLI、provenance |
| Afivo-streamer | 外部 backend | 三维非轴对称流注；本项目只保留配置、hook、合同和摘要 |
| COMSOL | 外部专有工具接口 | Stage-B 静电场/几何参照；v0.1.0 无通用 COMSOL 自动导入器 |
| Jefimenko | 项目内置 Python 分析 | 由 Stage-F 冻结源计算观测点原生场与可信频谱 |
| thermal/RLC | 项目内置 Python | 热通道参考模型、双向热-电路耦合、port transient |
| openEMS/CSXCAD | 外部 backend | H1/H2/H4 全波运行；项目内保留运行脚本、后处理和合同 |
| Stage-I validation | 项目内置 Python/合同 | VNA、示波器、重复性、距离、方向和不确定度处理 |

# 3. 模块总表

| 模块 | 对应 Stage | 核心功能 | 主要输入 | 主要输出 | 外部依赖 | 当前状态 | 对课题的作用 |
|---|---|---|---|---|---|---|---|
| COMSOL 接口 | B | 静电几何和场分布对照 | 电极几何、电压、场结果元数据 | 外部场分布、`Emax`、增强指标 | COMSOL | 只提供合同/文档 | 电极结构与高场区判定 |
| PETSc 2D streamer | C | 轴对称 Poisson、输运反应、光电离、电流诊断 | 内建针-板几何、电压、种子、时步参数 | `ne/np/nn/rho/phi/E`、head、Joule、电极电流 | PETSc/MPI | 已编译工作流，非统一 CLI | 放电发展、电流和冷-热交接 |
| Afivo 3D | D/E | 非轴对称三维放电和几何对照 | `.cfg`、外部 hook、transport/化学参数 | Silo、线剖面、紧凑 CSV/JSON | Afivo-streamer | 外部 backend | 箔片、错位针对等三维影响 |
| native RF | F | `J_RF`、电流矩、Jefimenko 场、频谱可信门 | Stage4/Stage5 冻结时空源 | 观测点 `E/B`、FFT/CWT/ESD、RFTrustReport | 无新外部依赖 | GHz 掩码内 `TRUSTED_PHYSICS` | 原生放电 RF 的带宽与形状 |
| thermal/RLC/port | G | 热通道、RLC 耦合、端口波形 | 热初值、电路参数、电流/电导信号 | `T`、`Rsp`、`Vport`、`Iport`、能量诊断 | SciPy | 数值参考/开发验证 | 击穿后频率机制与端口激励 |
| openEMS full-wave | H | S11/S21、结构传递、接收器响应 | 几何、mesh、port/native field、接收器轴 | Touchstone、transfer、received spectrum/waveform | openEMS/CSXCAD | 开发链通过；有已知债务 | 发射结构、传播、接收选频 |
| experimental validation | I | 采集合同、质量门、比较、重复性/不确定度 | Touchstone、waveform CSV、几何/校准/不确定度 | metrics、matrix、ledger、status | 实验仪器数据 | 工具通过，真实数据缺失 | 仿真-实验交叉验证 |
| CLI/config | RP-1 | case 校验、轻量运行、provenance | YAML | 统一 case 目录 | 无 | 可直接使用 | 日常入口和状态防误用 |
| report | RP-1 | 从 manifest/status 生成 Markdown | `CASE_ID` | `report/report.md` | 无 | 可直接使用 | 保留 case 摘要和科学状态 |
| literature/change | RP-1 | 登记、证据卡、影响分析、变更请求 | registry/evidence/CR YAML/Markdown | `NO_CHANGE/REVIEW/CHANGE_REQUEST_REQUIRED` | 无 | 框架可用 | 防止新文献直接改变冻结物理 |

# 4. 这个工具具体可以解决什么问题

## 4.1 微间隙电极静电场

- 可回答：局部电场分布、`Emax`、field enhancement、高场区和潜在起始位置。
- COMSOL 职责：在外部建模并计算实际三维静电场。
- microgap-rf 职责：提供 Stage-B 元数据/交接约定，以及 PETSc/Afivo 静电基准。
- v0.1.0 不支持从 COMSOL 工程文件一键导入。

## 4.2 二维流注放电

`stage_c2_dynamic` 计算电子/正负离子、空间电荷、电势、径向/轴向电场、光电离源和流注头传播。
可选诊断包括局部场近似、反应源分解、head tracking 和冷态 Joule handoff。
`stage_c3_current` 输出导电/位移/总电流、表面电荷、瞬时电导和电阻。

## 4.3 非轴对称三维流注

Afivo 工作流用于三角箔片、对齐/错位针对等不能用轴对称表示的工况。
Stage D 用于与 PETSc 2D 共同基准对照，Stage E 用于受控三维几何和对称破缺。
Afivo 求解器不内置；项目中的 `.cfg`、hook、比较脚本和摘要是接口层。

## 4.4 原生宽带 RF

Stage F 从冻结的总 `J_RF=-e Gamma_e`、电荷和电流矩出发，通过 Jefimenko 观测场计算 FFT/CWT/ESD 和可信频率掩码。
可回答：观测点原生 `E/B` 的时频特征、GHz 可信频点及机制诊断。
Stage4 可信带为 `2.941408508909–7.966314711629 GHz`，Stage5 为 `3.047273105187–10.233758844919 GHz`，仅 RFTrustReport 掩码内可用。
F-R4 机制百分比归因仍为 `NOT_RESOLVED`。**350 MHz native RF 没有被验证，状态是 `NOT_RESOLVED`，不代表辐射为零。**

## 4.5 击穿后 thermal/RLC

Stage G 将冷放电交接量、导电通道热响应和串联 RLC 电路分层处理。
G1 输出径向温度、压力、电导率、通道半径和 `Rsp`；G2 输出热-电路耦合时序；G3 输出 `V_port_V`、`I_port_A`、分支电流、动态阻抗和能量诊断。
冷-热真实初始化和生产标定仍未完成。

## 4.6 全波传播和接收

- H2/openEMS：开发 Tx/Rx 几何的 `S11`、`S21`、近/远场与 200–500 MHz 参考传递。
- H3：将 G3 `Vport` 通过冻结的结构+接收传递函数得到参考接收频谱/波形。
- H4：将 Stage-F 已传播到观测点的 GHz native `E` 投影到 2.5 mm 短偶极开发夹具，应用 `H_rx,E`。
- 不能将 H3 和 H4 当前结果相干相加。H3 加载失配较高；H4 绝对幅值仅是数值参考。

## 4.7 实验交叉验证

Stage I 定义 VNA、示波器、几何、校准、reference plane、重复次数和不确定度合同。
可比较：峰值频率、频谱质心、归一化频谱相关、S 参数误差、允许时的绝对幅值、波形 NRMSE/相关/到达时间、距离和方向趋势。
CLI `validate` 在 v0.1.0 **只直接执行合成 WP-I-B bundle dry-run**；真实数据应按 Stage-I 合同与脚本/API 重入。

## 4.8 距离/方向/重复性/不确定度

`validation/stage_i/wp_i_d/generate_wp_i_d.py` 及 `streamer_rf.validation` 可处理：
频率稳定性、幅值距离趋势、频谱质心、均值/SD/CV/95% CI、方向响应、极化对比和不确定度传播。
现有结果是 `SYNTHETIC_DRY_RUN`，只证明分析代码可运行。

## 4.9 科研模型维护

| 对象 | 位置 | 用法 |
|---|---|---|
| literature registry | `literature/literature_registry.yaml` | 登记论文元数据、关联 Stage、review/impact 状态 |
| evidence card | `literature/evidence_cards/` | 记录方程、参数、验证、适用范围和冲突 |
| decision record | `docs/decisions/` | 记录专家决策、重跑边界和测试要求 |
| change request | `literature/change_requests/` | 在确认需改模型后提出可追溯变更 |
| impact analysis | `microgap-rf literature impact PAPER_ID` | 根据 registry 输出 `NO_CHANGE/REVIEW/CHANGE_REQUEST_REQUIRED`，不自动改代码 |

固定流程：新论文 → 登记 → 证据卡 → 专家复核 → 影响分析 → 决策记录 → 可选 CR → 最小修改与回归。

# 5. 对博士课题的具体用途

课题：基于可调控微间隙击穿放电的多机制无线传感方法研究。

| 课题问题 | 使用模块 | 输入 | 得到的结果 | 可支撑的论文证据 |
|---|---|---|---|---|
| 电极几何影响 | COMSOL/PETSc C1/Afivo E1 | 几何、间隙、电压 | 场分布、`Emax`、高场区 | 结构-场增强关系 |
| 放电稳定性 | PETSc C2/C3、Afivo D/E | 电压、种子、光电离、三维对称破缺 | head、bridge、电流、三维差异 | 放电演化和重复性机制的数值证据 |
| 原生 RF 频谱 | Stage F/H4 | 冻结 `J_RF`、observer、receiver axis | 可信 GHz bins、接收选频后电压 | 原生 RF 的可信频域特征 |
| 不同放电阶段频率机制 | F 分期、G3/H3 | stage label、port transient | GHz native 与 200–500 MHz port/structure 两条独立路径 | “原生”与“击穿后端口”机制分离 |
| RLC 调制 | G1/G2/G3 | `L/C/Rs/Vs`、热通道 `Rsp` | `Vport/Iport`、能量和动态阻抗 | 电路与热通道如何形成端口频谱 |
| 天线/传播 | H2/H3 | Tx/Rx 几何、S 参数、G3 port | S11/S21、transfer、参考接收波形 | 结构选频与传播链路 |
| 接收链路 | H3/H4/H5 | `Vport` 或 native `E`、接收器响应 | receiver-terminal voltage | 两条可追溯信号链，不做未标定相干相加 |
| 距离影响 | Stage-I WP-I-B/D | 多距离独立重复波形 | 幅值比、幂律拟合、频率稳定性 | 实验距离趋势及不确定度 |
| 实验验证 | Stage I | VNA/scope、几何、校准、reference plane、uncertainty | 比较 metrics、ledger、validation status | 仿真-实验一致性和模型差异证据（待真实数据） |

## 我现在应该用哪条路径？

| 当前问题 | 直接进入 | 实际入口 |
|---|---|---|
| 比较不同针尖/电极的高场区 | COMSOL Stage B | 外部 COMSOL 静电工作流，按 Stage-B 合同导出 |
| 看轴对称流注发展和击穿过程 | PETSc Stage C | `build/bin/stage_c2_dynamic`、`build/bin/stage_c3_current` |
| 几何无法用轴对称表示 | Afivo Stage E | 外部 `streamer3d` + `solver3d/afivo_reference/stage_e/` |
| 看放电自身 GHz 辐射 | Stage F | `rf/**/generate_*.py` 和 RFTrustReport |
| 看击穿后 RLC/端口频谱 | Stage G | `thermal/g1/`、`thermal/g2/`、`thermal/g3_port/` 的生成脚本和 G3 port contract |
| 看 350 MHz 天线传播和接收 | Stage H3 | 冻结 H3 contract；新全波计算需 openEMS Stage 脚本 |
| 看 native GHz 经过接收器后的结果 | Stage H4 | H4 received spectrum/contract；新接收器计算需 openEMS |
| 已有 VNA/示波器真实数据 | Stage I | Stage-I 重入合同、parser/API 和 discrepancy ledger |
| 新论文可能影响当前模型 | literature impact | `microgap-rf literature impact PAPER_ID` |

# 6. 安装与环境

## Python

- 要求：Python `>=3.11`。
- v0.1.0 未发布到 PyPI；使用 GitHub Release wheel 或源码安装。

```bash
# 安装正式 wheel
python3 -m venv .venv
.venv/bin/pip install ./microgap_rf-0.1.0-py3-none-any.whl

# 从源码开发安装
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/pip install -e .
```

## 使用模式与发行内容边界

| 使用模式 | 可以直接做什么 | 需要 Git 仓库吗 |
|---|---|---|
| **仅安装 wheel** | `--help`、`--version`、`doctor`、自定义 smoke config、`report`、Python API、包内 `literature impact` 和 `release audit` | 否，但仓库级示例/数据不随 wheel 完整提供 |
| **Git 仓库 + 项目环境** | `examples/`、PETSc Stage、350 MHz 冻结 contract、Stage-I synthetic dry-run、Stage-F/G/H 脚本与冻结结果 | 是 |
| **仓库 + 外部 backend** | Afivo/openEMS/COMSOL 实际计算 | 是，并需单独安装/配置外部软件 |

**使用提示：** 本文的 `examples/`、Stage C–I、350 MHz 和 synthetic validation 示例均以 Git 仓库工作副本为前提。仅安装 wheel 时，这些仓库级脚本、合同和 fixture 并不全部随包提供；请自行创建 smoke YAML，或使用 Git 仓库。

## C++/PETSc/MPI

```bash
cmake -S . -B build -G Ninja
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

PETSc 发现可使用 `PETSC_DIR`/`PETSC_ARCH`；MPI 由 CMake 查找。Afivo 设置 `AFIVO_STREAMER_ROOT`，openEMS 设置 `OPENEMS_ROOT` 和需要时的 `OPENEMS_PYTHON`。COMSOL 在项目外安装和运行。

```bash
microgap-rf doctor
```

| 状态 | 含义 | 处理 |
|---|---|---|
| `AVAILABLE` | 项目可找到且能查询 | 可执行对应工作流 |
| `OPTIONAL_MISSING` | 外部 backend/可选环境变量缺失 | 核心 CLI 仍可用，相关大型工作流不可用 |
| `REQUIRED_MISSING` | **完整科研工具链不完整**：Python 核心依赖或 CMake/C++/MPI/PETSc 中至少一项缺失 | `doctor` 返回退出码 2；执行 PETSc/C++ 仿真前必须修复对应依赖 |

`REQUIRED_MISSING` 不代表所有 CLI 命令都不能用。若仅使用已安装 wheel 的轻量配置、报告或不依赖缺失项的 Python 后处理，相关命令仍可能正常工作。

# 7. 软件最基本使用方式

- **case**：一次不覆盖的运行实例，用 `case_id` 标识。
- **config**：描述 workflow、backend、输入和输出的 YAML。
- **results**：新 CLI case 的统一目录；Stage 历史结果仍保留在各自目录。

相对 `output.directory` 的基准目录取决于运行模式：

| 运行模式 | 相对输出的基准目录 |
|---|---|
| 源码 Git 仓库 | repository root |
| 仅安装 wheel，传入外部 YAML | 该 YAML 所在目录 |

`microgap-rf report CASE_ID` 默认从当前执行目录的 `results/` 查找安装包 case（在源码仓库内为仓库 `results/`）。如果 YAML 不在当前目录，应显式传入 `--results-root`。

```text
选择/编辑 config
       -> microgap-rf doctor
       -> microgap-rf run CONFIG   或   microgap-rf validate CONFIG
       -> 查看 results/<case_id>/status.json 和 data/
       -> microgap-rf report CASE_ID
```

```bash
microgap-rf --help
microgap-rf doctor
microgap-rf run examples/configs/smoke.yaml
microgap-rf validate examples/configs/validation_synthetic.yaml
microgap-rf report rp1_smoke
microgap-rf literature impact DANGOLA_2008_EQUILIBRIUM_AIR
microgap-rf release audit
```

`run`/`validate` 的目标目录已存在时会拒绝覆盖。新运行应使用新 `case_id` 或新输出目录。

Python API 仍可独立导入，CLI 不替代这些模块：

```python
import streamer_rf
import streamer_rf.fullwave
import streamer_rf.thermal
import streamer_rf.validation
```

`release audit` 是只读的包内发行门审计，不连接 GitHub 查询公开发行。因此 v0.1.0 的审计 JSON 仍可显示 `GITHUB_RELEASE_CREATED=false`/`tag_created=false`；这表示命令本身没有发布副作用，不是 GitHub v0.1.0 的实时状态。

# 8. 配置文件怎么写

v0.1.0 的统一 schema 只支持 `smoke`、`system_rf`、`validation`。通用层与 workflow-specific 层通过 YAML 段分开；各 Stage 的几何/源参数仍在它们的专用 CLI/配置中。

## 实例 1：smoke

```yaml
schema_version: "1.0"
contract_version: rp1.1
model_version: frozen-v2.0
case_id: rp1_smoke
workflow:
  type: smoke
backend:
  core: core
output:
  directory: results/rp1_smoke
scientific_claims:
  native_rf_trusted_at_target: false
```

## 实例 2：350 MHz system 冻结结果复用

```yaml
case_id: system_350mhz_development
workflow:
  type: system_rf
frequency:
  analysis_band_MHz: [200, 500]
  target_MHz: 350
backend:
  fullwave: openems
inputs:
  h3_contract: fullwave/h3/h3_result_contract.json
output:
  directory: results/system_350mhz_development
scientific_claims:
  native_rf_trusted_at_target: false
```

## 实例 3：synthetic validation

```yaml
case_id: validation_synthetic
workflow:
  type: validation
backend:
  validation: stage_i
inputs:
  dataset: validation/stage_i/system_350mhz/synthetic_input
output:
  directory: results/validation_synthetic
scientific_claims:
  requested_status: SYNTHETIC_DRY_RUN
  native_rf_trusted_at_target: false
```

| 字段 | v0.1.0 含义 |
|---|---|
| `case_id` | 1–80 个字符，允许字母、数字、`_ . -`；建议每次运行唯一 |
| `workflow.type` | 仅 `smoke/system_rf/validation` |
| `backend` | 逻辑名到 `core/petsc/openems/afivo/stage_i` 的映射 |
| `backend_options.execution_mode` | `FROZEN_RESULT_REUSE` 或 `LIVE`；LIVE openEMS/Afivo 要求环境变量 |
| `frequency` | 当前键为 `analysis_band_MHz` 和 `target_MHz`；`system_rf` 必须为 200–500 MHz |
| `inputs` | workflow-specific 文件/目录，如 H3 contract 或 synthetic dataset |
| `output.directory` | case 输出目录，不得已存在 |
| `scientific_claims` | 科学状态门；不允许把 350 MHz native RF 或 synthetic 数据升级 |

**代码差异提示：** v0.1.0 统一 schema 不定义通用 `geometry` 或 `source` 段。这些输入由 Stage/backend-specific 入口读取，不应自行加到统一 YAML 并期待 CLI 执行物理仿真。
`system_rf` 校验器接受 200–500 MHz 内的 target，但 v0.1.0 的轻量 runner 在 `data/result.json` 中固定记录 `target_Hz=350e6`。因此应使用冻结 350 MHz 示例；修改 `target_MHz` 不会重算 H3 或生成其他 target 的预测。

# 9. 不同仿真任务怎么开始

## 9.1 PETSc 二维流注

- **目的**：静电基准、流注动力学、电流/电阻诊断。
- **准备**：PETSc/MPI/C++17；需改的电压、种子密度/尺度/偏移、时步和光电离开关。几何与网格在 v0.1.0 C++ app 中固定。
- **运行**：

```bash
cmake --build build --parallel
build/bin/stage_c1_smoke
build/bin/stage_c2_dynamic results/stage_c2/my_case \
  --case my_case --voltage 3000 --steps 120 --sp3 0 \
  --head-tracking --joule-handoff
build/bin/stage_c3_current results/stage_c3/my_case \
  --case my_case --mode case-b --voltage 500 --steps 500
```

- **主要输出**：C2 的 `fields_initial/final.csv`、`axis_profile_*.csv`、`diagnostics.csv`、`summary.txt`，可选 `streamer_head_tracking.csv`、`cold_thermal_handoff.csv`；C3 的 `current_diagnostics.csv` 和 `summary.txt`。
- **怎么看**：先检查 `summary.txt` 状态，再用 `diagnostics.csv` 看时间演化，用 field/axis CSV 看空间分布。
- **边界**：v0.1.0 没有统一 CLI 直接启动此大型计算。

## 9.2 Afivo 三维流注

- **目的**：非轴对称几何、AMR 三维流注和 PETSc/Afivo 交叉对照。
- **准备**：外部 Afivo 的冻结 commit，`AFIVO_STREAMER_ROOT`，目标 `.cfg`，需要时的 `afivo_user` hook。
- **入口**：`solver3d/afivo_reference/common_benchmark/` 与 `solver3d/afivo_reference/stage_e/`。
- **运行形式**：

```bash
export AFIVO_STREAMER_ROOT=/path/to/afivo-streamer
OMP_NUM_THREADS=4 "$AFIVO_STREAMER_ROOT/programs/standard_3d/streamer3d" \
  solver3d/afivo_reference/stage_e/configs/e3_misaligned_needle_pair_500V.cfg
.venv/bin/python solver3d/afivo_reference/stage_e/scripts/e3_summarize_misaligned.py
```

- **主要输出**：Afivo Silo/日志/线剖面，项目 hook 的 field/source CSV，`results_summary/` 中的紧凑 JSON/CSV。
- **边界**：当前 `.cfg` 可包含开发机绝对路径；复制后必须改为本机路径。v0.1.0 没有统一 CLI 直接启动 Afivo。

## 9.3 COMSOL 电场结果接入

- **目的**：将实际电极几何的静电场作为 Stage-B 对照。
- **准备**：几何定义、单位、电压/边界、坐标系、mesh 描述、`E/phi` 导出和 `Emax` 元数据。
- **运行**：在 COMSOL 外部工作流中完成，按 Stage-B 文档保存输入/输出合同。
- **边界**：v0.1.0 当前不支持直接解析 COMSOL 工程或用统一 CLI 启动 COMSOL。

## 9.4 Stage-F native RF

- **目的**：从冻结 Stage4/Stage5 总源计算 Jefimenko 观测场、频谱和信任门。
- **准备**：Stage-F 冻结源文件/哈希、observer、时间网格、RFTrustReport。
- **轻量入口**：

```bash
.venv/bin/python rf/source/audit/generate_stage_f1_audits.py
.venv/bin/python rf/jefimenko/validation/generate_stage_f2_validation.py
.venv/bin/python rf/spectral/validation/generate_stage_f3_validation.py
.venv/bin/python rf/production/f4_attribution/generate_stage_f4_attribution.py
.venv/bin/python rf/report/generate_stage_f_report.py
```

- **主要输出**：Jefimenko waveform、stage RF summary、频谱产品、RFTrustReport 与 report。
- **边界**：不得在信任空缺内插值；200–500 MHz native 不可信。大型 Stage-F 源不由统一 CLI 重算。

## 9.5 thermal/RLC

- **目的**：生成热通道参考响应、双向热-电路耦合和 G3 port transient。
- **准备**：G0 handoff 状态、热模型系数、RLC 参数、电导/电流输入。
- **运行**：

```bash
.venv/bin/python thermal/g1/generate_g1_reference.py --out-dir thermal/g1
.venv/bin/python thermal/g2/generate_g2_reference.py
.venv/bin/python thermal/g3_port/generate_g3_port.py
```

- **主要输出**：`g1_reference_timeseries.csv`、`g1_profile_*.csv`、`g2_*_timeseries.csv`、`g2_reference_summary.json`、`g3_port_uniform.csv`、`g3_port_summary.json`。
- **边界**：这些是冻结参考/开发模型；运行会重写同目录参考文件，新研究应先分支输出或备份，不要覆盖 v0.1.0 冻结结果。

## 9.6 openEMS 全波

- **目的**：求 Tx/Rx 结构、S 参数、空间传递和 GHz 短偶极接收器响应。
- **准备**：`OPENEMS_ROOT`、隔离 Python 环境/版本、几何、mesh/PML、port、频带和输出目录。
- **实际入口**：`fullwave/h1-h4/run_*.py` 启动 openEMS；`generate_*.py` 处理冻结/raw 结果；`fullwave/h3/generate_h3_results.py` 和 `fullwave/h5/generate_h5_results.py` 执行链路集成。
- **示例**：

```bash
"${OPENEMS_PYTHON:-python}" fullwave/h2/run_h2_350mhz.py --output /path/to/raw_case
.venv/bin/python fullwave/h2/generate_h2_350mhz_results.py
.venv/bin/python fullwave/h3/generate_h3_results.py
```

- **边界**：不要把 `run_*` 当成轻量命令。v0.1.0 没有 `microgap-rf simulate fullwave`。

## 9.7 350 MHz 系统分析

- **目的**：检查冻结 G3/H3 击穿后 port/structure 路径。
- **准备**：`fullwave/h3/h3_result_contract.json`、正式 200–500 MHz 频带，约 350 MHz target。
- **运行**：

```bash
microgap-rf run examples/configs/system_350mhz_dry_run.yaml
microgap-rf report system_350mhz_development
```

- **主要输出**：统一 case 目录及 `data/result.json`，其中记录 H3 contract 路径/哈希、频带和状态。
- **重要**：该命令复用冻结 H3 contract，**不重跑 openEMS，也不生成新物理预测**。350 MHz native RF 仍为 `NOT_RESOLVED`。

## 9.8 实验数据 validation

- **目的**：执行 provenance、质量、校准/reference-plane、频谱、重复性和差异台账流程。
- **准备**：VNA Touchstone、scope CSV、measurement/geometry/calibration/uncertainty metadata、SHA-256 manifest。
- **当前 CLI dry-run**：

```bash
microgap-rf validate examples/configs/validation_synthetic.yaml
microgap-rf report validation_synthetic
```

- **主要输出**：`data/validation_result.json`、manifest/provenance/status 和报告。
- **真实数据替代方法**：按 `validation/stage_i/final/stage_i_real_data_reentry_contract.json` 和 Stage-I 生成脚本/API 执行。v0.1.0 的统一 `validate` 尚不是通用真实数据摄取器。

# 10. 做一项新仿真需要准备什么

| 输入 | 必须/可选 | 实际使用者 | 位置/形式 |
|---|---|---|---|
| `case_id`、workflow、backend、output | 统一 CLI 必须 | CLI | YAML，`examples/configs/` |
| 分析带与 target | `system_rf` 必须 | CLI/H3 | `analysis_band_MHz: [200,500]`、`target_MHz: 350` |
| H3 contract | 350 MHz 复用必须 | CLI | `fullwave/h3/h3_result_contract.json` |
| PETSc 电压、种子 `n0/sigma/offset`、steps/dt、SP3 | C2/C3 可调 | C++ app | CLI 参数 |
| PETSc 电极几何、网格、气体基线 | backend-specific，v0.1.0 app 内置 | C++ source | 新几何需明确代码/配置变更，不在统一 YAML |
| Afivo `.cfg`、hook、transport/化学、AMR | Afivo 必须 | 外部 Afivo | `solver3d/afivo_reference/` |
| COMSOL 几何、电压、网格、坐标/单位 | COMSOL 必须 | 外部 COMSOL | 外部模型+导出元数据 |
| Stage-F source hash、time grid、observer、RFTrustReport | native RF 必须 | Stage-F/H4 | 冻结 contract/manifest |
| RLC `L/Cext/Cgap/Rs/Vs`、thermal grid/初值 | G1/G2 必须 | thermal/circuit Python API | backend-specific Python 对象/脚本 |
| Tx/Rx 几何、port、mesh/PML、距离/方向、频带 | openEMS 必须 | H1/H2/H4 | Stage run script 和 CSXCAD 几何 |
| VNA/scope 原始数据 | 实验比较必须 | Stage I | `.s1p/.s2p`、CSV |
| 实际 Tx/Rx 几何、distance/orientation/polarization | 绝对/趋势比较必须 | Stage I | measurement/geometry JSON |
| reference plane、线缆/连接器、校准 | 绝对幅值比较必须 | Stage I | calibration metadata |
| 重复索引、仪器/几何/频率不确定度 | 统计比较必须 | Stage I | uncertainty JSON/contract |

# 11. 外部数据如何导入

| 数据源 | 文件类型 | 关键字段/单位 | 接入方式 |
|---|---|---|---|
| COMSOL | 外部导出表/场文件 | 坐标 m、`phi` V、`E` V/m、几何 ID、网格/边界 | v0.1.0 无通用 parser；按 Stage-B 合同整理后用于人工/脚本对照 |
| Afivo | `.cfg`、Silo、text/CSV | `x/y/z_m`、`phi_V`、`E*_Vpm`、`rho_Cpm3`、`ne_m3`、`J*_Apm2`、level | 由 `solver3d/afivo_reference/` hook 与摘要脚本处理 |
| openEMS | Touchstone/raw run directory | frequency Hz、复数 S 参数、Z0，几何/mesh 元数据 | H1/H2/H4 `generate_*.py`和 `streamer_rf.fullwave` |
| VNA | `RX_S11.s1p`、`TX_S11.s1p`、`TX_RX_S21.s2p` | Touchstone 频率单位、复数 S11/S21、参考阻抗；校准类型、IFBW、功率、平均 | `streamer_rf.fullwave.receiver.read_touchstone`；原始数据不平滑 |
| 示波器 | CSV | 正式合同列 `time_s,voltage_V,channel_id`；另需 sample rate/dt、阻抗、BW limit、trigger、pretrigger、probe、record length、repetition | Stage-I scope contract/API；原始波形不归一化 |
| 几何元数据 | JSON | Tx/Rx ID 和几何、`distance_m`、方向、极化、环境 | `stage_i_measurement_contract.json` 字段集 |
| 校准/reference plane | JSON | `SOURCE_PORT/TX_FEED/FREE_SPACE_REFERENCE/RX_FEED/INSTRUMENT_INPUT`，线缆/连接器状态 | Stage-I reference-plane 和 calibration 合同 |
| 不确定度 | JSON | repeatability、instrument、geometry、distance、orientation、calibration、sampling | `stage_i_uncertainty_contract.json`；未知值保留 `NOT_PROVIDED` |

注：WP-I-B 合成 scope CSV 的底层 dry-run reader 使用两列 `time_s,voltage_V`，`channel_id` 位于事件元数据中；未来正式导入应遵循三列 Stage-I 合同或显式转换，不应假设两种格式等价。

# 12. 实验数据怎么进入 Stage I

```text
真实 VNA/scope/几何/校准/不确定度数据
                        |
                        v
provenance gate: 来源、原始路径、SHA-256、重复索引
                        |
                        v
quality gate: 有限值、时间单调、dt/长度、剪切、基线、频带
                        |
                        v
calibration/reference plane: 校准、线缆/连接器、阻抗、比较平面
                        |
                        v
comparison: H3 200–500 MHz 或 Stage-F/H4 信任掩码交集
                        |
                        v
discrepancy ledger: 差异、不确定度、可能分类、待办
                        |
                        v
validation status: 仅在真实证据和状态门通过后更新
```

- 原始层 `RAW` 不覆盖；`CALIBRATED` 只允许可追溯校准/去嵌；`DERIVED` 保存 FFT、PSD/ESD、频带能量和统计。
- H3 正式比较仅为 200–500 MHz；不向 50–200 MHz 外推。
- native GHz 比较掩码 = 测量支持 ∩ RFTrustReport ∩ receiver-transfer 支持。
- synthetic 必须保留 `SYNTHETIC_DEVELOPMENT_INPUT`/`SYNTHETIC_DRY_RUN` 和 `scientific_validation_allowed=false`，因为它只测试软件，不是仪器证据。
- 重入合同：`validation/stage_i/final/stage_i_real_data_reentry_contract.json`。无需重跑所有 synthetic dry-run。

# 13. 输出目录怎么看

`results/<case_id>/` 是逻辑结构，不总是固定在当前目录。源码仓库模式下，相对路径位于仓库根目录；wheel 模式下，相对路径位于输入 YAML 所在目录。

```text
results/<case_id>/
|-- config_input.yaml
|-- config_resolved.yaml
|-- run_manifest.json
|-- provenance.json
|-- status.json
|-- logs/
|-- data/
|   |-- result.json              # run
|   `-- validation_result.json   # validate
|-- figures/
`-- report/
    `-- report.md
```

| 项目 | 用途 | 长期保存 | 可再生 |
|---|---|---|---|
| `config_input.yaml` | 用户实际提交的配置 | 是 | 否，是原始输入记录 |
| `config_resolved.yaml` | 加入默认值和科学状态后的配置 | 是 | 可，但应保存对应版本 |
| `run_manifest.json` | 命令、版本、hash、backend、runtime、退出状态 | 是 | 不建议重写 |
| `provenance.json` | Git 与输入/输出 hash 摘要 | 是 | 不建议重写 |
| `status.json` | 运行状态与冻结科学状态 | 是 | 可由同版本重建，不应手改 |
| `logs/` | 预留/实际运行日志 | 失败 case 应保存 | 部分 |
| `data/` | workflow 结果 | 是 | 依赖冻结输入 |
| `figures/` | 预留图表位置 | 按需 | 是 |
| `report/report.md` | 自动 case 摘要 | 是 | 是 |

Stage C–I 的冻结结果不会迁移到该目录，仍位于 `rf/`、`thermal/`、`fullwave/`、`validation/` 等 Stage 目录。

# 14. 关键结果怎么看

| 结果量 | 工具中的含义 | 主要位置 | 能回答什么 |
|---|---|---|---|
| 静电场 | `phi_V`、`Er/Ez/Eabs`、`Emax` | C1/C2 field CSV，Afivo field export，COMSOL 外部结果 | 哪里是高场区，几何是否增强电场 |
| streamer | `ne/np/nn/rho`、head position/status、bridge flag | `fields_*.csv`、`streamer_head_tracking.csv`、diagnostics | 流注何时/向哪里发展，是否桥接 |
| current | 高压/地电极的 conduction/displacement/total current | `current_diagnostics.csv` | 电流成分和时间尺度 |
| Joule | `P_cond`、`PJ_*`、channel energy/radius diagnostics | C2 `cold_thermal_handoff.csv`、C3/G1/G2 | 冷放电向热模型传递了什么 |
| native RF | Jefimenko `E/B`、FFT/CWT/ESD、trusted mask | `rf/production/`、RFTrustReport，H4 spectrum | 信任掩码内的原生 GHz 频谱形状 |
| RLC | `Rsp`、`V_port_V`、`I_port_A`、能量残差 | `thermal/g2/`、`thermal/g3_port/` | 击穿后热通道/电路如何形成端口激励 |
| S 参数 | 复数 `S11/S21`、Z0、reference plane | H2 Touchstone、Stage-I VNA | 匹配、传输和频率选择 |
| received waveform | H3 端口路径或 H4 native 路径的负载电压 | `fullwave/h3/`、`fullwave/h4/` | 参考接收电压的时间形状；H4 waveform 仅次级重建 |
| FFT/峰值频率 | 规定窗和归一化下的频谱/最大 bin | Stage F/H3/H4/Stage-I derived | 主要频率特征在哪里 |
| spectral centroid | 频率按频谱强度加权的中心 | Stage-I metrics、H4/H5 | 接收器/条件是否使频谱整体偏移 |
| SNR | 事件峰值对基线 noise RMS 的 dB 比 | WP-I-B/D event metrics | 当前数据是否足以支持频谱特征 |
| repeatability | 事件均值、SD、CV、95% CI | `wp_i_b_repeatability.csv`、`wp_i_d_repeatability_matrix.csv` | 同一工况的散布 |
| uncertainty | 重复性、几何、方向、幅值和频率项的传播 | `wp_i_d_uncertainty_*` | 哪些输入主导派生结果的不确定度 |
| validation status | `NOT_MEASURED/.../VALIDATED/MODEL_DISCREPANCY` | Stage-I status/final contract | 证据已到哪一层；有数据不等于已验证 |

# 15. 如何生成报告

```bash
microgap-rf report CASE_ID
# 非默认 results 根目录：
microgap-rf report CASE_ID --results-root /path/to/results
```

输入 case 目录必须已包含 `run_manifest.json` 和 `status.json`。
对于位于 `/path/to/my_smoke_001.yaml` 且配置为 `output.directory: results/my_smoke_001` 的 wheel case，应使用：

```bash
microgap-rf run /path/to/my_smoke_001.yaml
microgap-rf report my_smoke_001 --results-root /path/to/results
```

输出为 `<results-root>/<case_id>/report/report.md`，包含 backend、退出状态、runtime、科学状态、Git 提交和配置 hash。
该报告不自动解读全部 Stage 物理结果，也不自动产生论文结论；幅值、机制和差异原因仍需科研者解释。

# 16. Provenance 和可复现性

| 字段 | 作用 |
|---|---|
| `input_config_hash` | 确认用户原始 YAML 是否改变 |
| `resolved_config_hash` | 确认实际执行配置 |
| `git_commit` | 定位源码或安装包构建状态 |
| `backend_version` | 定位 PETSc/Afivo/openEMS 等 backend |
| `input_hashes`/`output_hashes` | 检查输入/输出完整性 |
| `model_version` | 指明物理模型版本 |
| `contract_version` | 指明模块交接语义 |
| `schema_version` | 指明配置/数据结构 |

历史结果不应就地修改；新模型应使用新版本字段和新 case，保留旧 hash 链。

# 17. 软件状态和科学状态

| 项目 | 软件状态 | 科学状态 |
|---|---|---|
| PETSc 2D solver | 求解器、诊断和测试可用 | 参考/开发工况；新实验几何需另行校准 |
| Stage4/Stage5 native GHz | Jefimenko 到接收器链路已开发验证 | RFTrustReport 稀疏 bins 内 `TRUSTED_PHYSICS`；Stage5 `FULL_MAXWELL_REFERENCE_PENDING` |
| 350 MHz system | G3/H3 开发路径可复用 | `SYSTEM_350MHZ_VALIDATION=NOT_MEASURED` |
| 350 MHz native RF | 状态门已实现 | `NATIVE_RF_350MHZ=NOT_RESOLVED`，不是零 |
| H3 full-wave | port→structure/receiver 后处理已开发验证 | `FULL_WAVE_LOADING_MISMATCH_HIGH`；`FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED` |
| H4 receiver | native field→receiver 链路已开发验证 | `H4_ABSOLUTE_AMPLITUDE_STATUS=NUMERICAL_REFERENCE_ONLY`；mesh sensitivity 存在 |
| Stage-I validation | 合同、parser、metrics、不确定度和 synthetic dry-run 可用 | `STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT` |
| 公开科学验证 | 软件 v0.1.0 已发布 | `PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE=false` |

# 18. 当前已知限制

1. Afivo、openEMS/CSXCAD 和 COMSOL 需分别安装，不在 wheel/sdist 中。
2. 统一 CLI 不直接启动 PETSc 生产运行、Afivo、openEMS、COMSOL、Stage-F 大型负载或 thermal/RLC 重算。
3. `system_rf` 是冻结 H3 contract 复用，不是新 full-wave 求解。
4. CLI `validate` 目前只直接支持 WP-I-B synthetic bundle dry-run；真实数据需使用 Stage-I 合同、脚本/API。
5. Stage-F 信任仅适用于冻结稀疏 GHz bins；350 MHz native 未解析。
6. H3 加载反馈未耦合，H4 绝对幅值仅为数值参考，Stage5 全 Maxwell 参考待完成。
7. 生产 Tx/Rx 几何和真实 VNA/scope/校准/不确定度尚未提供。
8. `solver3d/afivo_reference/` 的历史说明/配置含本机绝对路径；新机器上要用 `AFIVO_STREAMER_ROOT` 和本地路径替换。
9. `system_rf` 是固定 350 MHz 的 contract 复用入口；v0.1.0 不会根据其他 `target_MHz` 重算 H3。

# 19. 常用任务速查表

| 我要做什么 | 准备什么 | 用哪个模块 | 执行什么 | 看哪个结果 |
|---|---|---|---|---|
| 看电极场增强 | 几何/电压 | C1、Afivo E1 或 COMSOL | `build/bin/stage_c1_smoke` 或外部 workflow | `Emax`、field CSV/线剖面 |
| 跑 2D streamer | PETSc、电压/种子/时步 | C2 | `build/bin/stage_c2_dynamic ...` | `diagnostics.csv`、`fields_final.csv`、`summary.txt` |
| 看电极电流 | C3 参数 | C3 | `build/bin/stage_c3_current ...` | `current_diagnostics.csv` |
| 跑 3D streamer | Afivo、`.cfg`、hook | D/E | 外部 `streamer3d CONFIG` | Silo/raw CSV 和 `results_summary/` |
| 计算 native RF | 冻结 Stage-F source/contract | F | `rf/**/generate_*.py` | RFTrustReport、waveform/spectrum/report |
| 看 RLC/port | G1/G2/G3 冻结输入 | G | `thermal/g*/generate_*.py` | `g2_*`、`g3_port_uniform.csv` |
| 看接收端 | H2 raw/H3 contract/H4 native field | H3/H4 | `fullwave/h3/generate_h3_results.py` 或 H4 generator | received spectrum/timeseries/contract |
| 检查 350 MHz 路径 | H3 contract | CLI `system_rf` | `microgap-rf run examples/configs/system_350mhz_dry_run.yaml` | `data/result.json`、manifest |
| 导入 VNA | `.s1p/.s2p`+校准/reference plane | Stage I | Stage-I parser/API；synthetic 可用 CLI validate | S 参数 metrics/quality gate |
| 导入 scope | waveform CSV+采样/触发/重复元数据 | Stage I | Stage-I scope API/script | event metrics、FFT、SNR |
| 复现 synthetic 距离比较 | 冻结 WP-I-B/D synthetic events | WP-I-D | `validation/stage_i/wp_i_d/generate_wp_i_d.py` | distance statistics/fit/ledger |
| 复现 synthetic 方向比较 | 冻结 WP-I-D synthetic 多角度 events | WP-I-D | 同上 | orientation statistics/fit/contrast |
| 生成 case 报告 | 已有 manifest/status | report | `microgap-rf report CASE_ID` | `report/report.md` |
| 检查环境 | 已激活环境 | doctor | `microgap-rf doctor` | AVAILABLE/MISSING JSON |
| 分析新论文影响 | registry 中的 paper ID | literature | `microgap-rf literature impact PAPER_ID` | 建议 JSON，不修改代码 |

# 20. 推荐日常工作流程

```text
1. 为新计算/实验分配唯一 case_id，不覆盖冻结结果。
2. 记录几何、坐标、单位、电压、气体、距离、方向和参考平面。
3. 运行 microgap-rf doctor，确认本次所需 backend 可用。
4. 轻量 case 使用统一 YAML；大型物理使用对应 Stage 入口。
5. 保存原始输入、config、backend 版本、命令、hash 和原始输出。
6. 先通过数值/质量门，再解读频谱、幅值或机制。
7. 真实实验按 Stage-I 重入合同导入，不修改 RAW。
8. 仅在 reference plane、接收器、几何和校准兼容时比较绝对幅值。
9. 更新 discrepancy ledger 和 validation status，生成 case 报告。
10. 新论文先进入 registry/evidence/decision 流程；确认需变更后再建 CR。
```

真实距离/方向数据不应直接交给硬编码 synthetic 路径的 `generate_wp_i_d.py`；应复用 `streamer_rf.validation` 统计 API 并按 Stage-I 重入合同建立新摄取层。
