# Microgap RF Replica 科研软件工具包

## 项目目标与软件定位

本项目是面向微间隙放电、streamer、本征射频、热火花通道、电路、全波传播、
接收机响应和实验交叉验证的多物理科研软件。核心由 C++17/PETSc/MPI 求解器、
可导入的 `streamer_rf` Python 包、外部 Afivo/openEMS/COMSOL 接口和可追溯数据合同组成。
当前软件工具链已完成，但真实实验科学验证仍待进行。

当前候选版本为 `0.1.0rc1`。变更记录见 [CHANGELOG](CHANGELOG.md)，候选版本说明见
[0.1.0rc1 release notes](release/0.1.0rc1_release_notes.md)。该候选版本尚未公开发布。

## 主要能力与 A-J 架构

- A：配置、生命周期和数值基础；B：外部 COMSOL 静电/几何接口。
- C：二维 PETSc streamer、电极电流及冷态到热态交接量。
- D-E：外部 Afivo 三维交叉验证；F：Jefimenko 本征场与可信频谱。
- G：LTE 热通道、RLC 双向耦合和瞬态端口；H：openEMS 结构与接收机链路。
- I：测量合同、合成 dry-run 和真实数据重入；J：包装与可复现性。

权威交接关系见 [软件架构](docs/zh/软件架构.md)，当前状态见
[科学状态](docs/zh/科学状态.md)。

## CLI 快速开始

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/pip install -e .
microgap-rf --help
microgap-rf --version
microgap-rf doctor
microgap-rf run examples/configs/smoke.yaml
microgap-rf validate examples/configs/validation_synthetic.yaml
microgap-rf report rp1_smoke
```

完整命令映射见 [CLI 架构说明](docs/zh/CLI架构说明.md)。CLI 是薄层，科学逻辑仍位于
`streamer_rf` 模块及冻结阶段脚本中。

## Python API 与配置文件

设置 `PYTHONPATH=python` 后可直接 `import streamer_rf`。统一 YAML 配置包含公共的
`case_id`、`workflow`、`backend`、`output` 和版本字段，并允许工作流专用的
`frequency`、`inputs` 等段。配置在执行前检查单位、范围、路径、后端和科学状态冲突。

## 外部后端

Afivo-streamer 与 openEMS/CSXCAD 均保持外置，分别通过 `AFIVO_STREAMER_ROOT`、
`OPENEMS_ROOT` 和 `OPENEMS_PYTHON` 配置。COMSOL 是 Stage-B 专有外部工作流。
后端提交与接口见 `packaging/external_backends.json`，源码不被复制进本仓库。

## 结果目录与实验验证

新 CLI case 写入 `results/<case_id>/`，包含输入/解析配置、运行清单、来源追踪、日志、
数据、图表、报告和状态。历史冻结结果不迁移。Stage-I 合成夹具仅用于
`SYNTHETIC_DRY_RUN`；真实 VNA、示波器、几何、校准及不确定度数据按
[实验重入](docs/zh/实验数据重入.md)导入，不改变核心 API。

## 科学状态边界

- `STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT`。
- `PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE=false`。
- 350 MHz 本征等离子体 RF 为 `NOT_RESOLVED`，不能解释为零。
- H3 全波负载反馈未耦合；H4 绝对幅值仅是数值参考；Stage5 全 Maxwell 参考待完成。

## 文献与模型变更维护

新文献必须依次经过登记、证据卡、专家审阅、影响分析、决策记录及可选变更请求，
之后才能进行最小实现和回归测试。详见[文献维护](docs/zh/文献维护流程.md)和
[科学变更管理](docs/zh/科研模型变更流程.md)。新结果记录 `model_version`、
`contract_version`、`schema_version`，不得覆盖旧冻结结果。

## 数据、复现、许可证与引用

小型合同和参考数据保留；多 GB Afivo/Stage-F 数据外部归档或重建。详见
[数据政策](docs/zh/数据政策.md)。轻量复现命令为 `./scripts/reproduce_smoke.sh`。
本项目代码采用 [Apache-2.0](LICENSE) 许可证，软件作者署名为 `Zach`，引用信息见
[`CITATION.cff`](CITATION.cff)。当前仅达到 release candidate 准备条件，尚未创建公开发行。
`STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT`，软件发行准备状态不得解释为
物理模型已经通过真实实验验证。
