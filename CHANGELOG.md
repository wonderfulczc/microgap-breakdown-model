# 更新日志

本项目遵循语义化版本的开发、候选和正式发行序列。科学状态与软件版本相互独立。

## 0.1.0 - Initial release

首个正式软件版本提供：

- `microgap-rf` 统一 CLI 与配置驱动工作流；
- PETSc/MPI 二维 streamer 框架与 Afivo 三维外部后端接口；
- Jefimenko 本征 RF、热通道/RLC 与 openEMS 全波接口；
- Stage-I 验证架构及距离、方位、重复性和不确定度处理管线；
- 文献维护、科学变更治理与资源友好的可复现检查。

Scientific validation remains pending real experiments. 软件发行不代表物理模型已完成真实实验验证。

## 0.1.0rc1 - Initial release candidate

首个候选版本提供：

- PETSc/MPI 二维 streamer 求解框架；
- Afivo 三维外部后端接口；
- Jefimenko 本征 RF 分析；
- 热通道、RLC 和物理端口链路；
- openEMS 外部全波接口与接收机传递链；
- Stage-I 仿真/测量验证框架；
- `microgap-rf` 统一 CLI 与配置驱动工作流；
- 文献证据维护和科学模型变更治理。

本版本的真实实验科学验证仍待完成：
`STAGE_I_SCIENTIFIC_VALIDATION=PENDING_REAL_EXPERIMENT`。候选版本不表示生产系统或
350 MHz 本征放电 RF 已经得到实验验证。
