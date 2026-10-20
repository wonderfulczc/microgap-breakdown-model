# I1仿真工具：源码入口与发布范围

本次将此前本地GENERIC_DEV中的工具代码纳入Git，不上传测试文件、
原生库、执行程序、检查点、轨迹、实验CSV或测试结果。既有Git历史中的测试不删除。
提交代码并不改变历史数学/工程Gate，也不代表实际器件或RF实验已通过认证。

## 入口和复用关系

| 功能 | 入口 |
|---|---|
| 参数化装置及来源校验 | `adapters/i1_p1/device.py` |
| 接触分离TENG及冻结位置电荷源 | `adapters/i1_p1/sources.py`、`adapters/i1_teng_source/source.py` |
| 串联电路及端口事务 | `adapters/i1_teng_source/circuit.py`、`adapters/i1_fp3_port_circuit/` |
| 冷/给定热态/源网络统一入口 | `adapters/i1_cold_physics_v3_5/scenario.py` |
| E2E、ROM、传播和接收 | `adapters/i1_e2e/pipeline.py` |
| 原生冷放电、离子、表面发射及SP3 | `cpp/src/streamer/`、`cpp/src/sp3_geometry.cpp` |
| 持久接受态恢复 | `adapters/i1_cold_physics_v3_1/checkpoints.py` |
| 给定热态及热柱串联反馈 | `adapters/i1_fp4_handoff/given_thermal.py`、`thermal_series/provider.py` |
| 原生物种空间源、导体和回流 | `python/streamer_rf/rf/source/` |
| Toepler/RW形成阶段候选 | `adapters/i1_spark_formation/spark_formation_models.py` |
| 空间冻结电导及均匀反应诊断 | `adapters/i1_spatial_bridge/network.py`、`adapters/i1_conductive_channel/kinetics.py` |
| COMSOL/参数化B0几何接入 | `adapters/comsol_fp1/`、`adapters/comsol_fp2/` |

这些适配器按仓库根目录运行，目前不是独立发布的pip插件：

```bash
PYTHONPATH=python:. python -m adapters.i1_e2e.pipeline \
  path/to/scenario.json path/to/new_output_directory
```

输出目录必须是新目录。输入引用的文件须存在，SHA必须与实际内容一致。
`source_contract.file`等来源路径相对于仓库根目录；不要将说明中的历史
`validation/...`路径误当成此次Git提交已经包含的数据。
依赖沿用项目的`pyproject.toml`和requirements文件，不包含本地虚拟环境。

配置、状态和物理来源合同见：

- [E2E接口及运行说明](i1_implementation/i1_predictive_pipeline_architecture.md)
- [ROM推导与单位](i1_implementation/I1_E2E_V1_MODEL_DERIVATION.md)
- [Part VI端口与能量合同](i1_part_vi_port_circuit_energy_contract.md)
- [实际装置测绘指南](I1_EXPERIMENTAL_ASSEMBLY_MEASUREMENT_GUIDE.md)
- `config/i1/`中的声明式schema；运行时语义检查仍以提供器实现为准。

## 原生构建

沿用原项目MPI/PETSc依赖，在新的构建目录配置：

```bash
cmake -S . -B build-i1-local -DI1_BUILD_LOCAL_VALIDATION_TOOLS=OFF
cmake --build build-i1-local --target i1_cold_scenario i1_cold_v31_recovery \
  i1_cold_v32_smoke i1_cold_v34_smoke b0_profile_static
```

`*_smoke`是现有统一入口使用的版本化原生数据驱动程序名称，
保留名称以兼容既有场景；不是把本次排除的测试程序重新上传。
新增加的本地验证目标默认关闭。只有本地另有相应验证源文件时才应开启
`I1_BUILD_LOCAL_VALIDATION_TOOLS`。本次不上传也不覆盖任何已验收二进制。
新构建有新的源码/库哈希，不可自动恢复另一个构建的持久检查点。

## 文档和数据边界

`docs/i1_implementation/`保存重要实现说明的原文副本，
覆盖冷物理V3至V3.5、FP4交接、P1源/装置和E2E。
它们是历史执行记录，内部相对路径仍以原本地结果目录为上下文。
报告引用的CSV、保护清单、测试记录和执行源快照按要求留在本地，不随此次提交发布。
长期路线图中的历史结果链接同样可能指向未发布的本地研究档案。

本次发布整理仅增加：

1. 排除本地测试时仍可配置的CMake开关；
2. E2E等效源元数据运行模块，消除对制造测试脚本的导入依赖；
3. 源码入口、合同schema及说明文档索引。

这不重跑历史科学轨迹，不修改其结果；先前执行版SHA与此次发布版SHA可能不同。
来源绑定输入、有效域拒绝和历史文件保护逻辑保持，不绕过缺失参数。

## 当前能力限制

- E2E功能链已经运行，但合成冒烟的启动快RC模态欠分辨；
  不能将完整链路成功解释成全部时域波形精度通过。
- 一维ROM是条件性物理模型，不等价于二维FULL；短窗粒子差异尚无独立收敛界。
- 给定LTE热态是独立参考，不自动继承冷段材料能量。
- 自动冷热交接仍为`MISSING_PHYSICAL_CLOSURE`。
- 真实TENG源初态、离子/表面标定、全局端口、实际天线及接收RF仍未验证。
- 测试配置中的电压、有效截面和电容不是实际装置的默认标定值。

具体工程范围与限制见[实施报告](i1_implementation/I1_E2E_V1_IMPLEMENTATION_REPORT.md)。
