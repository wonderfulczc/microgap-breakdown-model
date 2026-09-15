# 统一 CLI 架构说明

`microgap-rf` 是现有模块的薄层，不包含新的物理算法。

| CLI 命令 | 现有模块或合同 | 后端 | 输入 | 输出 |
|---|---|---|---|---|
| `doctor` | Stage-J 后端清单及系统探测 | 核心/可选外部后端 | 环境 | 可用性 JSON |
| `run` smoke | Stage F/G/H/I 可导入接口及最终合同 | 核心 Python | case YAML | 统一 case 目录 |
| `run` system_rf | 冻结 H3 结果合同 | openEMS 结果复用 | case YAML/H3 合同 | 开发状态与来源清单 |
| `validate` | `streamer_rf.validation.wp_i_b` | Stage-I | 测量/合成数据合同 | 验证 dry-run 结果 |
| `report` | `run_manifest.json` 与 `status.json` | 无 | `CASE_ID` | Markdown 报告 |
| `literature impact` | `streamer_rf.literature` | 无 | registry paper ID | 影响建议 JSON |

未直接暴露为统一 CLI 的冻结入口仍保持原样：`cpp/apps/stage_c1_smoke.cpp`、
`stage_c2_dynamic.cpp`、`stage_c3_current.cpp` 对应 PETSc/MPI 核心；`rf/**/generate_*.py`、
`thermal/g1-g3`、`fullwave/h1-h5` 和 `validation/stage_i` 脚本对应阶段级合同生成与参考运行。
这些入口只有在其 backend-specific 配置与资源策略被正式映射后才会加入 `simulate` 子命令，
RP-1 不包装可能意外启动大型仿真的命令。

命令返回 0 表示成功，配置、路径或运行错误返回 2。可选外部后端缺失不会使 `doctor`
崩溃；核心依赖缺失会明确列入 `REQUIRED_MISSING`。配置规范见
[配置与结果合同](配置与结果合同.md)。
