# Final project structure

## cpp/

C++17 求解器源码，包括网格、Poisson、SP3、ISG-0、三粒子流注求解器、Stage 4/5 运行入口和 C++ 单元测试。

## python/

Python 配置调度、独立后处理、趋势分析、绘图和总结生成脚本。后续新课题工况应优先复用 `python/streamer_rf/` 中的诊断和辐射后处理模块。

## config/

正式运行配置和阶段工况参数。启动新课题工况时应复制现有 Stage 4/5 配置模板，再替换几何、电压、种子和输运参数。

## tests/

Python 测试与阶段回归测试。修改求解器或后处理后应先运行 pytest，再运行阶段 validator。

## tools/

阶段关闭验证器、证据审计、库存审计、清理脚本、清理完整性验证器和 HTML 总结验证器。

## docs/

阶段关闭报告、验证矩阵、决策记录、公式审计、课题迁移说明、清理报告、课题书 Markdown 总结和 HTML 总结。

## results/

GitHub 归档状态下，`results/` 已删除，不随仓库上传。该目录是运行输出目录，新仿真会重新生成 `results/`。历史 Stage 关闭报告中保留的 `results/...` 路径是历史证据路径说明，不代表当前归档仍携带原始结果。

## archive/

GitHub 归档状态下，`archive/audit_only/` 已删除，不随仓库上传。历史失效结果归档不再作为仓库内容保留。

## build/

CMake/Ninja 构建目录。可通过 `cmake --build build --parallel` 重新构建。

## .venv/

本地 Python 虚拟环境。`.gitignore` 已排除该目录；个人 GitHub 归档不应上传 `.venv/`。在新机器上应重新创建环境并安装所需依赖。

## 常用命令

重新构建：

```bash
cmake --build build --parallel
```

运行 C++ 测试：

```bash
ctest --test-dir build --output-on-failure
```

运行阶段验证（需要重新生成或恢复 `results/` 证据后才适用）：

```bash
.venv/bin/python tools/validate_stage1_closure.py
.venv/bin/python tools/validate_stage2_closure.py
.venv/bin/python tools/validate_stage3_resource_preparation.py
.venv/bin/python tools/validate_stage4_closure.py
.venv/bin/python tools/validate_stage5_closure.py
```

短时 smoke test：

```bash
.venv/bin/python tools/validate_cleanup_integrity.py
```

打开 HTML 总结：

```bash
xdg-open docs/simulation_reconstruction_summary/index.html
```

后续新课题工况建议从 Stage 4/5 的 YAML 配置和 Python 调度脚本派生，先替换电极几何、边界条件和电压波形，再接入实际输运参数和外部电路模型。
