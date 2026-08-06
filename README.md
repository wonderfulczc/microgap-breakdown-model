# 微间隙击穿仿真模型工具链

本仓库保存的是一个面向微间隙击穿研究的 C++/PETSc/MPI 流注流体仿真工具链。它最初通过对 Shi 2019 流注碰撞与宽频辐射仿真流程的资源感知复现建立，随后被整理为可上传 GitHub 的轻量化求解器、工具链和工作流基础。

本仓库包含：

- C++17 流注流体求解器核心组件；
- PETSc、MPI 与 CMake 构建集成；
- Python 配置、运行调度、后处理和绘图工具；
- C++ 与 Python 测试；
- 阶段验证器和证据审计工具；
- Stage 1 至 Stage 5 的方法文档、关闭报告和项目迁移说明；
- 面向微间隙电极放电仿真的后续开发参考。

为了便于个人 GitHub 归档，历史大规模仿真结果已经删除。因此，本仓库应被理解为“求解器、仿真环境、工具链和流程归档”，不是包含完整原始仿真结果的证据包。

## 项目目的

本项目的直接目的是保留一套已经跑通的数值框架，使其可以从简化自由空间流注算例逐步迁移到真实微间隙放电模型。目标研究流程为：

实际微间隙电极几何
→ 电极电势边界和实际电压波形
→ 电极附近种子电子/表面发射
→ 微间隙流注与击穿
→ 电极回路时变电流
→ 电流矩与原生宽频辐射
→ 外部 RLC 耦合
→ 接收端信号

当前代码只覆盖上述链路的一部分，包括流注流体输运、Poisson 空间电荷场、SP3 光电离、ISG-0 电子通量、checkpoint 风格运行逻辑、电流矩积分、FFT/ESD 后处理和验证基础设施。

## 推荐运行环境

推荐使用 Linux 或 Windows + WSL2 Ubuntu。

不推荐直接在 Windows 原生环境下编译，因为当前 CMake 和 Python 工作流直接使用 PETSc、MPI、`pkg-config` 和 `mpirun`。这些依赖在 Linux/WSL2 下配置更稳定。

## 依赖安装

Ubuntu 或 WSL2 Ubuntu 下安装系统依赖：

```bash
sudo apt update
sudo apt install -y build-essential cmake ninja-build pkg-config \
  openmpi-bin libopenmpi-dev libpetsc-real-dev python3-venv
```

创建 Python 虚拟环境并安装 Python 依赖：

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## 构建方法

配置并构建 C++ 程序：

```bash
cmake -S . -B build -G Ninja
cmake --build build --parallel
```

运行 C++ 单元测试：

```bash
ctest --test-dir build --output-on-failure
```

运行 Python 导入 smoke test：

```bash
PYTHONPATH=python PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -c \
  "import streamer_rf.stage4, streamer_rf.stage5; print('python smoke passed')"
```

## 仓库结构

仓库主要目录如下：

- `cpp/`：C++ 求解器源码、PETSc/MPI 数值模块和 C++ 测试；
- `python/`：运行脚本、分析脚本、电流矩和辐射后处理工具；
- `config/`：阶段算例和运行配置模板；
- `tests/`：Python 回归测试；
- `tools/`：验证器、证据审计、清理和总结验证工具；
- `docs/`：方法说明、阶段关闭报告、课题迁移说明和最终总结；
- `requirements.txt`：Python 依赖列表；
- `CMakeLists.txt`：C++ 构建入口。

## 使用方式

典型使用流程如下：

1. 在 `config/` 下修改或新增算例配置；
2. 使用 CMake/Ninja 构建 C++ 可执行程序；
3. 使用对应 Python driver 启动运行，例如 `python/stage4/run_stage4.py` 或 `python/stage5/finalize_resource_run.py`；
4. 使用对应 Python 分析脚本处理输出；
5. 在存在结果证据时，使用 `tools/` 中的 validator 进行一致性检查。

新仿真会在本地重新生成 `results/` 目录。该目录已被 `.gitignore` 忽略，不会上传到 GitHub。

## 当前科学状态

本项目在资源感知标准下完成了求解器、流注碰撞、辐射后处理链、趋势逻辑和课题迁移准备。它没有完成 Shi 2019 的 1:1 定量复刻，也没有获得作者原始输运表。

历史 Stage closure validator 依赖已经删除的 `results/` 目录。这些验证器作为历史审计工具保留；只有重新生成或恢复对应结果证据后，旧阶段关闭验证才会再次完整通过。

## 当前局限性

本仓库还不是完整的微间隙击穿模型。后续仍需要实现或替换：

- 真实电极几何；
- 电极电势边界；
- 实际电压波形输入；
- 电极附近种子电子或表面发射模型；
- 电极回路时变电流；
- 外部 RLC 耦合；
- 接收端或天线信号模型；
- 面向目标气体、气压和间隙尺度的输运参数校准。

此外，本项目不包含 Shi 2019 作者原始输运表，也不声称已经逐点复现 Shi 2019 的图 1 至图 4。

## 文档索引

- `docs/github_archive_assessment.md`：GitHub 归档范围和删除历史结果后的影响；
- `docs/project_simulation_transfer.md`：如何将当前工具链迁移到微间隙课题；
- `docs/simulation_toolchain_summary_for_project_book.md`：可写入课题书的中文总结；
- `docs/simulation_reconstruction_summary/index.html`：可离线打开的图文总结页面；
- `docs/final_project_structure.md`：最终项目目录结构说明。

## 上传前已完成的最小验证

上传前已经执行并通过：

```bash
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=python PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -c \
  "import streamer_rf.stage4, streamer_rf.stage5; print('python smoke passed')"
```

## 重要说明

本仓库用于保存仿真环境、求解器、工具和研究流程。若需要恢复完整历史结果，需要重新运行对应阶段算例或从外部备份恢复 `results/` 目录。
