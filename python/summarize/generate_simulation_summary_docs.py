#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
SUMMARY = DOCS / "simulation_reconstruction_summary"
ASSETS = SUMMARY / "assets"
DATA = SUMMARY / "data"
SUMMARY.mkdir(parents=True, exist_ok=True)
ASSETS.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)


def read_one(path: str) -> pd.Series:
    return pd.read_csv(ROOT / path).iloc[0]


def fmt_ns(v: float) -> str:
    return f"{v * 1e9:.3f} ns"


def fmt_pct(v: float) -> str:
    return f"{v * 100:.2f}%"


def fmt_sci(v: float, unit: str = "") -> str:
    return f"{v:.3e}{(' ' + unit) if unit else ''}"


def metrics() -> dict[str, object]:
    s4_event = read_one("results/stage4/collision/collision_event.csv")
    s4_pulse = read_one("results/stage4/current_moment/delta_current_pulse_metrics.csv")
    s4_deriv = pd.read_csv(ROOT / "results/stage4/radiation/delta_current_derivative.csv")
    s4_band = pd.read_csv(ROOT / "results/stage4/radiation/band_energy.csv")
    s5_comp = pd.read_csv(ROOT / "results/stage5/trends/case_comparison.csv")
    high = s5_comp[s5_comp["case_id"] == "F"].iloc[0]
    gap = s5_comp[s5_comp["factor"] == "seed_gap"].iloc[0]
    asym = s5_comp[s5_comp["factor"] == "seed_sigma_asymmetry"].iloc[0]
    s5_pulse = read_one("results/stage5/current_moment/highfield_delta_current_pulse_metrics.csv")
    s5_freq = read_one("results/stage5/radiation/highfield_frequency_trust.csv")
    s5_band = pd.read_csv(ROOT / "results/stage5/radiation/highfield_band_energy.csv")
    return {
        "s4_event": s4_event,
        "s4_pulse": s4_pulse,
        "s4_deriv_peak": float(s4_deriv["dDeltaI_dt_local_poly"].abs().max()),
        "s4_band": s4_band,
        "high": high,
        "gap": gap,
        "asym": asym,
        "s5_pulse": s5_pulse,
        "s5_freq": s5_freq,
        "s5_band": s5_band,
    }


def band_text(df: pd.DataFrame) -> str:
    return "；".join(f"{r.band}={r.integrated_energy_J:.3e} J（{r.trusted_status}）" for r in df.itertuples())


def write_project_book(m: dict[str, object]) -> None:
    s4 = m["s4_event"]
    s4p = m["s4_pulse"]
    high = m["high"]
    gap = m["gap"]
    asym = m["asym"]
    text = f"""# 第一部分 仿真工具链搭建总结

## 1. 仿真方法定位

本项目建立的核心仿真链路不是通过 COMSOL Plasma Module 完成，而是以自主开发的 C++17 求解器为主体。PETSc 负责稀疏矩阵、线性方程组和并行求解接口，MPI 负责多进程区域划分、并行归约和一致性验证。Python 负责配置生成、运行调度、独立后处理、绘图、证据审计和阶段关闭验证。CMake 与 Ninja 用于可重复构建，HDF5、CSV 和 JSON 用于保存全场、标量历史、频谱和证据链数据。

COMSOL 在后续课题中更适合作为复杂电极几何、静电场分布和边界条件的交叉验证工具，而不是当前流注碰撞与辐射链路的主求解平台。

## 2. 已建立的模型

项目已建立二维轴对称 r-z 流体模型，求解电子、正离子和负离子三类粒子。电子方程包含漂移、扩散、碰撞电离、二体附着、三体附着、电子-离子复合以及光电离源项；正离子方程包含电离产生、复合损失和光电离产生；负离子方程包含附着产生和正负离子复合损失。电势由 Poisson 方程给出，空间电荷密度为正离子、电子和负离子的电荷差。电子通量采用 ISG-0 格式，开放边界采用 OpenCharge 轴对称近似，光电离采用三组三阶 SP3 方程及耦合 Robin 边界。

核心变量包括电子数密度 n_e、正离子数密度 n_p、负离子数密度 n_n、电势 phi、电场 E、电子迁移率 mu_e、扩散系数 D_e、电离频率 nu_i、附着频率 nu_a、光电离源 S_ph 和电流矩 I_CM。主时间推进为一阶显式有限体积/有限差分更新，并配有自适应时间步、非负性检查、守恒诊断、checkpoint 与断点续算。

## 3. 各工具协同逻辑

整体流程为：

配置文件 → Python 运行调度 → C++ 主求解器 → PETSc 矩阵与 KSP → MPI 并行计算 → HDF5/CSV 结果 → Python 物理诊断 → 电流矩 → FFT/ESD → 图表 → validator 与证据链。

配置文件给出网格、气体条件、种子参数、背景场、输出策略和终止条件。C++ 求解器负责推进粒子方程、Poisson 场、SP3 光电离和电流矩积分。PETSc/MPI 提供可扩展的线性求解与归约能力。Python 后处理从正式输出中提取流注头、碰撞时间、桥接密度、电场塌缩、电流矩增量、导数、频谱和频带能量。每个阶段的 validator 检查 run_id、原始文件、SHA256、验证矩阵和关闭报告之间的一致性。

## 4. 已完成的通用功能

已完成单种子双向流注、双种子相向传播、流注碰撞、左右 isolated 对照、碰撞增量电流矩、局部多项式与五点差分导数、FFT、ESD、VHF/UHF/SHF 频带能量、MPI 一致性、输出采样敏感性、checkpoint 断点续算和运行证据追踪。

## 5. 面向博士课题的复用方法

可直接复用的部分包括 PETSc/MPI 求解框架、网格与场数据结构、Poisson 模块、SP3 光电离模块、ISG-0 电子通量、自适应时间步、checkpoint、电流矩积分、FFT/ESD 后处理、参数扫描框架、证据审计和关闭验证器。这些模块构成后续微间隙放电仿真的底层数值基础。

## 6. 面向微间隙击穿需要替换的内容

当前模型仍使用自由空间 Gaussian 种子、均匀背景场、Morrow–Lowke 解析输运、无实际电极边界、固定空气气氛和无外部电路耦合。这些设置适合文献流程复现和工具链验证，但不能直接代表真实微间隙击穿。

建议后续课题建模顺序为：

实际电极几何 → 电极边界与电压波形 → 种子/表面发射 → 微间隙流注和击穿 → 电极回路电流 → 电流矩与原生宽频辐射 → 外部 RLC 耦合 → 接收信号。

# 第二部分 本次文献仿真复刻结论说明

## 1. 复刻目标

本次复刻的目标是理解 Shi 2019 相关仿真模型，跑通自主求解器，跑通流注碰撞和辐射后处理链，建立可迁移到本课题的计算工具和研究逻辑。项目不追求 Shi 2019 的 1:1 数值复现，也不声称已经获得作者原始输运表或逐图定量复刻能力。

## 2. Stage 1 至 Stage 5 的逻辑

Stage 1 重建文献参数、公式、FFT 和 ESD 链路。Stage 2 验证 ISG-0、Poisson OpenCharge 和 SP3 基础模块。Stage 3 将三粒子流体方程耦合为单种子双向流注求解器。Stage 4 建立双流注碰撞、isolated subtraction、电流矩和辐射链。Stage 5 在资源感知条件下验证背景场、种子间距和种子不对称趋势，并形成面向博士课题的迁移说明。

## 3. 主要结论

三粒子求解器可以稳定运行，光电离对正流注传播具有可观测影响。双种子工况可以形成相向传播和碰撞，碰撞后桥接电子密度显著上升，碰撞区域电场明显下降。左右 isolated 对照可用于提取碰撞附加电流矩，电流矩导数可进一步构建宽频辐射频谱。Stage 5 表明增强背景场会明显提前碰撞并增强电流矩，增大种子间距会延长碰撞过程并可能在资源窗口内不碰撞，不对称种子会引入传播不对称和头部识别困难。

## 4. 代表性定量结果

Stage 4 基准工况 S4-PAPERLIKE-20UM 的碰撞时间为 {fmt_ns(float(s4.t_collision_s))}，碰撞区电场下降 {fmt_pct(float(s4.field_drop_fraction))}，桥接密度增长 {float(s4.bridge_growth_ratio):.2f} 倍，Delta I 峰值为 {fmt_sci(float(s4p.Delta_I_peak_abs), "A m")}，Delta I SNR 为 {float(s4p.Delta_I_snr):.2f}，局部多项式导数峰值约 {fmt_sci(float(m["s4_deriv_peak"]), "A m s^-1")}。Stage 4 频带能量为：{band_text(m["s4_band"])}。

Stage 5 高背景场工况 S5-HIGHFIELD-COLLISION 的碰撞时间为 {fmt_ns(float(high.t_collision))}，Delta I 峰值为 {fmt_sci(float(high.Delta_I_peak), "A m")}，谱质心约 {float(high.spectral_centroid) / 1e9:.3f} GHz，可信频率上限为 {float(m["s5_freq"].f_trust_Hz) / 1e9:.1f} GHz。高场频带能量为：{band_text(m["s5_band"])}。相对于 Stage 4 基准，高场工况碰撞明显提前，碰撞增量电流矩增大，频带能量整体提高。

长间距工况 S5-LARGERGAP-COLLISION 的状态为 {gap.collision_status}，用于说明更大初始间距会增加传播距离并可能超出资源窗口。不对称工况 S5-ASYMMETRIC-COLLISION 的状态为 {asym.collision_status}，速度比诊断约为 {float(asym.head_velocity_ratio):.3f}，说明种子尺寸不对称会改变传播和头部识别特征。

## 5. 复刻完成程度

已经完成的内容包括仿真模型结构、数值求解器、PETSc/MPI 并行、流注传播、碰撞、电流矩、FFT/ESD、趋势逻辑和课题迁移准备。

尚未完成的内容包括作者原始输运表、10/5 µm 完整网格收敛、Shi 2019 Figures 1–4 逐点复刻、完整 Maxwell 远场传播、天线与接收链路以及 1:1 定量预测。

## 6. 对课题的意义

当前成果已经形成课题仿真的底层求解器和完整工作流。下一步工作重点应从文献复刻转移到真实微间隙电极几何、电压驱动、种子形成、放电电流和外部 RLC 耦合。当前代码可以作为微间隙击穿与宽频辐射仿真的工程起点，但不能直接替代面向真实电极和外部回路的专用模型。
"""
    (DOCS / "simulation_toolchain_summary_for_project_book.md").write_text(text, encoding="utf-8")


def figure_card(num: int, src: str, title: str, run_id: str, source: str) -> str:
    return f"""
    <figure class="fig">
      <img src="assets/{src}" alt="图{num} {title}">
      <figcaption>图{num}. {title}<br><span>run_id: {run_id}；source: {source}</span></figcaption>
    </figure>
    """


def write_html(m: dict[str, object]) -> None:
    s4 = m["s4_event"]
    s4p = m["s4_pulse"]
    high = m["high"]
    gap = m["gap"]
    asym = m["asym"]
    figs = [
        figure_card(1, "figure_01_toolchain_architecture.svg", "工具链架构图", "DOCUMENTATION", "docs/stage5_closure_report.md"),
        figure_card(2, "figure_02_stage_timeline.svg", "Stage 1 至 Stage 5 时间线", "DOCUMENTATION", "docs/stage5_closure_report.md"),
        figure_card(3, "figure_03_stage3_streamer_evolution.png", "单种子双向流注演化", "S3-COARSE_ML_RECOVERY_RESUME1", "results/stage3/single_seed/coarse_ml_recovery_resume1/"),
        figure_card(4, "figure_04_sp3_on_off.png", "SP3 on/off 对比", "S3-COARSE_NO_SP3; S3-SP3_FORK_OFF_RESUME1", "results/stage3/sensitivity/photoionization_sensitivity.csv"),
        figure_card(5, "figure_05_collision_sequence.png", "双流注碰撞序列", "S4-PAPERLIKE-20UM", "results/stage4/figures/"),
        figure_card(6, "figure_06_collision_diagnostics.png", "头距、桥接密度与间隙电场", "S4-PAPERLIKE-20UM", "results/stage4/collision/collision_event.csv"),
        figure_card(7, "figure_07_current_moment.png", "collision/left/right 漂移电流矩", "S4-PAPERLIKE-20UM", "results/stage4/current_moment/"),
        figure_card(8, "figure_08_delta_current.png", "碰撞增量电流矩 Delta I", "S4-PAPERLIKE-20UM", "results/stage4/current_moment/delta_current_moment.csv"),
        figure_card(9, "figure_09_current_derivative.png", "dDelta I/dt", "S4-PAPERLIKE-20UM", "results/stage4/radiation/delta_current_derivative.csv"),
        figure_card(10, "figure_10_fft_esd.png", "FFT 与 ESD", "S4-PAPERLIKE-20UM", "results/stage4/radiation/"),
        figure_card(11, "figure_11_band_energy.png", "VHF/UHF/SHF 频带能量", "S4-PAPERLIKE-20UM; S5-HIGHFIELD-COLLISION", "results/stage4/radiation/band_energy.csv; results/stage5/radiation/highfield_band_energy.csv"),
        figure_card(12, "figure_12_stage5_trends.png", "Stage 5 趋势比较", "S5-HIGHFIELD-COLLISION; S5-LARGERGAP-COLLISION; S5-ASYMMETRIC-COLLISION", "results/stage5/trends/case_comparison.csv"),
        figure_card(13, "figure_13_project_transfer.svg", "面向微间隙击穿课题的迁移路线", "DOCUMENTATION", "docs/project_simulation_transfer.md"),
    ]
    html = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>基于C++/PETSc/MPI的流注碰撞与宽频辐射仿真工具链</title>
<style>
body{{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK SC","Microsoft YaHei",Arial,sans-serif;color:#172033;background:#f5f7fb;line-height:1.65}}
header{{background:#14213d;color:white;padding:48px 7vw 36px}}
h1{{font-size:34px;margin:0 0 8px}} h2{{margin-top:42px;border-left:5px solid #315f9b;padding-left:12px}} h3{{margin-top:28px}}
main{{max-width:1180px;margin:0 auto;padding:30px 4vw 80px;background:white}}
.subtitle{{font-size:19px;opacity:.88}} .cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin-top:26px}}
.card{{background:#eef4ff;color:#14213d;border:1px solid #bdd3f4;border-radius:12px;padding:16px;font-weight:600}}
.flow{{background:#f8fafc;border:1px solid #d9e2ef;border-radius:10px;padding:14px;margin:16px 0}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}}
.fig{{border:1px solid #d9e2ef;border-radius:10px;padding:10px;margin:0;background:#fff}}
.fig img{{width:100%;height:auto;display:block}}
figcaption{{font-size:14px;color:#344055;margin-top:8px}} figcaption span{{font-size:12px;color:#697386}}
table{{border-collapse:collapse;width:100%;margin:14px 0}} th,td{{border:1px solid #d8dee9;padding:8px;text-align:left}} th{{background:#eef4ff}}
.note{{background:#fff8e8;border-left:4px solid #d99b1a;padding:12px;margin:14px 0}} .limit{{background:#fef2f2;border-left:4px solid #c0392b;padding:12px;margin:14px 0}}
@media print{{body{{background:white}} main{{box-shadow:none}} .grid{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<header>
<h1>基于C++/PETSc/MPI的流注碰撞与宽频辐射仿真工具链</h1>
<div class="subtitle">Shi 2019仿真流程的资源感知复现及其向微间隙击穿无线传感课题的迁移</div>
<div class="cards">
<div class="card">求解器：三粒子流体 + Poisson + SP3</div>
<div class="card">数值框架：C++17 + PETSc + MPI</div>
<div class="card">后处理：Current moment + FFT + ESD</div>
<div class="card">项目状态：Workflow complete, not 1:1 reproduction</div>
</div>
</header>
<main>
<section><h2>第一章 仿真方法和工具链</h2>
<p>核心求解器由自主 C++17 代码实现，PETSc 提供线性代数与 KSP 求解，MPI 提供并行计算。Python 负责配置、调度、后处理、绘图和证据管理。COMSOL 后续仅作为复杂电极几何与静电场交叉验证工具。</p>
<div class="grid">{figs[0]}{figs[12]}</div>
<div class="flow">配置文件 → Python运行调度 → C++主求解器 → PETSc矩阵与KSP → MPI区域并行 → HDF5/CSV结果 → Python物理诊断 → 电流矩 → FFT/ESD → 图表 → validator与证据链</div>
</section>
<section><h2>第二章 文献复刻过程</h2>
<p>Stage 1 重建参数、公式和频谱链；Stage 2 验证 ISG-0、Poisson OpenCharge 和 SP3；Stage 3 验证三粒子单流注求解器；Stage 4 跑通碰撞、isolated subtraction、电流矩和辐射链；Stage 5 做趋势逻辑和课题迁移。</p>
{figs[1]}
</section>
<section><h2>第三章 单流注与光电离</h2>
<p>Stage 3 资源感知结果表明，单 Gaussian 种子可形成双向流注，SP3 关闭后正头传播不能维持同等响应。</p>
<div class="grid">{figs[2]}{figs[3]}</div>
</section>
<section><h2>第四章 双流注碰撞</h2>
<p>Stage 4 基准 S4-PAPERLIKE-20UM 的碰撞时间为 {fmt_ns(float(s4.t_collision_s))}，桥接密度增长 {float(s4.bridge_growth_ratio):.2f} 倍，碰撞区电场下降 {fmt_pct(float(s4.field_drop_fraction))}。</p>
<div class="grid">{figs[4]}{figs[5]}</div>
</section>
<section><h2>第五章 电流矩与辐射后处理</h2>
<p>主电流矩采用电子漂移电流积分，扩散电流作为敏感性。isolated subtraction 得到碰撞附加电流矩，Stage 4 Delta I 峰值为 {fmt_sci(float(s4p.Delta_I_peak_abs), "A m")}，SNR 为 {float(s4p.Delta_I_snr):.2f}。</p>
<div class="grid">{figs[6]}{figs[7]}{figs[8]}{figs[9]}{figs[10]}</div>
<div class="note">频谱中的 trusted 表示 sampling/Nyquist 可信，不代表空间网格、输运模型和辐射模型已经实现高频定量收敛。</div>
</section>
<section><h2>第六章 Stage 5趋势</h2>
<p>高背景场工况 S5-HIGHFIELD-COLLISION 的碰撞时间为 {fmt_ns(float(high.t_collision))}，Delta I 峰值为 {fmt_sci(float(high.Delta_I_peak), "A m")}，谱质心约 {float(high.spectral_centroid)/1e9:.3f} GHz。长间距工况状态为 {gap.collision_status}；不对称种子工况状态为 {asym.collision_status}，速度比诊断约 {float(asym.head_velocity_ratio):.3f}。</p>
{figs[11]}
<div class="note">G 和 A 工况没有完整 isolated 对照，因此只报告 event_local_proxy；不得与 strict_delta 绝对幅值直接排序。</div>
</section>
<section><h2>第七章 面向博士课题的迁移</h2>
<p>可复用模块包括 PETSc/MPI 框架、网格与场数据结构、Poisson、SP3、ISG-0、时间推进、checkpoint、电流矩积分、FFT/ESD 和证据审计。必须替换的部分包括自由空间 Gaussian 种子、均匀背景场、Morrow–Lowke 输运、无电极边界、固定空气气氛和无外部电路耦合。</p>
{figs[12]}
</section>
<section><h2>第八章 结论与限制</h2>
<p>项目已经完成仿真模型和流程主要复刻，跑通自主 C++/PETSc/MPI 求解工具链，并完成碰撞、电流矩、FFT 和 ESD 处理，可支撑课题仿真开发。</p>
<div class="limit">限制：未完成 Shi 2019 的 1:1 定量复刻；未获得作者原始输运表；未完成严格高分辨率收敛；当前结果不直接代表实际微间隙电极击穿。</div>
</section>
</main>
</body>
</html>
"""
    (SUMMARY / "index.html").write_text(html, encoding="utf-8")
    (SUMMARY / "README.md").write_text(
        "本目录为离线 HTML 项目总结。直接用浏览器打开 `index.html` 即可；所有图片在 `assets/`，来源登记在 `data/figure_sources.csv`。\n",
        encoding="utf-8",
    )


def write_final_structure() -> None:
    text = """# Final project structure

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

正式仿真结果、run registry、derived artifact manifest、标量历史、电流矩、频谱、图表和清理审计结果。正式结论只应引用 registry、validation matrix 或 manifest 中登记的结果。

## archive/

仅保存 audit_only 归档，例如失效 Stage 2 结果和手动中断/损坏运行。归档不得用于最终科学结论。

## build/

CMake/Ninja 构建目录。可通过 `cmake --build build --parallel` 重新构建。

## .venv/

项目 Python 虚拟环境。用于运行测试、验证器、后处理和 HTML 总结生成。

## 常用命令

重新构建：

```bash
cmake --build build --parallel
```

运行 C++ 测试：

```bash
ctest --test-dir build --output-on-failure
```

运行阶段验证：

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
"""
    (DOCS / "final_project_structure.md").write_text(text, encoding="utf-8")


def main() -> None:
    m = metrics()
    write_project_book(m)
    write_html(m)
    write_final_structure()
    print("Generated simulation summary markdown, HTML and final structure document.")


if __name__ == "__main__":
    main()
