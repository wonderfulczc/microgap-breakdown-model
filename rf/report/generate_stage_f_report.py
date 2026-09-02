#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.spectral.bands import PROJECT_BANDS, STANDARD_BANDS  # noqa: E402
from streamer_rf.rf.spectral.cwt import band_energy_vs_time, morlet_cwt  # noqa: E402

REPORT = ROOT / "rf/report"
FIG_DIR = REPORT / "figures"
F1 = ROOT / "rf/source/audit/stage_f1_synthetic_audit.json"
F2 = ROOT / "rf/jefimenko/validation/stage_f2_validation.json"
F3 = ROOT / "rf/spectral/validation/stage_f3_validation.json"
F4 = ROOT / "rf/production/f4_attribution/stage_f4_attribution_summary.json"
F4_STAGE = ROOT / "rf/production/f4_attribution/stage_rf_summary.csv"

plt.rcParams.update(
    {
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel(ax, label: str) -> None:
    ax.text(-0.13, 1.06, label, transform=ax.transAxes, fontweight="bold", va="top", ha="left")


def savefig(fig: plt.Figure, name: str) -> tuple[str, str]:
    png = FIG_DIR / f"{name}.png"
    pdf = FIG_DIR / f"{name}.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(png.relative_to(ROOT)), str(pdf.relative_to(ROOT))


def add_stage_spans(ax, labels: pd.DataFrame, *, time_offset_s: float = 0.0) -> None:
    colors = ["#e9eef7", "#eef7e9", "#fff4d9", "#fde9e7", "#ece7f5"]
    ymin, ymax = ax.get_ylim()
    for i, row in enumerate(labels.itertuples()):
        x0 = (row.start_time_s + time_offset_s) * 1e9
        x1 = (row.end_time_s + time_offset_s) * 1e9
        ax.axvspan(x0, x1, color=colors[i % len(colors)], alpha=0.6, lw=0)
        ax.text(0.5 * (x0 + x1), ymax, row.stage.replace("_", "\n"), ha="center", va="top", fontsize=6.3)
    ax.set_ylim(ymin, ymax)


def plot_f1(f1: dict, f3: dict) -> tuple[str, str]:
    fig = plt.figure(figsize=(7.3, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gs[0, 0])
    panel(ax, "(a)")
    steps = ["AMR\nsource", "J=-eGamma", "conservative\nremap", "Jefimenko", "RF audits", "trusted\nspectrum"]
    ax.plot(range(len(steps)), np.zeros(len(steps)), "-o", color="#1f5f85", lw=1.8)
    for i, s in enumerate(steps):
        ax.text(i, 0.08, s, ha="center", va="bottom")
    ax.set_xlim(-0.4, len(steps) - 0.6)
    ax.set_ylim(-0.15, 0.35)
    ax.axis("off")

    ax = fig.add_subplot(gs[0, 1])
    panel(ax, "(b)")
    labels = ["Q error", "M error"]
    cons = [f1["conservative_remap"]["relative_charge_error"], f1["conservative_remap"]["relative_M_error"]]
    bad = [f1["bad_nearest_remap"]["relative_charge_error"], f1["bad_nearest_remap"]["relative_M_error"]]
    x = np.arange(2)
    ax.bar(x - 0.18, cons, width=0.36, label="conservative", color="#3b7a57")
    ax.bar(x + 0.18, bad, width=0.36, label="bad nearest", color="#b04a4a")
    ax.set_xticks(x, labels)
    ax.set_ylabel("relative error")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.legend(frameon=False)

    ax = fig.add_subplot(gs[0, 2])
    panel(ax, "(c)")
    ax.bar(["conservative", "bad nearest"], [f1["dQdt_artifact_conservative"], f1["dQdt_artifact_bad_nearest"]], color=["#3b7a57", "#b04a4a"])
    ax.set_ylabel(r"$|dQ/dt|$ artifact")

    ax = fig.add_subplot(gs[1, 0])
    panel(ax, "(d)")
    ax.bar(["conservative", "bad nearest"], [f1["dMdt_artifact_conservative_norm"], f1["dMdt_artifact_bad_nearest_norm"]], color=["#3b7a57", "#b04a4a"])
    ax.set_ylabel(r"$|dM/dt|$ artifact")

    ax = fig.add_subplot(gs[1, 1:])
    panel(ax, "(e)")
    trust = f3["rf_trust_report"]
    names = ["derivative", "interpolation", "mesh", "sampling", "Nyquist"]
    vals = [
        trust["derivative_trust_frequency_Hz"],
        trust["interpolation_trust_frequency_Hz"],
        trust["mesh_trust_frequency_Hz"],
        trust["sampling_trust_frequency_Hz"],
        trust["nyquist_Hz"],
    ]
    ax.bar(names, np.asarray(vals) / 1e6, color="#557a95")
    ax.axhline(trust["trusted_frequency_high_Hz"] / 1e6, color="#b04a4a", ls="--", label="trusted upper bound")
    ax.set_ylabel("frequency (MHz)")
    ax.set_title("trusted upper frequency = minimum active audit bound")
    ax.legend(frameon=False)
    ax.tick_params(axis="x", rotation=20)
    return savefig(fig, "Fig_F1_rf_trustworthiness")


def plot_f2(f2: dict) -> tuple[str, str]:
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 5.2), constrained_layout=True)
    r = np.array([0.1, 0.2, 0.4, 0.8])
    near = r ** f2["static_charge"]["near_field_slope"]
    rad = r ** f2["radiation"]["radiation_1_over_R_slope"]
    ax = axs[0, 0]
    panel(ax, "(a)")
    ax.loglog(r, near / near[0], "o-", color="#1f5f85")
    ax.set(xlabel="distance (m)", ylabel="normalized near field", title=f"slope {f2['static_charge']['near_field_slope']:.4f}")
    ax.grid(alpha=0.25)
    ax = axs[0, 1]
    panel(ax, "(b)")
    ax.loglog(r, rad / rad[0], "o-", color="#7a5c1f")
    ax.set(xlabel="distance (m)", ylabel="normalized radiation field", title=f"slope {f2['radiation']['radiation_1_over_R_slope']:.4f}")
    ax.grid(alpha=0.25)
    ax = axs[0, 2]
    panel(ax, "(c)")
    ax.bar(["R1", "R2"], [100 * f2["radiation"]["current_moment_error_R1"], 100 * f2["radiation"]["current_moment_error_R2"]], color="#557a95")
    ax.set_ylabel("CM vs full error (%)")
    ax.set_title("far-field convergence")
    ax = axs[1, 0]
    panel(ax, "(d)")
    ax.bar(["E/B", "angular"], [100 * f2["radiation"]["E_over_B_max_relative_error"], f2["radiation"]["angular_sin_theta_max_abs_error"]], color=["#557a95", "#3b7a57"])
    ax.set_yscale("log")
    ax.set_ylabel("relative / absolute error")
    ax = axs[1, 1]
    panel(ax, "(e)")
    ax.bar(["delay"], [f2["radiation"]["propagation_delay_error_s"]], color="#7a5c1f")
    ax.set_yscale("log")
    ax.set_ylabel("arrival delay error (s)")
    ax = axs[1, 2]
    panel(ax, "(f)")
    names = ["static E", "steady B", "interp fine", "interp coarse"]
    vals = [
        f2["static_charge"]["relative_error"],
        f2["steady_current_magnetic"]["relative_error"],
        f2["radiation"]["interpolation_fine_relative_error"],
        f2["radiation"]["interpolation_coarse_relative_error"],
    ]
    ax.bar(names, vals, color="#5b6f8f")
    ax.set_yscale("log")
    ax.set_ylabel("validation error")
    ax.tick_params(axis="x", rotation=25)
    return savefig(fig, "Fig_F2_jefimenko_validation")


def cwt_for_waveform(wave: pd.DataFrame, trust: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = wave["time_s"].to_numpy(dtype=float)
    y = wave["Ez_dJ"].to_numpy(dtype=float)
    f0 = max(trust["trusted_frequency_low_Hz"], 1.0 / max(t[-1] - t[0], 1e-30))
    f1 = trust["trusted_frequency_high_Hz"]
    freq = np.geomspace(f0, f1, 48)
    cwt = morlet_cwt(t, y, freq, n_cycles=5.0)
    power = np.where(cwt.coi_mask, cwt.power, np.nan)
    return cwt.time_s, cwt.frequency_Hz, power


def plot_dataset_column(fig: plt.Figure, axes: list[plt.Axes], case_id: str, title: str) -> None:
    summary = load_json(ROOT / f"rf/production/f4_attribution/{case_id}_summary.json")
    trust = summary["rf_trust_report"]
    source_dir = ROOT / summary["source_dir"]
    scalar = pd.read_csv(source_dir / "scalar_history.csv")
    cm = pd.read_csv(source_dir / "current_moment.csv")
    wave = pd.read_csv(ROOT / f"rf/production/f4_attribution/{case_id}_jefimenko_waveform.csv")
    labels = pd.read_csv(ROOT / f"rf/production/f4_attribution/{case_id}_stage_labels.csv")
    t_ns = scalar["time_s"].to_numpy(dtype=float) * 1e9
    axes[0].plot(t_ns, scalar["total_electrons"], color="#1f5f85", label="Ne")
    axes[0].set_title(title)
    axes[0].set_ylabel("total electrons")
    axes[0].set_yscale("log")
    add_stage_spans(axes[0], labels)
    axes[1].plot(cm["time"] * 1e9, cm["I_CM_electron"], color="#3b7a57", label="M")
    dM = np.gradient(cm["I_CM_electron"].to_numpy(dtype=float), cm["time"].to_numpy(dtype=float))
    axes[1].plot(cm["time"] * 1e9, dM / max(np.max(np.abs(dM)), 1e-300) * np.max(np.abs(cm["I_CM_electron"])), color="#b04a4a", lw=1, label="scaled dM/dt")
    axes[1].set_ylabel("M (A m)")
    add_stage_spans(axes[1], labels)
    axes[1].legend(frameon=False, loc="upper left")
    axes[2].plot(wave["time_s"] * 1e9, wave["Ez_dJ"], color="#1f5f85")
    axes[2].set_ylabel(r"$E_{dJ,z}$ (V/m)")
    add_stage_spans(axes[2], labels, time_offset_s=0.20 / 299792458.0)
    ct, cf, cp = cwt_for_waveform(wave, trust)
    im = axes[3].pcolormesh(ct * 1e9, cf / 1e9, cp, shading="auto", cmap="magma")
    axes[3].axhspan(0, trust["trusted_frequency_low_Hz"] / 1e9, color="0.8", alpha=0.6)
    axes[3].axhspan(trust["trusted_frequency_high_Hz"] / 1e9, max(cf) / 1e9, color="0.8", alpha=0.6)
    axes[3].set_ylabel("frequency (GHz)")
    axes[3].set_yscale("log")
    axes[3].set_title(f"trusted {trust['trusted_frequency_low_Hz']/1e9:.2f}-{trust['trusted_frequency_high_Hz']/1e9:.2f} GHz")
    fig.colorbar(im, ax=axes[3], fraction=0.046, pad=0.02, label="CWT power")
    band_names = ["SHF", "3-10 GHz"]
    stage = pd.read_csv(ROOT / f"rf/production/f4_attribution/{case_id}_stage_rf_summary.csv")
    x = np.arange(len(stage))
    for name in band_names:
        col = f"{name}_energy"
        if col in stage:
            vals = pd.to_numeric(stage[col], errors="coerce").to_numpy(dtype=float)
            if np.any(np.isfinite(vals)):
                axes[4].plot(x, vals, "o-", label=name)
    axes[4].set_yscale("symlog", linthresh=1e-30)
    axes[4].set_xticks(x, [s.replace("_", "\n") for s in stage["stage"]], rotation=0)
    axes[4].set_ylabel("trusted/partial band CWT energy")
    axes[4].legend(frameon=False, loc="best")


def plot_f3() -> tuple[str, str]:
    fig, axs = plt.subplots(5, 2, figsize=(8.2, 8.8), constrained_layout=True)
    plot_dataset_column(fig, list(axs[:, 0]), "F4-P-stage4-left-isolated", "Stage4 isolated streamer")
    plot_dataset_column(fig, list(axs[:, 1]), "F4-C-stage5-highfield-collision", "Stage5 collision")
    labels = ["(a)", "(f)", "(b)", "(g)", "(c)", "(h)", "(d)", "(i)", "(e)", "(j)"]
    for ax, lab in zip(axs.ravel(), labels):
        panel(ax, lab)
    for ax in axs[-1, :]:
        ax.set_xlabel("stage")
    return savefig(fig, "Fig_F3_stage_time_frequency")


def plot_f4(stage: pd.DataFrame) -> tuple[str, str]:
    rows = [f"{r.dataset}:{r.stage}" for r in stage.itertuples()]
    cols = ["VHF", "UHF", "1-3 GHz", "3-10 GHz", "SHF"]
    status_col = {"1-3 GHz": "1_3GHz_status", "3-10 GHz": "3_10GHz_status"}
    unresolved_by_policy = {"VHF", "UHF", "1-3 GHz"}
    data = np.full((len(rows), len(cols)), np.nan)
    status = [["UNTRUSTED"] * len(cols) for _ in rows]
    for i, row in stage.iterrows():
        for j, band in enumerate(cols):
            status[i][j] = "UNTRUSTED" if band in unresolved_by_policy else row[status_col.get(band, f"{band}_status")]
            if status[i][j] != "UNTRUSTED":
                val = row.get(f"{band}_energy", np.nan)
                data[i, j] = float(val) if pd.notna(val) else np.nan
    finite = np.isfinite(data)
    norm = data.copy()
    if np.any(finite):
        norm[finite] = data[finite] / max(np.nanmax(data), 1e-300)
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    masked = np.ma.masked_invalid(norm)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("0.82")
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(cols)), cols)
    ax.set_yticks(np.arange(len(rows)), rows)
    ax.set_title("Stage-resolved trusted spectral fingerprints")
    for i in range(len(rows)):
        for j in range(len(cols)):
            text = "NOT\nRESOLVED" if status[i][j] == "UNTRUSTED" else status[i][j].replace("_", "\n")
            ax.text(j, i, text, ha="center", va="center", fontsize=6.5, color="black")
    fig.colorbar(im, ax=ax, label="normalized trusted/partial CWT energy")
    return savefig(fig, "Fig_F4_stage_frequency_fingerprint")


def write_html(figs: dict[str, tuple[str, str]], payload: dict) -> None:
    sections = [
        (
            "Fig. F1",
            "Conservation-aware RF trustworthiness framework for streamer-radiation simulations",
            "图F1  守恒约束的射频源可信度框架。该图展示从等离子体源项、连续性一致电流、保守重映射、Jefimenko积分到频谱可信区间判定的完整链路，并用受控算例说明非保守重映射会在时间导数中引入伪射频信号。",
            "保守重映射在受控AMR测试中保持Q和M相对误差为0，而bad nearest重映射产生约0.113的Q/M误差。对应的dQ/dt artifact从0.03339增至0.2062，dM/dt artifact从0.006678增至0.04124。",
            "非守恒后处理足以制造人工RF；因此后续频谱必须先通过源项守恒、连续性和可信频带审核。",
            "可以支撑RF后处理需要守恒门控和可信频带门控的论文方法结论。",
            "这些数值来自controlled/synthetic validation，不能被解释为真实放电频谱强度。",
            "方法贡献",
        ),
        (
            "Fig. F2",
            "Validation of the full Jefimenko electromagnetic-field solver",
            "图F2  完整Jefimenko电磁场求解器验证。图中给出静电近场、稳恒电流磁场、辐射场距离标度、电流矩近似与完整Jefimenko积分、E/B关系和传播延迟误差。",
            "近场距离斜率为-2.0000，辐射项斜率为-0.9467。角向误差为2.13e-5，E/B相对误差为5.71%，传播延迟误差为1.03e-25 s。电流矩近似误差从R1的19.4%降至R2的1.93%。",
            "在远场条件改善时，电流矩近似向完整Jefimenko解收敛；完整积分保留了空间源分布和延迟传播信息。",
            "可以支撑完整Jefimenko求解器和远场电流矩近似收敛性的数值验证结论。",
            "不能声称电流矩在近场或任意观测距离均可替代Jefimenko积分。",
            "数值验证",
        ),
        (
            "Fig. F3",
            "Temporal correspondence between discharge evolution and trusted radio-frequency emission",
            "图F3  放电演化与可信射频辐射的时间对应关系。左列为Stage4单流注传播，右列为Stage5高场碰撞窗口；每列依次给出放电诊断、电流矩、远场辐射波形、CWT时频图和可信频带能量。",
            "Stage4在avalanche、inception和early propagation期间均出现处于2.94-7.97 GHz可信区间内的可解析辐射活动。Stage5在encounter/collision/post-collision窗口中伴随快速current-moment变化，并在3.05-10.23 GHz可信区间内产生增强的GHz/low-SHF响应。",
            "快速电流矩变化控制远场主要辐射时序；碰撞/相互作用窗口比平滑早期传播更容易产生高频时间结构。",
            "可以支撑GHz/low-SHF子区间与早期传播和碰撞窗口存在时间对应关系。",
            "VHF/UHF由于当前时间窗不足，没有得到可解释的stage attribution；不能写成传播主要产生VHF，也不能写成collision产生全部GHz辐射。",
            "核心主图",
        ),
        (
            "Fig. F4",
            "Stage-resolved trusted spectral fingerprints of microgap discharge dynamics",
            "图F4  微间隙放电阶段分辨的可信频谱指纹。矩阵仅对通过RFTrustReport的频带或其重叠子区间显示能量；未可信频段以灰色NOT RESOLVED显示。",
            "当前可信结果集中在GHz/low-SHF子区间。VHF、UHF和1-3 GHz均为NOT_RESOLVED，这表示时间窗和可信频带不足以归因，而不是证明这些频段没有辐射。",
            "该图把stage mask、CWT能量和RFTrustReport合并，防止把不可解释频段误读为低能量或无辐射。",
            "可以支撑Stage F当前只完成GHz/low-SHF子区间pilot归因的边界清晰结论。",
            "不能建立普适的stage-frequency阈值，也不能从未解析频带推断物理缺失。",
            "重要机制支撑",
        ),
    ]
    html_sections = []
    for idx, (fig, title, caption, results, physics, supported, cannot, role) in enumerate(sections, 1):
        key = f"F{idx}"
        png, _ = figs[key]
        html_sections.append(
            f"""
<section>
<h2>{fig}</h2>
<img src="{html.escape(png.replace('rf/report/', ''))}" alt="{html.escape(fig)}">
<h3>拟投稿图题</h3><p>{html.escape(title)}</p>
<h3>中文图注</h3><p>{html.escape(caption)}</p>
<h3>Results式结果描述</h3><p>{html.escape(results)}</p>
<h3>物理解释</h3><p>{html.escape(physics)}</p>
<h3>可支撑的论文结论</h3><p>{html.escape(supported)}</p>
<h3>当前不能得出的结论</h3><p>{html.escape(cannot)}</p>
<h3>推荐论文角色</h3><p>{html.escape(role)}</p>
</section>
"""
        )
    matrix = pd.read_csv(REPORT / "summary_matrix.csv").to_html(index=False, classes="matrix")
    text = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>Stage F Trusted RF Simulation Smoke-Test and Paper-Figure Assessment</title>
<style>
body{{font-family:Arial, sans-serif; max-width:1100px; margin:32px auto; line-height:1.55; color:#222;}}
img{{max-width:100%; border:1px solid #ddd;}}
h1,h2,h3{{line-height:1.25;}}
h3{{margin-bottom:0.2rem;}}
table{{border-collapse:collapse; width:100%; font-size:0.9rem;}}
td,th{{border:1px solid #ccc; padding:0.35rem; text-align:left;}}
.note{{background:#f6f6f6; padding:0.8rem; border-left:4px solid #777;}}
</style></head><body>
<h1>Stage F Trusted RF Simulation Smoke-Test and Paper-Figure Assessment</h1>
<p>本报告汇总Stage F1-F4的真实输出，用于评估“原生宽频电磁辐射仿真”在未来投稿论文中的图件潜力。Stage F已经建立守恒可信源、完整Jefimenko场求解、频谱可信带审核和阶段分辨归因pilot，但当前仍属于development/smoke-level validation，不是最终实验校准仿真。</p>
<p class="note">关键限制：没有用zero padding、人工静态尾段或UNTRUSTED频带生成物理结论。Stage E三维源仍非正式continuity-consistent scientific RF source，因此Fig. F5标记为DATA_INSUFFICIENT。</p>
{''.join(html_sections)}
<section>
<h2>Fig. F5</h2>
<h3>拟投稿图题</h3><p>Three-dimensional polarization and radiation characteristics</p>
<h3>状态</h3><p><strong>FIG_F5_STATUS = DATA_INSUFFICIENT</strong></p>
<h3>说明</h3><p>F3已经验证polarization和radiation-direction算法能力，但Stage E三维source仍是drift-only历史输出，缺少continuity-consistent 3D scientific J_RF、source time series、RFTrustReport和far-field audit的完整闭环。因此本报告不生成Fig. F5物理结果图。</p>
<h3>推荐论文角色</h3><p>暂不建议使用。</p>
</section>
<section>
<h2>Summary Matrix</h2>
{matrix}
</section>
<section>
<h2>总体判断</h2>
<p><strong>当前Stage F是否已经证明“VHF/UHF/GHz分别来自哪个放电阶段”？PARTIAL。</strong> 已经支持GHz/low-SHF子区间与Stage4早期传播及Stage5 encounter/collision窗口的时间关联；VHF、UHF和1-3 GHz仍未解析。</p>
<p><strong>RF trustworthiness framework贡献评级：methodological contribution。</strong> 最强方法证据是保守重映射、Jefimenko解析验证、RFTrustReport和current-moment/Jefimenko一致性。最强物理证据是Stage5 collision window在可信GHz/low-SHF子区间内的增强响应。</p>
<p>正式论文前必须重跑：更长物理有效VHF/UHF窗口、正式mesh sensitivity、continuity-consistent 3D source与polarization/direction case、实验/Stage B约束几何和电压。可直接保留的方法验证：F1 remap trust、F2 Jefimenko analytic validation、F3 FFT/CWT/trust synthetic validation。</p>
</section>
</body></html>
"""
    (REPORT / "stage_f_smoke_paper_report.html").write_text(text)


def write_manifests(figs: dict[str, tuple[str, str]]) -> None:
    entries = {
        "FIG_F1": {
            "final_files": figs["F1"],
            "source_files": [str(F1.relative_to(ROOT)), str(F3.relative_to(ROOT))],
            "caption": "守恒约束的RF可信度框架和重映射artifact验证。",
            "conclusion": "non-conservative remap can generate artificial RF derivatives.",
            "limitation": "controlled synthetic validation, not discharge physics spectrum.",
        },
        "FIG_F2": {
            "final_files": figs["F2"],
            "source_files": [str(F2.relative_to(ROOT))],
            "caption": "完整Jefimenko求解器解析/半解析验证。",
            "conclusion": "full Jefimenko solver passes static, magnetic, radiation, delay, and far-field checks.",
            "limitation": "method validation only.",
        },
        "FIG_F3": {
            "final_files": figs["F3"],
            "source_files": [
                "rf/production/f4_attribution/F4-P-stage4-left-isolated_jefimenko_waveform.csv",
                "rf/production/f4_attribution/F4-C-stage5-highfield-collision_jefimenko_waveform.csv",
                "rf/production/f4_attribution/F4-P-stage4-left-isolated_stage_labels.csv",
                "rf/production/f4_attribution/F4-C-stage5-highfield-collision_stage_labels.csv",
            ],
            "caption": "Stage4/Stage5阶段分辨时频对应。",
            "conclusion": "trusted GHz/low-SHF activity is temporally resolved for early propagation and collision windows.",
            "limitation": "VHF/UHF not resolved.",
        },
        "FIG_F4": {
            "final_files": figs["F4"],
            "source_files": [str(F4_STAGE.relative_to(ROOT))],
            "caption": "阶段-频段可信频谱指纹矩阵。",
            "conclusion": "only trusted or partially trusted GHz/low-SHF cells receive energy values.",
            "limitation": "not a universal stage-frequency threshold.",
        },
        "FIG_F5": {
            "status": "DATA_INSUFFICIENT",
            "source_files": [],
            "caption": "3D polarization/direction not generated.",
            "conclusion": "deferred until continuity-consistent 3D scientific J source exists.",
            "limitation": "Stage E 3D data are drift-only for scientific RF purposes.",
        },
    }
    for item in entries.values():
        source_hashes = {}
        for src in item.get("source_files", []):
            p = ROOT / src
            if p.exists():
                source_hashes[src] = sha256(p)
        item["source_sha256"] = source_hashes
    (REPORT / "figure_manifest.yaml").write_text(yaml.safe_dump(entries, sort_keys=False, allow_unicode=True))
    rows = []
    for fig, item in entries.items():
        for src in item.get("source_files", []):
            p = ROOT / src
            rows.append(
                {
                    "figure": fig,
                    "source_file": src,
                    "exists": p.exists(),
                    "sha256": sha256(p) if p.exists() else "",
                    "trust_report": "rf/production/f4_attribution/stage_f4_attribution_summary.json" if fig in {"FIG_F3", "FIG_F4"} else "",
                    "limitation": item.get("limitation", ""),
                }
            )
    pd.DataFrame(rows).to_csv(REPORT / "data_manifest.csv", index=False)


def write_summary_matrix(stage: pd.DataFrame) -> None:
    cols = ["VHF", "UHF", "1-3 GHz", "3-10 GHz", "SHF"]
    status_col = {"1-3 GHz": "1_3GHz_status", "3-10 GHz": "3_10GHz_status"}
    unresolved_by_policy = {"VHF", "UHF", "1-3 GHz"}
    rows = []
    for _, row in stage.iterrows():
        out = {"dataset": row["dataset"], "stage": row["stage"]}
        for band in cols:
            status = "UNTRUSTED" if band in unresolved_by_policy else row[status_col.get(band, f"{band}_status")]
            if status == "UNTRUSTED":
                out[band] = "NOT_RESOLVED"
            elif row["stage_frequency_attribution"] == "SUPPORTED" and status == "PARTIALLY_TRUSTED":
                out[band] = "PARTIALLY_SUPPORTED"
            elif row["stage_frequency_attribution"] == "SUPPORTED":
                out[band] = "SUPPORTED"
            else:
                out[band] = "NOT_RESOLVED"
        rows.append(out)
    pd.DataFrame(rows).to_csv(REPORT / "summary_matrix.csv", index=False)


def main() -> None:
    REPORT.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    f1 = load_json(F1)
    f2 = load_json(F2)
    f3 = load_json(F3)
    f4 = load_json(F4)
    stage = pd.read_csv(F4_STAGE)
    figs = {
        "F1": plot_f1(f1, f3),
        "F2": plot_f2(f2),
        "F3": plot_f3(),
        "F4": plot_f4(stage),
    }
    write_summary_matrix(stage)
    write_html(figs, f4)
    write_manifests(figs)
    print(json.dumps({"figures": figs, "FIG_F5_STATUS": "DATA_INSUFFICIENT"}, indent=2))


if __name__ == "__main__":
    main()
