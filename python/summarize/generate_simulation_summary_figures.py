#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "docs/simulation_reconstruction_summary/assets"
DATA = ROOT / "docs/simulation_reconstruction_summary/data"
SCRIPT = "python/summarize/generate_simulation_summary_figures.py"
ASSETS.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

SOURCE_ROWS: list[dict[str, str]] = []


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def add_sources(figure: str, sources: list[tuple[Path, str, str]]) -> None:
    for source, run_id, metric in sources:
        SOURCE_ROWS.append(
            {
                "figure": figure,
                "source_file": rel(source),
                "run_id": run_id,
                "metric": metric,
                "sha256": sha256(source),
                "generation_script": SCRIPT,
            }
        )


def savefig(name: str, sources: list[tuple[Path, str, str]]) -> None:
    out = ASSETS / name
    plt.tight_layout()
    plt.savefig(out, dpi=180)
    plt.close()
    add_sources(name, sources)


def write_svg(name: str, title: str, boxes: list[str], arrows: list[tuple[int, int]], sources: list[tuple[Path, str, str]]) -> None:
    w, h = 1200, 220 + 95 * math.ceil(len(boxes) / 3)
    coords: list[tuple[int, int]] = []
    for idx, _ in enumerate(boxes):
        row, col = divmod(idx, 3)
        coords.append((70 + col * 370, 95 + row * 95))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        "<style>text{font-family:Arial,'Noto Sans CJK SC',sans-serif}.title{font-size:28px;font-weight:700}.box{fill:#f6f8fb;stroke:#2b4c7e;stroke-width:2;rx:12}.label{font-size:17px}.arrow{stroke:#4b5563;stroke-width:2;marker-end:url(#m)}</style>",
        '<defs><marker id="m" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#4b5563"/></marker></defs>',
        f'<text x="60" y="45" class="title">{title}</text>',
    ]
    for i, j in arrows:
        x1, y1 = coords[i]
        x2, y2 = coords[j]
        parts.append(f'<line class="arrow" x1="{x1+250}" y1="{y1+30}" x2="{x2}" y2="{y2+30}"/>')
    for (x, y), label in zip(coords, boxes):
        parts.append(f'<rect class="box" x="{x}" y="{y}" width="250" height="60"/>')
        parts.append(f'<text class="label" x="{x+125}" y="{y+37}" text-anchor="middle">{label}</text>')
    parts.append("</svg>")
    (ASSETS / name).write_text("\n".join(parts), encoding="utf-8")
    add_sources(name, sources)


def heatmap_from_field(ax, field_csv: Path, value: str, title: str, log: bool = True) -> None:
    df = pd.read_csv(field_csv)
    pivot = df.pivot(index="r_m", columns="z_m", values=value)
    arr = pivot.to_numpy()
    if log:
        arr = np.log10(np.maximum(arr, 1.0))
        label = f"log10 {value}"
    else:
        label = value
    im = ax.imshow(
        arr,
        origin="lower",
        aspect="auto",
        extent=[df.z_m.min() * 1e3, df.z_m.max() * 1e3, df.r_m.min() * 1e3, df.r_m.max() * 1e3],
        cmap="viridis",
    )
    ax.set_title(title)
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("r (mm)")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label)


def figure_01() -> None:
    write_svg(
        "figure_01_toolchain_architecture.svg",
        "仿真工具链架构",
        ["YAML配置", "Python调度", "C++17求解器", "PETSc/MPI", "HDF5/CSV", "Python后处理", "Current moment", "FFT/ESD", "图表与报告", "Evidence validator"],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)],
        [(ROOT / "docs/stage5_closure_report.md", "DOCUMENTATION", "toolchain architecture")],
    )


def figure_02() -> None:
    write_svg(
        "figure_02_stage_timeline.svg",
        "Stage 1–5 资源感知复现路线",
        ["Stage 1\n公式/FFT", "Stage 2\nPoisson/SP3", "Stage 3\n单流注求解器", "Stage 4\n碰撞/辐射链", "Stage 5\n趋势/课题迁移"],
        [(0, 1), (1, 2), (2, 3), (3, 4)],
        [(ROOT / "docs/stage5_closure_report.md", "DOCUMENTATION", "stage timeline")],
    )


def figure_03() -> None:
    hist = ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/scalar_history.csv"
    heads = ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/head_trajectory.csv"
    field = ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/fields_final.csv"
    h = pd.read_csv(hist)
    tr = pd.read_csv(heads)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    axs[0, 0].plot(h.time_s * 1e9, h.E_max_V_m / 1e6)
    axs[0, 0].set_ylabel("Emax (MV/m)")
    axs[0, 0].set_xlabel("time (ns)")
    axs[0, 0].set_title("Stage 3 single-seed field enhancement")
    axs[0, 1].plot(h.time_s * 1e9, h.channel_ne_m_3 / 1e20, label="channel ne")
    axs[0, 1].plot(h.time_s * 1e9, h.ne_max_m_3 / 1e20, label="max ne", alpha=0.8)
    axs[0, 1].legend()
    axs[0, 1].set_ylabel("density (1e20 m$^{-3}$)")
    axs[0, 1].set_xlabel("time (ns)")
    axs[0, 1].set_title("Conductive channel density")
    axs[1, 0].plot(tr.time_s * 1e9, tr.lower_head_z_m * 1e3, label="lower")
    axs[1, 0].plot(tr.time_s * 1e9, tr.upper_head_z_m * 1e3, label="upper")
    axs[1, 0].legend()
    axs[1, 0].set_ylabel("head z (mm)")
    axs[1, 0].set_xlabel("time (ns)")
    axs[1, 0].set_title("Opposite head propagation")
    heatmap_from_field(axs[1, 1], field, "ne_m_3", "Final electron density")
    savefig("figure_03_stage3_streamer_evolution.png", [(hist, "S3-COARSE_ML_RECOVERY_RESUME1", "scalar history"), (heads, "S3-COARSE_ML_RECOVERY_RESUME1", "head trajectory"), (field, "S3-COARSE_ML_RECOVERY_RESUME1", "final field")])


def figure_04() -> None:
    src = ROOT / "results/stage3/sensitivity/photoionization_sensitivity.csv"
    df = pd.read_csv(src)
    labels = df["run"].str.replace("_", "\n")
    x = np.arange(len(df))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].bar(x, df["upper_displacement_m"] * 1e3, color="#4477aa")
    axs[0].set_xticks(x, labels)
    axs[0].set_ylabel("upper displacement (mm)")
    axs[0].set_title("Positive-head propagation weakens without SP3")
    axs[1].bar(x, df["photoionization_ahead"], color="#66aa55")
    axs[1].set_yscale("symlog", linthresh=1)
    axs[1].set_xticks(x, labels)
    axs[1].set_ylabel("Sph ahead (m$^{-3}$ s$^{-1}$)")
    axs[1].set_title("SP3-on/off photoionization source")
    savefig("figure_04_sp3_on_off.png", [(src, "S3-COARSE_NO_SP3", "photoionization off sensitivity"), (src, "S3-SP3_FORK_OFF_RESUME1", "SP3 fork-off sensitivity")])


def figure_05() -> None:
    files = [
        ROOT / "results/stage4/figures/01_double_seed_initial_density.png",
        ROOT / "results/stage4/figures/02_opposing_propagation_density.png",
        ROOT / "results/stage4/figures/03_precollision_field_enhancement.png",
        ROOT / "results/stage4/figures/04_bridge_channel_density.png",
        ROOT / "results/stage4/figures/05_postcollision_field_drop.png",
    ]
    titles = ["double seed", "opposing propagation", "pre-collision E", "bridge channel", "post-collision drop"]
    fig, axs = plt.subplots(1, 5, figsize=(15, 4))
    for ax, f, title in zip(axs, files, titles):
        ax.imshow(plt.imread(f))
        ax.set_axis_off()
        ax.set_title(title, fontsize=10)
    savefig("figure_05_collision_sequence.png", [(f, "S4-PAPERLIKE-20UM", "collision sequence panel") for f in files])


def figure_06() -> None:
    ev = ROOT / "results/stage4/collision/collision_event.csv"
    met = ROOT / "results/stage4/runs/paperlike_20um/collision_metrics.csv"
    e = pd.read_csv(ev).iloc[0]
    m = pd.read_csv(met)
    fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    t = m["time"] * 1e9
    for ax in axs:
        ax.axvline(e.t_collision_s * 1e9, color="crimson", ls="--", label="t_collision")
    gap_col = "inner_gap_m" if "inner_gap_m" in m else "d_head"
    if gap_col in m:
        axs[0].plot(t, m[gap_col] * 1e3)
        axs[0].set_ylabel("head gap (mm)")
    bridge_min = "bridge_min_ne_m_3" if "bridge_min_ne_m_3" in m else "bridge_min_ne"
    bridge_mean = "bridge_mean_ne_m_3" if "bridge_mean_ne_m_3" in m else "bridge_mean_ne"
    if bridge_min in m:
        axs[1].plot(t, m[bridge_min] / 1e20, label="bridge min")
    if bridge_mean in m:
        axs[1].plot(t, m[bridge_mean] / 1e20, label="bridge mean")
    axs[1].set_ylabel("bridge ne (1e20 m$^{-3}$)")
    axs[1].legend()
    col = "E_gap_max_V_m" if "E_gap_max_V_m" in m else ("E_gap_max" if "E_gap_max" in m else "E_max_V_m")
    axs[2].plot(t, m[col] / 1e6)
    axs[2].set_ylabel("E gap (MV/m)")
    axs[2].set_xlabel("time (ns)")
    axs[0].set_title("Stage 4 collision diagnostics")
    savefig("figure_06_collision_diagnostics.png", [(ev, "S4-PAPERLIKE-20UM", "collision event"), (met, "S4-PAPERLIKE-20UM", "collision metrics")])


def figure_07() -> None:
    files = [
        (ROOT / "results/stage4/current_moment/collision.csv", "collision"),
        (ROOT / "results/stage4/current_moment/left_isolated.csv", "left isolated"),
        (ROOT / "results/stage4/current_moment/right_isolated.csv", "right isolated"),
    ]
    plt.figure(figsize=(9, 4.5))
    for f, label in files:
        df = pd.read_csv(f)
        plt.plot(df["time"] * 1e9, df["I_CM_drift"], label=label)
    plt.xlabel("time (ns)")
    plt.ylabel("I_CM,drift (A m)")
    plt.title("Drift current moment: collision and isolated controls")
    plt.legend()
    savefig("figure_07_current_moment.png", [(f, "S4-PAPERLIKE-20UM", label) for f, label in files])


def figure_08() -> None:
    f = ROOT / "results/stage4/current_moment/delta_current_moment.csv"
    p = ROOT / "results/stage4/current_moment/delta_current_pulse_metrics.csv"
    df = pd.read_csv(f)
    pulse = pd.read_csv(p).iloc[0]
    plt.figure(figsize=(9, 4.5))
    plt.plot(df["time"] * 1e9, df["Delta_I"], label="Delta I")
    plt.axvline(pulse.Delta_I_peak_time_s * 1e9, color="crimson", ls="--", label="peak")
    plt.xlabel("time (ns)")
    plt.ylabel("Delta I_CM (A m)")
    plt.title("Collision-added current moment from isolated subtraction")
    plt.legend()
    savefig("figure_08_delta_current.png", [(f, "S4-PAPERLIKE-20UM", "Delta I"), (p, "S4-PAPERLIKE-20UM", "pulse metrics")])


def figure_09() -> None:
    f = ROOT / "results/stage4/radiation/delta_current_derivative.csv"
    df = pd.read_csv(f)
    plt.figure(figsize=(9, 4.5))
    plt.plot(df["time"] * 1e9, df["dDeltaI_dt_local_poly"], label="local polynomial")
    plt.plot(df["time"] * 1e9, df["dDeltaI_dt_uniform_five_point"], label="uniform five-point", alpha=0.8)
    plt.xlabel("time (ns)")
    plt.ylabel("dDeltaI/dt (A m s$^{-1}$)")
    plt.title("Derivative of collision-added current moment")
    plt.legend()
    savefig("figure_09_current_derivative.png", [(f, "S4-PAPERLIKE-20UM", "Delta I derivative")])


def figure_10() -> None:
    spec = ROOT / "results/stage4/radiation/delta_current_spectrum.csv"
    esd = ROOT / "results/stage4/radiation/delta_current_esd.csv"
    s = pd.read_csv(spec)
    e = pd.read_csv(esd)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    axs[0].loglog(np.maximum(s.frequency_Hz, 1.0), np.maximum(s.abs_dDeltaI_dt_transform, 1e-30))
    axs[0].set_xlabel("frequency (Hz)")
    axs[0].set_ylabel("|FFT(dDeltaI/dt)|")
    axs[0].set_title("FFT amplitude")
    axs[1].loglog(np.maximum(e.frequency_Hz, 1.0), np.maximum(e.ESD_J_per_Hz, 1e-40))
    axs[1].set_xlabel("frequency (Hz)")
    axs[1].set_ylabel("ESD (J/Hz)")
    axs[1].set_title("Energy spectral density")
    savefig("figure_10_fft_esd.png", [(spec, "S4-PAPERLIKE-20UM", "FFT"), (esd, "S4-PAPERLIKE-20UM", "ESD")])


def figure_11() -> None:
    b0 = ROOT / "results/stage4/radiation/band_energy.csv"
    hf = ROOT / "results/stage5/radiation/highfield_band_energy.csv"
    a = pd.read_csv(b0).assign(case="baseline")
    b = pd.read_csv(hf).assign(case="high field")
    df = pd.concat([a, b], ignore_index=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    pivot = df.pivot(index="band", columns="case", values="integrated_energy_J")
    pivot.plot(kind="bar", ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("band energy proxy (J)")
    ax.set_title("VHF/UHF/SHF band energy (sampling/Nyquist trusted)")
    ax.tick_params(axis="x", rotation=0)
    savefig("figure_11_band_energy.png", [(b0, "S4-PAPERLIKE-20UM", "band energy"), (hf, "S5-HIGHFIELD-COLLISION", "band energy")])


def figure_12() -> None:
    comp = ROOT / "results/stage5/trends/case_comparison.csv"
    df = pd.read_csv(comp)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    labels = df["case_id"].tolist()
    axs[0, 0].bar(labels, df["t_collision"] * 1e9)
    axs[0, 0].set_ylabel("t_collision (ns)")
    axs[0, 0].set_title("Collision-time trend")
    axs[0, 1].bar(labels, df["E_gap_drop_fraction"] * 100)
    axs[0, 1].set_ylabel("field drop (%)")
    axs[0, 1].set_title("Gap field drop")
    axs[1, 0].bar(labels, df["bridge_growth"])
    axs[1, 0].set_yscale("log")
    axs[1, 0].set_ylabel("bridge growth")
    axs[1, 0].set_title("Bridge-density growth")
    colors = ["#4477aa" if m == "strict_delta" else "#cc8844" for m in df["metric_type"]]
    axs[1, 1].bar(labels, df["I_CM_peak"], color=colors)
    axs[1, 1].set_yscale("log")
    axs[1, 1].set_ylabel("I_CM / proxy peak (A m)")
    axs[1, 1].set_title("strict_delta vs event-local proxy")
    savefig("figure_12_stage5_trends.png", [(comp, "S5-HIGHFIELD-COLLISION", "case comparison"), (comp, "S5-LARGERGAP-COLLISION", "case comparison"), (comp, "S5-ASYMMETRIC-COLLISION", "case comparison")])


def figure_13() -> None:
    write_svg(
        "figure_13_project_transfer.svg",
        "微间隙击穿课题迁移路线",
        ["实际电极几何", "局部电场增强", "种子/表面发射", "流注与击穿", "放电电流", "电流矩", "原生宽频辐射", "RLC调制", "接收信号"],
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)],
        [(ROOT / "docs/project_simulation_transfer.md", "DOCUMENTATION", "project transfer route")],
    )


def main() -> None:
    for fn in [figure_01, figure_02, figure_03, figure_04, figure_05, figure_06, figure_07, figure_08, figure_09, figure_10, figure_11, figure_12, figure_13]:
        fn()
    with (DATA / "figure_sources.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["figure", "source_file", "run_id", "metric", "sha256", "generation_script"])
        writer.writeheader()
        writer.writerows(SOURCE_ROWS)
    print(f"Generated {len(list(ASSETS.glob('figure_*')))} summary figures.")


if __name__ == "__main__":
    main()
