#!/usr/bin/env python3
"""Generate the Stage E paper-style smoke-test report and figures."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-streamer-rf-replica")

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, Polygon, Rectangle


REPORT = Path(__file__).resolve().parent
STAGE_E = REPORT.parent
REPO = STAGE_E.parents[2]
COMMON = REPO / "solver3d/afivo_reference/common_benchmark"
D3 = COMMON / "d3"
RAW = STAGE_E / "results_raw"
SUMMARY = STAGE_E / "results_summary"
FIGURES = REPORT / "figures"

E_CHARGE = 1.602176634e-19
GAP_UM = 70.0


plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "font.size": 8,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_figure(fig: mpl.figure.Figure, name: str) -> dict[str, str]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    png = FIGURES / f"{name}.png"
    pdf = FIGURES / f"{name}.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return {
        "png": str(png.relative_to(REPORT)),
        "pdf": str(pdf.relative_to(REPORT)),
        "png_sha256": sha256(png),
        "pdf_sha256": sha256(pdf),
    }


def label(ax: mpl.axes.Axes, text: str) -> None:
    ax.text(
        0.015,
        0.965,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.0},
    )


def parse_afivo_log(path: Path, has_z_cm: bool = False) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    user0 = 33
    user = {
        "vol_E_gt_0p5_Emax_m3": user0,
        "vol_E_gt_0p8_Emax_m3": user0 + 1,
        "conductor_volume_m3": user0 + 2,
        "Emax_user_x_m": user0 + 3,
        "Emax_user_y_m": user0 + 4,
        "Emax_user_z_m": user0 + 5,
        "min_dx_m": user0 + 6,
        "x_cm_m": user0 + 7,
        "y_cm_m": user0 + 8,
        "r_cm_m": user0 + 9,
        "head_z_m": user0 + 10,
        "total_charge_C": user0 + 11,
        "min_ne_m3": user0 + 12,
        "ne_moment_asym": user0 + 13,
        "rho_moment_asym": user0 + 14,
        "total_electrons_user": user0 + 15,
    }
    if has_z_cm:
        user["z_cm_m"] = user0 + 16
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("it "):
            continue
        data = line.split()
        if len(data) < user0:
            continue
        row = {
            "it": float(data[0]),
            "time_s": float(data[1]),
            "dt_s": float(data[2]),
            "total_electrons": float(data[4]),
            "total_positive_ions": float(data[5]),
            "sum_charge_number": float(data[6]),
            "Emax_Vpm": float(data[8]),
            "ne_max_m3": float(data[12]),
            "voltage_V": float(data[16]),
            "wc_time_s": float(data[25]),
            "n_cells": float(data[26]),
            "min_dx_standard_m": float(data[27]),
            "highest_level": float(data[32]),
        }
        for key, idx in user.items():
            row[key] = float(data[idx]) if idx < len(data) else math.nan
        rows.append(row)
    return pd.DataFrame(rows)


def read_lineout(path: Path) -> pd.DataFrame:
    header = path.read_text().splitlines()[0].lstrip("#").split()
    return pd.read_csv(path, sep=r"\s+", comment="#", names=header, skiprows=1)


def source_path_from_summary(summary: dict[str, Any], case: str) -> Path:
    rel = summary["dynamic_cases"][case]["source_file"]
    return STAGE_E / rel


def e3_source_path(summary: dict[str, Any], case: str) -> Path:
    rel = summary["cases"][case]["source_file"]
    return STAGE_E / rel


def source_slice(
    path: Path,
    plane: str,
    variables: tuple[str, ...],
    lsf_gas_only: bool = True,
) -> pd.DataFrame:
    coord = "y_m" if plane == "xz" else "x_m"
    s = "x_m" if plane == "xz" else "y_m"
    usecols = ["x_m", "y_m", "z_m", "lsf_m", *variables]
    df = pd.read_csv(path, usecols=usecols)
    if lsf_gas_only:
        df = df[df["lsf_m"] > 0.0]
    min_abs = df[coord].abs().min()
    sl = df[df[coord].abs() <= min_abs + 1.0e-15].copy()
    sl["s_um"] = sl[s] * 1.0e6
    sl["z_um"] = sl["z_m"] * 1.0e6
    return sl


def field_slice(path: Path, plane: str) -> pd.DataFrame:
    coord = "y_m" if plane == "xz" else "x_m"
    s = "x_m" if plane == "xz" else "y_m"
    df = pd.read_csv(path, usecols=["x_m", "y_m", "z_m", "Eabs_Vpm", "lsf_m"])
    df = df[df["lsf_m"] > 0.0]
    min_abs = df[coord].abs().min()
    sl = df[df[coord].abs() <= min_abs + 1.0e-15].copy()
    sl["s_um"] = sl[s] * 1.0e6
    sl["z_um"] = sl["z_m"] * 1.0e6
    return sl


def scatter_map(
    ax: mpl.axes.Axes,
    data: pd.DataFrame,
    value: str,
    title: str,
    cmap: str = "magma",
    log: bool = False,
    cbar_label: str = "",
) -> None:
    vals = data[value].to_numpy()
    if log:
        positive = vals[np.isfinite(vals) & (vals > 0.0)]
        vmin = max(float(np.nanpercentile(positive, 2)) if positive.size else 1.0, 1.0)
        vmax = max(float(np.nanpercentile(positive, 99.8)) if positive.size else 10.0, vmin * 10.0)
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    sc = ax.scatter(
        data["s_um"],
        data["z_um"],
        c=vals,
        s=1.4,
        cmap=cmap,
        norm=norm,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(title)
    ax.set_xlabel("transverse coordinate (um)")
    ax.set_ylabel("z (um)")
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(cbar_label)


def plot_geometry_axisymmetric(ax: mpl.axes.Axes) -> None:
    ax.add_patch(Rectangle((-8, 70), 16, 50, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((0, 70), 5, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Rectangle((-40, -3), 80, 3, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.plot([0, 0], [0, 70], "k--", lw=0.8)
    ax.annotate("70 um gap", xy=(8, 35), xytext=(18, 35), arrowprops={"arrowstyle": "-|>", "lw": 0.8})
    ax.set_xlim(-45, 45)
    ax.set_ylim(-8, 125)
    ax.set_aspect("equal")
    ax.set_xlabel("r / x (um)")
    ax.set_ylabel("z (um)")
    ax.set_title("PETSc 2D section")


def plot_geometry_3d_aligned(ax: mpl.axes.Axes) -> None:
    ax.add_patch(Rectangle((-5, 70), 10, 55, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((0, 70), 5, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Rectangle((-5, -55), 10, 55, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((0, 0), 5, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.plot([0, 0], [0, 70], "k--", lw=0.8)
    ax.set_xlim(-35, 35)
    ax.set_ylim(-60, 130)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("z (um)")
    ax.set_title("Afivo 3D Cartesian")


def plot_triangular_schematic(ax: mpl.axes.Axes) -> None:
    tri = np.array([[0, 70], [-20, 120], [20, 120]])
    ax.add_patch(Polygon(tri, closed=True, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((0, 70), 3, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Rectangle((-35, -3), 70, 3, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.text(0, 124, "40 um width", ha="center", va="bottom", fontsize=7)
    ax.text(13, 70, "rounded tip", ha="left", va="center", fontsize=7)
    ax.set_xlim(-40, 40)
    ax.set_ylim(-8, 130)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("z (um)")
    ax.set_title("development triangular foil")


def plot_needle_pair(ax: mpl.axes.Axes, offset_um: float, title: str) -> None:
    ax.add_patch(Rectangle((-5, 70), 10, 55, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((0, 70), 5, facecolor="#d08c2c", edgecolor="black", lw=0.8))
    ax.add_patch(Rectangle((offset_um - 5, -55), 10, 55, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.add_patch(Circle((offset_um, 0), 5, facecolor="#7d7d7d", edgecolor="black", lw=0.8))
    ax.plot([0, 0], [0, 70], "k--", lw=0.6)
    if offset_um:
        ax.annotate("10 um", xy=(offset_um, 8), xytext=(0, 8), arrowprops={"arrowstyle": "<->", "lw": 0.8}, ha="center")
    ax.set_xlim(-30, 35)
    ax.set_ylim(-60, 130)
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("z (um)")
    ax.set_title(title)


def fig_e1() -> tuple[dict[str, str], dict[str, Any]]:
    d3 = pd.read_csv(D3 / "results/stage_d3_comparison.csv")
    metrics = load_json(D3 / "results/stage_d3_metrics.json")
    t = d3["time_s"] * 1.0e12
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8), constrained_layout=True)
    plot_geometry_axisymmetric(axs[0, 0])
    label(axs[0, 0], "(a)")
    plot_geometry_3d_aligned(axs[0, 1])
    label(axs[0, 1], "(b)")
    axs[0, 2].plot(t, d3["petsc_Emax"] / 1e6, "o-", label="PETSc 2D")
    axs[0, 2].plot(t, d3["afivo_Emax"] / 1e6, "s-", label="Afivo true 3D")
    axs[0, 2].set_xlabel("time (ps)")
    axs[0, 2].set_ylabel("Emax (MV/m)")
    axs[0, 2].set_title("Emax(t)")
    axs[0, 2].legend(frameon=False)
    label(axs[0, 2], "(c)")
    axs[1, 0].plot(t, d3["petsc_ne_max"] / 1e15, "o-", label="PETSc 2D")
    axs[1, 0].plot(t, d3["afivo_ne_max"] / 1e15, "s-", label="Afivo true 3D")
    axs[1, 0].set_xlabel("time (ps)")
    axs[1, 0].set_ylabel("ne,max (1e15 m$^{-3}$)")
    axs[1, 0].set_title("ne,max(t)")
    label(axs[1, 0], "(d)")
    axs[1, 1].plot(t, d3["petsc_total_e"], "o-", label="PETSc 2D")
    axs[1, 1].plot(t, d3["afivo_total_e"], "s-", label="Afivo true 3D")
    axs[1, 1].set_xlabel("time (ps)")
    axs[1, 1].set_ylabel("total electrons")
    axs[1, 1].set_title("N_e(t)")
    label(axs[1, 1], "(e)")
    err_names = ["Emax", "ne,max", "N_e", "head z"]
    err_vals = [
        metrics["max_Emax_rel_err"] * 100,
        metrics["max_ne_max_rel_err"] * 100,
        metrics["max_total_electrons_rel_err"] * 100,
        metrics["max_head_position_diff_m"] * 1e6,
    ]
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]
    axs[1, 2].bar(err_names[:3], err_vals[:3], color=colors[:3])
    axr = axs[1, 2].twinx()
    axr.bar([err_names[3]], [err_vals[3]], color=colors[3])
    axs[1, 2].set_ylabel("max relative error (%)")
    axr.set_ylabel("head z diff (um)")
    axs[1, 2].set_title("error summary")
    label(axs[1, 2], "(f)")
    files = save_figure(fig, "Fig_E1_solver_validation")
    return files, {"metrics": metrics}


def fig_e2() -> tuple[dict[str, str], dict[str, Any]]:
    e1 = load_json(SUMMARY / "e1_electrostatic_summary.json")
    e2 = load_json(SUMMARY / "e2_dynamic_summary.json")
    refined = e2["electrostatic_refined"]
    xz = field_slice(RAW / "e2_refined_x_field_cells_000000.csv", "xz")
    yz = field_slice(RAW / "e2_refined_y_field_cells_000000.csv", "yz")
    prof = pd.DataFrame(
        {
            "s_um": read_lineout(RAW / "e2_refined_x_line_000000.txt")["x"] * 1e6,
            "E_x_MVpm": read_lineout(RAW / "e2_refined_x_line_000000.txt")["electric_fld"] / 1e6,
            "E_y_MVpm": read_lineout(RAW / "e2_refined_y_line_000000.txt")["electric_fld"] / 1e6,
        }
    )
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 5.0), constrained_layout=True)
    plot_triangular_schematic(axs[0, 0])
    label(axs[0, 0], "(a)")
    scatter_map(axs[0, 1], xz, "Eabs_Vpm", "x-z center slice |E|", cbar_label="V/m")
    label(axs[0, 1], "(b)")
    scatter_map(axs[0, 2], yz, "Eabs_Vpm", "y-z center slice |E|", cbar_label="V/m")
    label(axs[0, 2], "(c)")
    axs[1, 0].plot(prof["s_um"], prof["E_x_MVpm"], label="x-z")
    axs[1, 0].plot(prof["s_um"], prof["E_y_MVpm"], label="y-z")
    axs[1, 0].set_xlabel("transverse coordinate (um)")
    axs[1, 0].set_ylabel("|E| (MV/m)")
    axs[1, 0].set_title("selected field profiles")
    axs[1, 0].legend(frameon=False)
    label(axs[1, 0], "(d)")
    axs[1, 1].bar(
        ["medium", "refined"],
        [e1["medium_profile"]["electrostatic_asymmetry_metric"], refined["refined_asymmetry_metric"]],
        color=["#4c78a8", "#f58518"],
    )
    axs[1, 1].set_ylabel("normalized L2 asymmetry")
    axs[1, 1].set_title("asymmetry stability")
    label(axs[1, 1], "(e)")
    axs[1, 2].bar(
        ["E>0.5Emax", "E>0.8Emax"],
        [
            refined["refined_high_field_volumes_m3"]["E_gt_0p5_Emax"],
            refined["refined_high_field_volumes_m3"]["E_gt_0p8_Emax"],
        ],
        color=["#54a24b", "#e45756"],
    )
    axs[1, 2].set_yscale("log")
    axs[1, 2].set_ylabel("volume (m$^3$)")
    axs[1, 2].set_title("high-field region")
    axs[1, 2].text(
        0.05,
        0.82,
        "Emax: 28.27 -> 57.34 MV/m\nlocal tip magnitude not converged",
        transform=axs[1, 2].transAxes,
        va="top",
        fontsize=7,
    )
    label(axs[1, 2], "(f)")
    files = save_figure(fig, "Fig_E2_triangular_electrostatic")
    return files, {"e1": e1, "e2_refined": refined}


def fig_e3() -> tuple[dict[str, str], dict[str, Any]]:
    e2 = load_json(SUMMARY / "e2_dynamic_summary.json")
    log_a = parse_afivo_log(RAW / "e2_case_a_pi_off_500V_log.txt")
    log_b = parse_afivo_log(RAW / "e2_case_b_pi_on_500V_log.txt")
    src0 = RAW / "e2_case_a_pi_off_500V_source_000000.csv"
    src5 = source_path_from_summary(e2, "case_a_pi_off_500V")
    t0 = source_slice(src0, "xz", ("ne_m3", "rho_Cpm3"))
    xz5 = source_slice(src5, "xz", ("ne_m3", "rho_Cpm3"))
    yz5 = source_slice(src5, "yz", ("ne_m3", "rho_Cpm3"))
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 5.0), constrained_layout=True)
    scatter_map(axs[0, 0], t0, "ne_m3", "t=0 ps ne, x-z", log=True, cbar_label="m$^{-3}$")
    label(axs[0, 0], "(a)")
    scatter_map(axs[0, 1], xz5, "ne_m3", "t=5 ps ne, x-z", log=True, cbar_label="m$^{-3}$")
    label(axs[0, 1], "(b)")
    scatter_map(axs[0, 2], yz5, "ne_m3", "t=5 ps ne, y-z", log=True, cbar_label="m$^{-3}$")
    label(axs[0, 2], "(c)")
    scatter_map(axs[1, 0], xz5, "rho_Cpm3", "t=5 ps rho, x-z", cbar_label="C/m$^3$")
    label(axs[1, 0], "(d)")
    axs[1, 1].plot(log_a["time_s"] * 1e12, log_a["total_electrons"], "o-", label="PI OFF")
    axs[1, 1].plot(log_b["time_s"] * 1e12, log_b["total_electrons"], "s--", label="PI ON")
    axs[1, 1].set_xlabel("time (ps)")
    axs[1, 1].set_ylabel("total electrons")
    axs[1, 1].set_title("N_e(t)")
    axs[1, 1].legend(frameon=False)
    label(axs[1, 1], "(e)")
    vals = [
        e2["dynamic_cases"]["case_a_pi_off_500V"]["ne_plane_asymmetry"],
        e2["dynamic_cases"]["case_a_pi_off_500V"]["rho_plane_asymmetry"],
        e2["dynamic_cases"]["case_b_pi_on_500V"]["ne_plane_asymmetry"],
        e2["dynamic_cases"]["case_b_pi_on_500V"]["rho_plane_asymmetry"],
    ]
    axs[1, 2].bar(["ne OFF", "rho OFF", "ne ON", "rho ON"], vals, color=["#4c78a8", "#4c78a8", "#f58518", "#f58518"])
    axs[1, 2].set_ylabel("normalized L2 asymmetry")
    axs[1, 2].set_title("dynamic asymmetry")
    axs[1, 2].tick_params(axis="x", rotation=25)
    label(axs[1, 2], "(f)")
    files = save_figure(fig, "Fig_E3_triangular_dynamic")
    return files, {"e2": e2}


def fig_e4() -> tuple[dict[str, str], dict[str, Any]]:
    e3 = load_json(SUMMARY / "e3_summary.json")
    aligned_prof = pd.read_csv(SUMMARY / "e3_profiles_aligned.csv")
    mis_prof = pd.read_csv(SUMMARY / "e3_profiles_misaligned.csv")
    traj = pd.read_csv(SUMMARY / "e3_centroid_trajectory.csv")
    comparison = pd.read_csv(SUMMARY / "e3_case_comparison.csv")
    fig, axs = plt.subplots(2, 4, figsize=(9.6, 4.9), constrained_layout=True)
    plot_needle_pair(axs[0, 0], 0.0, "aligned needle pair")
    label(axs[0, 0], "(a)")
    plot_needle_pair(axs[0, 1], 10.0, "10 um misaligned pair")
    label(axs[0, 1], "(b)")
    for ax, data, title, lab in [
        (axs[0, 2], aligned_prof, "aligned |E| profiles", "(c)"),
        (axs[0, 3], mis_prof, "misaligned |E| profiles", "(d)"),
    ]:
        for plane, style in [("xz", "-"), ("yz", "--")]:
            sub = data[(data["plane"] == plane) & (data["z_m"].between(0, 70e-6))]
            prof = sub.groupby("s_m", as_index=False)["Eabs_Vpm"].max()
            ax.plot(prof["s_m"] * 1e6, prof["Eabs_Vpm"] / 1e6, style, label=plane)
        ax.set_xlabel("transverse coordinate (um)")
        ax.set_ylabel("max |E| across gap (MV/m)")
        ax.set_title(title)
        ax.legend(frameon=False)
        label(ax, lab)
    axs[1, 0].bar(
        ["aligned", "misaligned"],
        [e3["cases"]["aligned"]["field_plane_asymmetry"], e3["cases"]["misaligned"]["field_plane_asymmetry"]],
        color=["#4c78a8", "#e45756"],
    )
    axs[1, 0].set_yscale("symlog", linthresh=1e-8)
    axs[1, 0].set_ylabel("field asymmetry")
    axs[1, 0].set_title("field symmetry breaking")
    label(axs[1, 0], "(e)")
    x = np.arange(2)
    width = 0.35
    axs[1, 1].bar(
        x - width / 2,
        [e3["cases"]["aligned"]["ne_plane_asymmetry"], e3["cases"]["misaligned"]["ne_plane_asymmetry"]],
        width,
        label="ne",
        color="#4c78a8",
    )
    axs[1, 1].bar(
        x + width / 2,
        [e3["cases"]["aligned"]["rho_plane_asymmetry"], e3["cases"]["misaligned"]["rho_plane_asymmetry"]],
        width,
        label="rho",
        color="#f58518",
    )
    axs[1, 1].set_xticks(x, ["aligned", "misaligned"])
    axs[1, 1].set_ylabel("normalized L2 asymmetry")
    axs[1, 1].set_title("plasma asymmetry")
    axs[1, 1].legend(frameon=False)
    label(axs[1, 1], "(f)")
    for case, style in [("aligned", "o-"), ("misaligned", "s-")]:
        sub = traj[traj["case"] == case]
        axs[1, 2].plot(sub["time_s"] * 1e12, sub["r_cm_m"] * 1e6, style, label=case)
    axs[1, 2].set_xlabel("time (ps)")
    axs[1, 2].set_ylabel("electron centroid r_cm (um)")
    axs[1, 2].set_title("transverse centroid")
    axs[1, 2].legend(frameon=False)
    label(axs[1, 2], "(g)")
    metrics = ["Emax_max_Vpm", "ne_max_final_m3", "total_electrons_final"]
    labels = ["Emax", "ne,max", "N_e"]
    vals_aligned = [e3["cases"]["aligned"][m] for m in metrics]
    vals_mis = [e3["cases"]["misaligned"][m] for m in metrics]
    axs[1, 3].bar(np.arange(3) - width / 2, [1.0, 1.0, 1.0], width, label="aligned", color="#4c78a8")
    axs[1, 3].bar(np.arange(3) + width / 2, [b / a for a, b in zip(vals_aligned, vals_mis)], width, label="misaligned/aligned", color="#e45756")
    axs[1, 3].set_xticks(np.arange(3), labels)
    axs[1, 3].set_ylabel("normalized value")
    axs[1, 3].set_ylim(0.98, 1.02)
    axs[1, 3].set_title("bulk intensity comparison")
    axs[1, 3].legend(frameon=False)
    label(axs[1, 3], "(h)")
    files = save_figure(fig, "Fig_E4_misaligned_needle")
    return files, {"e3": e3, "comparison": comparison.to_dict(orient="records")}


def fig_e5() -> tuple[dict[str, str], dict[str, Any]]:
    d3_metrics = load_json(D3 / "results/stage_d3_metrics.json")
    e1 = load_json(SUMMARY / "e1_electrostatic_summary.json")
    e2 = load_json(SUMMARY / "e2_dynamic_summary.json")
    e3 = load_json(SUMMARY / "e3_summary.json")
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.5), constrained_layout=True)
    names = ["Emax", "ne,max", "N_e"]
    vals = [
        d3_metrics["max_Emax_rel_err"] * 100,
        d3_metrics["max_ne_max_rel_err"] * 100,
        d3_metrics["max_total_electrons_rel_err"] * 100,
    ]
    axs[0].bar(names, vals, color=["#4c78a8", "#f58518", "#54a24b"])
    axs[0].set_ylabel("D3 max relative difference (%)")
    axs[0].set_title("axisymmetric 2D/3D validation")
    axs[0].text(
        0.02,
        0.86,
        f"head z max diff = {d3_metrics['max_head_position_diff_m']*1e6:.2f} um\nAfivo symmetry metric = {d3_metrics['afivo_max_symmetry_metric']:.2e}",
        transform=axs[0].transAxes,
        va="top",
        fontsize=7,
    )
    label(axs[0], "(a)")
    cases = ["D3 aligned", "E2 triangular", "E3 aligned", "E3 misaligned"]
    field_asym = [
        d3_metrics["afivo_max_symmetry_metric"],
        e2["electrostatic_refined"]["refined_asymmetry_metric"],
        e3["cases"]["aligned"]["field_plane_asymmetry"],
        e3["cases"]["misaligned"]["field_plane_asymmetry"],
    ]
    plasma_asym = [
        d3_metrics["afivo_max_symmetry_metric"],
        max(
            e2["dynamic_cases"]["case_a_pi_off_500V"]["ne_plane_asymmetry"],
            e2["dynamic_cases"]["case_a_pi_off_500V"]["rho_plane_asymmetry"],
        ),
        max(e3["cases"]["aligned"]["ne_plane_asymmetry"], e3["cases"]["aligned"]["rho_plane_asymmetry"]),
        max(e3["cases"]["misaligned"]["ne_plane_asymmetry"], e3["cases"]["misaligned"]["rho_plane_asymmetry"]),
    ]
    colors = ["#4c78a8", "#f58518", "#54a24b", "#e45756"]
    for name, x, y, color in zip(cases, field_asym, plasma_asym, colors):
        axs[1].scatter(x, y, s=55, label=name, color=color)
    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    axs[1].set_xlabel("field asymmetry metric")
    axs[1].set_ylabel("plasma asymmetry metric")
    axs[1].set_title("evidence map, not a threshold")
    axs[1].legend(frameon=False, loc="best")
    label(axs[1], "(b)")
    files = save_figure(fig, "Fig_E5_2d_3d_validity_summary")
    return files, {"d3": d3_metrics, "field_asym": field_asym, "plasma_asym": plasma_asym}


def html_section(
    fig_id: str,
    title: str,
    caption: str,
    results: str,
    physics: str,
    supports: list[str],
    cannot_claim: list[str],
    role: str,
    role_reason: str,
    fig_file: str,
) -> str:
    bullet_supports = "\n".join(f"<li>{html.escape(x)}</li>" for x in supports)
    bullet_cannot = "\n".join(f"<li>{html.escape(x)}</li>" for x in cannot_claim)
    return f"""
<section class="figure-block">
<h2>{html.escape(fig_id)}</h2>
<figure>
<img src="{html.escape(fig_file)}" alt="{html.escape(fig_id)}">
</figure>
<h3>图号与拟投稿标题</h3>
<p><strong>{html.escape(title)}</strong></p>
<h3>中文图注</h3>
<p>{caption}</p>
<h3>论文式结果描述</h3>
{results}
<h3>物理解释</h3>
<p>{physics}</p>
<h3>可支撑的结论</h3>
<ul>{bullet_supports}</ul>
<h3>当前不能得出的结论</h3>
<ul>{bullet_cannot}</ul>
<h3>对未来论文的作用</h3>
<p><strong>{html.escape(role)}</strong>。{html.escape(role_reason)}</p>
</section>
"""


def update_manifests(figures: dict[str, dict[str, Any]], captions: dict[str, str]) -> None:
    manifest_path = REPORT / "figure_manifest.yaml"
    existing = yaml.safe_load(manifest_path.read_text())
    for fig in existing["figures"]:
        fid = fig["id"]
        if fid in figures:
            fig["final_files"] = figures[fid]["files"]
            fig["caption_zh"] = captions[fid]
            fig["data_hash_metadata"] = figures[fid]["files"]
            fig["report_html"] = "stage_e_smoke_paper_report.html"
    manifest_path.write_text(yaml.safe_dump(existing, allow_unicode=True, sort_keys=False), encoding="utf-8")

    data_path = REPORT / "data_manifest.csv"
    rows = list(csv.DictReader(data_path.open()))
    existing_keys = {(r["stage"], r["figure"], r["panel"], r["source_path"]) for r in rows}
    for fid, item in figures.items():
        for fmt in ("png", "pdf"):
            rel = item["files"][fmt]
            row = {
                "stage": "StageClosure",
                "figure": fid,
                "panel": f"final_{fmt}",
                "source_path": rel,
                "description": f"Generated Stage E paper-style {fmt.upper()} figure",
                "variables": "see figure_manifest.yaml",
                "tracked_in_git": "yes",
                "notes": f"sha256={item['files'][fmt + '_sha256']}",
            }
            key = (row["stage"], row["figure"], row["panel"], row["source_path"])
            if key not in existing_keys:
                rows.append(row)
                existing_keys.add(key)
    html_row = {
        "stage": "StageClosure",
        "figure": "ALL",
        "panel": "html_report",
        "source_path": "stage_e_smoke_paper_report.html",
        "description": "Stage E 3D non-axisymmetric smoke-test and paper-figure assessment",
        "variables": "all Stage D3/E1/E2/E3 summary values",
        "tracked_in_git": "yes",
        "notes": "",
    }
    key = (html_row["stage"], html_row["figure"], html_row["panel"], html_row["source_path"])
    if key not in existing_keys:
        rows.append(html_row)
    with data_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_html(figures: dict[str, dict[str, Any]], captions: dict[str, str]) -> str:
    d3 = load_json(D3 / "results/stage_d3_metrics.json")
    e2 = load_json(SUMMARY / "e2_dynamic_summary.json")
    e3 = load_json(SUMMARY / "e3_summary.json")
    e2a = e2["dynamic_cases"]["case_a_pi_off_500V"]
    e2b = e2["dynamic_cases"]["case_b_pi_on_500V"]
    e3a = e3["cases"]["aligned"]
    e3m = e3["cases"]["misaligned"]
    sections = [
        html_section(
            "Fig. E1",
            "Validation of the 3D Afivo model against the axisymmetric PETSc model",
            captions["FIG_E1_solver_validation"],
            f"""
<p>在相同的轴对称物理设定下，Afivo true 3D Cartesian 计算给出的 Emax、峰值电子密度和总电子数均跟随 PETSc 2D 轴对称基准演化。最大相对差异分别为 {d3['max_Emax_rel_err']*100:.2f}%、{d3['max_ne_max_rel_err']*100:.2f}% 和 {d3['max_total_electrons_rel_err']*100:.2f}%，前沿位置最大差异为 {d3['max_head_position_diff_m']*1e6:.2f} um。</p>
<p>Afivo 解的横向质心偏移低于 {d3['afivo_max_cm_offset_m']:.2e} m，中心切面对称性指标为 {d3['afivo_max_symmetry_metric']:.2e}。由于 leading edge 基本 pinned，head velocity 的相对误差不作为主要判据。</p>
""",
            "该结果说明，当几何、电压、输运、化学和种子均保持轴对称时，3D Cartesian AMR 求解不会自发引入显著非轴对称响应；剩余差异主要来自电极离散和数值格式差异。",
            [
                "Afivo 3D backend 可作为真正 3D 非轴对称问题的可信求解器基础。",
                "PETSc 2D 与 Afivo 3D 在轴对称退化工况下具有工程一致性。",
            ],
            [
                "不能声称两个不同数值框架逐点场完全一致。",
                "不能用该图证明真实实验击穿或成熟火花阶段。",
            ],
            "方法验证图",
            "它是进入 3D 非轴对称结果之前的必要数值可信度证据，适合正文 Methods/validation 或补充材料主图。",
            figures["FIG_E1_solver_validation"]["files"]["png"],
        ),
        html_section(
            "Fig. E2",
            "Three-dimensional electrostatic field of the triangular-foil electrode",
            captions["FIG_E2_triangular_electrostatic"],
            """
<p>三角铜箔 development geometry 在 x-z 与 y-z 中心截面给出明显不同的电场分布。局部加密后，电场非对称性指标由 medium 的 0.1586 变为 refined 的 0.1031，非轴对称结论保持稳定。</p>
<p>同时，局部峰值场从 2.827e7 V/m 增至 5.734e7 V/m，变化约 102.8%。因此图中可展示尖端/边缘的高场集中位置，但不能把 57.3 MV/m 解释为已经收敛的实验物理峰值场。</p>
""",
            "有限宽度、有限厚度和有限边缘圆角破坏了任何绕 z 轴的旋转对称性，使高场区域同时受尖端曲率和侧边缘影响；这是 2D 轴对称投影无法保持的几何信息。",
            [
                "三角箔电极的局部电场是本质 3D 问题。",
                "高场区域稳定出现在 rounded tip/edge 附近。",
            ],
            [
                "不能声称 Emax=57.3 MV/m 是真实实验峰值场。",
                "不能声称 tip/edge Emax 已经网格收敛。",
            ],
            "重要支撑图",
            "它支撑 Stage E 进入 3D 的几何必要性；正式投稿前应使用 Stage B 校准几何重跑。",
            figures["FIG_E2_triangular_electrostatic"]["files"]["png"],
        ),
        html_section(
            "Fig. E3",
            "Transfer of geometric asymmetry to early discharge dynamics",
            captions["FIG_E3_triangular_dynamic"],
            f"""
<p>在三角箔 development geometry 中，PI OFF case 的总电子数由 {e2a['initial_total_electrons']:.6g} 增至 {e2a['final_total_electrons']:.7g}，PI ON case 由 {e2b['initial_total_electrons']:.6g} 增至 {e2b['final_total_electrons']:.7g}。两者均表现出早期 avalanche 增长，但未出现 streamer propagation 或 bridge。</p>
<p>5 ps 时，electron-density plane asymmetry 约为 {e2a['ne_plane_asymmetry']:.4f}，charge-density plane asymmetry 约为 {e2a['rho_plane_asymmetry']:.4f}。PI ON/OFF 在这个 5 ps development case 中差异很小，因此该图不把 photoionization 作为主结论。</p>
""",
            "初始种子仍位于尖端附近，但三角箔电场在两个中心截面上不同，电子漂移和反应源在早期即采样到这种几何非对称性，使 avalanche 的 ne/rho 分布继承电极的 3D 特征。",
            [
                "early discharge/avalanche inherits the 3D geometric asymmetry。",
                "rho、ne、E、J 的 3D source export pipeline 可用于后续 Stage F 数据输入。",
            ],
            [
                "不能写 triangular foil streamer has been demonstrated。",
                "不能写 photoionization generally negligible，只能说本 5 ps development case 中差异小。",
            ],
            "Supplementary候选",
            "它证明动态和 source export 管线工作，但当前只是早期 avalanche，尚不足以作为主物理图。",
            figures["FIG_E3_triangular_dynamic"]["files"]["png"],
        ),
        html_section(
            "Fig. E4",
            "Controlled symmetry breaking by lateral needle misalignment",
            captions["FIG_E4_misaligned_needle"],
            f"""
<p>双针 controlled benchmark 中，aligned case 的 field asymmetry 为 {e3a['field_plane_asymmetry']:.2e}，misaligned case 增至 {e3m['field_plane_asymmetry']:.4f}。对应的电子密度非对称性从 {e3a['ne_plane_asymmetry']:.2e} 增至 {e3m['ne_plane_asymmetry']:.4e}，电荷密度非对称性从 {e3a['rho_plane_asymmetry']:.4e} 增至 {e3m['rho_plane_asymmetry']:.4e}。</p>
<p>与此同时，两组整体放电强度量基本不变：aligned/misaligned 的 Emax 分别为 {e3a['Emax_max_Vpm']:.4e} 和 {e3m['Emax_max_Vpm']:.4e} V/m，ne,max 分别为 {e3a['ne_max_final_m3']:.4e} 和 {e3m['ne_max_final_m3']:.4e} m^-3，总电子数分别为 {e3a['total_electrons_final']:.4f} 和 {e3m['total_electrons_final']:.4f}。电子横向质心最大偏移从 {e3a['max_r_cm_m']:.2e} m 增至 {e3m['max_r_cm_m']:.2e} m。</p>
""",
            "在相同电压、输运、化学、seed 和分辨率下，10 um lateral offset 改变的是边界几何本身，而不是总体电离强度。整体 Emax、ne,max 和 Ne 接近，说明观测到的空间非对称性主要由 geometry symmetry breaking 触发。",
            [
                "真实几何错位能够在 3D Afivo 中产生强的场非对称性和可测的早期 plasma asymmetry。",
                "该响应不是 Cartesian numerical noise，因为 aligned control 保持近对称。",
            ],
            [
                "不能写 10 um offset causes streamer deflection。",
                "不能声称已经产生 bridge 或完整 streamer crossing。",
            ],
            "核心主图",
            "这是 Stage E 最强的 controlled 证据，可作为论文中说明 3D 几何必要性的机制支撑主图。",
            figures["FIG_E4_misaligned_needle"]["files"]["png"],
        ),
        html_section(
            "Fig. E5",
            "Evidence map for the applicability of axisymmetric and three-dimensional models",
            captions["FIG_E5_2d_3d_validity_summary"],
            """
<p>证据图分为两部分：一部分展示 D3 轴对称退化验证误差，另一部分比较 Stage E 的非轴对称场/等离子体响应指标。该图不定义 universal 2D/3D threshold，而是把不同 benchmark 的证据边界清楚分开。</p>
<p>D3 表明轴对称 geometry 可由 PETSc 2D 高效表示；E2/E3 表明三角箔和错位双针会产生 2D 轴对称模型无法表达的空间非对称性。</p>
""",
            "轴对称模型的有效性取决于几何、边界和初始条件是否保持旋转对称。一旦电极具有有限厚度/边缘或 deliberate lateral offset，空间电场和源项不再能由单一 r-z 截面代表。",
            [
                "axisymmetric geometries can be efficiently represented by PETSc 2D。",
                "intrinsically or deliberately symmetry-broken geometries require 3D resolution。",
            ],
            [
                "不能声称建立了 asymmetry > X 的通用阈值。",
                "不能用当前 smoke 数据替代 Stage B/实验校准后的生产仿真。",
            ],
            "重要支撑图",
            "它适合作为 Discussion/summary 图，帮助读者理解 2D 与 3D backend 的职责边界。",
            figures["FIG_E5_2d_3d_validity_summary"]["files"]["png"],
        ),
    ]
    summary_rows = [
        ("Fig E1", "3D solver 是否可信", "Strong", "early-time only", "Method validation", "No, production refinement only"),
        ("Fig E2", "triangular foil 是否真正 3D", "Strong qualitative", "tip Emax not converged", "Important support", "Yes, Stage B calibrated geometry"),
        ("Fig E3", "几何非对称是否进入 early dynamics", "Moderate", "avalanche only, no streamer", "Supplementary candidate", "Yes, if used as main physics"),
        ("Fig E4", "misalignment 是否产生 controlled symmetry breaking", "Strong", "early-time only, no bridge", "Mechanism-supporting main figure", "Yes, for final experimental geometry"),
        ("Fig E5", "2D/3D 适用性边界如何表达", "Strong as synthesis", "not a universal threshold", "Important support", "No for method summary; update values after production"),
    ]
    table = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in row) + "</tr>"
        for row in summary_rows
    )
    forbidden = [
        "triangular foil streamer has been demonstrated",
        "10 um offset causes streamer deflection",
        "Emax=57.3 MV/m is the physical electrode field",
        "photoionization has negligible effect generally",
        "a universal 2D/3D threshold was established",
    ]
    forbidden_html = "\n".join(f"<li>{html.escape(x)}</li>" for x in forbidden)
    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>Stage E 3D Non-axisymmetric Simulation Smoke-Test and Paper-Figure Assessment</title>
<style>
body {{ font-family: Arial, "Noto Sans CJK SC", sans-serif; margin: 32px auto; max-width: 1100px; line-height: 1.55; color: #222; }}
h1, h2, h3 {{ line-height: 1.25; }}
h1 {{ font-size: 24px; }}
h2 {{ margin-top: 36px; border-top: 1px solid #ddd; padding-top: 20px; }}
h3 {{ font-size: 16px; margin-bottom: 4px; }}
img {{ max-width: 100%; border: 1px solid #ddd; }}
figure {{ margin: 16px 0; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
td, th {{ border: 1px solid #ccc; padding: 6px; vertical-align: top; }}
th {{ background: #f3f3f3; }}
.note {{ background: #f7f7f7; border-left: 4px solid #888; padding: 10px 14px; }}
</style>
</head>
<body>
<h1>Stage E 3D Non-axisymmetric Simulation Smoke-Test and Paper-Figure Assessment</h1>
<p>本报告汇总 Stage D3、Stage E1、Stage E2 和 Stage E3 已完成的真实 smoke-test 数据，用于评估未来投稿论文中 3D 仿真部分的候选图。Stage D 的作用是验证 Afivo true 3D Cartesian solver 在轴对称条件下能够退化到 PETSc 2D 轴对称结果；Stage E 的作用是验证真实 non-axisymmetric geometry 会导致 3D field/plasma response。</p>
<p class="note">本报告属于 development/smoke-level validation，不是 final calibrated experimental simulation。所有三角箔和双针尺寸均为 development geometry；正式实验定量图应在 Stage B 几何和 measured voltage waveform 完成后重跑。</p>
{''.join(sections)}
<h2>总结表</h2>
<table>
<thead><tr><th>Figure</th><th>Scientific question</th><th>Evidence strength</th><th>Current limitation</th><th>Recommended paper role</th><th>Need rerun later?</th></tr></thead>
<tbody>{table}</tbody>
</table>
<h2>总体论文贡献判断</h2>
<p><strong>A. Stage D/E 是否值得进入论文正文？</strong> PARTIAL。D/E 值得进入论文，但不应被包装成主科学创新本身。它们更适合作为 numerical validation 与 mechanism-supporting evidence：Fig E1 证明 3D backend 可信，Fig E4 证明受控几何破缺确实产生 3D 响应。</p>
<p><strong>B. 推荐角色。</strong> Fig E4 可承担 mechanism-supporting evidence 并作为正文核心支撑图；Fig E1 主要属于 numerical validation；Fig E2 是 3D 几何必要性的支撑图；Fig E3 更适合作为 Supplementary 候选；Fig E5 可作为 discussion summary。</p>
<p><strong>C. 当前最强结果。</strong> Fig E4：仅改变 10 um lateral offset，field asymmetry 从近零增至 0.3719，而 Emax、ne,max 和 Ne 基本保持一致，这是 controlled symmetry-breaking 的最干净证据。</p>
<p><strong>D. 当前最弱结果。</strong> Fig E3：只显示 early avalanche 继承三角箔非对称性，没有 streamer propagation 或 bridge，PI ON/OFF 差异也很小。</p>
<p><strong>E. Stage B 和实验数据完成后需要正式重跑。</strong> Fig E2 和 Fig E4 的几何/电压应以 Stage B calibrated geometry 和 measured voltage waveform 重跑；若 Fig E3 要作为正文动态物理图，也必须重跑。</p>
<p><strong>F. 可保留为方法验证的现有 smoke 结果。</strong> Fig E1 作为 2D/3D backend validation 可保留；Fig E5 的方法边界逻辑可保留但数值可随生产 case 更新。</p>
<h2>禁止过度表述</h2>
<ul>{forbidden_html}</ul>
</body>
</html>
"""
    path = REPORT / "stage_e_smoke_paper_report.html"
    path.write_text(html_text, encoding="utf-8")
    return html_text


def validate_outputs(figures: dict[str, dict[str, Any]]) -> None:
    for item in figures.values():
        for fmt in ("png", "pdf"):
            path = REPORT / item["files"][fmt]
            if not path.exists() or path.stat().st_size <= 0:
                raise RuntimeError(f"missing figure output: {path}")
    report = REPORT / "stage_e_smoke_paper_report.html"
    html_text = report.read_text(encoding="utf-8")
    for item in figures.values():
        if item["files"]["png"] not in html_text:
            raise RuntimeError(f"HTML missing image link {item['files']['png']}")
    if "nan" in html_text.lower():
        raise RuntimeError("HTML contains NaN")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    generated: dict[str, dict[str, Any]] = {}
    for fid, func in [
        ("FIG_E1_solver_validation", fig_e1),
        ("FIG_E2_triangular_electrostatic", fig_e2),
        ("FIG_E3_triangular_dynamic", fig_e3),
        ("FIG_E4_misaligned_needle", fig_e4),
        ("FIG_E5_2d_3d_validity_summary", fig_e5),
    ]:
        files, meta = func()
        generated[fid] = {"files": files, "meta": meta}

    captions = {
        "FIG_E1_solver_validation": "图 E1. Afivo true 3D Cartesian 模型相对于 PETSc 2D 轴对称模型的退化验证。(a) PETSc 轴对称截面示意；(b) Afivo 真实 3D Cartesian 对齐电极示意；(c-e) Emax、ne,max 与总电子数随时间的对比；(f) 稳健误差指标汇总。最大差异为 Emax 5.94%、ne,max 12.99%、Ne 6.04%，head position 最大差异 0.66 um。",
        "FIG_E2_triangular_electrostatic": "图 E2. 三角铜箔 development geometry 的三维静电场。(a) 有限圆角三角箔示意；(b,c) refined x-z 与 y-z 中心截面 |E|；(d) 两个截面的 selected field profile；(e) medium 与 refined 非对称性指标；(f) 高场体积。refined asymmetry metric 为 0.1031，medium 为 0.1586；Emax 从 2.827e7 增至 5.734e7 V/m，局部峰值场尚未网格收敛。",
        "FIG_E3_triangular_dynamic": "图 E3. 三角箔几何非对称性向早期放电动力学的传递。(a) 初始电子密度切片；(b,c) 5 ps 时 x-z 与 y-z 电子密度切片；(d) 5 ps 电荷密度切片；(e) PI OFF/ON 总电子数演化；(f) ne/rho 非对称性。PI OFF 总电子数 4.249896→5.1771357，PI ON 为 4.249896→5.1776558；Avalanche=YES，Streamer propagation=NO，Bridge=NO。",
        "FIG_E4_misaligned_needle": "图 E4. 双针 lateral misalignment 导致的受控对称性破缺。(a,b) aligned 与 10 um misaligned 双针几何；(c,d) 两组中心切面场 profile；(e) field asymmetry；(f) ne/rho asymmetry；(g) 电子横向质心；(h) Emax、ne,max、Ne 归一化对比。field asymmetry 从 9.59e-10 增至 0.3719，ne asymmetry 从 2.80e-5 增至 3.52e-3，rho asymmetry 从 4.10e-3 增至 1.154e-2；整体强度指标保持接近。",
        "FIG_E5_2d_3d_validity_summary": "图 E5. 轴对称模型与三维模型适用性的证据图。(a) D3 轴对称 PETSc/Afivo 动态验证误差；(b) D3/E2/E3 的场非对称性与等离子体非对称性证据图。该图不定义通用阈值，而用于说明轴对称几何可由 2D 表示，内禀或受控破缺几何需要 3D 分辨。",
    }
    build_html(generated, captions)
    update_manifests(generated, captions)
    validate_outputs(generated)
    summary = {
        fid: {
            "png": item["files"]["png"],
            "pdf": item["files"]["pdf"],
            "png_sha256": item["files"]["png_sha256"],
            "pdf_sha256": item["files"]["pdf_sha256"],
        }
        for fid, item in generated.items()
    }
    (REPORT / "stage_e_smoke_paper_report_assets.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
