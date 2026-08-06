#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import SymLogNorm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/stage3/figures/current"


def _grid(field: pd.DataFrame, quantity: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = np.sort(field.z_m.unique())
    r = np.sort(field.r_m.unique())
    q = field.pivot(index="r_m", columns="z_m", values=quantity).loc[r, z].to_numpy()
    return z * 1e3, r * 1e3, q


def plot_field_triplet(run: str, field_file: str, title: str, output_name: str) -> None:
    f = pd.read_csv(ROOT / field_file)
    z, r, ne = _grid(f, "ne_m_3")
    _, _, e = _grid(f, "E_V_m")
    _, _, rho = _grid(f, "rho_C_m_3")
    ne_log = np.log10(np.maximum(ne, 1.0))

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
    panels = [
        (ne_log, r"$\log_{10}(n_e\,[m^{-3}])$", "viridis", None),
        (e / 1e6, r"$|E|$ [MV/m]", "magma", None),
        (rho, r"$\rho$ [C/m$^3$]", "coolwarm", SymLogNorm(linthresh=1e-3, vmin=-np.nanmax(abs(rho)), vmax=np.nanmax(abs(rho)))),
    ]
    for ax, (q, label, cmap, norm) in zip(axes, panels):
        im = ax.pcolormesh(z, r, q, shading="auto", cmap=cmap, norm=norm)
        ax.set_title(label)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(f"{title} ({run})")
    fig.savefig(OUT / output_name, dpi=180)
    plt.close(fig)


def plot_axis_profiles(run: str, field_file: str, title: str, output_name: str) -> None:
    f = pd.read_csv(ROOT / field_file)
    axis = f[f.i == f.i.min()].sort_values("z_m")
    z = axis.z_m.to_numpy() * 1e3
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    axes[0].plot(z, axis.E_V_m / 1e6)
    axes[0].set_ylabel("|E| [MV/m]")
    axes[1].semilogy(z, np.maximum(axis.ne_m_3, 1.0), label="ne")
    axes[1].semilogy(z, np.maximum(axis.np_m_3, 1.0), label="np", alpha=0.75)
    axes[1].semilogy(z, np.maximum(axis.nn_m_3, 1.0), label="nn", alpha=0.75)
    axes[1].set_ylabel("density [m$^{-3}$]")
    axes[1].legend(loc="best")
    axes[2].plot(z, axis.rho_C_m_3)
    axes[2].axhline(0.0, color="k", lw=0.8)
    axes[2].set_ylabel(r"$\rho$ [C/m$^3$]")
    axes[2].set_xlabel("z [mm]")
    fig.suptitle(f"{title} axis profiles ({run})")
    fig.savefig(OUT / output_name, dpi=180)
    plt.close(fig)


def plot_head_trajectory() -> None:
    traj = pd.read_csv(ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/head_trajectory.csv")
    t = traj.time_s.to_numpy() * 1e9
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(t, traj.lower_z_rho_m * 1e3, "o-", label="lower z_rho")
    ax.plot(t, traj.lower_z_E_m * 1e3, "o--", label="lower z_E", alpha=0.75)
    ax.plot(t, traj.lower_z_gradient_m * 1e3, "o:", label="lower z_grad", alpha=0.75)
    ax.plot(t, traj.upper_z_rho_m * 1e3, "s-", label="upper z_rho")
    ax.plot(t, traj.upper_z_E_m * 1e3, "s--", label="upper z_E", alpha=0.75)
    ax.plot(t, traj.upper_z_gradient_m * 1e3, "s:", label="upper z_grad", alpha=0.75)
    ax.axhline(5.0, color="k", lw=0.8, alpha=0.5)
    ax.set_xlabel("time [ns]")
    ax.set_ylabel("head position z [mm]")
    ax.set_title("20 um connected-front head trajectory")
    ax.legend(ncol=2, fontsize=8)
    fig.savefig(OUT / "coarse_ml_head_trajectory.png", dpi=180)
    plt.close(fig)


def plot_sp3_comparison() -> None:
    on = pd.read_csv(ROOT / "results/stage3/single_seed/coarse_ml_recovery_resume1/scalar_history.csv")
    off = pd.read_csv(ROOT / "results/stage3/single_seed/sp3_fork_off_resume1/scalar_history.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    axes[0].plot(on.time_s * 1e9, on.E_max_V_m / 1e6, label="SP3 on")
    axes[0].plot(off.time_s * 1e9, off.E_max_V_m / 1e6, label="SP3 off")
    axes[0].set_ylabel("Emax [MV/m]")
    axes[1].plot(on.time_s * 1e9, on.total_electrons, label="SP3 on")
    axes[1].plot(off.time_s * 1e9, off.total_electrons, label="SP3 off")
    axes[1].set_ylabel("total electrons")
    axes[2].semilogy(on.time_s * 1e9, np.maximum(on.sph_ahead_m_3_s_1, 1.0), label="SP3 on")
    axes[2].semilogy(off.time_s * 1e9, np.maximum(off.sph_ahead_m_3_s_1, 1.0), label="SP3 off")
    axes[2].set_ylabel(r"$S_{ph}$ ahead [m$^{-3}$ s$^{-1}$]")
    for ax in axes:
        ax.set_xlabel("time [ns]")
        ax.legend()
    fig.suptitle("SP3 on/off measured response")
    fig.savefig(OUT / "sp3_on_off_comparison.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plot_field_triplet(
        "coarse_ml_recovery_resume1",
        "results/stage3/single_seed/coarse_ml_recovery_resume1/fields_final.csv",
        "20 um ML final field at 2.0 ns",
        "coarse_ml_final_fields.png",
    )
    plot_axis_profiles(
        "coarse_ml_recovery_resume1",
        "results/stage3/single_seed/coarse_ml_recovery_resume1/fields_final.csv",
        "20 um ML final field at 2.0 ns",
        "coarse_ml_final_axis_profiles.png",
    )
    plot_field_triplet(
        "stage3_verification_baseline_10um_resume1",
        "results/stage3/single_seed/stage3_verification_baseline_10um_resume1/fields_343.csv",
        "10 um partial field at 0.400 ns",
        "partial_10um_0p400ns_fields.png",
    )
    plot_head_trajectory()
    plot_sp3_comparison()
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
