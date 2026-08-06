#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from streamer_rf.stage4 import (
    align_current_moments,
    collision_pulse_metrics,
    derivative_products,
    detect_collision_event,
    event_to_frame,
    radiation_products,
    output_sampling_sensitivity,
    read_csv_clean,
    sha256,
)


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()


def read_field(path: Path) -> pd.DataFrame:
    if path.suffix == ".h5":
        with h5py.File(path, "r") as h:
            return pd.DataFrame({k: h[k][()] for k in h.keys()})
    return pd.read_csv(path)


def closest_field(run_dir: Path, target_time: float, scalar: pd.DataFrame) -> Path:
    fields = sorted(run_dir.glob("fields_*.h5"))
    if not fields:
        fields = sorted(run_dir.glob("fields_*.csv"))
    if target_time <= 0:
        return run_dir / ("fields_0.h5" if (run_dir / "fields_0.h5").exists() else "fields_0.csv")
    best = None
    best_dt = float("inf")
    for p in fields:
        name = p.stem
        if name == "fields_final":
            t = float(scalar.time_s.iloc[-1])
        else:
            try:
                step = int(name.split("_")[1])
            except Exception:
                continue
            hit = scalar[scalar.step == step]
            if len(hit) == 0:
                continue
            t = float(hit.time_s.iloc[0])
        dt = abs(t - target_time)
        if dt < best_dt:
            best = p
            best_dt = dt
    if best is None:
        raise FileNotFoundError(f"no field output matched t={target_time}")
    return best


def plot_field(run_dir: Path, scalar: pd.DataFrame, target_time: float, quantity: str, out: Path, title: str) -> None:
    p = closest_field(run_dir, target_time, scalar)
    d = read_field(p)
    r = np.sort(d.r_m.unique())
    z = np.sort(d.z_m.unique())
    grid = d.pivot(index="z_m", columns="r_m", values=quantity).loc[z, r].to_numpy()
    plt.figure(figsize=(7, 4))
    plt.imshow(grid, origin="lower", aspect="auto", extent=[r[0] * 1e3, r[-1] * 1e3, z[0] * 1e3, z[-1] * 1e3])
    plt.colorbar(label=quantity)
    plt.xlabel("r (mm)")
    plt.ylabel("z (mm)")
    plt.title(f"{title}\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(out)


def write_with_hash(df: pd.DataFrame, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return sha256(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/stage4")
    args = ap.parse_args()
    root = ROOT / args.root

    collision_run = root / "runs/paperlike_20um"
    left_run = root / "runs/left_isolated_20um"
    right_run = root / "runs/right_isolated_20um"

    coll_metrics = read_csv_clean(collision_run / "collision_metrics.csv")
    coll_scalar = read_csv_clean(collision_run / "scalar_history.csv")
    event = detect_collision_event(coll_metrics)

    collision_dir = root / "collision"
    event_hash = write_with_hash(event_to_frame(event), collision_dir / "collision_event.csv")

    c = read_csv_clean(root / "current_moment/collision.csv")
    l = read_csv_clean(root / "current_moment/left_isolated.csv")
    r = read_csv_clean(root / "current_moment/right_isolated.csv")
    delta = align_current_moments(c, l, r, t_max=event.t_collision_s + 0.25e-9)
    pulse = collision_pulse_metrics(delta, event.t_collision_s)
    delta["I_collision"] = delta.I_collision_I_CM_drift
    delta["I_left"] = delta.I_left_I_CM_drift
    delta["I_right"] = delta.I_right_I_CM_drift
    delta["Delta_I"] = delta.Delta_I_CM_drift
    delta_out = delta[["time", "I_collision", "I_left", "I_right", "Delta_I", "interpolation_flag", "local_max_original_dt",
                       "I_collision_I_CM_electron", "I_left_I_CM_electron", "I_right_I_CM_electron", "Delta_I_CM_electron"]]
    write_with_hash(delta_out, root / "current_moment/delta_current_moment.csv")
    write_with_hash(pd.DataFrame([pulse]), root / "current_moment/delta_current_pulse_metrics.csv")

    deriv, deriv_summary = derivative_products(delta.rename(columns={"Delta_I_CM_drift": "Delta_I_CM_drift"}), event.t_collision_s)
    write_with_hash(deriv, root / "radiation/delta_current_derivative.csv")
    write_with_hash(pd.DataFrame([deriv_summary]), root / "radiation/derivative_sensitivity.csv")

    spec, esd, bands, rad_summary = radiation_products(delta, event.t_collision_s, window="rectangular")
    spec_tukey, esd_tukey, _, rad_summary_tukey = radiation_products(delta, event.t_collision_s, window="tukey")
    write_with_hash(spec, root / "radiation/delta_current_spectrum.csv")
    write_with_hash(esd, root / "radiation/delta_current_esd.csv")
    write_with_hash(bands, root / "radiation/band_energy.csv")
    ws = pd.DataFrame([
        {"window": "rectangular", **rad_summary, "peak_spectrum": float(spec.abs_dDeltaI_dt_transform.max())},
        {"window": "tukey", **rad_summary_tukey, "peak_spectrum": float(spec_tukey.abs_dDeltaI_dt_transform.max())},
    ])
    write_with_hash(ws, root / "radiation/window_sensitivity.csv")
    sampling = output_sampling_sensitivity(delta, event.t_collision_s)
    write_with_hash(sampling, root / "radiation/output_sampling_sensitivity.csv")

    figs = root / "figures"
    peak_time = event.E_gap_peak_time_s
    bridge_time = event.t_bridge_s
    after_time = event.E_gap_min_after_time_s
    plot_field(collision_run, coll_scalar, 0.0, "ne_m_3", figs / "01_double_seed_initial_density.png", "Double-seed initial electron density")
    plot_field(collision_run, coll_scalar, max(event.t_collision_s - 0.5e-9, 0.0), "ne_m_3", figs / "02_opposing_propagation_density.png", "Opposing streamer propagation")
    plot_field(collision_run, coll_scalar, peak_time, "E_V_m", figs / "03_precollision_field_enhancement.png", "Pre-collision field enhancement")
    plot_field(collision_run, coll_scalar, bridge_time, "ne_m_3", figs / "04_bridge_channel_density.png", "Bridge channel formation")
    plot_field(collision_run, coll_scalar, after_time, "E_V_m", figs / "05_postcollision_field_drop.png", "Post-collision field decrease")

    plt.figure(figsize=(7, 4))
    plt.plot(coll_metrics.time * 1e9, coll_metrics.d_head * 1e3, label="d_head (diagnostic)")
    plt.plot(coll_metrics.time * 1e9, coll_metrics.bridge_mean_ne / event.bridge_threshold_m_3, label="bridge / threshold")
    plt.plot(coll_metrics.time * 1e9, coll_metrics.E_gap_max / event.E_gap_peak_V_m, label="E_gap / peak")
    plt.axvline(event.t_collision_s * 1e9, color="k", ls="--", label="derived collision time")
    plt.xlabel("time (ns)")
    plt.legend()
    plt.title("Collision metrics\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "06_collision_metrics.png")

    plt.figure(figsize=(7, 4))
    plt.plot(c.time * 1e9, c.I_CM_drift, label="collision")
    plt.plot(l.time * 1e9, l.I_CM_drift, label="left isolated")
    plt.plot(r.time * 1e9, r.I_CM_drift, label="right isolated")
    plt.xlabel("time (ns)")
    plt.ylabel("I_CM drift (A m)")
    plt.legend()
    plt.title("Measured C++ drift current moment\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "07_current_moments.png")

    plt.figure(figsize=(7, 4))
    plt.plot(delta_out.time * 1e9, delta_out.Delta_I)
    plt.axvline(event.t_collision_s * 1e9, color="k", ls="--")
    plt.xlabel("time (ns)")
    plt.ylabel("Delta I_CM drift (A m)")
    plt.title("Derived collision additional current moment\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "08_delta_current_moment.png")

    plt.figure(figsize=(7, 4))
    plt.plot(deriv.time * 1e9, deriv.dDeltaI_dt_local_poly, label="local poly")
    plt.plot(deriv.time * 1e9, deriv.dDeltaI_dt_uniform_five_point, label="uniform five-point", alpha=0.75)
    plt.xlabel("time (ns)")
    plt.ylabel("d Delta I / dt (A m s^-1)")
    plt.legend()
    plt.title("Current-moment derivative methods\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "09_delta_current_derivative.png")

    plt.figure(figsize=(7, 4))
    plt.loglog(spec.frequency_Hz[1:], spec.abs_dDeltaI_dt_transform[1:])
    plt.axvline(rad_summary["f_trust_Hz"], color="k", ls="--", label="trust limit")
    plt.xlabel("frequency (Hz)")
    plt.ylabel("|FFT(dDeltaI/dt)|")
    plt.legend()
    plt.title("FFT amplitude with conservative trust limit\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "10_fft_amplitude.png")

    plt.figure(figsize=(7, 4))
    plt.loglog(esd.frequency_Hz[1:], esd.ESD_J_per_Hz[1:])
    plt.axvline(rad_summary["f_trust_Hz"], color="k", ls="--", label="trust limit")
    plt.xlabel("frequency (Hz)")
    plt.ylabel("ESD (J Hz^-1)")
    plt.legend()
    plt.title("Energy spectral density\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "11_esd.png")

    plt.figure(figsize=(7, 4))
    plt.loglog(spec.frequency_Hz[1:], spec.abs_dDeltaI_dt_transform[1:], label="rectangular")
    plt.loglog(spec_tukey.frequency_Hz[1:], spec_tukey.abs_dDeltaI_dt_transform[1:], label="Tukey", alpha=0.8)
    plt.axvline(rad_summary["f_trust_Hz"], color="k", ls="--", label="trust limit")
    plt.xlabel("frequency (Hz)")
    plt.ylabel("|FFT(dDeltaI/dt)|")
    plt.legend()
    plt.title("Window-function sensitivity\nResource-aware workflow reconstruction using the verified Morrow–Lowke-based solver.")
    savefig(figs / "12_window_sensitivity.png")

    manifest = []
    for p in sorted((root / "collision").glob("*.csv")) + sorted((root / "current_moment").glob("*.csv")) + sorted((root / "radiation").glob("*.csv")) + sorted(figs.glob("*.png")) + sorted(figs.glob("*.pdf")):
        manifest.append({"file": str(p.relative_to(ROOT)), "sha256": sha256(p)})
    write_with_hash(pd.DataFrame(manifest), root / "provenance/derived_artifact_manifest.csv")

    summary = {
        **event.__dict__,
        **pulse,
        **deriv_summary,
        **rad_summary,
        "event_csv_sha256": event_hash,
    }
    write_with_hash(pd.DataFrame([summary]), root / "stage4_analysis_summary.csv")
    print(pd.DataFrame([summary]).to_string(index=False))


if __name__ == "__main__":
    main()
