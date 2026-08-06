#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.stage4 import (
    align_current_moments,
    collision_pulse_metrics,
    derivative_products,
    detect_collision_event,
    event_to_frame,
    output_sampling_sensitivity,
    radiation_products,
    read_csv_clean,
)
from streamer_rf.stage5 import (
    band_metric,
    event_local_proxy_products,
    head_velocity_metrics,
    pulse_fwhm,
    robust_status_from_termination,
    spectral_centroid,
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ensure_dirs(root: Path) -> None:
    for rel in ["collision", "current_moment", "radiation", "trends", "mpi", "provenance"]:
        (root / rel).mkdir(parents=True, exist_ok=True)


def run_dir(root: Path, name: str) -> Path:
    return root / "runs" / name


def load_run(root: Path, name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series]:
    d = run_dir(root, name)
    return (
        read_csv_clean(d / "scalar_history.csv"),
        read_csv_clean(d / "collision_metrics.csv"),
        read_csv_clean(d / "current_moment.csv"),
        pd.read_csv(d / "termination.csv").iloc[0],
    )


def write_highfield(root: Path) -> dict[str, float]:
    _, m_c, cm_c, term_c = load_run(root, "highfield_collision_20um")
    _, _, cm_l, _ = load_run(root, "highfield_left_isolated_20um")
    _, _, cm_r, _ = load_run(root, "highfield_right_isolated_20um")
    ev = detect_collision_event(m_c)
    event_to_frame(ev).to_csv(root / "collision/highfield_collision_event.csv", index=False)
    cm_c.to_csv(root / "current_moment/highfield_collision.csv", index=False)
    cm_l.to_csv(root / "current_moment/highfield_left_isolated.csv", index=False)
    cm_r.to_csv(root / "current_moment/highfield_right_isolated.csv", index=False)

    delta = align_current_moments(cm_c, cm_l, cm_r, t_max=ev.t_collision_s + 0.20e-9)
    delta.to_csv(root / "current_moment/highfield_delta_current_moment.csv", index=False)
    pulse = collision_pulse_metrics(delta, ev.t_collision_s)
    pd.DataFrame([pulse]).to_csv(root / "current_moment/highfield_delta_current_pulse_metrics.csv", index=False)
    deriv, dsum = derivative_products(delta, ev.t_collision_s)
    deriv.to_csv(root / "radiation/highfield_delta_current_derivative.csv", index=False)
    pd.DataFrame([dsum]).to_csv(root / "radiation/highfield_derivative_sensitivity.csv", index=False)
    spec, esd, bands, rsum = radiation_products(delta, ev.t_collision_s)
    spec.to_csv(root / "radiation/highfield_delta_current_spectrum.csv", index=False)
    esd.to_csv(root / "radiation/highfield_delta_current_esd.csv", index=False)
    bands.to_csv(root / "radiation/highfield_band_energy.csv", index=False)
    pd.DataFrame([rsum]).to_csv(root / "radiation/highfield_frequency_trust.csv", index=False)
    sampling = stage5_sampling_sensitivity(delta, ev.t_collision_s)
    sampling.to_csv(root / "radiation/highfield_output_sampling_sensitivity.csv", index=False)
    spec_t, esd_t, bands_t, rsum_t = radiation_products(delta, ev.t_collision_s, window="tukey")
    pd.DataFrame([{
        "window": "tukey",
        "spectral_centroid_Hz": spectral_centroid(spec_t, esd_t),
        "f_trust_Hz": rsum_t["f_trust_Hz"],
        "VHF_energy": band_metric(bands_t, "VHF"),
        "UHF_energy": band_metric(bands_t, "UHF"),
        "SHF_energy": band_metric(bands_t, "SHF"),
    }]).to_csv(root / "radiation/highfield_window_sensitivity.csv", index=False)

    v_left, v_right, ratio, zc = head_velocity_metrics(m_c, ev.t_collision_s)
    return {
        "collision_status": ev.collision_status,
        "t_collision": ev.t_collision_s,
        "z_collision": zc,
        "positive_head_velocity": v_left,
        "negative_head_velocity": v_right,
        "head_velocity_ratio": ratio,
        "E_gap_peak": ev.E_gap_peak_V_m,
        "E_gap_drop_fraction": ev.field_drop_fraction,
        "bridge_growth": ev.bridge_growth_ratio,
        "I_CM_peak": float(np.nanmax(np.abs(cm_c.I_CM_drift))),
        "Delta_I_peak": pulse["Delta_I_peak_abs"],
        "derivative_peak": dsum["local_poly_peak"],
        "pulse_FWHM": pulse_fwhm(delta.time.to_numpy(float), delta.Delta_I_CM_drift.to_numpy(float)),
        "spectral_centroid": spectral_centroid(spec, esd),
        "trusted_frequency_limit": rsum["f_trust_Hz"],
        "VHF_metric": band_metric(bands, "VHF"),
        "UHF_metric": band_metric(bands, "UHF"),
        "SHF_metric": band_metric(bands, "SHF"),
        "wall_time": term_c.wall_time,
        "accepted_steps": term_c.accepted_steps,
    }


def smoothed_derivative_peak(t: np.ndarray, y: np.ndarray) -> float:
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    max_dt = float(np.max(np.diff(t)))
    tu = np.arange(t[0], t[-1] + 0.5 * max_dt, max_dt)
    yu = np.interp(tu, t, y)
    window = max(5, int(round(20e-12 / max_dt)) | 1)
    if window >= len(yu):
        window = max(5, (len(yu) // 2) * 2 - 1)
    ys = savgol_filter(yu, window, 3) if window >= 5 and len(yu) > window else yu
    return float(np.max(np.abs(np.gradient(ys, max_dt))))


def smoothed_lowband_centroid(t: np.ndarray, y: np.ndarray) -> float:
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    max_dt = float(np.max(np.diff(t)))
    tu = np.arange(t[0], t[-1] + 0.5 * max_dt, max_dt)
    yu = np.interp(tu, t, y)
    window = max(5, int(round(20e-12 / max_dt)) | 1)
    if window >= len(yu):
        window = max(5, (len(yu) // 2) * 2 - 1)
    ys = savgol_filter(yu, window, 3) if window >= 5 and len(yu) > window else yu
    dy = np.gradient(ys, max_dt)
    f = np.fft.rfftfreq(len(tu), max_dt)
    power = np.abs(np.fft.rfft(dy)) ** 2
    mask = (f >= 30e6) & (f <= 1.0e9)
    if not np.any(mask) or np.sum(power[mask]) <= 0:
        return np.nan
    return float(np.sum(f[mask] * power[mask]) / np.sum(power[mask]))


def stage5_sampling_sensitivity(delta: pd.DataFrame, t_collision: float) -> pd.DataFrame:
    rows = []
    base = delta.copy()
    base_peak = float(np.max(np.abs(base.Delta_I_CM_drift)))
    base_deriv = smoothed_derivative_peak(base.time.to_numpy(float), base.Delta_I_CM_drift.to_numpy(float))
    base_centroid = smoothed_lowband_centroid(base.time.to_numpy(float), base.Delta_I_CM_drift.to_numpy(float))
    for stride in [1, 2, 4]:
        d = delta.iloc[::stride].copy()
        if float(d.time.iloc[-1]) < float(delta.time.iloc[-1]):
            d = pd.concat([d, delta.iloc[[-1]]], ignore_index=True)
        peak = float(np.max(np.abs(d.Delta_I_CM_drift)))
        deriv = smoothed_derivative_peak(d.time.to_numpy(float), d.Delta_I_CM_drift.to_numpy(float))
        centroid = smoothed_lowband_centroid(d.time.to_numpy(float), d.Delta_I_CM_drift.to_numpy(float))
        rows.append({
            "sampling": f"every_{stride}_sample",
            "sample_count": len(d),
            "Delta_I_peak": peak,
            "Delta_I_peak_rel_diff": abs(peak - base_peak) / max(base_peak, np.finfo(float).tiny),
            "smoothed_derivative_peak": deriv,
            "smoothed_derivative_peak_rel_diff": abs(deriv - base_deriv) / max(base_deriv, np.finfo(float).tiny),
            "spectral_centroid_Hz": centroid,
            "spectral_centroid_rel_diff": abs(centroid - base_centroid) / max(base_centroid, np.finfo(float).tiny),
            "status": "PASS" if (
                abs(peak - base_peak) / max(base_peak, np.finfo(float).tiny) <= 0.05
                and abs(deriv - base_deriv) / max(base_deriv, np.finfo(float).tiny) <= 0.05
                and abs(centroid - base_centroid) / max(base_centroid, np.finfo(float).tiny) <= 0.05
            ) else "FAIL",
        })
    return pd.DataFrame(rows)


def baseline_metrics() -> dict[str, float]:
    ev = pd.read_csv(ROOT / "results/stage4/collision/collision_event.csv").iloc[0]
    pulse = pd.read_csv(ROOT / "results/stage4/current_moment/delta_current_pulse_metrics.csv").iloc[0]
    dsum = pd.read_csv(ROOT / "results/stage4/radiation/derivative_sensitivity.csv").iloc[0]
    bands = pd.read_csv(ROOT / "results/stage4/radiation/band_energy.csv")
    spec = pd.read_csv(ROOT / "results/stage4/radiation/delta_current_spectrum.csv")
    esd = pd.read_csv(ROOT / "results/stage4/radiation/delta_current_esd.csv")
    cm = pd.read_csv(ROOT / "results/stage4/current_moment/collision.csv")
    metrics = pd.read_csv(ROOT / "results/stage4/runs/paperlike_20um/collision_metrics.csv")
    vl, vr, ratio, zc = head_velocity_metrics(metrics, float(ev.t_collision_s))
    return {
        "collision_status": ev.collision_status,
        "t_collision": float(ev.t_collision_s),
        "z_collision": zc,
        "positive_head_velocity": vl,
        "negative_head_velocity": vr,
        "head_velocity_ratio": ratio,
        "E_gap_peak": float(ev.E_gap_peak_V_m),
        "E_gap_drop_fraction": float(ev.field_drop_fraction),
        "bridge_growth": float(ev.bridge_growth_ratio),
        "I_CM_peak": float(np.nanmax(np.abs(cm.I_CM_drift))),
        "Delta_I_peak": float(pulse.Delta_I_peak_abs),
        "derivative_peak": float(dsum.local_poly_peak),
        "pulse_FWHM": np.nan,
        "spectral_centroid": spectral_centroid(spec, esd),
        "trusted_frequency_limit": np.nan,
        "VHF_metric": band_metric(bands, "VHF"),
        "UHF_metric": band_metric(bands, "UHF"),
        "SHF_metric": band_metric(bands, "SHF"),
    }


def proxy_case(root: Path, name: str, run_id: str, factor: str, value: str) -> dict[str, float | str]:
    _, metrics, current, term = load_run(root, name)
    ev = detect_collision_event(metrics)
    event_to_frame(ev).to_csv(root / f"collision/{run_id.lower()}_event_diagnostic.csv", index=False)
    proxy, deriv, spec, bands, summary = event_local_proxy_products(current, metrics)
    safe = run_id.lower().replace("-", "_")
    proxy.to_csv(root / f"current_moment/{safe}_event_local_proxy_current.csv", index=False)
    deriv.to_csv(root / f"radiation/{safe}_event_local_proxy_derivative.csv", index=False)
    spec.to_csv(root / f"radiation/{safe}_event_local_proxy_spectrum.csv", index=False)
    bands.to_csv(root / f"radiation/{safe}_event_local_proxy_band_energy.csv", index=False)
    vl, vr, ratio, zc = head_velocity_metrics(metrics, ev.t_collision_s if np.isfinite(ev.t_collision_s) else None)
    return {
        "case_id": run_id,
        "factor": factor,
        "parameter_value": value,
        "collision_status": robust_status_from_termination(term, ev.collision_status),
        "t_collision": ev.t_collision_s if ev.collision_status == "PASS" and term.termination_reason != "NO_COLLISION_WITHIN_RESOURCE_WINDOW" else np.nan,
        "z_collision": zc,
        "positive_head_velocity": vl,
        "negative_head_velocity": vr,
        "head_velocity_ratio": ratio,
        "E_gap_peak": ev.E_gap_peak_V_m,
        "E_gap_drop_fraction": ev.field_drop_fraction,
        "bridge_growth": ev.bridge_growth_ratio,
        "I_CM_peak": summary["I_CM_peak"],
        "Delta_I_peak": np.nan,
        "derivative_peak": summary["proxy_derivative_peak"],
        "pulse_FWHM": summary["proxy_FWHM_s"],
        "spectral_centroid": summary["spectral_centroid_Hz"],
        "VHF_metric": band_metric(bands, "VHF"),
        "UHF_metric": band_metric(bands, "UHF"),
        "SHF_metric": band_metric(bands, "SHF"),
        "metric_type": "event_local_proxy",
        "run_id": run_id,
        "evidence_status": "derived_from_measured",
    }


def mpi_consistency(root: Path) -> pd.DataFrame:
    one = read_csv_clean(run_dir(root, "highfield_mpi_probe_1rank") / "scalar_history.csv")
    two = read_csv_clean(run_dir(root, "highfield_mpi_probe_2rank") / "scalar_history.csv")
    cols = ["E_max_V_m", "ne_max_m_3", "total_electrons", "total_charge_C"]
    rows = []
    t_end = min(float(one.time_s.max()), float(two.time_s.max()))
    for col in cols:
        a = float(np.interp(t_end, one.time_s, one[col]))
        b = float(np.interp(t_end, two.time_s, two[col]))
        absdiff = abs(a - b)
        rel = absdiff / max(abs(a), abs(b), np.finfo(float).tiny)
        abs_threshold = 1e-24 if col == "total_charge_C" else 0.0
        ok = rel <= 1e-9 or absdiff <= abs_threshold
        rows.append({"metric": col, "time_s": t_end, "rank1": a, "rank2": b, "absolute_difference": absdiff, "relative_difference": rel, "relative_threshold": 1e-9, "absolute_threshold": abs_threshold, "status": "PASS" if ok else "FAIL"})
    return pd.DataFrame(rows)


def manifest(root: Path) -> None:
    rows = []
    for p in sorted(root.glob("**/*")):
        if p.is_file() and any(part in {"collision", "current_moment", "radiation", "trends", "mpi"} for part in p.parts):
            rows.append({"file": str(p.relative_to(ROOT)), "sha256": sha(p), "evidence_status": "derived_from_measured"})
    pd.DataFrame(rows).to_csv(root / "provenance/derived_artifact_manifest.csv", index=False)


def main() -> None:
    root = ROOT / "results/stage5"
    ensure_dirs(root)
    high = write_highfield(root)
    base = baseline_metrics()
    rows = []
    rows.append({"case_id": "B0", "factor": "baseline", "parameter_value": "Stage4 S4-PAPERLIKE-20UM", **base, "metric_type": "strict_delta", "run_id": "S4-PAPERLIKE-20UM", "evidence_status": "derived_from_measured"})
    rows.append({"case_id": "F", "factor": "background_field", "parameter_value": "2.0x3.2e6 V/m", **high, "metric_type": "strict_delta", "run_id": "S5-HIGHFIELD-COLLISION", "evidence_status": "derived_from_measured"})
    rows.append(proxy_case(root, "largergap_collision_20um", "S5-LARGERGAP-COLLISION", "seed_gap", "5.0 mm"))
    rows.append(proxy_case(root, "asymmetric_collision_20um", "S5-ASYMMETRIC-COLLISION", "seed_sigma_asymmetry", "sigma1=0.1mm,sigma2=0.5mm"))
    comp = pd.DataFrame(rows)
    ordered = [
        "case_id","factor","parameter_value","collision_status","t_collision","z_collision",
        "positive_head_velocity","negative_head_velocity","head_velocity_ratio","E_gap_peak",
        "E_gap_drop_fraction","bridge_growth","I_CM_peak","Delta_I_peak","derivative_peak",
        "pulse_FWHM","spectral_centroid","VHF_metric","UHF_metric","SHF_metric","metric_type",
        "run_id","evidence_status",
    ]
    comp[ordered].to_csv(root / "trends/case_comparison.csv", index=False)
    mpi = mpi_consistency(root)
    mpi.to_csv(root / "mpi/rank_consistency.csv", index=False)
    manifest(root)
    print("Stage 5 analysis products written")


if __name__ == "__main__":
    main()
