from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.constants import c, epsilon_0

from streamer_rf.esd import integrate_band_energy


BANDS_HZ = {
    "VHF": (30e6, 300e6),
    "UHF": (0.3e9, 3e9),
    "SHF": (3e9, 30e9),
}


def esd_per_hz_from_derivative_spectrum(dI_dt_transform: np.ndarray) -> np.ndarray:
    s = np.asarray(dI_dt_transform)
    return 2 * np.abs(s) ** 2 / (3 * epsilon_0 * c**3)


@dataclass(frozen=True)
class CollisionEvent:
    t_distance_s: float
    t_bridge_s: float
    t_field_collapse_s: float
    t_collision_s: float
    bridge_threshold_m_3: float
    bridge_pre_median_m_3: float
    bridge_late_median_m_3: float
    bridge_growth_ratio: float
    E_gap_peak_V_m: float
    E_gap_peak_time_s: float
    E_gap_min_after_V_m: float
    E_gap_min_after_time_s: float
    field_drop_fraction: float
    d_head_initial_m: float
    d_head_min_m: float
    d_head_min_time_s: float
    distance_status: str
    bridge_status: str
    field_status: str
    collision_status: str


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv_clean(path: Path) -> pd.DataFrame:
    """Read a CSV and drop rows corrupted by interrupted concurrent writes."""
    df = pd.read_csv(path, on_bad_lines="skip")
    for c in df.columns:
        converted = pd.to_numeric(df[c], errors="coerce")
        if converted.notna().sum() >= max(1, int(0.5 * len(df))):
            df[c] = converted
    numeric = df.select_dtypes(include=[np.number]).columns
    if len(numeric):
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=numeric)
    return df.reset_index(drop=True)


def detect_collision_event(metrics: pd.DataFrame) -> CollisionEvent:
    required = {"time", "d_head", "bridge_mean_ne", "E_gap_max"}
    if not required.issubset(metrics.columns):
        raise ValueError(f"collision metrics missing {sorted(required - set(metrics.columns))}")
    m = metrics.sort_values("time").reset_index(drop=True)
    if len(m) < 16:
        raise ValueError("not enough collision metric samples")

    early = m[m.time <= m.time.quantile(0.20)]
    late = m[m.time >= m.time.quantile(0.90)]
    bridge_pre = float(early.bridge_mean_ne.median())
    bridge_late = float(late.bridge_mean_ne.median())
    bridge_threshold = bridge_pre + 0.5 * max(bridge_late - bridge_pre, 0.0)
    bridge_hits = m[m.bridge_mean_ne >= bridge_threshold]
    t_bridge = float(bridge_hits.time.iloc[0]) if len(bridge_hits) else np.nan
    bridge_growth = bridge_late / max(bridge_pre, np.finfo(float).tiny)
    bridge_status = "PASS" if np.isfinite(t_bridge) and bridge_growth >= 10.0 else "FAIL"

    peak_idx = int(m.E_gap_max.idxmax())
    peak = m.loc[peak_idx]
    after = m[m.time > peak.time]
    threshold_e = float(peak.E_gap_max) * 0.8
    collapse_hits = after[after.E_gap_max <= threshold_e]
    t_collapse = float(collapse_hits.time.iloc[0]) if len(collapse_hits) else np.nan
    if len(after):
      min_after = after.loc[after.E_gap_max.idxmin()]
    else:
      min_after = peak
    drop = (float(peak.E_gap_max) - float(min_after.E_gap_max)) / max(float(peak.E_gap_max), np.finfo(float).tiny)
    field_status = "PASS" if np.isfinite(t_collapse) and drop >= 0.20 else "FAIL"

    d_initial = float(m.d_head.iloc[0])
    min_row = m.loc[m.d_head.idxmin()]
    d_min = float(min_row.d_head)
    t_distance = float(min_row.time)
    dz_guess = 2.0e-5
    distance_strong = d_min <= max(2.0 * dz_guess, 2.0e-4)
    distance_closing = d_min < 0.75 * d_initial
    if distance_strong:
        distance_status = "PASS"
    elif distance_closing:
        distance_status = "NOISY_NONCONTRADICTORY"
    else:
        distance_status = "FAIL"

    strong_times = [x for x in [t_bridge, t_collapse, t_distance if distance_strong else np.nan] if np.isfinite(x)]
    t_collision = float(np.median(strong_times[:])) if strong_times else np.nan
    strong_count = int(bridge_status == "PASS") + int(field_status == "PASS") + int(distance_status == "PASS")
    noncontradictory = distance_status in {"PASS", "NOISY_NONCONTRADICTORY"}
    collision_status = "PASS" if strong_count >= 2 or (bridge_status == "PASS" and field_status == "PASS" and noncontradictory) else "FAIL"

    return CollisionEvent(
        t_distance_s=t_distance,
        t_bridge_s=t_bridge,
        t_field_collapse_s=t_collapse,
        t_collision_s=t_collision,
        bridge_threshold_m_3=float(bridge_threshold),
        bridge_pre_median_m_3=bridge_pre,
        bridge_late_median_m_3=bridge_late,
        bridge_growth_ratio=float(bridge_growth),
        E_gap_peak_V_m=float(peak.E_gap_max),
        E_gap_peak_time_s=float(peak.time),
        E_gap_min_after_V_m=float(min_after.E_gap_max),
        E_gap_min_after_time_s=float(min_after.time),
        field_drop_fraction=float(drop),
        d_head_initial_m=d_initial,
        d_head_min_m=d_min,
        d_head_min_time_s=t_distance,
        distance_status=distance_status,
        bridge_status=bridge_status,
        field_status=field_status,
        collision_status=collision_status,
    )


def event_to_frame(ev: CollisionEvent) -> pd.DataFrame:
    return pd.DataFrame([ev.__dict__])


def align_current_moments(collision: pd.DataFrame, left: pd.DataFrame, right: pd.DataFrame, t_max: float | None = None) -> pd.DataFrame:
    for name, df in {"collision": collision, "left": left, "right": right}.items():
        if not {"time", "I_CM_drift", "I_CM_electron"}.issubset(df.columns):
            raise ValueError(f"{name} current moment missing required columns")
    c = collision.sort_values("time").drop_duplicates("time")
    l = left.sort_values("time").drop_duplicates("time")
    r = right.sort_values("time").drop_duplicates("time")
    t_end = min(float(c.time.max()), float(l.time.max()), float(r.time.max()))
    if t_max is not None:
        t_end = min(t_end, float(t_max))
    grid = c[(c.time >= max(float(l.time.min()), float(r.time.min()), float(c.time.min()))) & (c.time <= t_end)].time.to_numpy()
    if len(grid) < 8:
        raise ValueError("not enough overlapping current samples")
    out = pd.DataFrame({"time": grid})
    for col in ["I_CM_drift", "I_CM_electron"]:
        out[f"I_collision_{col}"] = np.interp(grid, c.time, c[col])
        out[f"I_left_{col}"] = np.interp(grid, l.time, l[col])
        out[f"I_right_{col}"] = np.interp(grid, r.time, r[col])
        out[f"Delta_{col}"] = out[f"I_collision_{col}"] - out[f"I_left_{col}"] - out[f"I_right_{col}"]
    dt_sources = np.r_[np.diff(c.time.to_numpy()), np.diff(l.time.to_numpy()), np.diff(r.time.to_numpy())]
    out["interpolation_flag"] = True
    out["local_max_original_dt"] = float(np.nanmax(dt_sources))
    return out


def collision_pulse_metrics(delta: pd.DataFrame, t_collision: float) -> dict[str, float | str]:
    y = delta["Delta_I_CM_drift"].to_numpy(float)
    t = delta["time"].to_numpy(float)
    pre = y[t < t_collision - 0.25e-9]
    if len(pre) < 8:
        pre = y[: max(8, len(y) // 5)]
    med = float(np.median(pre))
    sigma = float(1.4826 * np.median(np.abs(pre - med))) if len(pre) else 0.0
    floor = max(1e-18, 10 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(y)))))
    noise = max(sigma, floor)
    idx = int(np.argmax(np.abs(y)))
    peak = float(y[idx])
    snr = abs(peak) / noise
    near = abs(float(t[idx]) - t_collision) <= 0.4e-9
    status = "PASS" if snr >= 5.0 and near else "FAIL"
    return {
        "pre_median": med,
        "sigma_pre": sigma,
        "noise_floor": floor,
        "noise_used": noise,
        "Delta_I_peak": peak,
        "Delta_I_peak_abs": abs(peak),
        "Delta_I_peak_time_s": float(t[idx]),
        "Delta_I_snr": float(snr),
        "causality_status": status,
    }


def local_poly_derivative(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    if t.ndim != 1 or y.shape != t.shape or len(t) < 5:
        raise ValueError("matching 1D arrays with at least five points required")
    out = np.empty_like(y)
    for i in range(len(t)):
        lo = max(0, i - 2)
        hi = min(len(t), lo + 5)
        lo = max(0, hi - 5)
        x = t[lo:hi] - t[i]
        deg = min(3, len(x) - 1)
        coeff = np.polyfit(x, y[lo:hi], deg)
        out[i] = np.polyder(np.poly1d(coeff))(0.0)
    return out


def uniform_grid_from_original(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    max_dt = float(np.max(np.diff(t)))
    dt = max_dt
    n = int(np.floor((t[-1] - t[0]) / dt)) + 1
    if n < 8:
        raise ValueError("not enough uniform samples")
    tu = t[0] + np.arange(n) * dt
    return tu, np.interp(tu, t, y), max_dt


def five_point_uniform_derivative(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    dt = float(t[1] - t[0])
    out = np.gradient(y, dt, edge_order=2)
    if len(y) >= 5:
        out[2:-2] = (y[:-4] - 8 * y[1:-3] + 8 * y[3:-1] - y[4:]) / (12 * dt)
    return out


def spectral_derivative(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    dt = float(t[1] - t[0])
    freq = np.fft.rfftfreq(len(t), dt)
    z = np.fft.rfft(y)
    return np.fft.irfft(1j * 2 * np.pi * freq * z, n=len(t))


def derivative_products(delta: pd.DataFrame, t_collision: float) -> tuple[pd.DataFrame, dict[str, float]]:
    t = delta.time.to_numpy(float)
    y = delta.Delta_I_CM_drift.to_numpy(float)
    local = local_poly_derivative(t, y)
    tu, yu, max_dt = uniform_grid_from_original(t, y)
    five = five_point_uniform_derivative(tu, yu)
    spec = spectral_derivative(tu, yu)
    five_i = np.interp(t, tu, five)
    spec_i = np.interp(t, tu, spec)
    trusted = np.abs(t - t_collision) <= 0.5e-9
    out = pd.DataFrame({
        "time": t,
        "dDeltaI_dt_local_poly": local,
        "dDeltaI_dt_uniform_five_point": five_i,
        "dDeltaI_dt_spectral": spec_i,
        "trusted": trusted,
        "reason": np.where(trusted, "inside_collision_window", "outside_collision_window"),
    })
    win = trusted
    p1 = float(np.max(np.abs(local[win]))) if np.any(win) else float(np.max(np.abs(local)))
    p2 = float(np.max(np.abs(five_i[win]))) if np.any(win) else float(np.max(np.abs(five_i)))
    rel = abs(p1 - p2) / max(p1, p2, np.finfo(float).tiny)
    return out, {"local_poly_peak": p1, "uniform_five_point_peak": p2, "peak_relative_difference": float(rel), "max_original_dt_s": max_dt}


def radiation_products(delta: pd.DataFrame, t_collision: float, window: str = "rectangular") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, float]]:
    sub = delta[(delta.time >= t_collision - 0.5e-9) & (delta.time <= t_collision + 0.3e-9)].copy()
    if len(sub) < 16:
        sub = delta.copy()
    t = sub.time.to_numpy(float)
    y = sub.Delta_I_CM_drift.to_numpy(float)
    # Remove a linear baseline fitted to the edge samples outside the central pulse.
    edge = (t < t_collision - 0.25e-9) | (t > t_collision + 0.20e-9)
    if np.count_nonzero(edge) >= 2:
        p = np.polyfit(t[edge] - t_collision, y[edge], 1)
        y = y - np.polyval(p, t - t_collision)
    tu, yu, max_dt = uniform_grid_from_original(t, y)
    if window == "tukey":
        try:
            from scipy.signal.windows import tukey
            w = tukey(len(yu), alpha=0.25)
        except Exception:
            w = np.hanning(len(yu))
    else:
        w = np.ones(len(yu))
    yw = yu * w
    dy = five_point_uniform_derivative(tu, yw)
    n_fft = 1 << int(np.ceil(np.log2(max(len(tu), 16))))
    dt = float(tu[1] - tu[0])
    f = np.fft.rfftfreq(n_fft, dt)
    z = dt * np.fft.rfft(dy, n=n_fft) * np.exp(-2j * np.pi * f * tu[0])
    amp = np.abs(z)
    esd = esd_per_hz_from_derivative_spectrum(z)
    f_ny = 1.0 / (2.0 * max_dt)
    f_trust = 0.5 * f_ny
    spec = pd.DataFrame({"frequency_Hz": f, "abs_dDeltaI_dt_transform": amp, "trusted": f <= f_trust})
    esd_df = pd.DataFrame({"frequency_Hz": f, "ESD_J_per_Hz": esd, "trusted": f <= f_trust})
    bands = []
    for name, (lo, hi) in BANDS_HZ.items():
        cover_hi = min(hi, f_trust, float(f[-1]))
        if f_trust < lo:
            status, energy = "untrusted", np.nan
        elif f_trust < hi:
            status = "partial"
            energy = integrate_band_energy(f, esd, lo, cover_hi) if cover_hi > lo else np.nan
        else:
            status = "trusted"
            energy = integrate_band_energy(f, esd, lo, hi)
        bands.append({"band": name, "f_low_Hz": lo, "f_high_Hz": hi, "trusted_status": status, "integrated_energy_J": energy})
    return spec, esd_df, pd.DataFrame(bands), {"max_original_dt_s": max_dt, "f_Nyquist_conservative_Hz": f_ny, "f_trust_Hz": f_trust}


def output_sampling_sensitivity(delta: pd.DataFrame, t_collision: float) -> pd.DataFrame:
    base_pulse = collision_pulse_metrics(delta, t_collision)
    base_deriv, base_deriv_summary = derivative_products(delta, t_collision)
    base_peak_d = base_deriv_summary["local_poly_peak"]
    rows = []
    for stride in [1, 2, 4]:
        d = delta.iloc[::stride].copy()
        if float(d.time.iloc[-1]) < float(delta.time.iloc[-1]):
            d = pd.concat([d, delta.iloc[[-1]]], ignore_index=True)
        pulse = collision_pulse_metrics(d, t_collision)
        deriv, deriv_summary = derivative_products(d, t_collision)
        rows.append({
            "sampling": f"every_{stride}_sample",
            "sample_count": len(d),
            "Delta_I_peak": pulse["Delta_I_peak"],
            "Delta_I_peak_time_s": pulse["Delta_I_peak_time_s"],
            "Delta_I_peak_rel_diff": abs(pulse["Delta_I_peak_abs"] - base_pulse["Delta_I_peak_abs"]) / max(base_pulse["Delta_I_peak_abs"], np.finfo(float).tiny),
            "derivative_peak": deriv_summary["local_poly_peak"],
            "derivative_peak_rel_diff": abs(deriv_summary["local_poly_peak"] - base_peak_d) / max(base_peak_d, np.finfo(float).tiny),
            "causality_status": pulse["causality_status"],
        })
    return pd.DataFrame(rows)
