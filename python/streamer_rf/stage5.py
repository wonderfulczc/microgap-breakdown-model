from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from streamer_rf.stage4 import (
    BANDS_HZ,
    align_current_moments,
    collision_pulse_metrics,
    derivative_products,
    detect_collision_event,
    output_sampling_sensitivity,
    radiation_products,
    read_csv_clean,
)


@dataclass(frozen=True)
class CaseMetrics:
    case_id: str
    factor: str
    parameter_value: str
    collision_status: str
    t_collision: float
    z_collision: float
    positive_head_velocity: float
    negative_head_velocity: float
    head_velocity_ratio: float
    E_gap_peak: float
    E_gap_drop_fraction: float
    bridge_growth: float
    I_CM_peak: float
    Delta_I_peak: float
    derivative_peak: float
    pulse_FWHM: float
    spectral_centroid: float
    VHF_metric: float
    UHF_metric: float
    SHF_metric: float
    metric_type: str
    run_id: str
    evidence_status: str


def robust_status_from_termination(term: pd.Series, event_status: str) -> str:
    reason = str(term.get("termination_reason", "unknown"))
    if reason == "NO_COLLISION_WITHIN_RESOURCE_WINDOW":
        return "NO_COLLISION_WITHIN_RESOURCE_WINDOW"
    return event_status


def head_velocity_metrics(metrics: pd.DataFrame, t_ref: float | None = None) -> tuple[float, float, float, float]:
    m = metrics.sort_values("time").drop_duplicates("time").reset_index(drop=True)
    if len(m) < 8:
        return np.nan, np.nan, np.nan, np.nan
    if t_ref is None or not np.isfinite(t_ref):
        sub = m.iloc[max(0, len(m) // 2):]
    else:
        sub = m[(m.time >= max(float(m.time.min()), t_ref - 0.4e-9)) & (m.time <= t_ref)]
        if len(sub) < 8:
            sub = m.iloc[max(0, len(m) // 2):]
    def slope(col: str) -> float:
        x = sub.time.to_numpy(float)
        y = sub[col].to_numpy(float)
        if len(np.unique(x)) < 2:
            return np.nan
        return float(np.polyfit(x, y, 1)[0])
    v_left = slope("left_inner_z")
    v_right = slope("right_inner_z")
    ratio = abs(v_left) / max(abs(v_right), np.finfo(float).tiny)
    z_collision = float(np.nanmedian([sub.left_inner_z.iloc[-1], sub.right_inner_z.iloc[-1]]))
    return v_left, v_right, ratio, z_collision


def pulse_fwhm(t: np.ndarray, y: np.ndarray) -> float:
    t = np.asarray(t, float)
    y = np.abs(np.asarray(y, float))
    if len(t) < 3 or np.nanmax(y) <= 0:
        return np.nan
    half = 0.5 * float(np.nanmax(y))
    hit = np.where(y >= half)[0]
    if len(hit) < 2:
        return np.nan
    return float(t[hit[-1]] - t[hit[0]])


def spectral_centroid(spec: pd.DataFrame, esd: pd.DataFrame | None = None) -> float:
    s = spec[spec.get("trusted", True)].copy()
    if len(s) < 2:
        s = spec.copy()
    f = s.frequency_Hz.to_numpy(float)
    if esd is not None and "ESD_J_per_Hz" in esd:
        y = esd.loc[s.index, "ESD_J_per_Hz"].to_numpy(float)
    elif "abs_dDeltaI_dt_transform" in s:
        y = s.abs_dDeltaI_dt_transform.to_numpy(float)
    else:
        return np.nan
    den = float(np.nansum(y))
    return float(np.nansum(f * y) / den) if den > 0 else np.nan


def event_local_proxy_products(current: pd.DataFrame, metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, float]]:
    ev = detect_collision_event(metrics)
    t_event = ev.t_collision_s if np.isfinite(ev.t_collision_s) else float(metrics.time.iloc[len(metrics) // 2])
    cm = current.sort_values("time").drop_duplicates("time")
    proxy = pd.DataFrame({"time": cm.time, "Delta_I_CM_drift": cm.I_CM_drift - np.nanmedian(cm.I_CM_drift.iloc[: max(4, len(cm) // 10)])})
    deriv, dsum = derivative_products(proxy, t_event)
    spec, esd, bands, rsum = radiation_products(proxy, t_event)
    summary = {
        "t_event_proxy_s": t_event,
        "I_CM_peak": float(np.nanmax(np.abs(cm.I_CM_drift))),
        "proxy_derivative_peak": dsum["local_poly_peak"],
        "proxy_FWHM_s": pulse_fwhm(proxy.time.to_numpy(float), proxy.Delta_I_CM_drift.to_numpy(float)),
        "spectral_centroid_Hz": spectral_centroid(spec, esd),
        "trusted_frequency_limit_Hz": rsum["f_trust_Hz"],
    }
    return proxy, deriv, spec, bands, summary


def band_metric(bands: pd.DataFrame, name: str) -> float:
    row = bands[bands.band == name]
    if len(row) == 0:
        return np.nan
    return float(row.integrated_energy_J.iloc[0]) if pd.notna(row.integrated_energy_J.iloc[0]) else np.nan
