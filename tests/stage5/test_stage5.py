from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.stage5 import head_velocity_metrics, pulse_fwhm, spectral_centroid


def test_case_comparison_separates_strict_delta_and_proxy():
    comp = pd.read_csv(ROOT / "results/stage5/trends/case_comparison.csv")
    assert {"strict_delta", "event_local_proxy"}.issubset(set(comp.metric_type))
    assert comp.loc[comp.case_id == "F", "metric_type"].iloc[0] == "strict_delta"
    assert comp.loc[comp.run_id == "S5-LARGERGAP-COLLISION", "metric_type"].iloc[0] == "event_local_proxy"
    proxy = comp[comp.metric_type == "event_local_proxy"]
    assert proxy.Delta_I_peak.isna().all()


def test_no_collision_resource_status_is_explicit():
    comp = pd.read_csv(ROOT / "results/stage5/trends/case_comparison.csv")
    statuses = set(comp.collision_status)
    assert "NO_COLLISION_WITHIN_RESOURCE_WINDOW" in statuses
    assert comp.loc[comp.run_id == "S5-LARGERGAP-COLLISION", "t_collision"].isna().all()


def test_pulse_fwhm_and_spectral_centroid_are_measured():
    t = np.linspace(0.0, 2e-9, 101)
    y = np.exp(-((t - 1e-9) / 0.12e-9) ** 2)
    assert pulse_fwhm(t, y) > 0.0
    spec = pd.DataFrame({"frequency_Hz": [0.0, 1.0, 2.0], "abs_dDeltaI_dt_transform": [0.0, 1.0, 1.0], "trusted": [True, True, True]})
    assert 1.0 < spectral_centroid(spec) < 2.1


def test_head_velocity_ratio_uses_trajectories():
    t = np.linspace(0, 1e-9, 16)
    m = pd.DataFrame({
        "time": t,
        "left_inner_z": 0.003 + 1e6 * t,
        "right_inner_z": 0.007 - 2e6 * t,
        "d_head": 0.004 - 3e6 * t,
        "bridge_mean_ne": np.linspace(1e16, 1e19, len(t)),
        "E_gap_max": np.linspace(5e6, 6e6, len(t)),
    })
    vl, vr, ratio, zc = head_velocity_metrics(m, 1e-9)
    assert vl > 0
    assert vr < 0
    assert 0.45 < ratio < 0.55
    assert np.isfinite(zc)


def test_stage5_output_sampling_sensitivity_passes_thresholds():
    s = pd.read_csv(ROOT / "results/stage5/radiation/highfield_output_sampling_sensitivity.csv")
    assert set(s.status) == {"PASS"}
    assert s.Delta_I_peak_rel_diff.max() <= 0.05
    assert s.smoothed_derivative_peak_rel_diff.max() <= 0.05
    assert s.spectral_centroid_rel_diff.max() <= 0.05
