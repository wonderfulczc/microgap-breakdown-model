from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.constants import C0  # noqa: E402
from streamer_rf.rf.jefimenko.nath import (  # noqa: E402
    NathAvalancheConfig,
    avalanche_state,
    evaluate_nath_at_observer_time,
    evaluate_nath_field,
    make_nath_source_series,
    nath_K,
    retarded_time_for_moving_head,
    source_charge_conservation_summary,
    waveform_metrics,
)
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_waveform  # noqa: E402


def _observer_time(observer: Observer, config: NathAvalancheConfig, tr: float) -> float:
    state = avalanche_state(config, tr)
    return float(tr + np.linalg.norm(observer.position - state.position_m) / C0)


def test_nath_K_nonrelativistic_limit_and_constant_source_zero_B() -> None:
    rhat = np.array([1.0, 0.0, 0.0])
    assert abs(nath_K(rhat, np.zeros(3)) - 1.0) < 1e-15
    cfg = NathAvalancheConfig(growth_rate_s=0.0, velocity_m_s=0.0)
    state = avalanche_state(cfg, 1e-9)
    fields = evaluate_nath_field(state, np.array([1e-3, 0.0, 0.0]), config=None)
    assert np.linalg.norm(fields.B_total) == 0.0
    assert np.linalg.norm(fields.E_charge_growth) == 0.0


def test_constant_velocity_charge_produces_transverse_near_B_and_zero_Bz() -> None:
    cfg = NathAvalancheConfig(growth_rate_s=0.0, velocity_m_s=1e5)
    obs = np.array([0.25e-3, 0.25e-3, 1e-3])
    fields = evaluate_nath_field(avalanche_state(cfg, 1e-9), obs, config=None)
    assert np.linalg.norm(fields.B_near) > 0.0
    assert abs(fields.B_total[2]) < np.linalg.norm(fields.B_total) * 1e-14
    assert np.linalg.norm(fields.B_charge_growth_radiation) == 0.0
    assert np.linalg.norm(fields.B_acceleration_radiation) == 0.0


def test_charge_growth_and_acceleration_terms_are_separable() -> None:
    obs = np.array([0.25e-3, 0.25e-3, 1e-3])
    growth = NathAvalancheConfig(growth_rate_s=5e8, velocity_m_s=1e5)
    fg = evaluate_nath_field(avalanche_state(growth, 1e-9), obs, config=None)
    assert np.linalg.norm(fg.B_charge_growth_radiation) > 0.0
    assert np.linalg.norm(fg.E_charge_growth) > 0.0

    accel_state = avalanche_state(NathAvalancheConfig(growth_rate_s=0.0, velocity_m_s=1e5), 1e-9, acceleration_m_s2=5e11)
    fa = evaluate_nath_field(accel_state, obs, config=None)
    assert np.linalg.norm(fa.B_acceleration_radiation) > 0.0


def test_retarded_time_root_stationary_constant_velocity_and_prearrival() -> None:
    stationary = NathAvalancheConfig(growth_rate_s=0.0, velocity_m_s=0.0)
    obs = Observer("x", 0.3, 0.0, 0.0)
    tr, ok, _ = retarded_time_for_moving_head(obs.position, stationary, 2e-9)
    assert ok
    assert abs((2e-9 - tr) - 0.3 / C0) < 1e-18

    moving = NathAvalancheConfig(velocity_m_s=1e5)
    tr0 = 1e-9
    tobs = _observer_time(obs, moving, tr0)
    tr, ok, _ = retarded_time_for_moving_head(obs.position, moving, tobs)
    assert ok
    assert abs(tr - tr0) < 1e-18

    _, _, tr_pre, valid = evaluate_nath_at_observer_time(moving, obs, 0.2e-9)
    assert tr_pre < moving.start_time_s
    assert not valid


def test_regularized_source_charge_conservation() -> None:
    cfg = NathAvalancheConfig()
    series = make_nath_source_series(cfg, np.linspace(0.0, 3e-9, 101))
    summary = source_charge_conservation_summary(series)
    assert summary["Q_span_C"] < cfg.q0_C * 1e-12
    assert series.records[0].n_cells > 10


def test_existing_jefimenko_matches_nath_near_nontrivial_components() -> None:
    cfg = NathAvalancheConfig()
    series = make_nath_source_series(cfg, np.linspace(0.0, 3e-9, 301))
    obs = Observer("near", 0.25e-3, 0.25e-3, 1e-3)
    tr_values = np.linspace(0.4e-9, 2.4e-9, 61)
    tobs = np.asarray([_observer_time(obs, cfg, tr) for tr in tr_values])
    E_ref, B_ref = [], []
    for tr in tr_values:
        fields = evaluate_nath_field(avalanche_state(cfg, float(tr)), obs.position, config=cfg)
        E_ref.append(fields.E_total)
        B_ref.append(fields.B_total)
    numerical = evaluate_waveform(series, [obs], tobs, source_manifest_id="F-R1-test")
    assert all(sample.retarded_time_valid for sample in numerical)
    E_num = np.asarray([sample.E_total for sample in numerical])
    B_num = np.asarray([sample.B_total for sample in numerical])
    E_ref = np.asarray(E_ref)
    B_ref = np.asarray(B_ref)
    for comp in range(3):
        assert waveform_metrics(E_ref[:, comp], E_num[:, comp], tobs)["normalized_L2_error"] < 0.025
    for comp in (0, 1):
        assert waveform_metrics(B_ref[:, comp], B_num[:, comp], tobs)["normalized_L2_error"] < 0.002
    assert np.max(np.abs(B_num[:, 2])) < np.max(np.abs(B_num[:, :2])) * 1e-14


def test_far_observer_main_components_and_distance_scaling() -> None:
    cfg = NathAvalancheConfig()
    series = make_nath_source_series(cfg, np.linspace(0.0, 3e-9, 301))
    obs = Observer("far", 1.0, 1.0, 0.001)
    tr_values = np.linspace(0.4e-9, 2.4e-9, 61)
    tobs = np.asarray([_observer_time(obs, cfg, tr) for tr in tr_values])
    E_ref, B_ref = [], []
    for tr in tr_values:
        fields = evaluate_nath_field(avalanche_state(cfg, float(tr)), obs.position, config=cfg)
        E_ref.append(fields.E_total)
        B_ref.append(fields.B_total)
    numerical = evaluate_waveform(series, [obs], tobs, source_manifest_id="F-R1-test")
    E_num = np.asarray([sample.E_total for sample in numerical])
    B_num = np.asarray([sample.B_total for sample in numerical])
    E_ref = np.asarray(E_ref)
    B_ref = np.asarray(B_ref)
    assert waveform_metrics(E_ref[:, 2], E_num[:, 2], tobs)["normalized_L2_error"] < 0.01
    for comp in (0, 1):
        assert waveform_metrics(B_ref[:, comp], B_num[:, comp], tobs)["normalized_L2_error"] < 0.001

    distances = np.asarray([0.02, 0.05, 0.1, 0.2])
    B_near = []
    B_rad = []
    for R in distances:
        fields = evaluate_nath_field(avalanche_state(cfg, 1.5e-9), np.array([R, R, 1e-3]), config=cfg)
        B_near.append(np.linalg.norm(fields.B_near))
        B_rad.append(np.linalg.norm(fields.B_charge_growth_radiation))
    near_slope = np.polyfit(np.log(distances), np.log(B_near), 1)[0]
    rad_slope = np.polyfit(np.log(distances), np.log(B_rad), 1)[0]
    assert abs(near_slope + 2.0) < 0.08
    assert abs(rad_slope + 1.0) < 0.08
