from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.stage4 import (
    align_current_moments,
    collision_pulse_metrics,
    derivative_products,
    detect_collision_event,
    esd_per_hz_from_derivative_spectrum,
    output_sampling_sensitivity,
    radiation_products,
)


def test_collision_event_detection_uses_relative_bridge_and_field_drop():
    t = np.linspace(0, 3e-9, 200)
    bridge = 1e16 + 2e20 / (1 + np.exp(-(t - 1.9e-9) / 0.05e-9))
    egap = 5e6 + 7e6 * np.exp(-((t - 0.9e-9) / 0.35e-9) ** 2)
    egap[t > 1.9e-9] *= 0.5
    d = pd.DataFrame({"time": t, "d_head": 3e-3 - 1.5e-3 * np.minimum(t / 2e-9, 1), "bridge_mean_ne": bridge, "E_gap_max": egap})
    ev = detect_collision_event(d)
    assert ev.bridge_status == "PASS"
    assert ev.field_status == "PASS"
    assert ev.collision_status == "PASS"
    assert ev.bridge_threshold_m_3 > ev.bridge_pre_median_m_3


def test_three_case_alignment_and_delta_current():
    t = np.linspace(0, 2e-9, 51)
    left = pd.DataFrame({"time": t, "I_CM_drift": 0.1 * t, "I_CM_electron": 0.2 * t})
    right = pd.DataFrame({"time": t + 1e-15, "I_CM_drift": 0.3 * t, "I_CM_electron": 0.4 * t})
    pulse = 1e-3 * np.exp(-((t - 1e-9) / 0.1e-9) ** 2)
    collision = pd.DataFrame({"time": t, "I_CM_drift": 0.4 * t + pulse, "I_CM_electron": 0.6 * t + 2 * pulse})
    out = align_current_moments(collision, left, right)
    assert np.max(out.Delta_I_CM_drift) > 9e-4
    assert np.max(out.Delta_I_CM_electron) > 1.8e-3


def test_delta_pulse_causality_and_derivatives():
    t = np.linspace(0, 2e-9, 101)
    y = 1e-3 * np.exp(-((t - 1e-9) / 0.08e-9) ** 2)
    d = pd.DataFrame({"time": t, "Delta_I_CM_drift": y})
    pulse = collision_pulse_metrics(d, 1e-9)
    deriv, summary = derivative_products(d, 1e-9)
    assert pulse["causality_status"] == "PASS"
    assert summary["peak_relative_difference"] < 0.2
    assert {"dDeltaI_dt_local_poly", "dDeltaI_dt_uniform_five_point", "dDeltaI_dt_spectral"}.issubset(deriv.columns)


def test_radiation_products_and_band_trust_labels():
    t = np.linspace(0, 2e-9, 128)
    y = 1e-3 * np.exp(-((t - 1e-9) / 0.1e-9) ** 2)
    d = pd.DataFrame({"time": t, "Delta_I_CM_drift": y})
    spec, esd, bands, summary = radiation_products(d, 1e-9)
    assert summary["f_trust_Hz"] > 0
    assert np.all(esd.ESD_J_per_Hz >= 0)
    assert set(bands.band) == {"VHF", "UHF", "SHF"}
    assert set(bands.trusted_status).issubset({"trusted", "partial", "untrusted"})


def test_stage1_esd_formula_reused_for_arbitrary_derivative_spectrum():
    z = np.array([0.0, 1.0, 2.0])
    e = esd_per_hz_from_derivative_spectrum(z)
    assert e[0] == 0.0
    assert np.isclose(e[2] / e[1], 4.0)


def test_output_sampling_sensitivity_is_measured_not_fixed_pass():
    t = np.linspace(0, 2e-9, 101)
    y = 1e-3 * np.exp(-((t - 1e-9) / 0.08e-9) ** 2)
    d = pd.DataFrame({"time": t, "Delta_I_CM_drift": y})
    out = output_sampling_sensitivity(d, 1e-9)
    assert len(out) == 3
    assert np.all(out.Delta_I_peak_rel_diff >= 0)
