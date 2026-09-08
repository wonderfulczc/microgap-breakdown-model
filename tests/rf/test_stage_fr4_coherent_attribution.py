from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.coherent_attribution import (  # noqa: E402
    coherent_attribution,
    cwt_linearity_error,
    mechanism_radiative_fields,
    uniform_resample,
)
from streamer_rf.rf.spectral.cwt import morlet_cwt  # noqa: E402


def _mask(shape: tuple[int, int]) -> np.ndarray:
    return np.ones(shape, dtype=bool)


def test_single_mechanism_eta_is_one() -> None:
    w = np.ones((3, 4), dtype=complex)
    out = coherent_attribution({"m1": w}, w, _mask(w.shape))
    assert out.status == "OK"
    assert abs(out.eta["m1"] - 1.0) < 1e-14
    assert out.eta_sum_error < 1e-14


def test_two_equal_in_phase_mechanisms_split_half() -> None:
    w1 = np.ones((2, 5), dtype=complex)
    w2 = np.ones((2, 5), dtype=complex)
    out = coherent_attribution({"m1": w1, "m2": w2}, w1 + w2, _mask(w1.shape))
    assert abs(out.eta["m1"] - 0.5) < 1e-14
    assert abs(out.eta["m2"] - 0.5) < 1e-14


def test_destructive_interference_is_not_clamped() -> None:
    w1 = np.ones((1, 8), dtype=complex)
    w2 = -0.5 * np.ones((1, 8), dtype=complex)
    out = coherent_attribution({"m1": w1, "m2": w2}, w1 + w2, _mask(w1.shape))
    assert abs(out.eta["m1"] - 2.0) < 1e-14
    assert abs(out.eta["m2"] + 1.0) < 1e-14
    assert out.interference_role["m2"] == "DESTRUCTIVE_NET"
    assert out.eta_sum_error < 1e-14


def test_quadrature_complex_mechanisms() -> None:
    w1 = np.ones((2, 3), dtype=complex)
    w2 = 1j * np.ones((2, 3), dtype=complex)
    rec = w1 + w2
    out = coherent_attribution({"m1": w1, "m2": w2}, rec, _mask(rec.shape))
    assert abs(out.eta["m1"] - 0.5) < 1e-14
    assert abs(out.eta["m2"] - 0.5) < 1e-14


def test_three_mechanisms_sum_to_one() -> None:
    w1 = 2.0 * np.ones((3, 3), dtype=complex)
    w2 = -0.25 * np.ones((3, 3), dtype=complex)
    w3 = (0.5 + 0.5j) * np.ones((3, 3), dtype=complex)
    out = coherent_attribution({"m1": w1, "m2": w2, "m3": w3}, w1 + w2 + w3, _mask(w1.shape))
    assert out.eta_sum_error < 1e-14
    assert out.eta["m2"] < 0.0


def test_stage_frequency_and_trust_masks_exclude_coefficients() -> None:
    w1 = np.ones((4, 6), dtype=complex)
    w2 = np.zeros_like(w1)
    mask = np.zeros((4, 6), dtype=bool)
    mask[1:3, 2:5] = True
    out = coherent_attribution({"m1": w1, "m2": w2}, w1, mask)
    assert out.trusted_count == 6
    assert out.trusted_fraction == 6 / 24
    assert abs(out.eta["m1"] - 1.0) < 1e-14


def test_zero_power_attribution_is_invalid() -> None:
    w = np.zeros((2, 4), dtype=complex)
    out = coherent_attribution({"m1": w}, w, _mask(w.shape))
    assert out.status == "ZERO_OR_NEAR_ZERO_RECONSTRUCTED_POWER"


def test_no_trusted_coefficients_is_invalid() -> None:
    w = np.ones((2, 4), dtype=complex)
    out = coherent_attribution({"m1": w}, w, np.zeros_like(w, dtype=bool))
    assert out.status == "NO_TRUSTED_COEFFICIENTS"
    assert out.trusted_count == 0


def test_near_zero_power_floor_is_invalid() -> None:
    w = 1e-30 * np.ones((2, 4), dtype=complex)
    out = coherent_attribution({"m1": w}, w, _mask(w.shape), p_floor=1e-20)
    assert out.status == "ZERO_OR_NEAR_ZERO_RECONSTRUCTED_POWER"


def test_eta_is_invariant_to_amplitude_scaling_and_common_phase() -> None:
    w1 = np.ones((2, 4), dtype=complex)
    w2 = 0.25j * np.ones((2, 4), dtype=complex)
    base = coherent_attribution({"m1": w1, "m2": w2}, w1 + w2, _mask(w1.shape))
    phase = 3.0 * np.exp(0.7j)
    shifted = coherent_attribution({"m1": phase * w1, "m2": phase * w2}, phase * (w1 + w2), _mask(w1.shape))
    assert abs(base.eta["m1"] - shifted.eta["m1"]) < 1e-14
    assert abs(base.eta["m2"] - shifted.eta["m2"]) < 1e-14


def test_cwt_linearity_for_existing_morlet() -> None:
    t = np.linspace(0.0, 10e-9, 512)
    f = np.geomspace(100e6, 2e9, 16)
    a = np.sin(2.0 * np.pi * 300e6 * t)
    b = 0.3 * np.sin(2.0 * np.pi * 900e6 * t + 0.4)
    ca = morlet_cwt(t, a, f, n_cycles=5.0)
    cb = morlet_cwt(t, b, f, n_cycles=5.0)
    cab = morlet_cwt(t, a + b, f, n_cycles=5.0)
    err = cwt_linearity_error(ca.coefficients + cb.coefficients, cab.coefficients, cab.coi_mask)
    assert err < 1e-12


def test_uniform_resample_and_mechanism_field_linearity() -> None:
    t = np.array([0.0, 1.0, 2.1, 3.0, 4.2, 5.0])
    y = np.column_stack((t, t * t))
    tu, yu = uniform_resample(t, y, np.ones(t.size, dtype=bool))
    assert tu.size >= 3
    np.testing.assert_allclose(yu[:, 0], tu, rtol=1e-12, atol=1e-12)

    obs = np.array([1.0, 0.0, 0.0])
    d = {
        "a": np.column_stack((np.zeros(4), np.zeros(4), np.arange(4.0))),
        "b": np.column_stack((np.zeros(4), np.zeros(4), 2.0 * np.arange(4.0))),
    }
    fields = mechanism_radiative_fields(d, obs)
    total = mechanism_radiative_fields({"total": d["a"] + d["b"]}, obs)["total"]
    np.testing.assert_allclose(fields["a"] + fields["b"], total, rtol=1e-14, atol=1e-30)
