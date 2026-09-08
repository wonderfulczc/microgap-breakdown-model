from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.mechanisms import (  # noqa: E402
    compute_mechanism_decomposition,
    local_polynomial_derivative,
)


def _vec_z(z: np.ndarray) -> np.ndarray:
    return np.column_stack((np.zeros_like(z), np.zeros_like(z), z))


def _valid_inner(mask: np.ndarray) -> np.ndarray:
    assert np.count_nonzero(mask) > 0
    return mask


def test_constant_q_constant_v_has_no_charge_or_acceleration_term() -> None:
    t = np.linspace(0.0, 1.0, 31)
    q = np.full_like(t, -2.0)
    z = 3.0 * t + 0.4
    M = _vec_z(q * 3.0)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z), derivative_radius=2)
    mask = _valid_inner(out["derivative_valid"])
    np.testing.assert_allclose(out["dM_charge_evolution_Am_s"][mask], 0.0, atol=1e-11)
    np.testing.assert_allclose(out["dM_head_acceleration_Am_s"][mask], 0.0, atol=1e-10)
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_changing_q_constant_v_recovers_head_charge_evolution() -> None:
    t = np.linspace(0.0, 1.0, 41)
    q = 1.0 + 2.0 * t
    v0 = -4.0
    z = v0 * t
    M = _vec_z(q * v0)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z), derivative_radius=2)
    mask = out["derivative_valid"]
    np.testing.assert_allclose(out["dM_charge_evolution_Am_s"][mask, 2], 2.0 * v0, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(out["dM_head_acceleration_Am_s"][mask, 2], 0.0, atol=1e-9)
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_constant_q_accelerating_head_recovers_qa() -> None:
    t = np.linspace(0.0, 1.0, 51)
    q = np.full_like(t, 3.0)
    a0 = 5.0
    z = 0.5 * a0 * t * t
    M = _vec_z(q * a0 * t)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z), derivative_radius=2)
    mask = out["derivative_valid"]
    np.testing.assert_allclose(out["dM_charge_evolution_Am_s"][mask, 2], 0.0, atol=1e-9)
    np.testing.assert_allclose(out["dM_head_acceleration_Am_s"][mask, 2], q[0] * a0, rtol=1e-10, atol=1e-9)
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_simultaneous_q_and_v_product_rule_residual_is_reported() -> None:
    t = np.linspace(-1.0, 1.0, 61)
    q = 1.0 + 0.2 * t
    z = 0.3 * t + 0.4 * t * t
    v = 0.3 + 0.8 * t
    M = _vec_z(q * v)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z), derivative_radius=2)
    mask = out["derivative_valid"]
    assert np.nanmax(np.linalg.norm(out["product_rule_residual_Am_s"][mask], axis=1)) < 1e-12
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_known_redistribution_term_is_recovered() -> None:
    t = np.linspace(0.0, 1.0, 51)
    q = np.full_like(t, 2.0)
    z = t
    Mredis = 0.7 * t * t
    M = _vec_z(q + Mredis)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z), derivative_radius=2)
    mask = out["derivative_valid"]
    np.testing.assert_allclose(out["dM_redistribution_Am_s"][mask, 2], 1.4 * t[mask], rtol=1e-10, atol=1e-10)
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_two_heads_sum_support() -> None:
    t = np.linspace(0.0, 1.0, 41)
    q = np.column_stack((np.full_like(t, 2.0), np.full_like(t, -1.0)))
    r = np.zeros((t.size, 2, 3))
    r[:, 0, 2] = 2.0 * t
    r[:, 1, 2] = -3.0 * t
    M = _vec_z(q[:, 0] * 2.0 + q[:, 1] * -3.0)
    out = compute_mechanism_decomposition(t, M, q, r, derivative_radius=2)
    mask = out["derivative_valid"]
    np.testing.assert_allclose(out["M_head_Am"][mask, 2], 7.0, rtol=1e-10, atol=1e-10)
    assert out["metrics"].normalized_rms_closure < 1e-10


def test_nonuniform_timestamp_derivative() -> None:
    t = np.sort(np.linspace(0.0, 1.0, 41) + 1e-3 * np.sin(np.linspace(0.0, 8.0, 41)))
    y = 2.0 + 3.0 * t + 4.0 * t * t
    dy, ok = local_polynomial_derivative(t, y, radius=2)
    np.testing.assert_allclose(dy[ok], 3.0 + 8.0 * t[ok], rtol=1e-10, atol=1e-10)


def test_invalid_head_gap_blocks_derivative_across_gap() -> None:
    t = np.linspace(0.0, 1.0, 21)
    q = np.ones_like(t)
    z = t
    valid = np.ones_like(t, dtype=bool)
    valid[10] = False
    out = compute_mechanism_decomposition(t, _vec_z(np.ones_like(t)), q, _vec_z(z), head_valid=valid)
    assert not out["derivative_valid"][10]
    assert not out["derivative_valid"][9]
    assert not out["derivative_valid"][11]


def test_translated_head_trajectory_preserves_velocity_and_mechanism() -> None:
    t = np.linspace(0.0, 1.0, 31)
    q = np.full_like(t, -1.0)
    r0 = _vec_z(2.0 * t)
    r1 = r0 + np.array([5.0, -2.0, 3.0])
    M = _vec_z(q * 2.0)
    a = compute_mechanism_decomposition(t, M, q, r0)
    b = compute_mechanism_decomposition(t, M, q, r1)
    mask = a["derivative_valid"] & b["derivative_valid"]
    np.testing.assert_allclose(a["velocity_m_s"][mask, 0, :], b["velocity_m_s"][mask, 0, :], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(
        a["dM_charge_evolution_Am_s"][mask],
        b["dM_charge_evolution_Am_s"][mask],
        rtol=1e-10,
        atol=1e-10,
    )


def test_signed_charge_and_velocity_sign_reversal() -> None:
    t = np.linspace(0.0, 1.0, 31)
    q = np.full_like(t, -2.0)
    z = -3.0 * t
    M = _vec_z(q * -3.0)
    out = compute_mechanism_decomposition(t, M, q, _vec_z(z))
    mask = out["derivative_valid"]
    np.testing.assert_allclose(out["M_head_Am"][mask, 2], 6.0, rtol=1e-10, atol=1e-10)
    assert math.isfinite(out["metrics"].normalized_rms_closure)
