from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.applicability import (  # noqa: E402
    current_moment_tau_M,
    current_weighted_source_scale,
    em_applicability_audit,
)
from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord  # noqa: E402


def _record(x: np.ndarray, y: np.ndarray, z: np.ndarray, Jz: np.ndarray, *, coord: str = "cartesian") -> SourceRecord:
    n = x.size
    meta = SourceMetadata(
        "test",
        "synthetic",
        "unit",
        0.0,
        coord,
        1.0,
        1.0,
        "none",
        "none",
        "off",
        "synthetic J",
    )
    return SourceRecord(
        meta,
        {
            "level": np.zeros(n, dtype=int),
            "x_center": x,
            "y_center": y,
            "z_center": z,
            "dx": np.ones(n),
            "dy": np.ones(n),
            "dz": np.ones(n),
            "cell_volume": np.ones(n),
            "rho": np.zeros(n),
            "Jx": np.zeros(n),
            "Jy": np.zeros(n),
            "Jz": Jz,
        },
    )


def test_L95_compact_source_translation_invariant() -> None:
    x = np.linspace(-1.0, 1.0, 101)
    J = np.exp(-0.5 * (x / 0.2) ** 2)
    s0 = current_weighted_source_scale(_record(x, np.zeros_like(x), np.zeros_like(x), J))
    s1 = current_weighted_source_scale(_record(x + 5.0, np.zeros_like(x), np.zeros_like(x), J))
    assert s0.status == "OK"
    assert abs(s0.L95_m - s1.L95_m) < 1e-12
    assert abs(s1.centroid_m[0] - 5.0) < 1e-12


def test_L95_weak_tail_is_robust_but_bounding_is_conservative() -> None:
    x = np.r_[np.linspace(-1.0, 1.0, 101), 20.0]
    J = np.r_[np.exp(-0.5 * (x[:-1] / 0.2) ** 2), 1e-8]
    scale = current_weighted_source_scale(_record(x, np.zeros_like(x), np.zeros_like(x), J))
    assert scale.L95_m < 1.0
    assert scale.Lbounding_m > 20.0


def test_zero_and_nonfinite_source_status() -> None:
    x = np.arange(3.0)
    zero = current_weighted_source_scale(_record(x, x * 0.0, x * 0.0, np.zeros(3)))
    assert zero.status == "ZERO_CURRENT_SOURCE"
    try:
        _record(x, x * 0.0, x * 0.0, np.array([1.0, math.nan, 2.0]))
    except ValueError as exc:
        assert "NaN or Inf" in str(exc)
    else:
        raise AssertionError("SourceRecord should reject nonfinite source arrays before F-R2 audit")


def test_axisymmetric_scale_uses_radial_distance_to_axis() -> None:
    r = np.array([0.0, 1.0, 2.0])
    z = np.zeros_like(r)
    scale = current_weighted_source_scale(_record(r, np.zeros_like(r), z, np.ones_like(r), coord="axisymmetric_rz"))
    assert scale.centroid_m[0] == 0.0
    assert scale.L95_m >= 2.0


def test_tau_M_sinusoid_constant_zero_and_nonuniform() -> None:
    omega = 2.0 * math.pi * 1.0e9
    t = np.linspace(0.0, 10.0e-9, 1001)
    M = np.sin(omega * t)
    audit = current_moment_tau_M(t, M)
    assert audit.status == "OK"
    assert abs(audit.tau_M_s - 1.0 / omega) / (1.0 / omega) < 5e-3
    assert abs(audit.omega_consistency_ratio - 1.0) < 0.02

    const = current_moment_tau_M(t, np.ones_like(t))
    assert const.status == "CONSTANT_SIGNAL"
    assert math.isinf(const.tau_M_s)
    zero = current_moment_tau_M(t, np.zeros_like(t))
    assert zero.status == "ZERO_SIGNAL"
    assert math.isinf(zero.tau_M_s)

    tn = np.sort(t + 0.2e-12 * np.sin(np.linspace(0, 20, t.size)))
    nonuniform = current_moment_tau_M(tn, np.sin(omega * tn))
    assert nonuniform.status == "OK"
    assert abs(nonuniform.tau_M_s - 1.0 / omega) / (1.0 / omega) < 8e-3


def test_tau_M_invalid_and_insufficient() -> None:
    assert current_moment_tau_M(np.array([0.0, 1.0]), np.array([1.0, 2.0])).status == "INSUFFICIENT_SAMPLES"
    assert current_moment_tau_M(np.array([0.0, 1.0, 0.5]), np.array([1.0, 2.0, 3.0])).status == "NONFINITE_OR_NONMONOTONE_INPUT"
    assert current_moment_tau_M(np.array([0.0, 1.0, 2.0]), np.array([1.0, math.nan, 3.0])).status == "NONFINITE_OR_NONMONOTONE_INPUT"


def test_epsilon_decision_bands_and_identity() -> None:
    omega = 2.0 * math.pi * 1.0e9
    t = np.linspace(0.0, 10e-9, 1001)
    tau = current_moment_tau_M(t, np.sin(omega * t))
    audit = em_applicability_audit(L95_p95_m=1e-3, L95_max_m=1.2e-3, Lbounding_p95_m=1.5e-3, tau=tau)
    assert audit.epsilon_EM > 0.0
    assert abs(audit.epsilon_EM - audit.two_pi_L_over_lambda) / audit.epsilon_EM < 1e-12
    assert audit.decision in {
        "LOW_EM_FEEDBACK_RISK",
        "REVIEW_SELECTED_FULL_MAXWELL_REFERENCE",
        "SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED",
    }
