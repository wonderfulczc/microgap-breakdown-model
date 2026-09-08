from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np


@dataclass(frozen=True)
class MechanismMetrics:
    sample_count: int
    valid_count: int
    rms_D_total: float
    rms_head_charge_evolution: float
    rms_head_acceleration: float
    rms_redistribution: float
    rms_closure_residual: float
    normalized_rms_closure: float
    peak_normalized_closure: float
    product_rule_normalized_rms: float
    corr_charge_evolution: float
    corr_head_acceleration: float
    corr_redistribution: float
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MechanismStatus:
    head_charge_evolution: str
    head_acceleration: str
    current_moment_redistribution: str
    basis: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def local_polynomial_derivative(
    times_s: np.ndarray,
    values: np.ndarray,
    *,
    valid: np.ndarray | None = None,
    radius: int = 2,
    degree: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Local nonuniform-time polynomial derivative at each input sample.

    The stencil is centered on the target sample. Samples whose whole local
    stencil is not valid are returned as NaN/False; no extrapolation across
    invalid gaps is performed.
    """

    t = np.asarray(times_s, dtype=float)
    y = np.asarray(values, dtype=float)
    scalar = y.ndim == 1
    if scalar:
        y = y[:, None]
    if t.ndim != 1 or y.shape[0] != t.size:
        raise ValueError("time and values dimensions do not match")
    if t.size < 2 * radius + 1:
        out = np.full_like(y, np.nan, dtype=float)
        return (out[:, 0] if scalar else out), np.zeros(t.size, dtype=bool)
    if radius < 1 or degree < 1 or 2 * radius + 1 < degree + 1:
        raise ValueError("invalid derivative stencil")
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0.0):
        raise ValueError("times must be finite and strictly increasing")
    base_valid = np.all(np.isfinite(y), axis=1)
    if valid is not None:
        v = np.asarray(valid, dtype=bool)
        if v.shape != (t.size,):
            raise ValueError("valid mask has wrong shape")
        base_valid &= v

    out = np.full(y.shape, np.nan, dtype=float)
    ok = np.zeros(t.size, dtype=bool)
    for i in range(radius, t.size - radius):
        idx = np.arange(i - radius, i + radius + 1)
        if not np.all(base_valid[idx]):
            continue
        x = t[idx] - t[i]
        x_scale = float(np.max(np.abs(x)))
        if not (x_scale > 0.0):
            continue
        xs = x / x_scale
        if np.linalg.matrix_rank(np.vander(xs, degree + 1, increasing=True)) < degree + 1:
            continue
        for comp in range(y.shape[1]):
            coeff = np.polynomial.polynomial.polyfit(xs, y[idx, comp], degree)
            out[i, comp] = coeff[1] / x_scale
        ok[i] = True
    return (out[:, 0] if scalar else out), ok


def _norm_rms(v: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return math.nan
    n = np.linalg.norm(v[mask], axis=1)
    return float(math.sqrt(np.mean(n * n)))


def _peak_normalized(v: np.ndarray, ref: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return math.nan
    denom = float(np.max(np.linalg.norm(ref[mask], axis=1)))
    if denom <= 0.0 or not math.isfinite(denom):
        return math.nan
    return float(np.max(np.linalg.norm(v[mask], axis=1)) / denom)


def _corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    if np.count_nonzero(mask) < 3:
        return math.nan
    aa = np.linalg.norm(a[mask], axis=1)
    bb = np.linalg.norm(b[mask], axis=1)
    if np.std(aa) == 0.0 or np.std(bb) == 0.0:
        return math.nan
    return float(np.corrcoef(aa, bb)[0, 1])


def compute_mechanism_decomposition(
    times_s: np.ndarray,
    M_total_Am: np.ndarray,
    q_head_C: np.ndarray,
    head_position_m: np.ndarray,
    *,
    head_valid: np.ndarray | None = None,
    derivative_radius: int = 2,
    derivative_degree: int = 2,
) -> dict[str, np.ndarray | MechanismMetrics]:
    t = np.asarray(times_s, dtype=float)
    M = np.asarray(M_total_Am, dtype=float)
    q = np.asarray(q_head_C, dtype=float)
    r = np.asarray(head_position_m, dtype=float)
    if M.ndim == 1:
        M = M[:, None]
    if M.shape[1] == 1:
        M = np.column_stack((np.zeros(t.size), np.zeros(t.size), M[:, 0]))
    if q.ndim == 1:
        q = q[:, None]
    if r.ndim == 2:
        r = r[:, None, :]
    if t.ndim != 1 or M.shape != (t.size, 3) or q.shape[0] != t.size or r.shape[:2] != q.shape or r.shape[2] != 3:
        raise ValueError("mechanism input dimensions do not match")
    if head_valid is None:
        hv = np.all(np.isfinite(q), axis=1) & np.all(np.isfinite(r), axis=(1, 2))
    else:
        hv = np.asarray(head_valid, dtype=bool)
        if hv.shape != (t.size,):
            raise ValueError("head_valid mask has wrong shape")
        hv &= np.all(np.isfinite(q), axis=1) & np.all(np.isfinite(r), axis=(1, 2))

    v = np.full_like(r, np.nan, dtype=float)
    a = np.full_like(r, np.nan, dtype=float)
    qdot = np.full_like(q, np.nan, dtype=float)
    velocity_ok = np.zeros((t.size, q.shape[1]), dtype=bool)
    accel_ok = np.zeros((t.size, q.shape[1]), dtype=bool)
    qdot_ok = np.zeros((t.size, q.shape[1]), dtype=bool)

    for k in range(q.shape[1]):
        v[:, k, :], velocity_ok[:, k] = local_polynomial_derivative(
            t, r[:, k, :], valid=hv, radius=derivative_radius, degree=derivative_degree
        )
        a[:, k, :], accel_ok[:, k] = local_polynomial_derivative(
            t, v[:, k, :], valid=velocity_ok[:, k], radius=derivative_radius, degree=derivative_degree
        )
        qdot[:, k], qdot_ok[:, k] = local_polynomial_derivative(
            t, q[:, k], valid=hv, radius=derivative_radius, degree=derivative_degree
        )

    M_head = np.nansum(q[:, :, None] * v, axis=1)
    M_head[~np.any(velocity_ok, axis=1)] = np.nan
    M_redis = M - M_head
    D_total, D_total_ok = local_polynomial_derivative(
        t, M, valid=np.all(np.isfinite(M), axis=1), radius=derivative_radius, degree=derivative_degree
    )
    D_head, D_head_ok = local_polynomial_derivative(
        t, M_head, valid=np.all(np.isfinite(M_head), axis=1), radius=derivative_radius, degree=derivative_degree
    )
    R, R_ok = local_polynomial_derivative(
        t, M_redis, valid=np.all(np.isfinite(M_redis), axis=1), radius=derivative_radius, degree=derivative_degree
    )

    G = np.nansum(qdot[:, :, None] * v, axis=1)
    A = np.nansum(q[:, :, None] * a, axis=1)
    ga_ok = np.all(qdot_ok & velocity_ok & accel_ok, axis=1)
    valid = D_total_ok & D_head_ok & R_ok & ga_ok & np.all(np.isfinite(G), axis=1) & np.all(np.isfinite(A), axis=1)
    reconstructed = G + A + R
    closure = D_total - reconstructed
    product_rule = D_head - (G + A)
    denom = _norm_rms(D_total, valid)
    closure_rms = _norm_rms(closure, valid)
    product_rms = _norm_rms(product_rule, valid)
    term_scale = max(
        _norm_rms(G, valid) if np.any(valid) else math.nan,
        _norm_rms(A, valid) if np.any(valid) else math.nan,
        _norm_rms(R, valid) if np.any(valid) else math.nan,
        denom if math.isfinite(denom) else math.nan,
        1.0,
    )
    if denom and math.isfinite(denom) and denom > 1e-12 * term_scale:
        norm_closure = closure_rms / denom
        norm_product = product_rms / denom
    elif closure_rms <= 1e-10 * term_scale:
        norm_closure = 0.0
        norm_product = 0.0 if product_rms <= 1e-10 * term_scale else math.inf
    else:
        norm_closure = math.inf
        norm_product = math.inf
    metrics = MechanismMetrics(
        sample_count=int(t.size),
        valid_count=int(np.count_nonzero(valid)),
        rms_D_total=denom,
        rms_head_charge_evolution=_norm_rms(G, valid),
        rms_head_acceleration=_norm_rms(A, valid),
        rms_redistribution=_norm_rms(R, valid),
        rms_closure_residual=closure_rms,
        normalized_rms_closure=norm_closure,
        peak_normalized_closure=_peak_normalized(closure, D_total, valid),
        product_rule_normalized_rms=norm_product,
        corr_charge_evolution=_corr(G, D_total, valid),
        corr_head_acceleration=_corr(A, D_total, valid),
        corr_redistribution=_corr(R, D_total, valid),
        status="OK" if np.count_nonzero(valid) >= 3 else "INSUFFICIENT_VALID_SAMPLES",
    )
    return {
        "M_total_Am": M,
        "velocity_m_s": v,
        "acceleration_m_s2": a,
        "qdot_C_s": qdot,
        "M_head_Am": M_head,
        "M_redis_Am": M_redis,
        "dM_total_Am_s": D_total,
        "dM_charge_evolution_Am_s": G,
        "dM_head_acceleration_Am_s": A,
        "dM_redistribution_Am_s": R,
        "closure_residual_Am_s": closure,
        "product_rule_residual_Am_s": product_rule,
        "derivative_valid": valid,
        "metrics": metrics,
    }


def classify_mechanism_status(
    metrics: MechanismMetrics,
    *,
    derivative_relative_changes: dict[str, float] | None = None,
    segmentation_relative_changes: dict[str, float] | None = None,
    field_correlation: float | None = None,
) -> MechanismStatus:
    if metrics.status != "OK" or metrics.valid_count < 3 or not math.isfinite(metrics.normalized_rms_closure):
        return MechanismStatus("NOT_RESOLVED", "NOT_RESOLVED", "NOT_RESOLVED", "insufficient valid closure samples")
    if metrics.normalized_rms_closure > 0.5:
        return MechanismStatus(
            "NOT_RESOLVED",
            "NOT_RESOLVED",
            "NOT_RESOLVED",
            f"closure={metrics.normalized_rms_closure:.3g} exceeds reduced-order interpretation range",
        )
    base = "RESOLVED" if metrics.normalized_rms_closure <= 0.05 else "INTERPRET_WITH_CAUTION"
    basis = f"closure={metrics.normalized_rms_closure:.3g}"
    if field_correlation is not None and math.isfinite(field_correlation) and field_correlation < 0.95:
        base = "INTERPRET_WITH_CAUTION"
        basis += f"; current-moment/Jefimenko correlation={field_correlation:.3g}"
    accel_status = base
    charge_status = base
    redis_status = base
    if derivative_relative_changes:
        accel_change = derivative_relative_changes.get("rms_head_acceleration", math.nan)
        if math.isfinite(accel_change) and accel_change > 10.0:
            accel_status = "NOT_RESOLVED"
            basis += "; HEAD_ACCELERATION derivative-sensitive"
        max_change = max((v for v in derivative_relative_changes.values() if math.isfinite(v)), default=0.0)
        if max_change > 0.5 and base == "RESOLVED":
            charge_status = redis_status = "INTERPRET_WITH_CAUTION"
            basis += "; derivative sensitivity changes amplitudes"
    if segmentation_relative_changes:
        max_seg = max((v for v in segmentation_relative_changes.values() if math.isfinite(v)), default=0.0)
        if max_seg > 0.5 and base == "RESOLVED":
            charge_status = accel_status = redis_status = "INTERPRET_WITH_CAUTION"
            basis += "; segmentation sensitivity changes amplitudes"
    return MechanismStatus(charge_status, accel_status, redis_status, basis)


def relative_change(reference: float, candidate: float) -> float:
    if not math.isfinite(reference) or not math.isfinite(candidate):
        return math.nan
    denom = max(abs(reference), abs(candidate), 1e-300)
    return abs(candidate - reference) / denom
