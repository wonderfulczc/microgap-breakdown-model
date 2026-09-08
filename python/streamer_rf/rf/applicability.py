from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from streamer_rf.rf.jefimenko.constants import C0
from streamer_rf.rf.source.schema import SourceRecord
from streamer_rf.rf.spectral.fft import compute_one_sided_spectrum, parseval_spectral_energy


@dataclass(frozen=True)
class SourceScale:
    time_s: float
    total_current_weight_Am: float
    centroid_m: tuple[float, float, float]
    L95_m: float
    Lbounding_m: float
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TimeScaleAudit:
    tau_M_s: float
    inverse_tau_rad_s: float
    f_equiv_Hz: float
    omega_rms_rad_s: float
    omega_consistency_ratio: float
    status: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class EMAudit:
    L_source_m: float
    tau_M_s: float
    epsilon_EM: float
    epsilon_EM_conservative: float
    f_equiv_Hz: float
    lambda_equiv_m: float
    L_over_lambda: float
    two_pi_L_over_lambda: float
    decision: str
    decision_basis: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _relative_percentile(values: np.ndarray, weights: np.ndarray, fraction: float) -> float:
    if not 0.0 < fraction <= 1.0:
        raise ValueError("fraction must be in (0, 1]")
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    total = float(np.sum(w, dtype=np.float64))
    if total <= 0.0:
        raise ValueError("weights must have positive sum")
    cdf = np.cumsum(w, dtype=np.float64) / total
    return float(v[min(int(np.searchsorted(cdf, fraction, side="left")), v.size - 1)])


def current_weighted_source_scale(record: SourceRecord, *, fraction: float = 0.95) -> SourceScale:
    J = np.column_stack((record.columns["Jx"], record.columns["Jy"], record.columns["Jz"]))
    if not np.all(np.isfinite(J)):
        return SourceScale(record.metadata.time_s, math.nan, (math.nan, math.nan, math.nan), math.nan, math.nan, "NONFINITE_SOURCE")
    volume = record.columns["cell_volume"]
    weight = np.linalg.norm(J, axis=1) * volume
    active = weight > 0.0
    total = float(np.sum(weight[active], dtype=np.float64))
    if total <= 0.0 or not math.isfinite(total):
        return SourceScale(record.metadata.time_s, 0.0, (math.nan, math.nan, math.nan), math.nan, math.nan, "ZERO_CURRENT_SOURCE")

    x = record.columns["x_center"]
    y = record.columns["y_center"]
    z = record.columns["z_center"]
    if record.metadata.coordinate_system.startswith("axisymmetric"):
        zc = float(np.sum(weight * z, dtype=np.float64) / total)
        centroid = np.asarray([0.0, 0.0, zc], dtype=float)
        distances = np.sqrt(x * x + (z - zc) ** 2)
        r_active = x[active]
        z_active = z[active]
        radial_extent = float(np.max(r_active))
        axial_extent = float(np.max(z_active) - np.min(z_active))
        Lbounding = float(2.0 * math.sqrt(radial_extent**2 + (0.5 * axial_extent) ** 2))
    else:
        xyz = np.column_stack((x, y, z))
        centroid = np.sum(xyz * weight[:, None], axis=0, dtype=np.float64) / total
        distances = np.linalg.norm(xyz - centroid[None, :], axis=1)
        lo = np.min(xyz[active], axis=0)
        hi = np.max(xyz[active], axis=0)
        Lbounding = float(np.linalg.norm(hi - lo))
    r95 = _relative_percentile(distances[active], weight[active], fraction)
    return SourceScale(
        time_s=float(record.metadata.time_s),
        total_current_weight_Am=total,
        centroid_m=(float(centroid[0]), float(centroid[1]), float(centroid[2])),
        L95_m=2.0 * r95,
        Lbounding_m=Lbounding,
        status="OK",
    )


def derivative_nonuniform(times_s: np.ndarray, values: np.ndarray) -> np.ndarray:
    t = np.asarray(times_s, dtype=float)
    y = np.asarray(values, dtype=float)
    if t.ndim != 1 or y.shape[0] != t.size:
        raise ValueError("time and values dimensions do not match")
    if t.size < 3:
        raise ValueError("at least three samples are required")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("times must be strictly increasing")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)):
        raise ValueError("time and values must be finite")
    return np.gradient(y, t, axis=0, edge_order=2)


def current_moment_tau_M(times_s: np.ndarray, M: np.ndarray) -> TimeScaleAudit:
    t = np.asarray(times_s, dtype=float)
    m = np.asarray(M, dtype=float)
    if m.ndim == 1:
        m = m[:, None]
    if t.size < 3:
        return TimeScaleAudit(math.nan, math.nan, math.nan, math.nan, math.nan, "INSUFFICIENT_SAMPLES")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(m)) or np.any(np.diff(t) <= 0.0):
        return TimeScaleAudit(math.nan, math.nan, math.nan, math.nan, math.nan, "NONFINITE_OR_NONMONOTONE_INPUT")
    m2 = np.sum(m * m, axis=1)
    if float(np.max(m2)) == 0.0:
        return TimeScaleAudit(math.inf, 0.0, 0.0, 0.0, math.nan, "ZERO_SIGNAL")
    dmdt = derivative_nonuniform(t, m)
    dm2 = np.sum(dmdt * dmdt, axis=1)
    numerator = float(np.trapezoid(m2, t))
    denominator = float(np.trapezoid(dm2, t))
    span = float(t[-1] - t[0])
    if denominator <= max(numerator / max(span, 1e-300) ** 2, 1e-300) * 1e-24:
        return TimeScaleAudit(math.inf, 0.0, 0.0, 0.0, math.nan, "CONSTANT_SIGNAL")
    tau = math.sqrt(numerator / denominator)
    omega_rms = current_moment_omega_rms(t, m)
    inv = 1.0 / tau
    return TimeScaleAudit(
        tau_M_s=tau,
        inverse_tau_rad_s=inv,
        f_equiv_Hz=inv / (2.0 * math.pi),
        omega_rms_rad_s=omega_rms,
        omega_consistency_ratio=omega_rms / inv if inv > 0.0 and math.isfinite(omega_rms) else math.nan,
        status="OK",
    )


def current_moment_omega_rms(times_s: np.ndarray, M: np.ndarray) -> float:
    t = np.asarray(times_s, dtype=float)
    m = np.asarray(M, dtype=float)
    if m.ndim == 1:
        m = m[:, None]
    dt = np.diff(t)
    if not np.allclose(dt, dt[0], rtol=1e-6, atol=max(abs(dt[0]) * 1e-9, 1e-30)):
        uniform_t = np.linspace(float(t[0]), float(t[-1]), t.size)
        vals = np.column_stack([np.interp(uniform_t, t, m[:, i]) for i in range(m.shape[1])])
        t = uniform_t
        m = vals
    y = m
    if np.max(np.abs(y)) == 0.0:
        return 0.0
    dydt = derivative_nonuniform(t, y)
    spectrum_m = compute_one_sided_spectrum(t, y, window="rectangular", apply_window_energy_correction=False)
    spectrum_dm = compute_one_sided_spectrum(t, dydt, window="rectangular", apply_window_energy_correction=False)
    denom = float(np.sum(parseval_spectral_energy(spectrum_m)))
    numer = float(np.sum(parseval_spectral_energy(spectrum_dm)))
    if denom <= 0.0 or numer <= 0.0:
        return 0.0
    return math.sqrt(numer / denom)


def em_applicability_audit(
    *,
    L95_p95_m: float,
    L95_max_m: float,
    Lbounding_p95_m: float,
    tau: TimeScaleAudit,
    conservative_tau_s: float | None = None,
) -> EMAudit:
    if tau.status != "OK" or not math.isfinite(tau.tau_M_s) or tau.tau_M_s <= 0.0:
        return EMAudit(
            L_source_m=float(L95_p95_m),
            tau_M_s=tau.tau_M_s,
            epsilon_EM=math.nan,
            epsilon_EM_conservative=math.nan,
            f_equiv_Hz=math.nan,
            lambda_equiv_m=math.nan,
            L_over_lambda=math.nan,
            two_pi_L_over_lambda=math.nan,
            decision="REVIEW_SELECTED_FULL_MAXWELL_REFERENCE",
            decision_basis=f"tau_M status={tau.status}",
        )
    L_source = float(L95_p95_m)
    eps = L_source / (C0 * tau.tau_M_s)
    tau_cons = conservative_tau_s if conservative_tau_s is not None and conservative_tau_s > 0.0 else tau.tau_M_s
    L_cons = max(float(L95_max_m), float(Lbounding_p95_m))
    eps_cons = L_cons / (C0 * tau_cons)
    freq = tau.f_equiv_Hz
    wavelength = C0 / freq if freq > 0.0 else math.inf
    L_over_lambda = L_source / wavelength if math.isfinite(wavelength) else 0.0
    worst = max(eps, eps_cons)
    if worst <= 0.1:
        decision = "LOW_EM_FEEDBACK_RISK"
    elif worst <= 0.3:
        decision = "REVIEW_SELECTED_FULL_MAXWELL_REFERENCE"
    else:
        decision = "SELECTED_FULL_MAXWELL_REFERENCE_REQUIRED"
    return EMAudit(
        L_source_m=L_source,
        tau_M_s=tau.tau_M_s,
        epsilon_EM=float(eps),
        epsilon_EM_conservative=float(eps_cons),
        f_equiv_Hz=float(freq),
        lambda_equiv_m=float(wavelength),
        L_over_lambda=float(L_over_lambda),
        two_pi_L_over_lambda=float(2.0 * math.pi * L_over_lambda),
        decision=decision,
        decision_basis="decision uses max(primary epsilon_EM, conservative epsilon_EM)",
    )
