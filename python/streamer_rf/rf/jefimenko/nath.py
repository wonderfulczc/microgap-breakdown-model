from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from streamer_rf.rf.jefimenko.constants import C0, EPS0, K_E
from streamer_rf.rf.jefimenko.observer import Observer
from streamer_rf.rf.source.integrals import current_moment, total_charge
from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord, SourceSeries


@dataclass(frozen=True)
class NathAvalancheConfig:
    case_id: str = "F-R1-Nath-straight-avalanche"
    q0_C: float = 1.0e-15
    growth_rate_s: float = 5.0e8
    velocity_m_s: float = 1.0e5
    start_z_m: float = 0.0
    start_time_s: float = 0.0
    source_width_m: float = 8.0e-6
    z_min_m: float = -4.0e-5
    z_max_m: float = 4.5e-4
    dz_m: float = 2.0e-6
    transverse_size_m: float = 2.0e-6


@dataclass(frozen=True)
class NathAvalancheState:
    t_s: float
    q_C: float
    dq_dt_A: float
    position_m: np.ndarray
    velocity_m_s: np.ndarray
    acceleration_m_s2: np.ndarray

    def as_serializable(self) -> dict[str, object]:
        out = asdict(self)
        out["position_m"] = [float(x) for x in self.position_m]
        out["velocity_m_s"] = [float(x) for x in self.velocity_m_s]
        out["acceleration_m_s2"] = [float(x) for x in self.acceleration_m_s2]
        return out


@dataclass(frozen=True)
class NathFieldSample:
    E_moving_charge: np.ndarray
    E_charge_growth: np.ndarray
    E_positive_ion_trail: np.ndarray
    B_near: np.ndarray
    B_acceleration_radiation: np.ndarray
    B_charge_growth_radiation: np.ndarray

    @property
    def E_total(self) -> np.ndarray:
        return self.E_moving_charge + self.E_charge_growth + self.E_positive_ion_trail

    @property
    def B_total(self) -> np.ndarray:
        return self.B_near + self.B_acceleration_radiation + self.B_charge_growth_radiation


def avalanche_state(config: NathAvalancheConfig, t_s: float, *, acceleration_m_s2: float = 0.0) -> NathAvalancheState:
    dt = float(t_s - config.start_time_s)
    if acceleration_m_s2 == 0.0:
        z = config.start_z_m + config.velocity_m_s * dt
        vz = config.velocity_m_s
    else:
        z = config.start_z_m + config.velocity_m_s * dt + 0.5 * acceleration_m_s2 * dt**2
        vz = config.velocity_m_s + acceleration_m_s2 * dt
    q_electron = -config.q0_C * math.exp(config.growth_rate_s * max(dt, 0.0))
    dqdt = config.growth_rate_s * q_electron if dt >= 0.0 else 0.0
    return NathAvalancheState(
        t_s=float(t_s),
        q_C=float(q_electron),
        dq_dt_A=float(dqdt),
        position_m=np.asarray([0.0, 0.0, z], dtype=float),
        velocity_m_s=np.asarray([0.0, 0.0, vz], dtype=float),
        acceleration_m_s2=np.asarray([0.0, 0.0, acceleration_m_s2], dtype=float),
    )


def nath_K(rhat: np.ndarray, velocity_m_s: np.ndarray) -> float:
    return float(1.0 - np.dot(rhat, velocity_m_s) / C0)


def retarded_time_for_moving_head(
    observer_position_m: np.ndarray,
    config: NathAvalancheConfig,
    t_obs_s: float,
    *,
    acceleration_m_s2: float = 0.0,
    max_iter: int = 64,
    tol_s: float = 1e-18,
) -> tuple[float, bool, int]:
    t_obs_s = float(t_obs_s)
    observer_position_m = np.asarray(observer_position_m, dtype=float)
    tr = t_obs_s - np.linalg.norm(observer_position_m - np.asarray([0.0, 0.0, config.start_z_m])) / C0
    for it in range(max_iter):
        state = avalanche_state(config, tr, acceleration_m_s2=acceleration_m_s2)
        rvec = observer_position_m - state.position_m
        R = float(np.linalg.norm(rvec))
        if R <= 0.0:
            return tr, False, it + 1
        rhat = rvec / R
        f = tr + R / C0 - t_obs_s
        if abs(f) <= tol_s:
            return float(tr), True, it + 1
        K = nath_K(rhat, state.velocity_m_s)
        if K <= 0.0:
            return float(tr), False, it + 1
        tr -= f / K
    return float(tr), False, max_iter


def evaluate_nath_field(
    state: NathAvalancheState,
    observer_position_m: np.ndarray,
    *,
    config: NathAvalancheConfig | None = None,
    ion_quadrature: int = 200,
) -> NathFieldSample:
    observer_position_m = np.asarray(observer_position_m, dtype=float)
    rvec = observer_position_m - state.position_m
    R = float(np.linalg.norm(rvec))
    if R <= 0.0:
        raise ValueError("observer is colocated with Nath moving head")
    rhat = rvec / R
    v = state.velocity_m_s
    a = state.acceleration_m_s2
    beta = v / C0
    K = nath_K(rhat, v)
    if K <= 0.0:
        raise ValueError("retarded factor K must be positive")

    q = state.q_C
    dqdt = state.dq_dt_A
    beta2 = float(np.dot(beta, beta))
    charge_velocity = (rhat - beta) * (1.0 - beta2) / R**2
    charge_acceleration = np.cross(rhat, np.cross(rhat - beta, a / C0)) / (C0 * R)
    E_moving = K_E * q * (charge_velocity + charge_acceleration) / K**3
    E_qdot = K_E * dqdt * (rhat * float(np.dot(rhat, v)) - v) / (C0**2 * K**2 * R)

    E_ion = np.zeros(3)
    if config is not None and state.t_s >= config.start_time_s:
        E_ion = positive_ion_trail_E(config, state.t_s, observer_position_m, ion_quadrature=ion_quadrature)

    B_near = K_E / C0**2 * (q * C0 / (K**3 * R**2)) * np.cross(rhat, rhat - beta) * (1.0 - beta2)
    B_acc = K_E / C0**2 * (q / (K**3 * R)) * np.cross(rhat, np.cross(rhat - beta, a / C0))
    B_qdot = K_E / C0**2 * (dqdt / (K**2 * R)) * np.cross(rhat, rhat - beta)
    return NathFieldSample(
        E_moving_charge=E_moving,
        E_charge_growth=E_qdot,
        E_positive_ion_trail=E_ion,
        B_near=B_near,
        B_acceleration_radiation=B_acc,
        B_charge_growth_radiation=B_qdot,
    )


def evaluate_nath_at_observer_time(
    config: NathAvalancheConfig,
    observer: Observer,
    t_obs_s: float,
    *,
    acceleration_m_s2: float = 0.0,
    include_ion_trail: bool = True,
) -> tuple[NathFieldSample, NathAvalancheState, float, bool]:
    tr, ok, _ = retarded_time_for_moving_head(observer.position, config, t_obs_s, acceleration_m_s2=acceleration_m_s2)
    state = avalanche_state(config, tr, acceleration_m_s2=acceleration_m_s2)
    sample = evaluate_nath_field(state, observer.position, config=config if include_ion_trail else None)
    return sample, state, tr, ok and tr >= config.start_time_s


def positive_ion_line_density_Cpm(config: NathAvalancheConfig, z_m: np.ndarray | float) -> np.ndarray:
    z = np.asarray(z_m, dtype=float)
    if config.velocity_m_s <= 0.0:
        raise ValueError("positive ion trail requires positive z velocity")
    t_dep = config.start_time_s + (z - config.start_z_m) / config.velocity_m_s
    q_abs = config.q0_C * np.exp(config.growth_rate_s * np.maximum(t_dep - config.start_time_s, 0.0))
    lam = config.growth_rate_s * q_abs / config.velocity_m_s
    return np.asarray(lam, dtype=float)


def positive_ion_trail_E(
    config: NathAvalancheConfig,
    t_s: float,
    observer_position_m: np.ndarray,
    *,
    ion_quadrature: int = 200,
) -> np.ndarray:
    z_head = avalanche_state(config, t_s).position_m[2]
    if z_head <= config.start_z_m:
        # Initial balancing ion charge is represented as a narrow but finite line segment.
        zs = np.asarray([config.start_z_m], dtype=float)
        weights = np.asarray([1.0], dtype=float)
        charges = weights * config.q0_C
    else:
        zs = np.linspace(config.start_z_m, z_head, ion_quadrature)
        dz = float(zs[1] - zs[0]) if zs.size > 1 else 0.0
        weights = np.full_like(zs, dz)
        weights[0] *= 0.5
        weights[-1] *= 0.5
        charges = positive_ion_line_density_Cpm(config, zs) * weights
        charges[0] += config.q0_C
    src = np.column_stack((np.zeros_like(zs), np.zeros_like(zs), zs))
    rvec = np.asarray(observer_position_m, dtype=float)[None, :] - src
    R = np.linalg.norm(rvec, axis=1)
    rhat = rvec / R[:, None]
    return K_E * np.sum(charges[:, None] * rhat / R[:, None] ** 2, axis=0)


def make_nath_source_series(config: NathAvalancheConfig, times_s: np.ndarray) -> SourceSeries:
    times_s = np.asarray(times_s, dtype=float)
    if np.any(np.diff(times_s) <= 0.0):
        raise ValueError("times_s must be strictly increasing")
    z = np.arange(config.z_min_m, config.z_max_m + 0.5 * config.dz_m, config.dz_m)
    n = z.size
    x = np.zeros(n)
    y = np.zeros(n)
    area = config.transverse_size_m**2
    volume = area * config.dz_m
    base_cols = {
        "level": np.zeros(n, dtype=int),
        "x_center": x,
        "y_center": y,
        "z_center": z,
        "dx": np.full(n, config.transverse_size_m),
        "dy": np.full(n, config.transverse_size_m),
        "dz": np.full(n, config.dz_m),
        "cell_volume": np.full(n, volume),
    }
    records = []
    for t in times_s:
        state = avalanche_state(config, float(t))
        q_abs = abs(state.q_C)
        g = np.exp(-0.5 * ((z - state.position_m[2]) / config.source_width_m) ** 2)
        g /= np.sum(g * volume)
        rho_head = state.q_C * g
        J = np.zeros((n, 3), dtype=float)
        J[:, 2] = rho_head * state.velocity_m_s[2]

        ion_charge = np.zeros(n, dtype=float)
        active = (z >= config.start_z_m) & (z <= state.position_m[2])
        if np.any(active):
            weights = positive_ion_line_density_Cpm(config, z[active]) * config.dz_m
            growth_charge = q_abs - config.q0_C
            if growth_charge > 0.0 and np.sum(weights) > 0.0:
                ion_charge[active] = weights * growth_charge / np.sum(weights)
        origin_index = int(np.argmin(np.abs(z - config.start_z_m)))
        ion_charge[origin_index] += config.q0_C
        rho_ion = ion_charge / volume
        rho = rho_head + rho_ion

        meta = SourceMetadata(
            case_id=config.case_id,
            solver="Nath-analytical-regularized",
            solver_version="F-R1",
            time_s=float(t),
            coordinate_system="cartesian",
            pressure_Pa=101325.0,
            temperature_K=300.0,
            geometry_id="straight_z_avalanche",
            voltage_state="analytical constant-velocity avalanche",
            photoionization="off",
            source_definition="moving finite electron head plus stationary positive ion trail",
            units={
                "rho": "C m^-3",
                "J": "A m^-2",
                "cell_volume": "m^3",
            },
            extra={
                "reference": "Nath 2022 IEEE TEMC / TechRxiv v1",
                "source_width_m": config.source_width_m,
                "transverse_size_m": config.transverse_size_m,
                "dz_m": config.dz_m,
            },
        )
        cols = dict(base_cols)
        cols.update({"rho": rho, "Jx": J[:, 0], "Jy": J[:, 1], "Jz": J[:, 2]})
        records.append(SourceRecord(meta, cols))
    return SourceSeries(tuple(records))


def waveform_metrics(reference: np.ndarray, candidate: np.ndarray, times_s: np.ndarray) -> dict[str, float]:
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    times_s = np.asarray(times_s, dtype=float)
    denom = float(np.linalg.norm(reference))
    l2 = float(np.linalg.norm(candidate - reference) / max(denom, 1e-300))
    ref_peak = float(reference[np.argmax(np.abs(reference))])
    cand_peak = float(candidate[np.argmax(np.abs(reference))])
    peak_rel = abs(cand_peak - ref_peak) / max(abs(ref_peak), 1e-300)
    if np.std(reference) == 0.0 or np.std(candidate) == 0.0:
        corr = 1.0 if np.allclose(reference, candidate) else float("nan")
    else:
        corr = float(np.corrcoef(reference, candidate)[0, 1])
    t_ref = float(times_s[np.argmax(np.abs(reference))])
    t_cand = float(times_s[np.argmax(np.abs(candidate))])
    return {
        "normalized_L2_error": l2,
        "relative_peak_error": float(peak_rel),
        "correlation": corr,
        "arrival_time_difference_s": abs(t_cand - t_ref),
        "reference_peak": ref_peak,
        "candidate_peak": cand_peak,
    }


def source_charge_conservation_summary(series: SourceSeries) -> dict[str, float]:
    charges = np.asarray([total_charge(record) for record in series.records], dtype=float)
    moments = np.asarray([current_moment(record) for record in series.records], dtype=float)
    return {
        "Q_min_C": float(np.min(charges)),
        "Q_max_C": float(np.max(charges)),
        "Q_span_C": float(np.max(charges) - np.min(charges)),
        "Q_rms_C": float(np.sqrt(np.mean(charges**2))),
        "M_peak_Am": float(np.max(np.linalg.norm(moments, axis=1))),
    }
