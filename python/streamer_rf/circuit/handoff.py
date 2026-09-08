from __future__ import annotations

import csv
import heapq
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from streamer_rf.rf.jefimenko.constants import EPS0
from streamer_rf.streamer.morrow_lowke_reference import evaluate as morrow_lowke_evaluate

QE = 1.602176634e-19
KB = 1.380649e-23


@dataclass(frozen=True)
class PercolationConfig:
    sigma_relative_threshold: float = 0.1
    neighbor_mode: str = "von_neumann"


@dataclass(frozen=True)
class PercolationResult:
    percolation_valid: bool
    percolation_status: str
    threshold_sigma_S_m: float
    threshold_fraction: float
    path_length_m: float = math.nan
    bottleneck_sigma_S_m: float = math.nan
    percolated_cell_count: int = 0


@dataclass(frozen=True)
class PersistenceSummary:
    first_percolation_time_s: float
    number_of_percolated_samples: int
    longest_consecutive_percolated_samples: int


@dataclass(frozen=True)
class HandoffCalibration:
    source_id: str
    source_reference: str
    Q_required_J: float | None = None
    Pi_H_threshold: float | None = None
    Xi_sigma_threshold: float | None = None
    bridge_requirement: bool | None = None

    def validate(self) -> None:
        if not self.source_id or not self.source_reference:
            raise ValueError("calibration provenance metadata is required")
        for name in ("Q_required_J", "Pi_H_threshold", "Xi_sigma_threshold"):
            value = getattr(self, name)
            if value is not None and (not math.isfinite(value) or value <= 0.0):
                raise ValueError(f"{name} must be finite and positive when supplied")


@dataclass(frozen=True)
class ColdThermalHandoffState:
    time_s: float
    bridge_legacy: bool
    percolation_valid: bool
    percolation_status: str
    PJ_channel_W: float
    QJ_channel_J: float
    channel_volume_m3: float
    channel_length_m: float
    channel_effective_radius_m: float
    sigma_eff_S_m: float
    E_channel_mean_V_m: float
    ne_channel_mean_m3: float
    ne_channel_max_m3: float
    Gb_S: float
    Rb_ohm: float
    dGb_dt_S_s: float
    tau_sigma_s: float
    tau_evolution_s: float
    Xi_sigma: float
    Q_required_J: float
    Pi_H: float
    calibration_status: str
    handoff_status: str
    thermal_profile_status: str = "NOT_AVAILABLE"


def _finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


def tau_sigma(sigma_eff_S_m: float) -> tuple[float, str]:
    if not _finite_positive(float(sigma_eff_S_m)):
        return math.nan, "INVALID_SIGMA_EFF"
    return EPS0 / float(sigma_eff_S_m), "VALID"


def tau_evolution(Gb_S: float, dGb_dt_S_s: float) -> tuple[float, str]:
    Gb = float(Gb_S)
    dGb = float(dGb_dt_S_s)
    if not _finite_positive(Gb) or not math.isfinite(dGb) or dGb == 0.0:
        return math.nan, "ZERO_GB_OR_DGBDT"
    return abs(Gb / dGb), "VALID"


def xi_sigma(sigma_eff_S_m: float, Gb_S: float, dGb_dt_S_s: float) -> tuple[float, str]:
    ts, ts_status = tau_sigma(sigma_eff_S_m)
    te, te_status = tau_evolution(Gb_S, dGb_dt_S_s)
    if ts_status != "VALID" or te_status != "VALID":
        return math.nan, "UNAVAILABLE"
    return ts / te, "VALID"


def pi_h(QJ_channel_J: float, Q_required_J: float | None) -> tuple[float, str]:
    if Q_required_J is None or not _finite_positive(float(Q_required_J)):
        return math.nan, "THERMAL_ENERGY_REFERENCE_NOT_AVAILABLE"
    qj = float(QJ_channel_J)
    if not math.isfinite(qj) or qj < 0.0:
        return math.nan, "INVALID_QJ_CHANNEL"
    return qj / float(Q_required_J), "VALID"


def apply_handoff_logic(
    *,
    bridge_or_percolated: bool,
    QJ_channel_J: float,
    Xi_sigma: float,
    calibration: HandoffCalibration | None = None,
) -> tuple[str, float, str]:
    if calibration is None:
        return "HANDOFF_CALIBRATION_PENDING", math.nan, "NOT_AVAILABLE"
    calibration.validate()
    Pi_H, pi_status = pi_h(QJ_channel_J, calibration.Q_required_J)
    if pi_status != "VALID":
        return "HANDOFF_CALIBRATION_INCOMPLETE", Pi_H, "INVALID_OR_MISSING"
    checks: list[bool] = []
    if calibration.bridge_requirement is not None:
        checks.append(bool(bridge_or_percolated) == calibration.bridge_requirement)
    if calibration.Pi_H_threshold is not None:
        checks.append(Pi_H >= calibration.Pi_H_threshold)
    if calibration.Xi_sigma_threshold is not None:
        checks.append(math.isfinite(Xi_sigma) and Xi_sigma <= calibration.Xi_sigma_threshold)
    if not checks:
        return "HANDOFF_CALIBRATION_INCOMPLETE", Pi_H, "VALID_NO_CRITERIA"
    return ("HANDOFF_CRITERIA_MET" if all(checks) else "HANDOFF_CRITERIA_NOT_MET"), Pi_H, "VALID"


def conductive_percolation(
    sigma_S_m: np.ndarray,
    gas_mask: np.ndarray,
    hv_adjacent_mask: np.ndarray,
    ground_adjacent_mask: np.ndarray,
    dr_m: float,
    dz_m: float,
    config: PercolationConfig = PercolationConfig(),
) -> PercolationResult:
    sigma = np.asarray(sigma_S_m, dtype=float)
    gas = np.asarray(gas_mask, dtype=bool)
    hv = np.asarray(hv_adjacent_mask, dtype=bool)
    ground = np.asarray(ground_adjacent_mask, dtype=bool)
    if sigma.shape != gas.shape or sigma.shape != hv.shape or sigma.shape != ground.shape:
        raise ValueError("sigma and masks must have identical shape")
    if config.neighbor_mode != "von_neumann":
        raise ValueError("only deterministic von_neumann connectivity is supported")
    if not math.isfinite(config.sigma_relative_threshold) or config.sigma_relative_threshold < 0.0:
        return PercolationResult(False, "INVALID_THRESHOLD", math.nan, config.sigma_relative_threshold)
    if dr_m <= 0.0 or dz_m <= 0.0 or not math.isfinite(dr_m) or not math.isfinite(dz_m):
        return PercolationResult(False, "INVALID_GRID_SPACING", math.nan, config.sigma_relative_threshold)
    finite_gas = gas & np.isfinite(sigma) & (sigma >= 0.0)
    if not np.any(finite_gas):
        return PercolationResult(False, "NO_FINITE_GAS_CONDUCTIVITY", math.nan, config.sigma_relative_threshold)
    sigma_max = float(np.max(sigma[finite_gas]))
    if sigma_max <= 0.0:
        return PercolationResult(False, "NO_CONDUCTIVE_CELLS", 0.0, config.sigma_relative_threshold)
    threshold = config.sigma_relative_threshold * sigma_max
    active = finite_gas & (sigma >= threshold)
    starts = np.argwhere(active & hv)
    targets = active & ground
    if starts.size == 0 or not np.any(targets):
        return PercolationResult(False, "NO_ELECTRODE_ADJACENT_SEEDS", threshold, config.sigma_relative_threshold)

    distances = np.full(sigma.shape, np.inf)
    previous: dict[tuple[int, int], tuple[int, int]] = {}
    pq: list[tuple[float, int, int]] = []
    for i, j in sorted(map(tuple, starts.tolist())):
        distances[i, j] = 0.0
        heapq.heappush(pq, (0.0, i, j))

    target: tuple[int, int] | None = None
    neighbors = ((-1, 0, dr_m), (1, 0, dr_m), (0, -1, dz_m), (0, 1, dz_m))
    while pq:
        dist, i, j = heapq.heappop(pq)
        if dist != distances[i, j]:
            continue
        if targets[i, j]:
            target = (i, j)
            break
        for di, dj, step_len in neighbors:
            ni, nj = i + di, j + dj
            if ni < 0 or nj < 0 or ni >= sigma.shape[0] or nj >= sigma.shape[1]:
                continue
            if not active[ni, nj]:
                continue
            nd = dist + step_len
            if nd < distances[ni, nj]:
                distances[ni, nj] = nd
                previous[(ni, nj)] = (i, j)
                heapq.heappush(pq, (nd, ni, nj))

    if target is None:
        return PercolationResult(False, "NO_CONDUCTIVE_PATH", threshold, config.sigma_relative_threshold)

    path = [target]
    while path[-1] in previous:
        path.append(previous[path[-1]])
    bottleneck = min(float(sigma[i, j]) for i, j in path)
    return PercolationResult(
        True,
        "PERCOLATED",
        threshold,
        config.sigma_relative_threshold,
        float(distances[target]),
        bottleneck,
        len(path),
    )


def percolation_persistence(time_s: Iterable[float], percolated: Iterable[bool]) -> PersistenceSummary:
    times = list(map(float, time_s))
    flags = list(map(bool, percolated))
    if len(times) != len(flags):
        raise ValueError("time and percolation arrays must have equal length")
    first = math.nan
    count = 0
    longest = 0
    current = 0
    for t, flag in zip(times, flags):
        if flag:
            if math.isnan(first):
                first = t
            count += 1
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return PersistenceSummary(first, count, longest)


def state_from_cr4_row(
    row: dict[str, str],
    percolation: PercolationResult | None = None,
    calibration: HandoffCalibration | None = None,
) -> ColdThermalHandoffState:
    def f(name: str) -> float:
        try:
            return float(row.get(name, "nan"))
        except ValueError:
            return math.nan

    bridge = bool(int(float(row.get("bridge_flag", "0"))))
    perc_valid = percolation.percolation_valid if percolation is not None else False
    perc_status = percolation.percolation_status if percolation is not None else "NO_FIELD_SNAPSHOT"
    xi = f("Xi_sigma")
    status, Pi_H, cal_status = apply_handoff_logic(
        bridge_or_percolated=bridge or perc_valid,
        QJ_channel_J=f("QJ_channel_J"),
        Xi_sigma=xi,
        calibration=calibration,
    )
    q_required = calibration.Q_required_J if calibration is not None and calibration.Q_required_J is not None else math.nan
    return ColdThermalHandoffState(
        time_s=f("time_s"),
        bridge_legacy=bridge,
        percolation_valid=perc_valid,
        percolation_status=perc_status,
        PJ_channel_W=f("PJ_channel_W"),
        QJ_channel_J=f("QJ_channel_J"),
        channel_volume_m3=f("channel_volume_m3"),
        channel_length_m=f("channel_length_m"),
        channel_effective_radius_m=f("channel_effective_radius_m"),
        sigma_eff_S_m=f("sigma_eff_S_m"),
        E_channel_mean_V_m=f("E_channel_mean_V_m"),
        ne_channel_mean_m3=f("ne_channel_mean_m3"),
        ne_channel_max_m3=f("ne_channel_max_m3"),
        Gb_S=f("Gb_S"),
        Rb_ohm=f("Rb_ohm"),
        dGb_dt_S_s=f("dGb_dt_S_s"),
        tau_sigma_s=f("tau_sigma_s"),
        tau_evolution_s=f("tau_evolution_s"),
        Xi_sigma=xi,
        Q_required_J=q_required,
        Pi_H=Pi_H,
        calibration_status=cal_status,
        handoff_status=status,
    )


def g1_initial_state_contract(state: ColdThermalHandoffState | None) -> dict[str, float | str | bool]:
    if state is None:
        return {"handoff_time_candidate": math.nan, "thermal_profile_status": "NOT_AVAILABLE", "state_valid": False}
    return {
        "handoff_time_candidate": state.time_s,
        "channel_length_m": state.channel_length_m,
        "channel_effective_radius_m": state.channel_effective_radius_m,
        "channel_volume_m3": state.channel_volume_m3,
        "sigma_eff_S_m": state.sigma_eff_S_m,
        "QJ_channel_J": state.QJ_channel_J,
        "PJ_channel_W": state.PJ_channel_W,
        "E_channel_mean_V_m": state.E_channel_mean_V_m,
        "ne_channel_mean_m3": state.ne_channel_mean_m3,
        "ne_channel_max_m3": state.ne_channel_max_m3,
        "Gb_S": state.Gb_S,
        "Rb_ohm": state.Rb_ohm,
        "thermal_profile_status": "NOT_AVAILABLE",
        "state_valid": True,
    }


def _read_key_values(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def _load_fields_snapshot(path: Path, pressure_Pa: float, temperature_K: float):
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if not rows:
        raise ValueError(f"{path} contains no field rows")
    max_i = max(int(r["i"]) for r in rows)
    max_j = max(int(r["j"]) for r in rows)
    shape = (max_i + 1, max_j + 1)
    sigma = np.zeros(shape)
    cell_type = np.empty(shape, dtype=object)
    r = np.zeros(shape)
    z = np.zeros(shape)
    N = pressure_Pa / (KB * temperature_K)
    for row in rows:
        i, j = int(row["i"]), int(row["j"])
        cell_type[i, j] = row["cell_type"]
        r[i, j] = float(row["r_m"])
        z[i, j] = float(row["z_m"])
        if row["cell_type"] == "gas":
            E = float(row["E_V_m"])
            ne = float(row["ne_m_3"])
            sigma[i, j] = QE * morrow_lowke_evaluate(E, N)["mobility"] * ne
    return sigma, cell_type, r, z


def electrode_adjacent_masks(cell_type: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gas = cell_type == "gas"
    hv = np.zeros(cell_type.shape, dtype=bool)
    ground = np.zeros(cell_type.shape, dtype=bool)
    for i in range(cell_type.shape[0]):
        for j in range(cell_type.shape[1]):
            if not gas[i, j]:
                continue
            for ni, nj in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                if ni < 0 or nj < 0 or ni >= cell_type.shape[0] or nj >= cell_type.shape[1]:
                    continue
                hv[i, j] |= cell_type[ni, nj] == "hv"
                ground[i, j] |= cell_type[ni, nj] == "ground"
    return gas, hv, ground


def analyze_cr4_handoff_run(
    run_dir: Path,
    output_dir: Path,
    config: PercolationConfig = PercolationConfig(),
    sensitivity: tuple[float, ...] = (0.05, 0.1, 0.2),
) -> dict[str, object]:
    run_dir = Path(run_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meta = _read_key_values(run_dir / "case_metadata.txt")
    summary_in = _read_key_values(run_dir / "summary.txt")
    cr4_path = run_dir / "cold_thermal_handoff.csv"
    if not cr4_path.exists():
        raise FileNotFoundError(f"missing C-R4 handoff CSV: {cr4_path}")

    rows = list(csv.DictReader(cr4_path.open(encoding="utf-8")))
    by_step = {int(float(row["step"])): row for row in rows}
    pressure = float(meta.get("pressure_Pa", 101325.0))
    temperature = float(meta.get("temperature_K", 300.0))

    snapshot_results: dict[int, PercolationResult] = {}
    sensitivity_rows: list[dict[str, object]] = []
    for path in sorted(run_dir.glob("fields_*.csv")):
        if path.name == "fields_initial.csv":
            step = 0
        elif path.name == "fields_final.csv":
            step = int(float(summary_in.get("total_steps", max(by_step) if by_step else 0)))
        elif path.name.startswith("fields_snapshot_"):
            step = int(path.stem.split("_")[-1])
        else:
            continue
        sigma, cell_type, r, z = _load_fields_snapshot(path, pressure, temperature)
        gas, hv, ground = electrode_adjacent_masks(cell_type)
        dr = float(np.min(np.diff(np.unique(r[:, 0]))))
        dz = float(np.min(np.diff(np.unique(z[0, :]))))
        snapshot_results[step] = conductive_percolation(sigma, gas, hv, ground, dr, dz, config)
        for frac in sensitivity:
            result = conductive_percolation(sigma, gas, hv, ground, dr, dz, PercolationConfig(frac))
            sensitivity_rows.append({"step": step, **asdict(result)})

    states: list[ColdThermalHandoffState] = []
    diagnostics_rows: list[dict[str, object]] = []
    for row in rows:
        step = int(float(row["step"]))
        percolation = snapshot_results.get(step)
        state = state_from_cr4_row(row, percolation=percolation, calibration=None)
        states.append(state)
        diagnostics_rows.append(asdict(state))

    persistence = percolation_persistence([s.time_s for s in states], [s.percolation_valid for s in states])
    bridge_persistence = percolation_persistence([s.time_s for s in states], [s.bridge_legacy for s in states])
    bridged_or_percolated = [s for s in states if s.bridge_legacy or s.percolation_valid]
    candidate = bridged_or_percolated[0] if bridged_or_percolated else None
    latest_cold_state = states[-1] if states else None

    diagnostics_path = output_dir / "g0_handoff_diagnostics.csv"
    if diagnostics_rows:
        with diagnostics_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(diagnostics_rows[0].keys()))
            writer.writeheader()
            writer.writerows(diagnostics_rows)

    sensitivity_path = output_dir / "g0_percolation_sensitivity.csv"
    if sensitivity_rows:
        with sensitivity_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(sensitivity_rows[0].keys()))
            writer.writeheader()
            writer.writerows(sensitivity_rows)

    qj = np.array([s.QJ_channel_J for s in states], dtype=float)
    xi = np.array([s.Xi_sigma for s in states], dtype=float)
    sigma_eff = np.array([s.sigma_eff_S_m for s in states], dtype=float)
    Gb = np.array([s.Gb_S for s in states], dtype=float)
    Rb = np.array([s.Rb_ohm for s in states], dtype=float)
    finite = lambda a: a[np.isfinite(a)]
    summary = {
        "framework_status": "FRAMEWORK_VALIDATED",
        "production_case_status": "PRODUCTION_HANDOFF_NOT_REACHED"
        if not bridged_or_percolated
        else "PRODUCTION_HANDOFF_CANDIDATE_RECORDED",
        "calibration_status": "NOT_AVAILABLE",
        "thermal_energy_reference_status": "NOT_AVAILABLE",
        "handoff_status": "HANDOFF_CALIBRATION_PENDING",
        "case_id": summary_in.get("case_id", meta.get("case_id", "UNKNOWN")),
        "time_interval_s": [float(rows[0]["time_s"]), float(rows[-1]["time_s"])] if rows else [math.nan, math.nan],
        "accepted_samples": len(rows),
        "bridge_legacy_classification": "AXIAL_THRESHOLD_PROXY",
        "bridge_legacy_first_time_s": bridge_persistence.first_percolation_time_s,
        "bridge_legacy_samples": bridge_persistence.number_of_percolated_samples,
        "bridge_legacy_longest_consecutive_samples": bridge_persistence.longest_consecutive_percolated_samples,
        "percolation_first_time_s": persistence.first_percolation_time_s,
        "percolation_samples_with_field_snapshots": len(snapshot_results),
        "percolation_true_samples": persistence.number_of_percolated_samples,
        "percolation_longest_consecutive_samples": persistence.longest_consecutive_percolated_samples,
        "QJ_channel_final_J": float(qj[np.isfinite(qj)][-1]) if np.any(np.isfinite(qj)) else math.nan,
        "QJ_channel_monotonic_nonnegative": bool(np.all(np.diff(finite(qj)) >= -1e-30)) if finite(qj).size > 1 else True,
        "sigma_eff_range_S_m": [float(np.min(finite(sigma_eff))), float(np.max(finite(sigma_eff)))] if finite(sigma_eff).size else [math.nan, math.nan],
        "Gb_range_S": [float(np.min(finite(Gb))), float(np.max(finite(Gb)))] if finite(Gb).size else [math.nan, math.nan],
        "Rb_range_ohm": [float(np.min(finite(Rb))), float(np.max(finite(Rb)))] if finite(Rb).size else [math.nan, math.nan],
        "Xi_sigma_range": [float(np.min(finite(xi))), float(np.max(finite(xi)))] if finite(xi).size else [math.nan, math.nan],
        "G1_handoff_candidate_contract": g1_initial_state_contract(candidate),
        "latest_cold_state_contract": g1_initial_state_contract(latest_cold_state),
        "diagnostics_csv": str(diagnostics_path),
        "percolation_sensitivity_csv": str(sensitivity_path),
    }
    summary_path = output_dir / "g0_handoff_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8")
    return summary
