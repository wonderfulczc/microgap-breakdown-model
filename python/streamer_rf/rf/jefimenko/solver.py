from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from streamer_rf.rf.jefimenko.constants import C0, K_B, K_E
from streamer_rf.rf.jefimenko.observer import Observer
from streamer_rf.rf.jefimenko.temporal import RetardedSourceInterpolator
from streamer_rf.rf.source.schema import SourceSeries


@dataclass(frozen=True)
class FieldSample:
    time_s: float
    observer_id: str
    E_rho_near: np.ndarray
    E_drho_induction: np.ndarray
    E_dJ_radiation: np.ndarray
    B_J_near: np.ndarray
    B_dJ_radiation: np.ndarray
    retarded_time_valid: bool
    retarded_time_valid_fraction: float
    source_manifest_id: str

    @property
    def E_total(self) -> np.ndarray:
        return self.E_rho_near + self.E_drho_induction + self.E_dJ_radiation

    @property
    def B_total(self) -> np.ndarray:
        return self.B_J_near + self.B_dJ_radiation

    def to_row(self) -> dict[str, float | str | bool]:
        row: dict[str, float | str | bool] = {
            "time_s": self.time_s,
            "observer_id": self.observer_id,
            "retarded_time_valid": self.retarded_time_valid,
            "retarded_time_valid_fraction": self.retarded_time_valid_fraction,
            "source_manifest_id": self.source_manifest_id,
        }
        for prefix, vec in (
            ("total", self.E_total),
            ("rho", self.E_rho_near),
            ("drho", self.E_drho_induction),
            ("dJ", self.E_dJ_radiation),
        ):
            row[f"Ex_{prefix}"] = float(vec[0])
            row[f"Ey_{prefix}"] = float(vec[1])
            row[f"Ez_{prefix}"] = float(vec[2])
        for prefix, vec in (
            ("total", self.B_total),
            ("J", self.B_J_near),
            ("dJ", self.B_dJ_radiation),
        ):
            row[f"Bx_{prefix}"] = float(vec[0])
            row[f"By_{prefix}"] = float(vec[1])
            row[f"Bz_{prefix}"] = float(vec[2])
        return row

    def to_dict(self) -> dict[str, object]:
        out = asdict(self)
        for key in ("E_rho_near", "E_drho_induction", "E_dJ_radiation", "B_J_near", "B_dJ_radiation"):
            out[key] = [float(x) for x in out[key]]
        out["E_total"] = [float(x) for x in self.E_total]
        out["B_total"] = [float(x) for x in self.B_total]
        return out


def evaluate_observer(
    series: SourceSeries,
    observer: Observer,
    time_s: float,
    *,
    source_manifest_id: str = "unknown",
    chunk_size: int = 200_000,
) -> FieldSample:
    record0 = series.records[0]
    observer.validate_outside_source(record0)
    interp = RetardedSourceInterpolator(series)

    xyz = np.column_stack(
        (record0.columns["x_center"], record0.columns["y_center"], record0.columns["z_center"])
    )
    volume = record0.columns["cell_volume"]
    obs = observer.position

    E_rho = np.zeros(3)
    E_drho = np.zeros(3)
    E_dJ = np.zeros(3)
    B_J = np.zeros(3)
    B_dJ = np.zeros(3)
    valid_count = 0
    n_cells = record0.n_cells

    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        r_vec = obs[None, :] - xyz[start:stop]
        r_mag = np.linalg.norm(r_vec, axis=1)
        if np.any(r_mag <= 0.0):
            raise ValueError("observer coincides with a source cell center")
        r_hat = r_vec / r_mag[:, None]
        tr = time_s - r_mag / C0
        src = interp.evaluate(tr)
        valid = np.asarray(src["valid_mask"], dtype=bool)
        valid_count += int(np.count_nonzero(valid))
        if not np.any(valid):
            continue
        vv = volume[start:stop][valid]
        rr = r_mag[valid]
        rh = r_hat[valid]
        rho = np.asarray(src["rho"])[valid]
        J = np.asarray(src["J"])[valid]
        drho = np.asarray(src["drho_dt"])[valid]
        dJ = np.asarray(src["dJ_dt"])[valid]

        E_rho += K_E * np.sum((rho * vv / rr**2)[:, None] * rh, axis=0)
        E_drho += K_E * np.sum((drho * vv / (C0 * rr))[:, None] * rh, axis=0)
        E_dJ += -K_E * np.sum((vv / (C0**2 * rr))[:, None] * dJ, axis=0)
        B_J += K_B * np.sum(np.cross(J, rh) * (vv / rr**2)[:, None], axis=0)
        B_dJ += K_B * np.sum(np.cross(dJ, rh) * (vv / (C0 * rr))[:, None], axis=0)

    valid_fraction = float(valid_count / n_cells)
    return FieldSample(
        time_s=float(time_s),
        observer_id=observer.observer_id,
        E_rho_near=E_rho,
        E_drho_induction=E_drho,
        E_dJ_radiation=E_dJ,
        B_J_near=B_J,
        B_dJ_radiation=B_dJ,
        retarded_time_valid=valid_count == n_cells,
        retarded_time_valid_fraction=valid_fraction,
        source_manifest_id=source_manifest_id,
    )


def evaluate_waveform(
    series: SourceSeries,
    observers: list[Observer],
    times_s: np.ndarray,
    *,
    source_manifest_id: str = "unknown",
    chunk_size: int = 200_000,
) -> list[FieldSample]:
    samples: list[FieldSample] = []
    for observer in observers:
        for time_s in np.asarray(times_s, dtype=float):
            samples.append(
                evaluate_observer(
                    series,
                    observer,
                    float(time_s),
                    source_manifest_id=source_manifest_id,
                    chunk_size=chunk_size,
                )
            )
    return samples
