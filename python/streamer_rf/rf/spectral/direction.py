from __future__ import annotations

import numpy as np

from .fluence import Z0


def radial_fluence_pattern(time_s: np.ndarray, observer_positions: np.ndarray, E_vectors: np.ndarray) -> dict[str, object]:
    t = np.asarray(time_s, dtype=float)
    pos = np.asarray(observer_positions, dtype=float)
    E = np.asarray(E_vectors, dtype=float)
    if E.shape[:2] != (pos.shape[0], t.size):
        raise ValueError("E_vectors must have shape (n_observer, n_time, 3)")
    dt = float(np.diff(t)[0])
    rhat = pos / np.linalg.norm(pos, axis=1)[:, None]
    radial = np.einsum("otc,oc->ot", E, rhat)
    Eperp = E - radial[:, :, None] * rhat[:, None, :]
    fluence = np.sum(np.sum(Eperp * Eperp, axis=2), axis=1) * dt / Z0
    idx = int(np.argmax(fluence))
    return {
        "fluence_Jpm2": fluence,
        "normalized_pattern": fluence / max(float(np.max(fluence)), 1e-300),
        "peak_index": idx,
        "peak_direction": rhat[idx],
    }


def dipole_sin2_pattern_error(observer_positions: np.ndarray, fluence: np.ndarray, dipole_axis: np.ndarray = np.array([0.0, 0.0, 1.0])) -> float:
    pos = np.asarray(observer_positions, dtype=float)
    rhat = pos / np.linalg.norm(pos, axis=1)[:, None]
    axis = dipole_axis / np.linalg.norm(dipole_axis)
    sin2 = 1.0 - (rhat @ axis) ** 2
    sim = np.asarray(fluence, dtype=float) / max(float(np.max(fluence)), 1e-300)
    ref = sin2 / max(float(np.max(sin2)), 1e-300)
    return float(np.max(np.abs(sim - ref)))

