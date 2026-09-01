from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from streamer_rf.rf.source.schema import SourceRecord


@dataclass(frozen=True)
class Observer:
    observer_id: str
    x_m: float
    y_m: float
    z_m: float

    @property
    def position(self) -> np.ndarray:
        return np.asarray([self.x_m, self.y_m, self.z_m], dtype=float)

    def validate_outside_source(self, record: SourceRecord) -> None:
        lo, hi = record.bounds()
        p = self.position
        inside = np.all((lo <= p) & (p <= hi), axis=1)
        if np.any(inside):
            raise ValueError(f"observer {self.observer_id} is inside a source cell")

