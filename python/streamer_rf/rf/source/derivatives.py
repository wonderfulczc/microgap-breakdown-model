from __future__ import annotations

import numpy as np


def central_difference_three_point(times: np.ndarray, values: np.ndarray, index: int = 1) -> np.ndarray:
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if len(times) != 3:
        raise ValueError("central_difference_three_point requires exactly three times")
    if index != 1:
        raise ValueError("only the center index is supported")
    if np.any(np.diff(times) <= 0):
        raise ValueError("times must be strictly increasing")

    t0, t1, t2 = times
    y0, y1, y2 = values[0], values[1], values[2]
    a = (t1 - t2) / ((t0 - t1) * (t0 - t2))
    b = (2.0 * t1 - t0 - t2) / ((t1 - t0) * (t1 - t2))
    c = (t1 - t0) / ((t2 - t0) * (t2 - t1))
    return a * y0 + b * y1 + c * y2

