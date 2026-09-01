"""Validated Jefimenko field solver."""

from .observer import Observer
from .solver import FieldSample, evaluate_observer, evaluate_waveform

__all__ = ["Observer", "FieldSample", "evaluate_observer", "evaluate_waveform"]

