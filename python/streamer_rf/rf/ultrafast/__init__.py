"""Ultrafast current-moment diagnostics for the Stage-F v2.1 extension."""

from .processing import (
    CrossingResult,
    ProcessingResult,
    PulseDiagnostics,
    ResamplingResult,
    SpectrumResult,
    analyze_pulses,
    central_difference_second_order,
    compute_esd_spectrum,
    find_relative_db_crossing,
    fit_spectral_slope,
    process_current_moment,
    resample_uniform,
)
from .same_event import (
    TimebaseAudit,
    assess_event_window,
    audit_native_timebases,
    central_derivative_native,
    combined_temporal_trust,
    compare_proxy_to_derivative,
    interior_extremum_time,
    pulse_resolution_status,
)

__all__ = [
    "CrossingResult",
    "ProcessingResult",
    "PulseDiagnostics",
    "ResamplingResult",
    "SpectrumResult",
    "analyze_pulses",
    "central_difference_second_order",
    "compute_esd_spectrum",
    "find_relative_db_crossing",
    "fit_spectral_slope",
    "process_current_moment",
    "resample_uniform",
    "TimebaseAudit",
    "assess_event_window",
    "audit_native_timebases",
    "central_derivative_native",
    "combined_temporal_trust",
    "compare_proxy_to_derivative",
    "interior_extremum_time",
    "pulse_resolution_status",
]
