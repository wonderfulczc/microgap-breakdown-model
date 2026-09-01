from __future__ import annotations

import numpy as np

from streamer_rf.rf.jefimenko.constants import C0, MU0
from .fft import Spectrum, compute_one_sided_spectrum, parseval_time_energy

Z0 = MU0 * C0


def vector_power_esd(spectrum: Spectrum) -> np.ndarray:
    return np.sum(spectrum.one_sided_esd, axis=1)


def spectral_fluence_from_Erad(
    time_s: np.ndarray,
    E_rad_perp: np.ndarray,
    *,
    window: str = "hann",
    far_field_valid: bool = True,
) -> dict[str, object]:
    spectrum = compute_one_sided_spectrum(
        time_s, E_rad_perp, window=window, component_names=("Ex", "Ey", "Ez")
    )
    dFdf = vector_power_esd(spectrum) / Z0
    time_fluence = float(np.sum(parseval_time_energy(time_s, E_rad_perp)) / Z0)
    spectral_fluence = float(np.sum(dFdf) * spectrum.df_Hz)
    return {
        "spectrum": spectrum,
        "dFdf_Jpm2Hz": dFdf,
        "time_fluence_Jpm2": time_fluence,
        "spectral_fluence_Jpm2": spectral_fluence,
        "parseval_relative_error": abs(spectral_fluence - time_fluence) / max(abs(spectral_fluence), abs(time_fluence), 1e-300),
        "spectral_fluence_valid": bool(far_field_valid),
    }

