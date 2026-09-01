from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.spectral.bands import Band, PROJECT_BANDS, classify_band, integrate_band, integrate_bands  # noqa: E402
from streamer_rf.rf.spectral.cwt import band_energy_vs_time, cwt_ridge, morlet_cwt  # noqa: E402
from streamer_rf.rf.spectral.direction import dipole_sin2_pattern_error, radial_fluence_pattern  # noqa: E402
from streamer_rf.rf.spectral.fft import (  # noqa: E402
    compute_one_sided_spectrum,
    make_window,
    parseval_spectral_energy,
    parseval_time_energy,
)
from streamer_rf.rf.spectral.fluence import spectral_fluence_from_Erad  # noqa: E402
from streamer_rf.rf.spectral.polarization import polarization_ellipse  # noqa: E402
from streamer_rf.rf.spectral.sensitivity import compare_rf_mesh_levels, spectral_relative_error, trust_frequency_from_error  # noqa: E402
from streamer_rf.rf.spectral.synthetic import chirp_signal, multitone, staged_bursts  # noqa: E402
from streamer_rf.rf.spectral.trust import build_trust_report  # noqa: E402


def test_fft_frequency_and_convention() -> None:
    fs = 20e9
    dt = 1.0 / fs
    t = np.arange(4096) * dt
    f0 = 400e6
    y = np.sin(2.0 * np.pi * f0 * t)
    spec = compute_one_sided_spectrum(t, y, window="rectangular")
    f_peak = spec.frequency_Hz[int(np.argmax(spec.one_sided_esd[:, 0]))]
    assert abs(f_peak - f0) <= spec.df_Hz
    assert spec.df_Hz == fs / len(t)
    assert spec.nyquist_Hz == fs / 2.0


def test_parseval_and_esd_normalization() -> None:
    fs = 10e9
    t = np.arange(8192) / fs
    y = multitone(t, {80e6: 1.0, 400e6: 0.5, 2e9: 0.25})
    spec = compute_one_sided_spectrum(t, y, window="rectangular")
    time_energy = parseval_time_energy(t, y)[0]
    spec_energy = parseval_spectral_energy(spec)[0]
    assert abs(spec_energy - time_energy) / time_energy < 1e-12


def test_window_normalization_for_energy_and_amplitude() -> None:
    n = 4096
    fs = 4096e6
    t = np.arange(n) / fs
    f0 = 256e6
    y = 2.0 * np.sin(2.0 * np.pi * f0 * t)
    _, info = make_window(n, "hann")
    assert 0.49 < info.coherent_gain < 0.51
    assert 2.65 < info.energy_correction < 2.68
    spec = compute_one_sided_spectrum(t, y, window="hann")
    amp_peak = float(np.max(spec.amplitude[:, 0])) / (t[-1] - t[0])
    assert abs(amp_peak - 2.0) < 0.03


def test_spectral_fluence_parseval_consistency() -> None:
    fs = 20e9
    t = np.arange(4096) / fs
    E = np.column_stack((np.sin(2 * np.pi * 400e6 * t), np.zeros_like(t), np.zeros_like(t)))
    result = spectral_fluence_from_Erad(t, E, window="rectangular", far_field_valid=True)
    assert result["parseval_relative_error"] < 1e-12
    assert result["spectral_fluence_valid"] is True


def test_band_integration_and_project_band_config() -> None:
    f = np.linspace(0.0, 2e9, 2001)
    density = np.ones_like(f)
    out = integrate_bands(f, density, (Band("test", 100e6, 200e6),))
    assert abs(out["test"] - 100e6) < 1e-6
    assert any(b.name == "60-90 MHz" for b in PROJECT_BANDS)


def test_project_band_status_partial_trusted() -> None:
    assert classify_band(Band("partial", 100e6, 200e6), 150e6, 500e6) == "PARTIALLY_TRUSTED"


def test_low_frequency_cycle_bound_is_configurable() -> None:
    report = build_trust_report(
        source_valid=True,
        current_provenance="synthetic",
        continuity_status="PASS",
        remap_status="CONSERVATIVE",
        dt_s=1e-10,
        duration_s=1e-6,
        derivative_trust_frequency_Hz=1e9,
        interpolation_trust_frequency_Hz=1e9,
        mesh_trust_frequency_Hz=1e9,
        sampling_trust_frequency_Hz=1e9,
        far_field_valid=True,
        minimum_cycles_for_interpretation=5,
    )
    assert report.low_frequency_interpretation_bound_Hz == 5e6


def test_fluence_far_field_gate() -> None:
    t = np.arange(1024) / 10e9
    E = np.column_stack((np.sin(2 * np.pi * 400e6 * t), np.zeros_like(t), np.zeros_like(t)))
    assert spectral_fluence_from_Erad(t, E, window="rectangular", far_field_valid=False)["spectral_fluence_valid"] is False


def test_aliasing_detection_and_nyquist_rejection() -> None:
    fs = 4e9
    t = np.arange(4096) / fs
    y = multitone(t, {80e6: 1.0, 400e6: 1.0, 2e9: 1.0, 5e9: 1.0})
    spec = compute_one_sided_spectrum(t, y, window="hann")
    assert np.isclose(spec.nyquist_Hz, 2e9)
    aliases = [abs(spec.frequency_Hz[np.argmax(spec.one_sided_esd[:, 0])] - 1e9) < 2 * spec.df_Hz]
    assert any(aliases)
    report = build_trust_report(
        source_valid=True,
        current_provenance="synthetic",
        continuity_status="PASS",
        remap_status="CONSERVATIVE",
        dt_s=1 / fs,
        duration_s=t[-1] - t[0],
        derivative_trust_frequency_Hz=10e9,
        interpolation_trust_frequency_Hz=10e9,
        mesh_trust_frequency_Hz=10e9,
        sampling_trust_frequency_Hz=spec.nyquist_Hz,
        far_field_valid=True,
    )
    assert report.trusted_frequency_high_Hz == spec.nyquist_Hz


def test_derivative_and_interpolation_trust_frequency_interfaces() -> None:
    fs_ref = 40e9
    t_ref = np.arange(8192) / fs_ref
    sig_ref = multitone(t_ref, {400e6: 1.0, 2e9: 0.5})
    t_cand = t_ref[::4]
    sig_cand = sig_ref[::4]
    f, err = spectral_relative_error(t_ref, sig_ref, t_cand, sig_cand)
    trust = trust_frequency_from_error(f, err, tolerance=0.10, rolling_bins=5)
    assert 0.0 < trust < 0.5 / (t_cand[1] - t_cand[0])
    mesh = compare_rf_mesh_levels(t_cand, sig_cand, t_ref, sig_ref, tolerance=0.10, rolling_bins=5)
    assert mesh["mesh_trust_frequency_Hz"] == trust


def test_cwt_chirp_ridge_and_coi_mask() -> None:
    fs = 8e9
    t = np.arange(2048) / fs
    y, inst = chirp_signal(t, 100e6, 1.5e9)
    freqs = np.geomspace(80e6, 1.8e9, 80)
    cwt = morlet_cwt(t, y, freqs, n_cycles=5.0)
    ridge, valid = cwt_ridge(cwt)
    interior = valid & (t > t[0] + 60 / fs) & (t < t[-1] - 60 / fs)
    rel = np.nanmedian(np.abs(ridge[interior] - inst[interior]) / inst[interior])
    assert rel < 0.12
    assert not np.all(cwt.coi_mask[:, 0])
    assert not np.all(cwt.coi_mask[:, -1])


def test_cwt_log_frequency_grid_and_power_shape() -> None:
    t = np.arange(512) / 5e9
    y = np.sin(2 * np.pi * 300e6 * t)
    freqs = np.geomspace(100e6, 1e9, 24)
    cwt = morlet_cwt(t, y, freqs)
    assert cwt.power.shape == (24, 512)
    assert np.all(np.diff(np.log(cwt.frequency_Hz)) > 0)


def test_staged_burst_temporal_recovery() -> None:
    fs = 10e9
    t = np.arange(8192) / fs
    y, bands, centers = staged_bursts(t)
    freqs = np.geomspace(50e6, 3e9, 90)
    cwt = morlet_cwt(t, y, freqs, n_cycles=6.0)
    energies = band_energy_vs_time(cwt, bands)
    recovered = {name: t[int(np.argmax(vals))] for name, vals in energies.items()}
    assert recovered["VHF"] < recovered["UHF"] < recovered["GHz"]
    for name in centers:
        assert abs(recovered[name] - centers[name]) < 8e-9


def test_band_energy_vs_time_nonnegative() -> None:
    t = np.arange(1024) / 10e9
    y, bands, _ = staged_bursts(t)
    cwt = morlet_cwt(t, y, np.geomspace(50e6, 3e9, 48))
    energies = band_energy_vs_time(cwt, bands)
    assert all(np.all(v >= 0) for v in energies.values())


def test_polarization_linear_circular_elliptical_and_nearfield_gate() -> None:
    fs = 10e9
    t = np.arange(4096) / fs
    w = 2.0 * np.pi * 400e6
    linear = polarization_ellipse(np.cos(w * t), np.zeros_like(t))
    assert linear["polarization_valid"] is True
    assert abs(linear["ellipticity"]) < 1e-3
    assert linear["degree_of_linear_polarization"] > 0.99
    circ = polarization_ellipse(np.cos(w * t), np.sin(w * t))
    assert abs(abs(circ["ellipticity"]) - 1.0) < 0.03
    ell = polarization_ellipse(np.cos(w * t), 0.5 * np.sin(w * t))
    assert 0.45 < abs(ell["ellipticity"]) < 0.55
    assert polarization_ellipse(np.cos(w * t), np.sin(w * t), far_field_valid=False)["polarization_valid"] is False


def test_dipole_radiation_direction_pattern() -> None:
    t = np.arange(2048) / 10e9
    angles = np.deg2rad([0, 30, 60, 90, 120, 150, 180])
    pos = np.column_stack((np.sin(angles), np.zeros_like(angles), np.cos(angles)))
    y = np.sin(2.0 * np.pi * 400e6 * t)
    E = np.zeros((len(angles), len(t), 3))
    for i, rhat in enumerate(pos):
        theta_hat = np.array([np.cos(angles[i]), 0.0, -np.sin(angles[i])])
        E[i] = np.sin(angles[i]) * y[:, None] * theta_hat[None, :]
    pat = radial_fluence_pattern(t, pos, E)
    assert pat["peak_index"] == 3
    assert dipole_sin2_pattern_error(pos, pat["fluence_Jpm2"]) < 1e-12


def test_rf_trust_report_intersection_and_band_status() -> None:
    report = build_trust_report(
        source_valid=True,
        current_provenance="J_RF",
        continuity_status="PASS",
        remap_status="CONSERVATIVE",
        dt_s=1e-11,
        duration_s=1e-6,
        derivative_trust_frequency_Hz=3e9,
        interpolation_trust_frequency_Hz=2e9,
        mesh_trust_frequency_Hz=5e9,
        sampling_trust_frequency_Hz=10e9,
        far_field_valid=True,
        minimum_cycles_for_interpretation=3,
    )
    assert report.trusted_frequency_low_Hz == 3e6
    assert report.trusted_frequency_high_Hz == 2e9
    assert report.limiting_factor == "interpolation"
    assert classify_band(Band("x", 60e6, 90e6), report.trusted_frequency_low_Hz, report.trusted_frequency_high_Hz) == "TRUSTED"
    assert classify_band(Band("x", 3e9, 10e9), report.trusted_frequency_low_Hz, report.trusted_frequency_high_Hz) == "UNTRUSTED"


def test_nonconservative_remap_rejects_scientific_rf() -> None:
    report = build_trust_report(
        source_valid=True,
        current_provenance="J_RF",
        continuity_status="PASS",
        remap_status="NEAREST",
        dt_s=1e-11,
        duration_s=1e-6,
        derivative_trust_frequency_Hz=1e9,
        interpolation_trust_frequency_Hz=1e9,
        mesh_trust_frequency_Hz=1e9,
        sampling_trust_frequency_Hz=1e9,
        far_field_valid=True,
    )
    assert report.scientific_rf_valid is False


def test_window_rectangular_metadata() -> None:
    w, info = make_window(16, "rectangular")
    assert np.all(w == 1.0)
    assert info.coherent_gain == 1.0
    assert info.energy_correction == 1.0


def test_untrusted_source_rejection_and_stage_e_pipeline_smoke() -> None:
    manifest = json.loads((ROOT / "rf/source/audit/rf_source_manifest.json").read_text())
    report = build_trust_report(
        source_valid=False,
        current_provenance=manifest["CURRENT_SOURCE_PROVENANCE"],
        continuity_status="PIPELINE_ONLY",
        remap_status="CONSERVATIVE",
        dt_s=5e-12,
        duration_s=5e-12,
        derivative_trust_frequency_Hz=0.0,
        interpolation_trust_frequency_Hz=0.0,
        mesh_trust_frequency_Hz=0.0,
        sampling_trust_frequency_Hz=100e9,
        far_field_valid=False,
    )
    assert report.scientific_rf_valid is False
    csv_path = ROOT / "rf/jefimenko/validation/real_afivo_pipeline_smoke.csv"
    if csv_path.exists():
        with csv_path.open(newline="") as handle:
            row = next(csv.DictReader(handle))
        t = np.array([0.0, 5e-12])
        vals = np.array([0.0, float(row["Ex_dJ"])])
        spec = compute_one_sided_spectrum(t, vals, window="rectangular")
        assert np.all(np.isfinite(spec.one_sided_esd))
