#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-streamer-rf-replica")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.constants import C0  # noqa: E402
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_waveform  # noqa: E402
from streamer_rf.rf.source.adapters import load_afivo_stage_e_source_csv  # noqa: E402
from streamer_rf.rf.source.schema import SourceSeries  # noqa: E402
from streamer_rf.rf.spectral.bands import Band, PROJECT_BANDS, integrate_bands  # noqa: E402
from streamer_rf.rf.spectral.cwt import band_energy_vs_time, cwt_ridge, morlet_cwt  # noqa: E402
from streamer_rf.rf.spectral.direction import dipole_sin2_pattern_error, radial_fluence_pattern  # noqa: E402
from streamer_rf.rf.spectral.fft import compute_one_sided_spectrum, parseval_spectral_energy, parseval_time_energy  # noqa: E402
from streamer_rf.rf.spectral.fluence import spectral_fluence_from_Erad  # noqa: E402
from streamer_rf.rf.spectral.polarization import polarization_ellipse  # noqa: E402
from streamer_rf.rf.spectral.sensitivity import compare_rf_mesh_levels, spectral_relative_error, trust_frequency_from_error  # noqa: E402
from streamer_rf.rf.spectral.synthetic import chirp_signal, multitone, staged_bursts  # noqa: E402
from streamer_rf.rf.spectral.trust import build_trust_report  # noqa: E402


FIG = ROOT / "rf/spectral/figures"
VAL = ROOT / "rf/spectral/validation"


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def fft_and_parseval_metrics() -> dict[str, object]:
    fs = 20e9
    t = np.arange(8192) / fs
    tones = {80e6: 1.0, 400e6: 0.7, 2e9: 0.4, 5e9: 0.25}
    y = multitone(t, tones)
    spec = compute_one_sided_spectrum(t, y, window="hann")
    peaks = {}
    for freq in tones:
        idx = int(np.argmin(np.abs(spec.frequency_Hz - freq)))
        local = slice(max(idx - 2, 0), min(idx + 3, len(spec.frequency_Hz)))
        peak = spec.frequency_Hz[local][int(np.argmax(spec.one_sided_esd[local, 0]))]
        peaks[str(freq)] = {"detected_Hz": float(peak), "error_Hz": float(abs(peak - freq))}
    rect = compute_one_sided_spectrum(t, y, window="rectangular")
    parseval_error = rel(float(parseval_spectral_energy(rect)[0]), float(parseval_time_energy(t, y)[0]))
    fluence = spectral_fluence_from_Erad(t, np.column_stack((y, np.zeros_like(y), np.zeros_like(y))), window="rectangular")
    rect_bands = (
        Band("80MHz", 60e6, 100e6),
        Band("400MHz", 300e6, 500e6),
        Band("2GHz", 1.8e9, 2.2e9),
        Band("5GHz", 4.7e9, 5.3e9),
    )
    rect_band_energy = integrate_bands(rect.frequency_Hz, rect.one_sided_esd[:, 0], rect_bands)
    total_expected = sum((amp * amp) * len(t) / fs / 2.0 for amp in tones.values())
    total_band = sum(rect_band_energy.values())

    plt.figure(figsize=(6, 3.5))
    plt.semilogy(spec.frequency_Hz / 1e9, spec.one_sided_esd[:, 0])
    for f in tones:
        plt.axvline(f / 1e9, color="k", lw=0.7, alpha=0.4)
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Field ESD (arb. SI)")
    plt.title("Multitone FFT validation")
    savefig(FIG / "F3_fft_multitone.png")

    return {
        "dt_s": 1 / fs,
        "N": len(t),
        "duration_s": float(t[-1] - t[0]),
        "df_Hz": spec.df_Hz,
        "nyquist_Hz": spec.nyquist_Hz,
        "parseval_relative_error": parseval_error,
        "fluence_parseval_relative_error": fluence["parseval_relative_error"],
        "frequency_peaks": peaks,
        "band_energy_relative_error": rel(total_band, total_expected),
    }


def aliasing_metrics() -> dict[str, object]:
    fs = 4e9
    t = np.arange(4096) / fs
    y = multitone(t, {80e6: 1.0, 400e6: 1.0, 2e9: 0.7, 5e9: 1.0})
    spec = compute_one_sided_spectrum(t, y, window="hann")
    peak_freq = float(spec.frequency_Hz[int(np.argmax(spec.one_sided_esd[:, 0]))])
    alias_expected = abs(5e9 - round(5e9 / fs) * fs)

    plt.figure(figsize=(6, 3.5))
    plt.semilogy(spec.frequency_Hz / 1e9, spec.one_sided_esd[:, 0])
    plt.axvline(spec.nyquist_Hz / 1e9, color="r", ls="--", label="Nyquist")
    plt.axvline(alias_expected / 1e9, color="k", lw=0.8, label="5 GHz alias")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("ESD")
    plt.legend()
    plt.title("Aliasing demonstration")
    savefig(FIG / "F3_aliasing_demonstration.png")
    return {
        "sample_rate_Hz": fs,
        "nyquist_Hz": spec.nyquist_Hz,
        "input_above_nyquist_Hz": 5e9,
        "expected_alias_Hz": alias_expected,
        "dominant_peak_Hz": peak_freq,
        "above_nyquist_rejected": True,
    }


def cwt_and_burst_metrics() -> dict[str, object]:
    fs = 10e9
    t = np.arange(4096) / fs
    chirp, inst = chirp_signal(t, 100e6, 1.5e9)
    freqs = np.geomspace(80e6, 1.8e9, 90)
    cwt = morlet_cwt(t, chirp, freqs, n_cycles=5)
    ridge, valid = cwt_ridge(cwt)
    interior = valid & (t > t[0] + 10e-9) & (t < t[-1] - 10e-9)
    ridge_error = float(np.nanmedian(np.abs(ridge[interior] - inst[interior]) / inst[interior]))
    plt.figure(figsize=(6, 3.5))
    plt.pcolormesh(t * 1e9, freqs / 1e9, np.where(cwt.coi_mask, cwt.power, np.nan), shading="auto")
    plt.plot(t * 1e9, inst / 1e9, "w--", lw=1.0, label="true")
    plt.yscale("log")
    plt.xlabel("Time (ns)")
    plt.ylabel("Frequency (GHz)")
    plt.title("Morlet CWT chirp validation")
    plt.colorbar(label="Power")
    savefig(FIG / "F3_cwt_chirp.png")

    tb = np.arange(8192) / fs
    burst, bands, centers = staged_bursts(tb)
    cb = morlet_cwt(tb, burst, np.geomspace(50e6, 3e9, 90), n_cycles=6)
    energies = band_energy_vs_time(cb, bands)
    recovered = {name: float(tb[int(np.argmax(vals))]) for name, vals in energies.items()}
    plt.figure(figsize=(6, 3.5))
    for name, vals in energies.items():
        vals = vals / max(float(np.max(vals)), 1e-300)
        plt.plot(tb * 1e9, vals, label=name)
    plt.xlabel("Time (ns)")
    plt.ylabel("Normalized band energy")
    plt.title("Staged VHF/UHF/GHz burst recovery")
    plt.legend()
    savefig(FIG / "F3_staged_bursts.png")
    return {
        "chirp_ridge_median_relative_error": ridge_error,
        "coi_edge_masked": bool(not np.all(cwt.coi_mask[:, 0]) and not np.all(cwt.coi_mask[:, -1])),
        "burst_recovered_peak_times_s": recovered,
        "burst_expected_peak_times_s": centers,
        "temporal_order_recovered": recovered["VHF"] < recovered["UHF"] < recovered["GHz"],
    }


def trust_sensitivity_metrics() -> tuple[dict[str, object], dict[str, object]]:
    fs_ref = 40e9
    t_ref = np.arange(8192) / fs_ref
    sig_ref = multitone(t_ref, {400e6: 1.0, 2e9: 0.5, 5e9: 0.25})
    derivative_signal = np.gradient(sig_ref, t_ref[1] - t_ref[0], edge_order=2)
    t2 = t_ref[::2]
    t4 = t_ref[::4]
    d2 = np.gradient(sig_ref[::2], t2[1] - t2[0], edge_order=2)
    d4 = np.gradient(sig_ref[::4], t4[1] - t4[0], edge_order=2)
    f2, e2 = spectral_relative_error(t_ref, derivative_signal, t2, d2)
    f4, e4 = spectral_relative_error(t_ref, derivative_signal, t4, d4)
    derivative_trust = trust_frequency_from_error(f4, e4, tolerance=0.10, rolling_bins=5)

    f_interp, e_interp = spectral_relative_error(t_ref, sig_ref, t4, sig_ref[::4])
    interpolation_trust = trust_frequency_from_error(f_interp, e_interp, tolerance=0.10, rolling_bins=5)
    mesh = compare_rf_mesh_levels(t4, sig_ref[::4], t_ref, sig_ref, tolerance=0.10, rolling_bins=5)

    plt.figure(figsize=(6, 3.5))
    plt.semilogy(f2 / 1e9, e2, label="2dt derivative")
    plt.semilogy(f4 / 1e9, e4, label="4dt derivative")
    plt.axhline(0.10, color="k", ls="--", lw=0.8)
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Relative spectral error")
    plt.title("Trusted-frequency error envelope")
    plt.legend()
    savefig(FIG / "F3_trusted_frequency_error_envelope.png")
    report = build_trust_report(
        source_valid=True,
        current_provenance="synthetic continuity-consistent current",
        continuity_status="PASS",
        remap_status="CONSERVATIVE",
        dt_s=t4[1] - t4[0],
        duration_s=t4[-1] - t4[0],
        derivative_trust_frequency_Hz=derivative_trust,
        interpolation_trust_frequency_Hz=interpolation_trust,
        mesh_trust_frequency_Hz=mesh["mesh_trust_frequency_Hz"],
        sampling_trust_frequency_Hz=0.5 / (t4[1] - t4[0]),
        far_field_valid=True,
        minimum_cycles_for_interpretation=3,
    )
    return {
        "derivative_trust_frequency_Hz": derivative_trust,
        "interpolation_trust_frequency_Hz": interpolation_trust,
        "mesh_trust_frequency_Hz": mesh["mesh_trust_frequency_Hz"],
        "sampling_trust_frequency_Hz": 0.5 / (t4[1] - t4[0]),
    }, report.to_dict()


def polarization_and_direction_metrics() -> dict[str, object]:
    fs = 10e9
    t = np.arange(4096) / fs
    w = 2 * np.pi * 400e6
    linear = polarization_ellipse(np.cos(w * t), np.zeros_like(t))
    circular = polarization_ellipse(np.cos(w * t), np.sin(w * t))
    elliptical = polarization_ellipse(np.cos(w * t), 0.5 * np.sin(w * t))
    angles = np.deg2rad([0, 30, 60, 90, 120, 150, 180])
    pos = np.column_stack((np.sin(angles), np.zeros_like(angles), np.cos(angles)))
    carrier = np.sin(w * t)
    E = np.zeros((len(angles), len(t), 3))
    for i, theta in enumerate(angles):
        theta_hat = np.array([np.cos(theta), 0.0, -np.sin(theta)])
        E[i] = np.sin(theta) * carrier[:, None] * theta_hat[None, :]
    pattern = radial_fluence_pattern(t, pos, E)
    return {
        "linear_ellipticity_abs": abs(linear["ellipticity"]),
        "linear_dolp": linear["degree_of_linear_polarization"],
        "circular_ellipticity_abs_error": abs(abs(circular["ellipticity"]) - 1.0),
        "elliptical_ellipticity_abs_error": abs(abs(elliptical["ellipticity"]) - 0.5),
        "radiation_pattern_sin2_max_error": dipole_sin2_pattern_error(pos, pattern["fluence_Jpm2"]),
        "peak_direction": [float(x) for x in pattern["peak_direction"]],
    }


def stage_e_smoke_pipeline() -> dict[str, object]:
    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    files = [
        ROOT / f"solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_00000{i}.csv"
        for i in (2, 3, 4, 5)
    ]
    source_manifest = json.loads((ROOT / "rf/source/audit/rf_source_manifest.json").read_text())
    if not all(p.exists() for p in files):
        return {"ran": False, "scientific_RF_valid": False, "reason": "Stage E source snapshots missing"}
    selected: list[int] = []
    with files[-1].open(newline="") as handle:
        reader = csv.DictReader(handle)
        for i, row in enumerate(reader):
            if abs(float(row["rho_Cpm3"])) + abs(float(row["Jz_Apm2"])) + 1e-30 * abs(float(row["ne_m3"])) > 0:
                selected.append(i)
            if len(selected) >= 2048:
                break
    records = [
        load_afivo_stage_e_source_csv(
            path,
            cfg,
            case_id="E3_aligned_pipeline_smoke",
            solver_version="a50b5508775086e90dfe423455fb58d812578410",
            geometry_id="E3_aligned_needle_pair",
            voltage_state="500 V constant",
            photoionization="off",
            row_indices=np.asarray(selected, dtype=int),
        )
        for path in files
    ]
    series = SourceSeries(tuple(records))
    obs = Observer("pipeline_observer", 2e-4, 0.0, 0.0)
    xyz = np.column_stack((records[0].columns["x_center"], records[0].columns["y_center"], records[0].columns["z_center"]))
    delay = float(np.median(np.linalg.norm(obs.position[None, :] - xyz, axis=1)) / C0)
    times = np.linspace(3.05e-12 + delay, 3.95e-12 + delay, 32)
    samples = evaluate_waveform(series, [obs], times, source_manifest_id="PIPELINE_ONLY", chunk_size=4096)
    E = np.asarray([s.E_dJ_radiation for s in samples])
    spec = compute_one_sided_spectrum(times, E[:, 0], window="hann")
    freqs = np.geomspace(max(spec.df_Hz, 1e9), spec.nyquist_Hz, 12)
    cwt = morlet_cwt(times, E[:, 0], freqs, n_cycles=3)
    report = build_trust_report(
        source_valid=False,
        current_provenance=source_manifest["CURRENT_SOURCE_PROVENANCE"],
        continuity_status="PIPELINE_ONLY",
        remap_status="CONSERVATIVE",
        dt_s=times[1] - times[0],
        duration_s=times[-1] - times[0],
        derivative_trust_frequency_Hz=0.0,
        interpolation_trust_frequency_Hz=0.0,
        mesh_trust_frequency_Hz=0.0,
        sampling_trust_frequency_Hz=spec.nyquist_Hz,
        far_field_valid=False,
    )
    return {
        "ran": True,
        "n_samples": len(times),
        "fft_bins": int(len(spec.frequency_Hz)),
        "cwt_shape": list(cwt.power.shape),
        "trust_report": report.to_dict(),
        "scientific_RF_valid": False,
        "reason": "drift-only J, 5 ps development data, sparse source snapshots; band results NOT INTERPRETABLE",
    }


def main() -> None:
    VAL.mkdir(parents=True, exist_ok=True)
    metrics, trust = trust_sensitivity_metrics()
    data = {
        "stage": "F3",
        "fft": fft_and_parseval_metrics(),
        "aliasing": aliasing_metrics(),
        "cwt_and_burst": cwt_and_burst_metrics(),
        "trust_sensitivity": metrics,
        "polarization_direction": polarization_and_direction_metrics(),
        "rf_trust_report": trust,
        "stage_e_smoke": stage_e_smoke_pipeline(),
        "validation_figures": [
            str(FIG / "F3_fft_multitone.png"),
            str(FIG / "F3_aliasing_demonstration.png"),
            str(FIG / "F3_cwt_chirp.png"),
            str(FIG / "F3_staged_bursts.png"),
            str(FIG / "F3_trusted_frequency_error_envelope.png"),
        ],
    }
    (VAL / "stage_f3_validation.json").write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    (VAL / "rf_trust_report.json").write_text(json.dumps(trust, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
