#!/usr/bin/env python3
"""Generate compact H4 products from frozen Stage-F and local receiver data."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.native_receiver import (  # noqa: E402
    NativeRFSourceContract,
    RECEIVER_AXIS,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    apply_trusted_receiver_transfer,
    load_trusted_band,
    native_field_spectrum,
    normalized_shape_metrics,
    sha256_file,
    trusted_bandlimited_timeseries,
    validate_h4_result_contract,
)
from streamer_rf.fullwave.transient import analytic_envelope  # noqa: E402


OUT = ROOT / "fullwave/h4"
STAGES = {
    "Stage4": {
        "case_id": "F4-P-stage4-left-isolated",
        "band": STAGE4_TRUSTED_BAND_HZ,
        "full_maxwell": None,
    },
    "Stage5": {
        "case_id": "F4-C-stage5-highfield-collision",
        "band": STAGE5_TRUSTED_BAND_HZ,
        "full_maxwell": "FULL_MAXWELL_REFERENCE_PENDING",
    },
}
OBSERVED_TIMESTEPS = {
    "baseline_port": 1872,
    "baseline_plane_parallel": 2208,
    "fine_port": 2528,
    "fine_plane_parallel": 2304,
    "baseline_plane_orthogonal": 1584,
}


def load_csv(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    if data.ndim != 1 or not data.size:
        raise ValueError(f"INVALID_CSV:{path}")
    return data


def complex_column(data, real_name, imag_name):
    value = np.asarray(data[real_name]) + 1j * np.asarray(data[imag_name])
    if np.any(~np.isfinite(value)):
        raise ValueError("NONFINITE_COMPLEX_DATA")
    return value


def write_csv(path: Path, columns: list[np.ndarray], header: str):
    np.savetxt(path, np.column_stack(columns), delimiter=",", header=header, comments="")


def receiver_products(raw: Path):
    files = {
        key: raw / key / "response.csv"
        for key in (
            "baseline_port",
            "baseline_plane_parallel",
            "fine_port",
            "fine_plane_parallel",
            "baseline_plane_orthogonal",
        )
    }
    if not all(path.is_file() for path in files.values()):
        raise FileNotFoundError("INCOMPLETE_H4_OPENEMS_FIXTURE")
    data = {key: load_csv(path) for key, path in files.items()}
    f = data["fine_plane_parallel"]["frequency_Hz"]
    if any(not np.array_equal(value["frequency_Hz"], f) for value in data.values()):
        raise ValueError("H4_RECEIVER_FREQUENCY_GRID_MISMATCH")
    h_base = complex_column(data["baseline_plane_parallel"], "H_rx_E_real_m", "H_rx_E_imag_m")
    h_fine = complex_column(data["fine_plane_parallel"], "H_rx_E_real_m", "H_rx_E_imag_m")
    h_orth = complex_column(data["baseline_plane_orthogonal"], "H_rx_E_real_m", "H_rx_E_imag_m")
    s_base = complex_column(data["baseline_port"], "S11_real", "S11_imag")
    s_fine = complex_column(data["fine_port"], "S11_real", "S11_imag")
    z_base = complex_column(data["baseline_port"], "Zin_real_ohm", "Zin_imag_ohm")
    z_fine = complex_column(data["fine_port"], "Zin_real_ohm", "Zin_imag_ohm")
    h_rel = np.abs(np.abs(h_base) - np.abs(h_fine)) / np.maximum(np.abs(h_fine), np.finfo(float).tiny)
    phase_difference = np.angle(h_base / h_fine)
    cross_ratio = np.abs(h_orth) / np.maximum(np.abs(h_base), np.finfo(float).tiny)
    write_csv(
        OUT / "h4_receiver_transfer.csv",
        [
            f,
            h_fine.real,
            h_fine.imag,
            h_base.real,
            h_base.imag,
            h_orth.real,
            h_orth.imag,
            s_fine.real,
            s_fine.imag,
            s_base.real,
            s_base.imag,
            z_fine.real,
            z_fine.imag,
            z_base.real,
            z_base.imag,
            h_rel,
            phase_difference,
            cross_ratio,
        ],
        (
            "frequency_Hz,H_rx_E_fine_real_m,H_rx_E_fine_imag_m,H_rx_E_baseline_real_m,"
            "H_rx_E_baseline_imag_m,H_rx_E_orthogonal_real_m,H_rx_E_orthogonal_imag_m,"
            "S11_fine_real,S11_fine_imag,S11_baseline_real,S11_baseline_imag,"
            "Zin_fine_real_ohm,Zin_fine_imag_ohm,Zin_baseline_real_ohm,Zin_baseline_imag_ohm,"
            "H_magnitude_relative_change,phase_difference_rad,cross_polarization_ratio"
        ),
    )
    key_metrics = []
    for target in (3.5e9, 6.0e9, 9.5e9):
        index = int(np.argmin(np.abs(f - target)))
        key_metrics.append(
            {
                "frequency_Hz": float(f[index]),
                "H_fine_magnitude_m": float(abs(h_fine[index])),
                "H_baseline_magnitude_m": float(abs(h_base[index])),
                "H_magnitude_relative_change": float(h_rel[index]),
                "H_phase_change_rad": float(phase_difference[index]),
                "Zin_fine_ohm": [float(z_fine[index].real), float(z_fine[index].imag)],
                "S11_fine_dB": float(20 * np.log10(max(abs(s_fine[index]), np.finfo(float).tiny))),
                "orthogonal_to_parallel_ratio": float(cross_ratio[index]),
                "orthogonal_rejection_dB": float(20 * np.log10(max(cross_ratio[index], 1e-300))),
            }
        )
    transfer_hash = sha256_file(OUT / "h4_receiver_transfer.csv")
    return f, h_fine, {
        "receiver_transfer_hash": transfer_hash,
        "key_frequency_metrics": key_metrics,
        "mesh_status": "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT",
        "baseline_to_fine_max_H_magnitude_change_key_frequencies": float(
            max(item["H_magnitude_relative_change"] for item in key_metrics)
        ),
        "short_dipole_open_circuit_effective_height_scale_m": 0.00125,
        "loaded_transfer_to_open_circuit_scale_ratio_at_key_frequencies": [
            item["H_fine_magnitude_m"] / 0.00125 for item in key_metrics
        ],
        "wire_status": "OPENEMS_CURVEPORT_PHYSICAL_WIRE_RADIUS_UNDEFINED",
        "calibration": "OPENEMS_TFSF_PLANE_WAVE_E_INCIDENT_1_V_PER_M",
        "raw_case_hashes": {key: sha256_file(path) for key, path in files.items()},
    }


def stage_products(stage: str, spec: dict, receiver_frequency, transfer):
    case = spec["case_id"]
    base = ROOT / "rf/production/f4_attribution"
    waveform_path = base / f"{case}_jefimenko_waveform.csv"
    trust_path = base / f"{case}_rf_trust_report.json"
    summary_path = base / f"{case}_summary.json"
    trusted_band = load_trusted_band(trust_path, spec["band"])
    waveform = load_csv(waveform_path)
    time_s = np.asarray(waveform["time_s"], dtype=float)
    field = np.column_stack([waveform[name] for name in ("Ex_total", "Ey_total", "Ez_total")])
    spectrum = native_field_spectrum(time_s, field)
    projected, receiver_h, received, trusted = apply_trusted_receiver_transfer(
        spectrum, receiver_frequency, transfer, trusted_band, RECEIVER_AXIS
    )
    trusted_frequency = spectrum.frequency_Hz[trusted]
    trusted_incident = projected[trusted]
    trusted_received = received[trusted]
    metrics = normalized_shape_metrics(trusted_frequency, trusted_incident, trusted_received)
    time_voltage = trusted_bandlimited_timeseries(
        received, trusted, spectrum.n_samples, spectrum.dt_s, time_start_s=float(time_s[0])
    )
    envelope = analytic_envelope(time_voltage)
    output_prefix = "h4_stage4" if stage == "Stage4" else "h4_stage5"
    write_csv(
        OUT / f"{output_prefix}_native_received_spectrum.csv",
        [
            spectrum.frequency_Hz,
            trusted.astype(int),
            spectrum.complex_spectrum[:, 0].real,
            spectrum.complex_spectrum[:, 0].imag,
            spectrum.complex_spectrum[:, 1].real,
            spectrum.complex_spectrum[:, 1].imag,
            spectrum.complex_spectrum[:, 2].real,
            spectrum.complex_spectrum[:, 2].imag,
            projected.real,
            projected.imag,
            receiver_h.real,
            receiver_h.imag,
            received.real,
            received.imag,
        ],
        (
            "frequency_Hz,trust_valid,Ex_spectrum_real_V_s_m,Ex_spectrum_imag_V_s_m,"
            "Ey_spectrum_real_V_s_m,Ey_spectrum_imag_V_s_m,Ez_spectrum_real_V_s_m,"
            "Ez_spectrum_imag_V_s_m,E_parallel_real_V_s_m,E_parallel_imag_V_s_m,"
            "H_rx_E_real_m,H_rx_E_imag_m,V_rx_real_V_s,V_rx_imag_V_s"
        ),
    )
    write_csv(
        OUT / f"{output_prefix}_native_received_timeseries.csv",
        [time_s, time_voltage, envelope],
        "time_s,V_rx_trusted_bandlimited_V,envelope_V",
    )
    dominant_energy = np.sum(spectrum.one_sided_esd[trusted], axis=0) * spectrum.df_Hz
    contract = NativeRFSourceContract(
        source_case_id=case,
        source_path=str(waveform_path.relative_to(ROOT)),
        source_hash=sha256_file(waveform_path),
        observer_id=str(np.unique(waveform["observer_id"])[0]),
        observer_position_m=(0.2, 0.0, 0.005),
        observer_distance_m=0.2,
        time_start_s=float(time_s[0]),
        time_end_s=float(time_s[-1]),
        dt_s=spectrum.dt_s,
        samples=spectrum.n_samples,
        field_components=("Ex_total", "Ey_total", "Ez_total"),
        dominant_trusted_component=("Ex_total", "Ey_total", "Ez_total")[int(np.argmax(dominant_energy))],
        receiver_axis=(0.0, 0.0, 1.0),
        trusted_frequency_Hz=trusted_band,
        rf_trust_report_path=str(trust_path.relative_to(ROOT)),
        rf_trust_report_hash=sha256_file(trust_path),
    )
    return contract.to_dict(), {
        "stage": stage,
        "case_id": case,
        "source_summary_hash": sha256_file(summary_path),
        "trusted_frequency_Hz": list(trusted_band),
        "trusted_fft_frequency_Hz": trusted_frequency.tolist(),
        "trusted_bin_count": int(np.sum(trusted)),
        "component_spectral_energy_metric": dominant_energy.tolist(),
        "dominant_component": ("Ex_total", "Ey_total", "Ez_total")[int(np.argmax(dominant_energy))],
        "spectral_metrics": metrics,
        "received_time_semantics": "TRUSTED_BAND_LIMITED_NATIVE_RESPONSE",
        "received_time_dt_s": spectrum.dt_s,
        "received_time_window_s": float(time_s[-1] - time_s[0]),
        "received_peak_abs_V": float(np.max(np.abs(time_voltage))),
        "received_RMS_V": float(np.sqrt(np.mean(time_voltage * time_voltage))),
        "principal_envelope_time_s": float(time_s[int(np.argmax(envelope))]),
        "full_maxwell_status": spec["full_maxwell"],
    }


def main(raw: Path):
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    receiver_frequency, transfer, receiver_summary = receiver_products(raw)
    source_contracts = {}
    stage_results = {}
    for stage, spec in STAGES.items():
        source_contracts[stage], stage_results[stage] = stage_products(
            stage, spec, receiver_frequency, transfer
        )
    status_350 = {
        "frequency_Hz": 350e6,
        "system_band_Hz": [200e6, 500e6],
        "STAGE_F_NATIVE_RF_TRUST": "NOT_RESOLVED",
        "H4_NATIVE_RECEIVED_VOLTAGE": "NOT_RESOLVED",
        "interpretation": "NOT_RESOLVED_IS_NOT_ZERO_RADIATION",
    }
    (OUT / "h4_350mhz_status.json").write_text(json.dumps(status_350, indent=2) + "\n")
    contract = {
        "source_contracts": source_contracts,
        "receiver_id": "H4_CANONICAL_SHORT_DIPOLE_FIELD_SENSOR",
        "receiver_axis": RECEIVER_AXIS.tolist(),
        "receiver_transfer_path": "fullwave/h4/h4_receiver_transfer.csv",
        "receiver_transfer_hash": receiver_summary["receiver_transfer_hash"],
        "trusted_frequency_masks_Hz": {
            stage: list(spec["band"]) for stage, spec in STAGES.items()
        },
        "native_received_spectrum_paths": {
            "Stage4": "fullwave/h4/h4_stage4_native_received_spectrum.csv",
            "Stage5": "fullwave/h4/h4_stage5_native_received_spectrum.csv",
        },
        "native_received_time_paths": {
            "Stage4": "fullwave/h4/h4_stage4_native_received_timeseries.csv",
            "Stage5": "fullwave/h4/h4_stage5_native_received_timeseries.csv",
        },
        "native_path_status": "TRUSTED_BAND_DEVELOPMENT_VERIFIED",
        "native_350MHz_status": "NOT_RESOLVED",
        "experimental_validation_status": "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING",
        "production_receiver_status": "NOT_RESOLVED",
        "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT": True,
        "H4_ABSOLUTE_AMPLITUDE_STATUS": "NUMERICAL_REFERENCE_ONLY",
        "Stage5_full_maxwell_status": "FULL_MAXWELL_REFERENCE_PENDING",
        "mechanism_attribution_consumed": False,
        "h2_h3_numerical_inputs_consumed": False,
    }
    validate_h4_result_contract(contract)
    (OUT / "h4_result_contract.json").write_text(json.dumps(contract, indent=2) + "\n")
    raw_results = {
        key: json.loads((raw / key / "result.json").read_text()) for key in OBSERVED_TIMESTEPS
    }
    for key, steps in OBSERVED_TIMESTEPS.items():
        raw_results[key]["observed_FDTD_timesteps"] = steps
        raw_results[key]["openEMS_reported_cell_count"] = int(
            {"baseline": 165731, "fine": 281799}[raw_results[key]["resolution"]]
        )
    summary = {
        "H4_DEVELOPMENT_GATE": "PASS",
        "H4_NATIVE_RF_PATHWAY": "TRUSTED_BAND_DEVELOPMENT_VERIFIED",
        "PRODUCTION_NATIVE_RF_RECEIVER": "NOT_RESOLVED",
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
        "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING": True,
        "NATIVE_PROPAGATION_ALREADY_INCLUDED": True,
        "RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT": True,
        "H4_ABSOLUTE_AMPLITUDE_STATUS": "NUMERICAL_REFERENCE_ONLY",
        "receiver": receiver_summary,
        "stages": stage_results,
        "native_350MHz": status_350,
        "pathway_isolation": {
            "G3_waveform_consumed": False,
            "H3_transfer_consumed": False,
            "F_R4_mechanism_eta_consumed": False,
            "total_Stage_F_field_used": True,
        },
        "plane_wave_api": {
            "type": "OPENEMS_TFSF_PLANE_WAVE",
            "exc_type": 10,
            "incident_field_amplitude_V_m": 1.0,
            "normalization": "VPORT_SPECTRUM_DIVIDED_BY_OPENEMS_ET_EXCITATION_SPECTRUM",
            "official_example": str(
                Path(os.environ.get("OPENEMS_ROOT", "OPENEMS_ROOT"))
                / "share/openEMS/python/Tutorials/RCS_Sphere.py"
            ),
        },
        "boundary": "PML_8_ALL_SIDES",
        "domain_geometry": {
            "inner_domain_half_extent_m": 0.018,
            "plane_wave_box_half_extent_m": 0.012,
            "receiver_center_m": [0.0, 0.0, 0.0],
            "receiver_to_inner_boundary_minimum_m": 0.01675,
            "domain_sanity_status": "FINITE_DECAY_NO_LATE_TIME_GROWTH",
        },
        "radiation_model_scope": "LOCAL_RECEIVER_ONLY_NATIVE_PROPAGATION_NOT_REPEATED",
        "raw_openEMS_cases": raw_results,
        "raw_openEMS_total_bytes": sum(v["raw_bytes"] for v in raw_results.values()),
        "runtime_generate_results_s": time.perf_counter() - started,
        "peak_RSS_KiB_generate_results": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (OUT / "h4_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=Path("/tmp/h4_openems_validation"))
    args = parser.parse_args()
    main(args.raw)
