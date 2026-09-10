"""Generate the H2 350 MHz theory/openEMS comparison without running openEMS."""
import hashlib
import json
from pathlib import Path
import re
import resource
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.receiver import (  # noqa: E402
    read_touchstone,
    receiver_validation_status,
    validate_exact_frequency_grid,
)

RAW = Path("/tmp/h2_openems_validation/development_350mhz")
THEORY_NAMES = (
    "H2_THEORY_350MHz_RX_S11.s1p",
    "H2_THEORY_350MHz_TX_S11.s1p",
    "H2_THEORY_350MHz_TX_RX_S21.s2p",
    "H2_THEORY_350MHz_reference.csv",
    "H2_THEORY_350MHz_metadata.json",
)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def locate(name):
    for candidate in (ROOT / name, ROOT / "fullwave/h2/theory_inputs" / name, Path("/tmp") / name):
        if candidate.exists():
            return candidate
    return None


def runtime_stats():
    text = Path("/tmp/h2_openems_validation/development_350mhz.log").read_text()
    cells = re.search(r"Time for\s+(\d+) iterations with\s+([0-9.]+) cells", text)
    dt = re.search(r"FDTD timestep is:\s+([0-9.e+-]+) s", text)
    if not cells or not dt:
        raise ValueError("INCOMPLETE_OPENEMS_LOG")
    return {"actual_timesteps": int(cells.group(1)), "actual_FDTD_cells": int(float(cells.group(2))),
            "actual_FDTD_timestep_s": float(dt.group(1))}


def finite_complex(values, label):
    if np.any(~np.isfinite(values)):
        raise ValueError(f"NONFINITE_{label}")


def impedance_from_s11(s11, z0):
    denominator = 1.0 - s11
    if np.any(np.abs(denominator) <= np.finfo(float).eps):
        raise ValueError("S11_IMPEDANCE_SINGULAR")
    return z0 * (1.0 + s11) / denominator


def db(values):
    return 20.0 * np.log10(np.maximum(np.abs(values), np.finfo(float).tiny))


def validate_metadata(metadata):
    expected = {
        "status": "THEORETICAL_SURROGATE_ONLY",
        "frequency_start_Hz": 200e6,
        "frequency_stop_Hz": 500e6,
        "points": 601,
        "frequency_step_Hz": 0.5e6,
        "center_frequency_Hz": 350e6,
        "reference_impedance_ohm": 50.0,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"THEORY_METADATA_MISMATCH_{key}")
    if receiver_validation_status(metadata["status"]) != "SIMULATION_ONLY":
        raise ValueError("THEORY_PROVENANCE_ESCALATION")


def validate_reference_csv(path, frequency, s11, s21, z0):
    reference = np.genfromtxt(path, delimiter=",", names=True)
    required = {
        "frequency_Hz", "Zin_real_ohm", "Zin_imag_ohm", "S11_real", "S11_imag",
        "S11_dB", "S21_real", "S21_imag", "S21_dB", "S21_phase_deg",
    }
    if not required.issubset(reference.dtype.names or ()):
        raise ValueError("INCOMPLETE_THEORY_REFERENCE_CSV")
    validate_exact_frequency_grid(reference["frequency_Hz"], 200e6, 500e6, 601)
    for name in required:
        if np.any(~np.isfinite(reference[name])):
            raise ValueError(f"NONFINITE_THEORY_REFERENCE_{name}")
    csv_s11 = reference["S11_real"] + 1j * reference["S11_imag"]
    csv_s21 = reference["S21_real"] + 1j * reference["S21_imag"]
    csv_zin = reference["Zin_real_ohm"] + 1j * reference["Zin_imag_ohm"]
    if not np.allclose(reference["frequency_Hz"], frequency, rtol=0, atol=1e-6):
        raise ValueError("THEORY_CSV_FREQUENCY_MISMATCH")
    if not np.allclose(csv_s11, s11, rtol=2e-8, atol=2e-9):
        raise ValueError("THEORY_CSV_S11_MISMATCH")
    if not np.allclose(csv_s21, s21, rtol=2e-8, atol=2e-9):
        raise ValueError("THEORY_CSV_S21_MISMATCH")
    if not np.allclose(csv_zin, impedance_from_s11(s11, z0), rtol=2e-8, atol=2e-6):
        raise ValueError("THEORY_CSV_ZIN_MISMATCH")
    return reference


def main():
    started = time.perf_counter()
    out = Path(__file__).parent
    result = json.loads((RAW / "result.json").read_text())
    result.update(runtime_stats())
    found = {name: locate(name) for name in THEORY_NAMES}
    missing = [name for name, path in found.items() if path is None]
    theory_status = "MISSING_REQUIRED_FILES" if missing else "AVAILABLE_PENDING_VALIDATION"
    theory = {"status": theory_status, "required_files": list(THEORY_NAMES), "missing_files": missing,
              "required_provenance": "THEORETICAL_SURROGATE_ONLY", "VNA_gate_satisfied": False}
    comparisons = None
    if not missing:
        metadata = json.loads(found["H2_THEORY_350MHz_metadata.json"].read_text())
        validate_metadata(metadata)
        parsed = {name: read_touchstone(path) for name, path in found.items() if path.suffix.lower() in (".s1p", ".s2p")}
        for data in parsed.values():
            validate_exact_frequency_grid(data.frequency_Hz, 200e6, 500e6, 601)
            if data.Z0_ohm != 50:
                raise ValueError("THEORY_REFERENCE_IMPEDANCE_MISMATCH")
            for name, values in data.parameters.items():
                finite_complex(values, f"THEORY_{name}")
        theory.update({"status": "VALID_THEORETICAL_SURROGATE", "metadata": metadata,
                       "validation_status": receiver_validation_status(metadata["status"]),
                       "VNA_gate_satisfied": False})
        numeric = np.genfromtxt(RAW / "sparameters.csv", delimiter=",", names=True)
        validate_exact_frequency_grid(numeric["frequency_Hz"], 200e6, 500e6, 601)
        rx = parsed["H2_THEORY_350MHz_RX_S11.s1p"]
        tx = parsed["H2_THEORY_350MHz_TX_S11.s1p"]
        pair = parsed["H2_THEORY_350MHz_TX_RX_S21.s2p"]
        if not np.allclose(rx.parameters["S11"], tx.parameters["S11"], rtol=0, atol=1e-14):
            raise ValueError("THEORY_TX_RX_S11_MISMATCH")
        if not np.allclose(pair.parameters["S11"], tx.parameters["S11"], rtol=0, atol=1e-14):
            raise ValueError("THEORY_S2P_S11_ORDERING_MISMATCH")
        if not np.allclose(pair.parameters["S22"], rx.parameters["S11"], rtol=0, atol=1e-14):
            raise ValueError("THEORY_S2P_S22_ORDERING_MISMATCH")
        s11_openems = numeric["S11_real"] + 1j * numeric["S11_imag"]
        s21_openems = numeric["S21_real"] + 1j * numeric["S21_imag"]
        s11_theory = tx.parameters["S11"]
        s21_theory = pair.parameters["S21"]
        finite_complex(s11_openems, "OPENEMS_S11")
        finite_complex(s21_openems, "OPENEMS_S21")
        validate_reference_csv(found["H2_THEORY_350MHz_reference.csv"], tx.frequency_Hz,
                               s11_theory, s21_theory, tx.Z0_ohm)
        idx = int(np.flatnonzero(tx.frequency_Hz == 350e6)[0])
        theory_feature_idx = int(np.argmin(abs(s11_theory)))
        openems_feature_idx = int(np.argmin(abs(s11_openems)))
        zin_theory = impedance_from_s11(s11_theory, tx.Z0_ohm)
        zin_openems = impedance_from_s11(s11_openems, tx.Z0_ohm)
        feature_mismatch_Hz = float(numeric["frequency_Hz"][openems_feature_idx] - tx.frequency_Hz[theory_feature_idx])
        corr11 = float(np.corrcoef(abs(s11_theory), abs(s11_openems))[0, 1])
        corr21 = float(np.corrcoef(abs(s21_theory), abs(s21_openems))[0, 1])
        comparisons = {
            "theory_S11_feature_frequency_Hz": float(tx.frequency_Hz[theory_feature_idx]),
            "openems_S11_feature_frequency_Hz": float(numeric["frequency_Hz"][openems_feature_idx]),
            "feature_frequency_signed_mismatch_Hz": feature_mismatch_Hz,
            "feature_frequency_absolute_mismatch_Hz": abs(feature_mismatch_Hz),
            "feature_frequency_relative_mismatch": abs(feature_mismatch_Hz) / float(tx.frequency_Hz[theory_feature_idx]),
            "theory_minimum_S11_dB": float(db(s11_theory)[theory_feature_idx]),
            "openems_minimum_S11_dB": float(db(s11_openems)[openems_feature_idx]),
            "theory_S11_at_350MHz": [float(s11_theory[idx].real), float(s11_theory[idx].imag)],
            "openems_S11_at_350MHz": [float(s11_openems[idx].real), float(s11_openems[idx].imag)],
            "theory_S11_dB_at_350MHz": float(db(s11_theory)[idx]),
            "openems_S11_dB_at_350MHz": float(db(s11_openems)[idx]),
            "theory_Zin_ohm_at_350MHz": [float(zin_theory[idx].real), float(zin_theory[idx].imag)],
            "openems_Zin_ohm_at_350MHz": [float(zin_openems[idx].real), float(zin_openems[idx].imag)],
            "Zin_complex_difference_openems_minus_theory_ohm": [
                float((zin_openems[idx] - zin_theory[idx]).real),
                float((zin_openems[idx] - zin_theory[idx]).imag),
            ],
            "Zin_absolute_difference_ohm": float(abs(zin_openems[idx] - zin_theory[idx])),
            "S11_magnitude_correlation": corr11,
            "S11_dB_signed_difference_openems_minus_theory_at_350MHz": float(db(s11_openems)[idx] - db(s11_theory)[idx]),
            "S11_dB_absolute_difference_at_350MHz": float(abs(db(s11_openems)[idx] - db(s11_theory)[idx])),
            "theory_S21_at_350MHz": [float(s21_theory[idx].real), float(s21_theory[idx].imag)],
            "openems_S21_at_350MHz": [float(s21_openems[idx].real), float(s21_openems[idx].imag)],
            "theory_S21_dB_at_350MHz": float(db(s21_theory)[idx]),
            "openems_S21_dB_at_350MHz": float(db(s21_openems)[idx]),
            "S21_magnitude_correlation": corr21,
            "S21_dB_signed_difference_openems_minus_theory_at_350MHz": float(db(s21_openems)[idx] - db(s21_theory)[idx]),
            "S21_dB_absolute_difference_at_350MHz": float(abs(db(s21_openems)[idx] - db(s21_theory)[idx])),
            "principal_spectral_trend": "NUMERIC_MAGNITUDE_CORRELATION_REPORTED_WITHOUT_AMPLITUDE_SCALING",
            "phase_status": "THEORY_PHASE_MODEL_LIMITED_TO_FREE_SPACE_PROPAGATION",
            "wire_radius_status": "OPENEMS_THIN_WIRE_RADIUS_NOT_IDENTICAL_TO_THEORY",
        }
        table = np.column_stack((
            tx.frequency_Hz,
            s11_theory.real, s11_theory.imag, db(s11_theory),
            s11_openems.real, s11_openems.imag, db(s11_openems),
            zin_theory.real, zin_theory.imag, zin_openems.real, zin_openems.imag,
            s21_theory.real, s21_theory.imag, db(s21_theory),
            s21_openems.real, s21_openems.imag, db(s21_openems),
        ))
        np.savetxt(out / "h2_theory_openems_comparison.csv", table, delimiter=",", fmt="%.12e",
                   header=("frequency_Hz,theory_S11_real,theory_S11_imag,theory_S11_dB,"
                           "openems_S11_real,openems_S11_imag,openems_S11_dB,"
                           "theory_Zin_real_ohm,theory_Zin_imag_ohm,openems_Zin_real_ohm,openems_Zin_imag_ohm,"
                           "theory_S21_real,theory_S21_imag,theory_S21_dB,"
                           "openems_S21_real,openems_S21_imag,openems_S21_dB"), comments="")
        write_json(out / "h2_theory_surrogate_summary.json", {
            "data_provenance": metadata["status"],
            "theory_input_status": theory["status"],
            "VNA_gate_satisfied": False,
            "comparison": comparisons,
            "H2_DEVELOPMENT_GATE": "PASS",
            "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
            "H3_TOOL_DEVELOPMENT_ALLOWED": True,
            "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING": True,
            "production_receiver_status": "NOT_RESOLVED",
        })
    write_json(out / "h2_350mhz_theory_manifest.json", theory)

    shutil.copyfile(RAW / "sparameters.csv", out / "h2_350mhz_openems_sparameters.csv")
    shutil.copyfile(RAW / "nf2ff_cut.csv", out / "h2_350mhz_openems_nf2ff_cut.csv")
    geometry = {key: result[key] for key in (
        "case", "scientific_role", "dipole_total_length_m", "dipole_arm_length_m",
        "theory_wire_radius_m", "openems_wire_model", "wire_radius_status",
        "tx_rx_center_separation_m", "orientation", "environment", "Z0_ohm",
        "frequency_start_Hz", "frequency_stop_Hz", "frequency_spacing_Hz", "frequency_points",
    )}
    write_json(out / "h2_350mhz_geometry.json", geometry)

    stage4 = [2941408508.9091916, 7966314711.629051]
    stage5 = [3047273105.1868486, 10233758844.919172]
    correction = {
        "SYSTEM_TRANSMISSION_TARGET_CENTER_Hz": 350e6,
        "SYSTEM_TRANSMISSION_MAX_Hz": 500e6,
        "H3_SYSTEM_FULLWAVE_DEVELOPMENT_BAND": {
            "low_Hz": 200e6, "high_Hz": 500e6,
            "status": "PROJECT_TARGET_BAND_FOR_TOOL_DEVELOPMENT",
        },
        "legacy_reference": {
            "role": "H2_LEGACY_INFRASTRUCTURE_REFERENCE",
            "frequency_Hz": [0.7e9, 1.3e9], "dipole_length_m": 0.15,
            "runner_sha256": "7b0131bea5787cc1e1e793406982feb2b9e94081e278d1feddde16c9c4b40453",
            "raw_results_modified": False,
        },
        "Stage_F_native_bands": {"Stage4": stage4, "Stage5": stage5,
                                  "role": "H4_NATIVE_RF_PHYSICS"},
        "G3_spectrum_policy": "FROZEN_NOT_REINTERPRETED_BY_500_MHZ_SYSTEM_TARGET",
    }
    write_json(out / "h2_frequency_correction.json", correction)
    write_json(out / "h2_350mhz_future_vna_contract.json", {
        "status": "VNA_MEASUREMENT_PENDING",
        "required_files": ["RX_S11.s1p", "TX_S11.s1p", "TX_RX_S21.s2p"],
        "frequency_start_Hz": 200e6,
        "frequency_stop_Hz": 500e6,
        "frequency_points": 601,
        "frequency_spacing_Hz": 0.5e6,
        "Z0_ohm": 50.0,
        "parser": "streamer_rf.fullwave.receiver.read_touchstone",
        "comparison_API_unchanged": True,
        "replaceable_fields": ["data_values", "geometry_metadata", "provenance", "validation_status"],
        "calibration_metadata_required": True,
    })
    summary = {
        "H2_DEVELOPMENT_GATE": "MINIMAL_FIX_REQUIRED_MISSING_THEORY_FILES" if missing else "PASS",
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
        "H3_TOOL_DEVELOPMENT_ALLOWED": not missing,
        "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING": True,
        "production_receiver_status": "NOT_RESOLVED",
        "theory_input_status": theory["status"],
        "theory_openems_comparison": comparisons,
        "openems": result,
        "resource": {"runtime_s": result["runtime_s"], "peak_RSS_KiB": result["peak_RSS_KiB"],
                     "raw_bytes": result["raw_bytes"], "generation_runtime_s": time.perf_counter() - started,
                     "generation_peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "G3_waveform_injected": False,
    }
    write_json(out / "h2_350mhz_summary.json", summary)


if __name__ == "__main__":
    main()
