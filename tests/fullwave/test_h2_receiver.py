import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.receiver import *


def geometry():
    return {
        "geometry_kind": "H2_CANONICAL_REFERENCE_RX",
        "hardware_status": "GEOMETRY_AVAILABLE_ONLY",
        "microgap_geometry_included": False,
        "total_length_m": 0.15,
        "tx_rx_distance_m": 0.3,
        "Z0_ohm": 50.0,
        "port_reference": PORT_REFERENCE,
        "voltage_reference": "POSITIVE_Z_ARM_MINUS_NEGATIVE_Z_ARM",
        "positive_current_direction": "PORT_INTO_RECEIVER",
    }


def test_receiver_geometry_contract_and_reference_plane():
    g = geometry()
    assert validate_receiver_geometry(g)
    assert len(geometry_hash(g)) == 64
    g["port_reference"] = "UNKNOWN"
    with pytest.raises(ValueError, match="REFERENCE_PLANE"):
        validate_receiver_geometry(g)


def test_no_stage_h_microgap_refinement():
    g = geometry()
    g["microgap_geometry_included"] = True
    with pytest.raises(ValueError, match="MICROGAP"):
        validate_receiver_geometry(g)


def test_loaded_and_open_circuit_voltage_are_distinct():
    assert loaded_to_open_circuit(1.0, 50.0, 50.0) == pytest.approx(2.0)
    assert loaded_to_open_circuit(1.0, 100.0, 50.0) == pytest.approx(3.0)
    with pytest.raises(ValueError):
        loaded_to_open_circuit(1.0, 50.0, 0.0)


def test_complex_band_mask():
    f = np.array([1.0, 2.0, 3.0])
    assert np.array_equal(band_mask(f, 1.5, 2.5), [False, True, False])
    with pytest.raises(ValueError):
        band_mask(f[::-1], 1.5, 2.5)


def test_near_far_field_classification():
    far = near_far_classification(0.15, 0.3, 1e9)
    assert far["status"] == "FAR_FIELD"
    assert far["far_field_limit_2D2_over_lambda_m"] == pytest.approx(0.150103842)
    assert near_far_classification(0.15, 0.01, 1e9)["status"] == "REACTIVE_NEAR_FIELD"


def test_touchstone_ri_s1p(tmp_path):
    path = tmp_path / "rx.s1p"
    path.write_text("# GHz S RI R 50\n1 0.1 -0.2\n2 0.2 -0.3\n")
    data = read_touchstone(path)
    assert data.Z0_ohm == 50
    assert data.frequency_Hz.tolist() == [1e9, 2e9]
    assert data.parameters["S11"][0] == pytest.approx(0.1 - 0.2j)


@pytest.mark.parametrize("fmt,row,expected", [
    ("MA", "1 0.5 90 0.25 180 0.25 0 0.5 -90", 0.5j),
    ("DB", "1 -6.020599913 0 -12.041199826 0 -12.041199826 0 -6.020599913 0", 0.5 + 0j),
])
def test_touchstone_s2p_complex_formats(tmp_path, fmt, row, expected):
    path = tmp_path / "pair.s2p"
    path.write_text(f"# GHz S {fmt} R 50\n{row}\n2 {row.split(maxsplit=1)[1]}\n")
    data = read_touchstone(path)
    assert data.parameters["S11"][0] == pytest.approx(expected)
    assert set(data.parameters) == {"S11", "S21", "S12", "S22"}


@pytest.mark.parametrize("text", [
    "1 0 0\n", "# GHz S RI R 50\n2 0 0\n1 0 0\n",
    "# GHz S RI R 50\n1 nan 0\n", "# GHz Y RI R 50\n1 0 0\n",
])
def test_invalid_touchstone(tmp_path, text):
    path = tmp_path / "bad.s1p"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_touchstone(path)


def test_interpolation_is_complex_and_bounded():
    result = bounded_complex_interpolate([1, 2], [1 + 2j, 3 + 4j], [1.5])
    assert result[0] == pytest.approx(2 + 3j)
    with pytest.raises(ValueError, match="EXTRAPOLATION"):
        bounded_complex_interpolate([1, 2], [1 + 2j, 3 + 4j], [0.5])


def test_spectral_metrics():
    value = spectral_metrics([1 + 1j, 2 + 0j, 4 - 1j], [1 + 1j, 2 + 0j, 4 - 1j])
    assert value["normalized_magnitude_L2"] == pytest.approx(0)
    assert value["magnitude_correlation"] == pytest.approx(1)


def test_g3_spectrum_and_frozen_band_provenance():
    data = np.genfromtxt(ROOT / "thermal/g3_port/g3_port_uniform.csv", delimiter=",", names=True)
    plan = g3_spectrum_plan(data["time_s"], data["V_port_V"])
    assert plan["dt_s"] == pytest.approx(1.25e-11)
    assert plan["upper_frequency_at_cumulative_fraction_Hz"] <= 5e8
    for filename in (
        "F4-P-stage4-left-isolated_rf_trust_report.json",
        "F4-C-stage5-highfield-collision_rf_trust_report.json",
    ):
        report = json.loads((ROOT / "rf/production/f4_attribution" / filename).read_text())
        assert report["trusted_frequency_low_Hz"] < report["trusted_frequency_high_Hz"]


def test_missing_vna_contract_stays_simulation_only():
    contract = {
        "receiver_id": "H2_CANONICAL_DIPOLE_RX",
        "hardware_status": "PLANNED_NOT_AVAILABLE",
        "geometry_hash": "a" * 64,
        "port_reference": PORT_REFERENCE,
        "Z0_ohm": 50.0,
        "measurement_id": None,
        "VNA_validation_status": "SIMULATION_ONLY",
        "H_rx_E_status": "RX_INTRINSIC_FIELD_TRANSFER_PENDING",
    }
    assert validate_transfer_contract(contract)
    contract["VNA_validation_status"] = "VALIDATED_WITH_VNA"
    with pytest.raises(ValueError, match="MEASUREMENT"):
        validate_transfer_contract(contract)


def test_receiver_validity_mask_combines_finite_band_and_source():
    frequency = np.array([0.5e9, 0.8e9, 1.0e9, 1.4e9])
    values = np.array([1, 1, np.nan, 1], dtype=complex)
    mask = band_mask(frequency, 0.7e9, 1.3e9) & np.isfinite(values)
    assert np.array_equal(mask, [False, True, False, False])


def test_generated_geometry_and_transfer_contract():
    g = json.loads((ROOT / "fullwave/h2/h2_geometry.json").read_text())
    stored_hash = g.pop("geometry_hash")
    assert validate_receiver_geometry(g)
    assert geometry_hash(g) == stored_hash
    contract = json.loads((ROOT / "fullwave/h2/h2_receiver_transfer_contract.json").read_text())
    assert validate_transfer_contract(contract)
    assert contract["geometry_hash"] == stored_hash


def test_generated_complex_sparameters_and_valid_band():
    data = np.genfromtxt(ROOT / "fullwave/h2/h2_reference_sparameters.csv", delimiter=",", names=True)
    s11 = data["S11_real"] + 1j * data["S11_imag"]
    s21 = data["S21_real"] + 1j * data["S21_imag"]
    mask = band_mask(data["frequency_Hz"], 0.7e9, 1.3e9) & np.isfinite(s11) & np.isfinite(s21)
    assert np.all(mask)
    assert np.allclose(data["S21_dB"], 20 * np.log10(np.abs(s21)))


def test_generated_missing_vna_gate_and_separate_band_contracts():
    summary = json.loads((ROOT / "fullwave/h2/h2_summary.json").read_text())
    bands = json.loads((ROOT / "fullwave/h2/h2_frequency_bands.json").read_text())
    assert summary["formal_node_status"] == "MINIMAL_FIX_REQUIRED"
    assert summary["blocker"] == "VNA_MEASUREMENT_REQUIRED"
    assert bands["H2_RECEIVER_VALIDATED_BANDS"] == []
    assert bands["H3_CIRCUIT_STRUCTURE_BANDS"] != bands["H4_NATIVE_RF_TRUSTED_BANDS"]


def test_350mhz_theory_metadata_cannot_close_vna_gate():
    metadata = json.loads((ROOT / "fullwave/h2/theory_inputs/H2_THEORY_350MHz_metadata.json").read_text())
    assert metadata["status"] == "THEORETICAL_SURROGATE_ONLY"
    assert receiver_validation_status(metadata["status"]) == "SIMULATION_ONLY"


def test_supplied_350mhz_touchstone_files_have_exact_grid_and_ordering():
    inputs = ROOT / "fullwave/h2/theory_inputs"
    rx = read_touchstone(inputs / "H2_THEORY_350MHz_RX_S11.s1p")
    tx = read_touchstone(inputs / "H2_THEORY_350MHz_TX_S11.s1p")
    pair = read_touchstone(inputs / "H2_THEORY_350MHz_TX_RX_S21.s2p")
    for data in (rx, tx, pair):
        assert validate_exact_frequency_grid(data.frequency_Hz, 200e6, 500e6, 601)
        assert data.Z0_ohm == 50.0
        assert all(np.all(np.isfinite(values)) for values in data.parameters.values())
    assert set(pair.parameters) == {"S11", "S21", "S12", "S22"}
    assert np.allclose(pair.parameters["S11"], tx.parameters["S11"], rtol=0, atol=1e-14)
    assert np.allclose(pair.parameters["S22"], rx.parameters["S11"], rtol=0, atol=1e-14)


def test_theory_comparison_closes_development_gate_only():
    summary = json.loads((ROOT / "fullwave/h2/h2_theory_surrogate_summary.json").read_text())
    assert summary["data_provenance"] == "THEORETICAL_SURROGATE_ONLY"
    assert summary["VNA_gate_satisfied"] is False
    assert summary["H2_DEVELOPMENT_GATE"] == "PASS"
    assert summary["H2_SCIENTIFIC_VALIDATION"] == "VNA_MEASUREMENT_PENDING"
    assert summary["H3_TOOL_DEVELOPMENT_ALLOWED"] is True
    assert summary["STAGE_H_EXPERIMENTAL_VALIDATION_PENDING"] is True
    assert summary["production_receiver_status"] == "NOT_RESOLVED"
    assert summary["comparison"]["wire_radius_status"] == "OPENEMS_THIN_WIRE_RADIUS_NOT_IDENTICAL_TO_THEORY"


def test_exact_350mhz_touchstone_frequency_grid(tmp_path):
    frequency = np.linspace(200e6, 500e6, 601)
    path = tmp_path / "future_rx.s1p"
    path.write_text("# MHz S RI R 50\n" + "\n".join(
        f"{f/1e6:.1f} 0 0" for f in frequency
    ) + "\n")
    data = read_touchstone(path)
    assert validate_exact_frequency_grid(data.frequency_Hz, 200e6, 500e6, 601)
    assert data.Z0_ohm == 50


def test_future_vna_replacement_uses_same_parser_and_contract(tmp_path):
    path = tmp_path / "replacement.s1p"
    path.write_text("# MHz S RI R 50\n200 0 0\n500 0 0\n")
    assert read_touchstone(path).source_format == "RI"
    metadata = {
        "instrument_model": "TRACEABLE_VNA",
        "calibration_type": "SOLT",
        "calibration_reference_plane": "CABLE_ENDS",
        "Z0_ohm": 50.0,
    }
    assert receiver_validation_status("VNA_MEASUREMENT", metadata) == "VALIDATED_WITH_VNA"


def test_frequency_correction_preserves_legacy_and_frozen_native_bands():
    correction = json.loads((ROOT / "fullwave/h2/h2_frequency_correction.json").read_text())
    assert correction["legacy_reference"]["runner_sha256"] == "7b0131bea5787cc1e1e793406982feb2b9e94081e278d1feddde16c9c4b40453"
    assert correction["H3_SYSTEM_FULLWAVE_DEVELOPMENT_BAND"]["high_Hz"] == 500e6
    assert correction["Stage_F_native_bands"]["Stage4"] == [2941408508.9091916, 7966314711.629051]
    assert correction["Stage_F_native_bands"]["Stage5"] == [3047273105.1868486, 10233758844.919172]
