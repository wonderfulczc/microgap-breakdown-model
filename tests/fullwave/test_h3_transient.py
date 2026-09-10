import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.transient import *


def test_g3_waveform_reader_and_frozen_contract():
    time, voltage, current = read_g3_waveform(ROOT / "thermal/g3_port/g3_port_uniform.csv")
    assert len(time) == 801
    assert np.median(np.diff(time)) == pytest.approx(1.25e-11)
    assert (time[0], time[-1]) == pytest.approx((0.0, 1e-8))
    assert np.all(np.isfinite(voltage)) and np.all(np.isfinite(current))


def test_h3_band_is_exactly_200_to_500_mhz():
    frequency = np.array([199e6, 200e6, 350e6, 500e6, 501e6])
    assert np.array_equal(exact_band_mask(frequency), [False, True, True, True, False])
    assert H3_BAND_HZ == (200e6, 500e6)


def test_stage_f_native_trust_bands_are_unchanged():
    correction = json.loads((ROOT / "fullwave/h2/h2_frequency_correction.json").read_text())
    assert correction["Stage_F_native_bands"]["Stage4"] == [2941408508.9091916, 7966314711.629051]
    assert correction["Stage_F_native_bands"]["Stage5"] == [3047273105.1868486, 10233758844.919172]


def test_total_voltage_transfer_uses_total_tx_port_voltage():
    s11 = np.array([0.25 + 0.1j])
    s21 = np.array([0.5 - 0.2j])
    assert total_port_voltage_transfer(s11, s21)[0] == pytest.approx(s21[0] / (1 + s11[0]))


def test_complex_transfer_interpolation_and_no_extrapolation():
    source_f = np.array([200e6, 300e6, 500e6])
    source_h = 1 + 1j * source_f / 1e8
    target = np.array([250e6, 400e6])
    assert np.allclose(interpolate_transfer(source_f, source_h, target), 1 + 1j * target / 1e8)
    with pytest.raises(ValueError, match="EXTRAPOLATION"):
        interpolate_transfer(source_f, source_h, [199e6])


def test_analytic_lti_transfer_and_linear_scaling():
    source = np.array([1 + 2j, -0.5 + 1j])
    transfer = np.array([2 - 1j, 2 - 1j])
    received = apply_transfer(source, transfer)
    assert np.allclose(received, source * transfer)
    assert np.allclose(apply_transfer(3 * source, transfer), 3 * received)


def test_time_shift_has_expected_phase():
    frequency = np.array([200e6, 350e6, 500e6])
    delay = 3e-9
    shifted = np.exp(-2j * np.pi * frequency * delay)
    assert np.allclose(np.angle(shifted / shifted[0]),
                       np.angle(np.exp(-2j * np.pi * (frequency - frequency[0]) * delay)))


def test_zero_padding_preserves_original_dft_bins():
    time = np.arange(9) * 1e-9
    signal = np.sin(2 * np.pi * np.arange(9) / 9) + 0.2
    _, original, _, _ = zero_padded_rfft(time, signal, multiplier=1)
    _, padded, _, _ = zero_padded_rfft(time, signal, multiplier=4)
    assert np.allclose(padded[::4][: original.size], original)


def test_causal_linear_convolution_has_no_circular_wrap():
    source = np.array([1.0, 2.0])
    impulse = np.array([0.0, 0.5, 0.25])
    result = causal_linear_convolution(source, impulse)
    assert result.tolist() == pytest.approx([0.0, 0.5, 1.25, 0.5])
    assert len(result) == len(source) + len(impulse) - 1


def test_analytic_envelope_of_single_tone_is_constant():
    phase = 2 * np.pi * np.arange(101) / 101
    assert np.allclose(analytic_envelope(np.cos(phase)), 1.0, atol=1e-12)


def test_loading_current_and_status_thresholds():
    voltage = np.array([2 + 0j, 1 + 1j])
    admittance = np.array([0.5 + 0j, 0.25 - 0.25j])
    assert np.allclose(implied_loading_current(admittance, voltage), admittance * voltage)
    assert loading_status(0.25) == "FULL_WAVE_LOADING_MISMATCH_LOW"
    assert loading_status(0.5) == "FULL_WAVE_LOADING_MISMATCH_MODERATE"
    assert loading_status(0.8) == "FULL_WAVE_LOADING_MISMATCH_HIGH"


def test_low_current_impedance_is_invalid_not_infinite():
    impedance, valid = waveform_impedance([1 + 0j, 2 + 0j], [0 + 0j, 0 + 0j])
    assert not np.any(valid)
    assert np.all(np.isnan(impedance.real))
    impedance, valid = waveform_impedance([1 + 0j, 2 + 0j], [1 + 0j, 1e-12 + 0j], relative_gate=1e-6)
    assert np.array_equal(valid, [True, False])


def test_baseline_has_non_ryu_provenance_and_matched_peak():
    time = np.linspace(0, 1e-8, 101)
    reference = 100 * np.exp(-time / 1e-9)
    baseline, duration = matched_rise_smooth_step(time, reference)
    assert BASELINE_PROVENANCE == "NUMERICAL_SOURCE_SHAPE_BASELINE_ONLY"
    assert baseline.max() == pytest.approx(100)
    assert duration > 0


def test_generated_spectrum_identity_and_status_propagation():
    data = np.genfromtxt(ROOT / "fullwave/h3/h3_received_spectrum.csv", delimiter=",", names=True)
    source = data["V_G3_real_V"] + 1j * data["V_G3_imag_V"]
    transfer = data["H_V_real"] + 1j * data["H_V_imag"]
    received = data["V_rx_real_V"] + 1j * data["V_rx_imag_V"]
    assert np.allclose(received, source * transfer, rtol=2e-11, atol=1e-12)
    summary = json.loads((ROOT / "fullwave/h3/h3_summary.json").read_text())
    assert summary["H2_SCIENTIFIC_VALIDATION"] == "VNA_MEASUREMENT_PENDING"
    assert summary["full_wave_feedback_status"] == "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED"
    assert summary["openEMS_rerun"] is False


def test_generated_contract_is_development_only():
    contract = json.loads((ROOT / "fullwave/h3/h3_result_contract.json").read_text())
    assert contract["frequency_band_Hz"] == [200e6, 500e6]
    assert contract["result_voltage_semantics"] == "REFERENCE_RECEIVED_VOLTAGE"
    assert contract["receiver_status"] == "THEORY_AND_FULLWAVE_DEVELOPMENT_VERIFIED"
    assert contract["experimental_validation_status"] == "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING"


def test_invalid_g3_waveform_is_rejected(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("time_s,V_port_V,I_port_A\n0,1,1\n1,nan,1\n2,1,1\n3,1,1\n")
    with pytest.raises(ValueError, match="INVALID_G3"):
        read_g3_waveform(path)
