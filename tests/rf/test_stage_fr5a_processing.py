from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import c, epsilon_0

from streamer_rf.rf.ultrafast import (
    analyze_pulses,
    central_difference_second_order,
    compute_esd_spectrum,
    find_relative_db_crossing,
    fit_spectral_slope,
    process_current_moment,
    resample_uniform,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "rf" / "f_r5a"


def load(name: str):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tanh_signal(tau: float = 10e-12):
    dt = 0.1e-12
    time = np.arange(-250e-12, 250e-12 + dt / 2, dt)
    moment = 1e-8 * (1.0 + np.tanh(time / tau))
    return time, moment


def test_uniform_resampling_and_metadata():
    time = np.array([0.0, 1.0, 2.2, 3.0, 4.0]) * 1e-12
    result = resample_uniform(time, time**2, strategy="USER_SPECIFIED_DT", user_dt_s=0.5e-12)
    assert np.allclose(np.diff(result.time_s), 0.5e-12)
    assert result.strategy == "USER_SPECIFIED_DT"
    assert result.resolution_status == "INTERPOLATED_NOT_PHYSICALLY_RESOLVED"
    assert result.upsampling_factor == pytest.approx(2.0)


def test_resampling_strategies_are_explicit():
    time = np.array([0.0, 1.0, 2.2, 3.2, 4.2]) * 1e-12
    values = np.arange(5.0)
    median = resample_uniform(time, values)
    reference = resample_uniform(time, values, strategy="REFERENCE_PROFILE_DT", reference_dt_s=1e-12)
    assert median.uniform_dt_s == pytest.approx(np.median(np.diff(time)))
    assert reference.uniform_dt_s == pytest.approx(1e-12)


@pytest.mark.parametrize("time", [np.array([0, 1, 1, 2, 3.0]), np.array([0, 2, 1, 3, 4.0])])
def test_time_must_be_strictly_monotonic(time):
    with pytest.raises(ValueError, match="strictly monotonically increasing"):
        resample_uniform(time, np.arange(5.0))


def test_central_derivative_second_order_accuracy_and_boundary_policy():
    time = np.linspace(-2.0, 2.0, 101)
    derivative_time, derivative = central_difference_second_order(time, time**3)
    assert derivative_time.size == time.size - 2
    assert np.max(np.abs(derivative - 3 * derivative_time**2)) < 0.002


def test_koile_moving_average_duration_metadata():
    time, moment = tanh_signal()
    result = process_current_moment(time, moment, profile="KOILE_COMPATIBLE")
    assert result.smoothing_metadata["smoothing_samples"] == 10
    assert result.smoothing_metadata["smoothing_duration_s"] == pytest.approx(1e-12)
    assert result.smoothing_metadata["window_length_source"] == "PAPER_REPORTED"


def test_tukey_alpha_is_project_choice():
    time, moment = tanh_signal()
    result = process_current_moment(time, moment, profile="NO_SMOOTHING", tukey_alpha=0.75)
    assert result.spectrum.tukey_alpha == 0.75
    assert result.processing_provenance["tukey_alpha"] == "PROJECT_OPERATIONAL_CHOICE"


def test_zero_padding_semantics_and_physical_resolution_are_distinct():
    time, moment = tanh_signal()
    result = process_current_moment(time, moment, profile="NO_SMOOTHING")
    spectrum = result.spectrum
    assert spectrum.zero_padded_duration_s >= 12e-9
    assert spectrum.display_fft_bin_spacing_Hz < spectrum.physical_frequency_resolution_Hz
    assert spectrum.physical_frequency_resolution_Hz == pytest.approx(1 / spectrum.physical_record_duration_s)


def test_fft_dt_normalization_and_parseval():
    time = np.arange(0.0, 100e-12, 0.1e-12)
    signal = np.sin(2 * np.pi * 10e9 * time)
    spectrum, _ = compute_esd_spectrum(time, signal, tukey_alpha=0.0, zero_pad_duration_s=100e-12)
    expected = (time[1] - time[0]) * np.fft.rfft(signal)
    assert np.allclose(spectrum.transform_A_m, expected)
    assert spectrum.parseval_relative_error < 1e-12


def test_esd_is_finite_positive_and_uses_reported_formula():
    time = np.arange(0.0, 100e-12, 0.1e-12)
    signal = np.sin(2 * np.pi * 10e9 * time)
    spectrum, _ = compute_esd_spectrum(time, signal, tukey_alpha=0.0, zero_pad_duration_s=100e-12)
    expected = np.abs(spectrum.transform_A_m) ** 2 / (6 * np.pi * epsilon_0 * c**3)
    assert np.all(np.isfinite(spectrum.ESD))
    assert np.all(spectrum.ESD >= 0)
    assert np.allclose(spectrum.ESD, expected)


def test_fwhm_and_fw1e_match_sech_squared_pulse():
    tau = 10e-12
    time, moment = tanh_signal(tau)
    result = process_current_moment(time, moment, profile="NO_SMOOTHING")
    assert result.pulse.PW_FWHM_s == pytest.approx(2 * tau * math.acosh(math.sqrt(2)), rel=2e-4)
    assert result.pulse.PW_FW1E_s == pytest.approx(2 * tau * math.acosh(math.sqrt(math.e)), rel=2e-4)
    assert result.pulse.interpretation == "AMBIGUOUS_PENDING_WAVEFORM_REPRODUCTION"
    assert result.temporal_resolution_metadata["status"] == "RESOLVED_AT_10_SAMPLE_TARGET"
    assert result.temporal_resolution_metadata["samples_per_FWHM"] > 100


def test_multi_peak_requires_explicit_primary_event():
    time = np.linspace(0, 1, 1001)
    signal = np.exp(-((time - 0.3) / 0.02) ** 2) + 0.8 * np.exp(-((time - 0.7) / 0.02) ** 2)
    result = analyze_pulses(time, signal)
    assert result.status == "NOT_RESOLVED_MULTI_PEAK"
    assert len(result.peaks) == 2
    assert result.PW_FWHM_s is None


def test_one_ghz_reference_and_db_crossings():
    frequency = np.logspace(8, 11, 2001)
    esd = (frequency / 1e9) ** -2
    f3 = find_relative_db_crossing(frequency, esd, level_dB=-3)
    f10 = find_relative_db_crossing(frequency, esd, level_dB=-10)
    assert f3.reference_frequency_Hz == 1e9
    assert f3.status == "PASS"
    assert f10.status == "PASS"
    assert f3.validated_crossing_Hz == pytest.approx(1e9 * 10 ** (3 / 20), rel=5e-4)
    assert f10.validated_crossing_Hz == pytest.approx(1e9 * 10 ** 0.5, rel=5e-4)


def test_crossing_out_of_band_and_ambiguous():
    frequency = np.linspace(0.1e9, 2e9, 100)
    assert find_relative_db_crossing(frequency, np.ones(100), level_dB=-3).status == "OUT_OF_BAND"
    oscillatory = 1.0 + 0.9 * np.sin(np.linspace(0, 14 * np.pi, 100))
    assert find_relative_db_crossing(frequency, oscillatory, level_dB=-3).status in {"AMBIGUOUS_CROSSING", "OUT_OF_BAND"}


def test_spectral_slope_uses_explicit_band():
    frequency = np.logspace(8, 11, 1000)
    esd = (frequency / 1e9) ** -2
    result = fit_spectral_slope(frequency, esd, (1e9, 10e9))
    assert result["fit_band_Hz"] == [1e9, 10e9]
    assert result["slope_dB_per_decade"] == pytest.approx(-20.0, rel=1e-10)
    assert result["R2"] == pytest.approx(1.0)


def test_koile_table_constants_and_group_trends():
    reference = load("benchmark_reference.json")
    assert reference["group1"][0] == {"distance_cm": 1.2, "Eamb_over_Ek": 0.65, "PW_ps": 21.2, "f3dB_GHz": 2.98, "f10dB_GHz": 12.6}
    assert reference["group2"][-1]["PW_ps"] == 9.6
    assert reference["group2"][-1]["f10dB_GHz"] == 30.0
    regression = load("koile_table_regression.json")
    assert regression["status"] == "PASS"
    assert regression["group1"]["status"] == "PASS"
    assert regression["group2"]["status"] == "PASS"


def test_f_r5a_status_preserves_scientific_boundary():
    status = load("f_r5a_status.json")
    assert status["F_R5A_STATUS"] == "PASS"
    assert status["KOILE_WAVEFORM_REPRODUCTION"] == "PENDING_EXTERNAL_DATA"
    assert status["PULSE_WIDTH_INTERPRETATION"] == "AMBIGUOUS"
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert status["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"
    assert status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert status["LARGE_SIMULATION_RERUN"] is False


def test_stage_f_dry_run_is_read_only_and_does_not_upgrade_trust():
    dry_run = load("stage_f_dry_run.json")
    assert dry_run["mode"] == "READ_ONLY_POSTPROCESSING_DRY_RUN"
    assert dry_run["consumption_status"] == "PASS"
    assert dry_run["trust_decision"] == "NO_TRUST_UPGRADE_STAGE_F_MASK_HAS_PRIORITY"
    assert dry_run["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"


def test_r0_frozen_baseline_hashes_are_unchanged():
    r0 = json.loads((ROOT / "rf" / "ultrafast" / "v2_1_r0_status.json").read_text(encoding="utf-8"))
    for relative, expected in r0["frozen_baseline_hashes"].items():
        assert sha256(ROOT / relative) == expected


def test_a_to_j_stages_remain_unchanged():
    matrix = json.loads((ROOT / "rf" / "ultrafast" / "v2_1_stage_status_matrix.json").read_text(encoding="utf-8"))
    assert list(matrix["stages"]) == list("ABCDEFGHIJ")
    assert "K" not in matrix["stages"]


def test_processing_contract_parameter_provenance():
    contract = load("processing_contract.json")
    assert contract["derivative"]["source"] == "PROJECT_OPERATIONAL_CHOICE"
    assert contract["window"]["type_source"] == "PAPER_REPORTED"
    assert contract["window"]["alpha_source"] == "PROJECT_BASELINE_CANDIDATE"
    assert contract["zero_padding"]["semantic"] == "DISPLAY_GRID_ONLY_NOT_PHYSICAL_RESOLUTION"
    assert contract["trust"]["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
