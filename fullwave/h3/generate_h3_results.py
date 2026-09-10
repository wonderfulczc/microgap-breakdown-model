"""Generate H3 by applying frozen H2 transfer data to the frozen G3 transient."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.foundation import C0  # noqa: E402
from streamer_rf.fullwave.transient import (  # noqa: E402
    BASELINE_PROVENANCE,
    H3_BAND_HZ,
    LOADING_THRESHOLDS,
    apply_transfer,
    analytic_envelope,
    causal_band_impulse,
    causal_linear_convolution,
    exact_band_mask,
    implied_loading_current,
    input_admittance,
    interpolate_transfer,
    loading_metrics,
    loading_status,
    matched_rise_smooth_step,
    read_g3_waveform,
    retained_spectral_fraction,
    total_port_voltage_transfer,
    waveform_impedance,
    zero_padded_rfft,
)


G3_PATH = ROOT / "thermal/g3_port/g3_port_uniform.csv"
G3_SUMMARY = ROOT / "thermal/g3_port/g3_port_summary.json"
H2_PATH = ROOT / "fullwave/h2/h2_350mhz_openems_sparameters.csv"
H2_SUMMARY = ROOT / "fullwave/h2/h2_350mhz_summary.json"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def save_csv(path, names, columns):
    np.savetxt(path, np.column_stack(columns), delimiter=",", fmt="%.12e", header=",".join(names), comments="")


def complex_pair(value):
    return [float(value.real), float(value.imag)]


def spectral_summary(frequency, spectrum):
    power = np.abs(spectrum) ** 2
    total = float(np.sum(power))
    return {
        "peak_frequency_Hz": float(frequency[np.argmax(np.abs(spectrum))]),
        "centroid_frequency_Hz": float(np.sum(frequency * power) / total),
        "integrated_spectral_metric": float(np.trapezoid(power, frequency)),
    }


def waveform_summary(time_s, waveform):
    peak_index = int(np.argmax(np.abs(waveform)))
    envelope = analytic_envelope(waveform)
    envelope_peak_index = int(np.argmax(envelope))
    threshold = 0.01 * float(np.max(np.abs(waveform)))
    onset = np.flatnonzero(np.abs(waveform) >= threshold)
    return {
        "peak_V": float(np.max(np.abs(waveform))),
        "RMS_V": float(np.sqrt(np.mean(waveform**2))),
        "absolute_peak_time_s": float(time_s[peak_index]),
        "principal_envelope_peak_V": float(envelope[envelope_peak_index]),
        "principal_envelope_peak_time_s": float(time_s[envelope_peak_index]),
        "one_percent_onset_time_s": float(time_s[onset[0]]),
    }


def main():
    started = time.perf_counter()
    output = Path(__file__).parent
    g3_summary = json.loads(G3_SUMMARY.read_text())
    contract = g3_summary["contract"]
    expected = {
        "sample_count": 801,
        "dt_port_s": 1.25e-11,
        "reference_plane_id": "EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP",
        "voltage_reference": "GAP_NODE_MINUS_GROUND",
        "positive_current_direction": "EXTERNAL_TO_GAP",
    }
    for key, value in expected.items():
        if contract.get(key) != value:
            raise ValueError(f"G3_CONTRACT_MISMATCH_{key}")
    if contract.get("valid_time_interval_s") != [0.0, 1e-08]:
        raise ValueError("G3_TIME_INTERVAL_MISMATCH")
    time_s, voltage, current = read_g3_waveform(G3_PATH)

    h2 = np.genfromtxt(H2_PATH, delimiter=",", names=True)
    f_h2 = h2["frequency_Hz"]
    if len(f_h2) != 601 or f_h2[0] != 200e6 or f_h2[-1] != 500e6:
        raise ValueError("H2_TRANSFER_SUPPORT_MISMATCH")
    s11 = h2["S11_real"] + 1j * h2["S11_imag"]
    s21 = h2["S21_real"] + 1j * h2["S21_imag"]
    zin = h2["Zin_real_ohm"] + 1j * h2["Zin_imag_ohm"]
    h_voltage = total_port_voltage_transfer(s11, s21)
    yin = input_admittance(zin)
    save_csv(output / "h3_transfer_function.csv",
             ["frequency_Hz", "H_V_real", "H_V_imag", "H_V_magnitude", "H_V_phase_rad",
              "Zin_real_ohm", "Zin_imag_ohm", "Yin_real_S", "Yin_imag_S"],
             [f_h2, h_voltage.real, h_voltage.imag, abs(h_voltage), np.unwrap(np.angle(h_voltage)),
              zin.real, zin.imag, yin.real, yin.imag])

    f_fft, v_spectrum, dt, nfft = zero_padded_rfft(time_s, voltage, multiplier=16)
    _, i_spectrum, _, _ = zero_padded_rfft(time_s, current, multiplier=16)
    mask = exact_band_mask(f_fft)
    f_band = f_fft[mask]
    h_band = interpolate_transfer(f_h2, h_voltage, f_band)
    z_band = interpolate_transfer(f_h2, zin, f_band)
    y_band = input_admittance(z_band)
    v_band = v_spectrum[mask]
    i_g3_band = i_spectrum[mask]
    v_rx_band = apply_transfer(v_band, h_band)
    i_implied_band = implied_loading_current(y_band, v_band)
    z_waveform, z_valid = waveform_impedance(v_band, i_g3_band)
    retained = retained_spectral_fraction(v_spectrum, mask)

    save_csv(output / "h3_g3_source_spectrum.csv",
             ["frequency_Hz", "V_G3_real_V", "V_G3_imag_V", "V_G3_magnitude_V",
              "I_G3_real_A", "I_G3_imag_A", "I_G3_magnitude_A"],
             [f_band, v_band.real, v_band.imag, abs(v_band),
              i_g3_band.real, i_g3_band.imag, abs(i_g3_band)])
    save_csv(output / "h3_received_spectrum.csv",
             ["frequency_Hz", "V_G3_real_V", "V_G3_imag_V", "H_V_real", "H_V_imag",
              "V_rx_real_V", "V_rx_imag_V", "V_G3_magnitude_V", "H_V_magnitude",
              "V_rx_magnitude_V", "phase_source_rad", "phase_transfer_rad", "phase_received_rad"],
             [f_band, v_band.real, v_band.imag, h_band.real, h_band.imag,
              v_rx_band.real, v_rx_band.imag, abs(v_band), abs(h_band), abs(v_rx_band),
              np.angle(v_band), np.angle(h_band), np.angle(v_rx_band)])
    save_csv(output / "h3_loading_diagnostic.csv",
             ["frequency_Hz", "I_H_implied_real_A", "I_H_implied_imag_A",
              "I_G3_real_A", "I_G3_imag_A", "Z_waveform_G3_real_ohm",
              "Z_waveform_G3_imag_ohm", "Z_waveform_G3_valid"],
             [f_band, i_implied_band.real, i_implied_band.imag, i_g3_band.real, i_g3_band.imag,
              z_waveform.real, z_waveform.imag, z_valid.astype(float)])

    full_h = np.zeros_like(f_fft, dtype=complex)
    full_h[mask] = h_band
    free_space_delay = 1.0 / C0
    impulse = causal_band_impulse(f_fft, full_h, nfft, dt, free_space_delay)
    received_time = causal_linear_convolution(voltage, impulse)
    received_t = np.arange(received_time.size) * dt + time_s[0]

    baseline, baseline_rise = matched_rise_smooth_step(time_s, voltage)
    baseline_received = causal_linear_convolution(baseline, impulse)
    baseline_t = np.arange(baseline_received.size) * dt + time_s[0]
    if not np.array_equal(received_t, baseline_t):
        raise ValueError("BASELINE_OUTPUT_GRID_MISMATCH")
    source_g3_padded = np.zeros(received_t.shape)
    source_baseline_padded = np.zeros(received_t.shape)
    source_g3_padded[: voltage.size] = voltage
    source_baseline_padded[: baseline.size] = baseline
    save_csv(output / "h3_received_timeseries.csv",
             ["time_s", "V_rx_reference_V"], [received_t, received_time])
    save_csv(output / "h3_source_baseline_comparison.csv",
             ["time_s", "V_source_G3_full_V", "V_source_baseline_V",
              "V_rx_G3_V", "V_rx_baseline_V"],
             [received_t, source_g3_padded, source_baseline_padded,
              received_time, baseline_received])

    _, baseline_spectrum, _, _ = zero_padded_rfft(time_s, baseline, multiplier=16)
    baseline_band = baseline_spectrum[mask]
    baseline_rx_band = apply_transfer(baseline_band, h_band)
    g3_wave = waveform_summary(received_t, received_time)
    base_wave = waveform_summary(received_t, baseline_received)
    correlation = float(np.corrcoef(received_time, baseline_received)[0, 1])
    source_spectral = spectral_summary(f_band, v_band)
    received_spectral = spectral_summary(f_band, v_rx_band)
    baseline_spectral = spectral_summary(f_band, baseline_rx_band)
    loading = loading_metrics(i_implied_band, i_g3_band)
    load_status = loading_status(loading["normalized_L2_current_mismatch"])

    nearest350 = int(np.argmin(abs(f_band - 350e6)))
    nearest_h2_350 = int(np.argmin(abs(f_h2 - 350e6)))
    nf2ff = json.loads(H2_SUMMARY.read_text())["openems"]["nf2ff"]
    baseline_summary = {
        "provenance": BASELINE_PROVENANCE,
        "Ryu_equation_available": False,
        "construction": "CUBIC_SMOOTH_STEP_WITH_G3_90_TO_10_PERCENT_TRANSITION_DURATION",
        "peak_voltage_V": float(np.max(abs(baseline))),
        "derived_rise_duration_s": baseline_rise,
        "received": base_wave,
        "received_spectrum": baseline_spectral,
        "received_350MHz_magnitude_V": float(abs(baseline_rx_band[nearest350])),
    }
    summary = {
        "method_status": "PHYSICS_DERIVED_G3_TO_FULLWAVE_RECEIVER_PATH_VALIDATED",
        "source_case_id": contract["source_case_id"],
        "reference_structure_id": "H3_350MHZ_REFERENCE_RADIATING_STRUCTURE",
        "PRODUCTION_TX_STRUCTURE": "PENDING_ACTUAL_GEOMETRY",
        "receiver_id": "H2_350MHZ_DEVELOPMENT_REFERENCE_RX",
        "receiver_status": "THEORY_AND_FULLWAVE_DEVELOPMENT_VERIFIED",
        "H2_SCIENTIFIC_VALIDATION": "VNA_MEASUREMENT_PENDING",
        "experimental_validation_status": "STAGE_H_EXPERIMENTAL_VALIDATION_PENDING",
        "full_wave_feedback_status": "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED",
        "openEMS_rerun": False,
        "transfer_definition": "H_V=S21/(1+S11)=V_rx_50ohm/V_tx_total_port",
        "transfer_source": "FROZEN_H2_50_OHM_S_PARAMETERS",
        "time_transform": {
            "normalization": "COMPLEX_SINGLE_SIDED_PEAK_AMPLITUDE_2_OVER_N",
            "original_samples": int(len(time_s)), "original_dt_s": dt,
            "original_duration_s": float(time_s[-1] - time_s[0]),
            "intrinsic_Fourier_scale_Hz": float(1.0 / (time_s[-1] - time_s[0])),
            "nfft": nfft, "zero_padding_status": "ZERO_PADDING_FOR_LINEAR_CONVOLUTION_ONLY",
            "frequency_spacing_Hz": float(f_fft[1] - f_fft[0]),
        },
        "band": {"low_Hz": H3_BAND_HZ[0], "high_Hz": H3_BAND_HZ[1],
                 "frequency_samples": int(np.sum(mask)), "retained_source_spectral_fraction": retained},
        "source_representations": {
            "FULL_G3_WAVEFORM": "thermal/g3_port/g3_port_uniform.csv",
            "H3_BANDLIMITED_G3_WAVEFORM": "h3_g3_source_spectrum.csv",
            "frozen_G3_modified": False,
        },
        "source_spectrum": source_spectral,
        "received_spectrum": received_spectral,
        "received_waveform": {**g3_wave, "dt_output_s": dt,
                              "sample_count": int(received_t.size),
                              "duration_s": float(received_t[-1] - received_t[0]),
                              "causal_impulse_samples": int(impulse.size)},
        "propagation_delay": {"R_m": 1.0, "R_over_c_s": free_space_delay,
                              "onset_minus_R_over_c_s": g3_wave["one_percent_onset_time_s"] - free_space_delay},
        "at_nearest_350MHz": {
            "frequency_Hz": float(f_band[nearest350]),
            "V_G3_magnitude_V": float(abs(v_band[nearest350])),
            "H_V_magnitude": float(abs(h_band[nearest350])),
            "H_V_phase_rad": float(np.angle(h_band[nearest350])),
            "V_rx_magnitude_V": float(abs(v_rx_band[nearest350])),
            "Zin_ohm": complex_pair(z_band[nearest350]),
            "Yin_S": complex_pair(y_band[nearest350]),
            "H2_grid_350MHz_H_V": complex_pair(h_voltage[nearest_h2_350]),
        },
        "loading": {**loading, "thresholds": {"low_max": LOADING_THRESHOLDS[0],
                                               "moderate_max": LOADING_THRESHOLDS[1]},
                    "status": load_status, "feedback_applied": False},
        "baseline": baseline_summary,
        "G3_vs_baseline": {
            "received_peak_ratio_baseline_over_G3": base_wave["peak_V"] / g3_wave["peak_V"],
            "received_RMS_ratio_baseline_over_G3": base_wave["RMS_V"] / g3_wave["RMS_V"],
            "received_spectral_centroid_difference_Hz": (
                baseline_spectral["centroid_frequency_Hz"] - received_spectral["centroid_frequency_Hz"]),
            "received_350MHz_magnitude_ratio_baseline_over_G3": (
                abs(baseline_rx_band[nearest350]) / abs(v_rx_band[nearest350])),
            "time_domain_waveform_correlation": correlation,
        },
        "nf2ff_reused": {**nf2ff, "source": "FROZEN_H2_350MHZ_OPENEMS_RESULT",
                         "principal_polarization": "Z_DIRECTED_DIPOLE_BROADSIDE"},
        "energy_interpretation": "ONE_WAY_VOLTAGE_SOURCE_HANDOFF_NOT_G3_PORT_ENERGY_PRESERVATION",
        "runtime_s": time.perf_counter() - started,
        "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    write_json(output / "h3_summary.json", summary)
    result_contract = {
        "source_case_id": contract["source_case_id"],
        "G3_input_hash": sha256(G3_PATH),
        "H2_transfer_hash": sha256(H2_PATH),
        "reference_structure_id": summary["reference_structure_id"],
        "receiver_id": summary["receiver_id"],
        "frequency_band_Hz": list(H3_BAND_HZ),
        "V_source_file": "h3_g3_source_spectrum.csv",
        "V_received_time_file": "h3_received_timeseries.csv",
        "V_received_spectrum_file": "h3_received_spectrum.csv",
        "H_V_file": "h3_transfer_function.csv",
        "Z_in_file": "h3_transfer_function.csv",
        "loading_status": load_status,
        "receiver_status": summary["receiver_status"],
        "experimental_validation_status": summary["experimental_validation_status"],
        "full_wave_feedback_status": summary["full_wave_feedback_status"],
        "result_voltage_semantics": "REFERENCE_RECEIVED_VOLTAGE",
    }
    write_json(output / "h3_result_contract.json", result_contract)


if __name__ == "__main__":
    main()
