import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.native_receiver import (  # noqa: E402
    NativeRFSourceContract,
    STAGE4_TRUSTED_BAND_HZ,
    STAGE5_TRUSTED_BAND_HZ,
    apply_trusted_receiver_transfer,
    exact_trusted_mask,
    load_trusted_band,
    native_field_spectrum,
    project_field,
    trusted_bandlimited_timeseries,
    validate_h4_result_contract,
)


F4 = ROOT / "rf/production/f4_attribution"
H4 = ROOT / "fullwave/h4"


def test_stage_f_source_contract_loading_and_native_propagation_flag():
    result = json.loads((H4 / "h4_result_contract.json").read_text())
    for contract in result["source_contracts"].values():
        assert contract["field_components"] == ["Ex_total", "Ey_total", "Ez_total"]
        assert contract["observer_position_m"] == [0.2, 0.0, 0.005]
        assert contract["NATIVE_PROPAGATION_ALREADY_INCLUDED"] is True
        assert (ROOT / contract["source_path"]).is_file()


def test_exact_frozen_stage4_and_stage5_trust_masks():
    assert load_trusted_band(
        F4 / "F4-P-stage4-left-isolated_rf_trust_report.json", STAGE4_TRUSTED_BAND_HZ
    ) == STAGE4_TRUSTED_BAND_HZ
    assert load_trusted_band(
        F4 / "F4-C-stage5-highfield-collision_rf_trust_report.json", STAGE5_TRUSTED_BAND_HZ
    ) == STAGE5_TRUSTED_BAND_HZ
    frequency = np.array([2.9e9, STAGE4_TRUSTED_BAND_HZ[0], 6e9, STAGE4_TRUSTED_BAND_HZ[1], 8e9])
    assert np.array_equal(exact_trusted_mask(frequency, STAGE4_TRUSTED_BAND_HZ), [0, 1, 1, 1, 0])


def test_350mhz_native_voltage_remains_not_resolved():
    status = json.loads((H4 / "h4_350mhz_status.json").read_text())
    assert status["STAGE_F_NATIVE_RF_TRUST"] == "NOT_RESOLVED"
    assert status["H4_NATIVE_RECEIVED_VOLTAGE"] == "NOT_RESOLVED"
    assert status["interpretation"] == "NOT_RESOLVED_IS_NOT_ZERO_RADIATION"


def test_field_vector_projection_and_orthogonal_projection():
    field = np.array([[1 + 2j, 3 - 1j, 4 + 5j], [-1j, 2j, 6 - 3j]])
    assert np.array_equal(project_field(field, [0, 0, 1]), field[:, 2])
    assert np.array_equal(project_field(field, [0, 1, 0]), field[:, 1])


def test_complex_transfer_application_and_no_extrapolation():
    time = np.arange(8) * 1e-10
    field = np.column_stack([np.zeros(8), np.zeros(8), np.cos(2 * np.pi * np.arange(8) / 8)])
    spectrum = native_field_spectrum(time, field)
    receiver_f = np.linspace(0, spectrum.frequency_Hz[-1], 21)
    transfer = 0.002 + 1j * receiver_f * 1e-13
    projected, interpolated, voltage, trusted = apply_trusted_receiver_transfer(
        spectrum, receiver_f, transfer, (spectrum.frequency_Hz[1], spectrum.frequency_Hz[-2])
    )
    assert np.allclose(voltage[trusted], projected[trusted] * interpolated[trusted])
    assert np.all(np.isnan(voltage[~trusted].real))
    with pytest.raises(ValueError, match="EXTRAPOLATION"):
        apply_trusted_receiver_transfer(
            spectrum,
            np.array([2e9, 3e9]),
            np.array([0.002 + 0.1j, 0.002 + 0.2j]),
            (spectrum.frequency_Hz[1], spectrum.frequency_Hz[-2]),
        )


def test_trusted_bandlimited_inverse_preserves_physical_time_phase():
    n = 9
    dt = 1e-10
    start = 3.2e-9
    time = start + np.arange(n) * dt
    signal = np.cos(2 * np.pi * np.arange(n) / n)
    spectrum = native_field_spectrum(time, np.column_stack([signal, signal, signal]))
    voltage = spectrum.complex_spectrum[:, 2]
    recovered = trusted_bandlimited_timeseries(voltage, np.ones(voltage.shape, bool), n, dt, start)
    assert np.allclose(recovered, signal * np.hanning(n), atol=1e-12)


def test_generated_untrusted_bins_are_unavailable_not_zero_and_stages_separate():
    stage4 = np.genfromtxt(H4 / "h4_stage4_native_received_spectrum.csv", delimiter=",", names=True)
    stage5 = np.genfromtxt(H4 / "h4_stage5_native_received_spectrum.csv", delimiter=",", names=True)
    for data, count in ((stage4, 5), (stage5, 7)):
        trusted = data["trust_valid"].astype(bool)
        voltage = data["V_rx_real_V_s"] + 1j * data["V_rx_imag_V_s"]
        assert trusted.sum() == count
        assert np.all(np.isfinite(voltage[trusted]))
        assert np.all(np.isnan(voltage[~trusted].real))
    assert not np.array_equal(stage4["frequency_Hz"], stage5["frequency_Hz"])


def test_stage5_full_maxwell_pending_is_propagated():
    summary = json.loads((H4 / "h4_summary.json").read_text())
    assert summary["stages"]["Stage5"]["full_maxwell_status"] == "FULL_MAXWELL_REFERENCE_PENDING"


def test_f_r4_and_h2_h3_numerical_inputs_are_isolated():
    contract = json.loads((H4 / "h4_result_contract.json").read_text())
    assert contract["mechanism_attribution_consumed"] is False
    assert contract["h2_h3_numerical_inputs_consumed"] is False
    assert validate_h4_result_contract(contract)
    assert contract["RECEIVER_TRANSFER_MESH_SENSITIVITY_PRESENT"] is True
    assert contract["H4_ABSOLUTE_AMPLITUDE_STATUS"] == "NUMERICAL_REFERENCE_ONLY"


def test_openems_orthogonal_polarization_is_suppressed():
    transfer = np.genfromtxt(H4 / "h4_receiver_transfer.csv", delimiter=",", names=True)
    for target in (3.5e9, 6e9, 9.5e9):
        index = int(np.argmin(abs(transfer["frequency_Hz"] - target)))
        assert transfer["cross_polarization_ratio"][index] < 1e-5


def test_h5_result_contract_validation_rejects_native_trust_promotion():
    contract = json.loads((H4 / "h4_result_contract.json").read_text())
    assert validate_h4_result_contract(contract)
    changed = dict(contract, native_350MHz_status="TRUSTED")
    with pytest.raises(ValueError, match="350MHZ"):
        validate_h4_result_contract(changed)


def test_native_source_contract_rejects_non_total_field():
    contract = NativeRFSourceContract(
        source_case_id="synthetic",
        source_path="x.csv",
        source_hash="a" * 64,
        observer_id="observer",
        observer_position_m=(1.0, 0.0, 0.0),
        observer_distance_m=1.0,
        time_start_s=0.0,
        time_end_s=1e-9,
        dt_s=1e-10,
        samples=11,
        field_components=("Ex_dJ", "Ey_dJ", "Ez_dJ"),
        dominant_trusted_component="Ez_dJ",
        receiver_axis=(0.0, 0.0, 1.0),
        trusted_frequency_Hz=(3e9, 8e9),
        rf_trust_report_path="trust.json",
        rf_trust_report_hash="b" * 64,
    )
    with pytest.raises(ValueError, match="TOTAL_NATIVE_FIELD"):
        contract.validate()
