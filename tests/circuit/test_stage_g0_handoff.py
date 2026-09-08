from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.circuit import (  # noqa: E402
    HandoffCalibration,
    PercolationConfig,
    apply_handoff_logic,
    conductive_percolation,
    g1_initial_state_contract,
    percolation_persistence,
    pi_h,
    tau_evolution,
    tau_sigma,
    xi_sigma,
)
from streamer_rf.circuit.handoff import state_from_cr4_row  # noqa: E402
from streamer_rf.rf.jefimenko.constants import EPS0  # noqa: E402


def _masks(shape=(5, 5)):
    gas = np.ones(shape, dtype=bool)
    hv = np.zeros(shape, dtype=bool)
    ground = np.zeros(shape, dtype=bool)
    hv[:, -1] = True
    ground[:, 0] = True
    return gas, hv, ground


def test_complete_conductive_path_percolates_with_path_length_and_bottleneck() -> None:
    sigma = np.zeros((5, 5))
    sigma[2, :] = [2.0, 3.0, 4.0, 5.0, 6.0]
    gas, hv, ground = _masks()
    result = conductive_percolation(sigma, gas, hv, ground, 1.0, 2.0, PercolationConfig(0.1))
    assert result.percolation_valid
    assert result.percolation_status == "PERCOLATED"
    assert result.path_length_m == pytest.approx(8.0)
    assert result.bottleneck_sigma_S_m == pytest.approx(2.0)
    assert result.percolated_cell_count == 5


def test_broken_path_and_two_half_paths_are_false() -> None:
    gas, hv, ground = _masks()
    sigma = np.zeros((5, 5))
    sigma[2, [0, 1, 3, 4]] = 1.0
    assert conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0).percolation_status == "NO_CONDUCTIVE_PATH"
    sigma2 = np.zeros((5, 5))
    sigma2[1, [0, 1]] = 1.0
    sigma2[3, [3, 4]] = 1.0
    assert conductive_percolation(sigma2, gas, hv, ground, 1.0, 1.0).percolation_status == "NO_CONDUCTIVE_PATH"


def test_isolated_conductive_island_does_not_bridge() -> None:
    gas, hv, ground = _masks()
    sigma = np.zeros((5, 5))
    sigma[2, 2] = 10.0
    result = conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0)
    assert not result.percolation_valid
    assert result.percolation_status == "NO_ELECTRODE_ADJACENT_SEEDS"


def test_conductor_gas_mask_blocks_traversal() -> None:
    gas, hv, ground = _masks()
    gas[2, 2] = False
    sigma = np.zeros((5, 5))
    sigma[2, :] = 1.0
    result = conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0)
    assert result.percolation_status == "NO_CONDUCTIVE_PATH"


def test_electrode_adjacent_start_end_detection() -> None:
    gas, hv, ground = _masks()
    sigma = np.ones((5, 5))
    hv[:] = False
    result = conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0)
    assert result.percolation_status == "NO_ELECTRODE_ADJACENT_SEEDS"


def test_connectivity_is_deterministic_under_repeated_calls() -> None:
    gas, hv, ground = _masks()
    sigma = np.ones((5, 5))
    a = conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0)
    b = conductive_percolation(sigma, gas, hv, ground, 1.0, 1.0)
    assert a == b


def test_threshold_parameter_handling_and_zero_conductivity() -> None:
    gas, hv, ground = _masks()
    sigma = np.ones((3, 3))
    assert conductive_percolation(sigma, gas[:3, :3], hv[:3, :3], ground[:3, :3], 1.0, 1.0, PercolationConfig(-1.0)).percolation_status == "INVALID_THRESHOLD"
    zero = conductive_percolation(np.zeros((3, 3)), gas[:3, :3], hv[:3, :3], ground[:3, :3], 1.0, 1.0)
    assert zero.percolation_status == "NO_CONDUCTIVE_CELLS"


def test_tau_sigma_tau_evolution_and_xi_sigma() -> None:
    ts, ts_status = tau_sigma(2.0)
    assert ts_status == "VALID"
    assert ts == pytest.approx(EPS0 / 2.0)
    te, te_status = tau_evolution(4.0, -2.0)
    assert te_status == "VALID"
    assert te == pytest.approx(2.0)
    xi, xi_status = xi_sigma(2.0, 4.0, -2.0)
    assert xi_status == "VALID"
    assert xi == pytest.approx(EPS0 / 4.0)


def test_invalid_timescale_inputs_remain_unavailable() -> None:
    assert tau_sigma(0.0)[1] == "INVALID_SIGMA_EFF"
    assert tau_evolution(0.0, 1.0)[1] == "ZERO_GB_OR_DGBDT"
    assert tau_evolution(1.0, 0.0)[1] == "ZERO_GB_OR_DGBDT"
    assert xi_sigma(math.nan, 1.0, 1.0)[1] == "UNAVAILABLE"


def test_pi_h_missing_and_valid_reference() -> None:
    value, status = pi_h(1.0, None)
    assert math.isnan(value)
    assert status == "THERMAL_ENERGY_REFERENCE_NOT_AVAILABLE"
    value, status = pi_h(2.0, 4.0)
    assert value == 0.5
    assert status == "VALID"


def test_no_calibration_never_outputs_production_handoff_ready() -> None:
    status, Pi_H, cal_status = apply_handoff_logic(
        bridge_or_percolated=True,
        QJ_channel_J=1.0,
        Xi_sigma=0.1,
        calibration=None,
    )
    assert status == "HANDOFF_CALIBRATION_PENDING"
    assert math.isnan(Pi_H)
    assert cal_status == "NOT_AVAILABLE"


def test_valid_synthetic_calibration_is_deterministic() -> None:
    calibration = HandoffCalibration(
        source_id="synthetic",
        source_reference="unit_test",
        Q_required_J=10.0,
        Pi_H_threshold=0.2,
        Xi_sigma_threshold=0.5,
        bridge_requirement=True,
    )
    assert apply_handoff_logic(bridge_or_percolated=True, QJ_channel_J=3.0, Xi_sigma=0.1, calibration=calibration)[0] == "HANDOFF_CRITERIA_MET"
    assert apply_handoff_logic(bridge_or_percolated=False, QJ_channel_J=3.0, Xi_sigma=0.1, calibration=calibration)[0] == "HANDOFF_CRITERIA_NOT_MET"


def test_malformed_calibration_rejected() -> None:
    with pytest.raises(ValueError):
        HandoffCalibration("", "ref", Q_required_J=1.0).validate()
    with pytest.raises(ValueError):
        HandoffCalibration("id", "ref", Q_required_J=-1.0).validate()


def test_initial_sample_and_missing_fields_state() -> None:
    row = {
        "time_s": "0.0",
        "bridge_flag": "0",
        "PJ_channel_W": "nan",
        "QJ_channel_J": "nan",
        "channel_volume_m3": "nan",
        "channel_length_m": "nan",
        "channel_effective_radius_m": "nan",
        "sigma_eff_S_m": "nan",
        "E_channel_mean_V_m": "nan",
        "ne_channel_mean_m3": "nan",
        "ne_channel_max_m3": "nan",
        "Gb_S": "0.0",
        "Rb_ohm": "inf",
        "dGb_dt_S_s": "nan",
        "tau_sigma_s": "nan",
        "tau_evolution_s": "nan",
        "Xi_sigma": "nan",
    }
    state = state_from_cr4_row(row)
    assert state.handoff_status == "HANDOFF_CALIBRATION_PENDING"
    assert state.percolation_status == "NO_FIELD_SNAPSHOT"
    assert math.isnan(state.Pi_H)


def test_percolation_persistence_summary() -> None:
    p = percolation_persistence([0.0, 1.0, 2.0, 3.0], [False, True, True, False])
    assert p.first_percolation_time_s == 1.0
    assert p.number_of_percolated_samples == 2
    assert p.longest_consecutive_percolated_samples == 2


def test_g1_initial_state_contract_has_no_thermal_profile() -> None:
    assert g1_initial_state_contract(None)["thermal_profile_status"] == "NOT_AVAILABLE"
    row = {
        "time_s": "1.0",
        "bridge_flag": "1",
        "PJ_channel_W": "2.0",
        "QJ_channel_J": "3.0",
        "channel_volume_m3": "4.0",
        "channel_length_m": "5.0",
        "channel_effective_radius_m": "6.0",
        "sigma_eff_S_m": "7.0",
        "E_channel_mean_V_m": "8.0",
        "ne_channel_mean_m3": "9.0",
        "ne_channel_max_m3": "10.0",
        "Gb_S": "11.0",
        "Rb_ohm": "12.0",
        "dGb_dt_S_s": "13.0",
        "tau_sigma_s": "14.0",
        "tau_evolution_s": "15.0",
        "Xi_sigma": "16.0",
    }
    state = state_from_cr4_row(row)
    contract = g1_initial_state_contract(state)
    assert contract["state_valid"] is True
    assert contract["handoff_time_candidate"] == 1.0
    assert contract["thermal_profile_status"] == "NOT_AVAILABLE"
