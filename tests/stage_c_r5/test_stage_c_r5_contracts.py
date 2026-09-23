from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "rf" / "c_r5"
CONTRACTS = BASE / "contracts"
DEV = BASE / "development"


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_stage_c_interface_audit_reuses_existing_quantities():
    audit = load(CONTRACTS / "c_r5_stage_c_interface_audit.json")
    assert audit["directly_available"]["impact_ionization_frequency"] == "evaluate_morrow_lowke.ionization_frequency"
    assert audit["directly_available"]["current_moment"] == "ElectronTransportCurrentSource"
    assert audit["physics_formula_duplicated"] is False


def test_event_topology_enum_and_conservative_rules():
    contract = load(CONTRACTS / "c_r5_event_topology_contract.json")
    assert contract["enum"] == [
        "PROPAGATION",
        "HEAD_INTERACTION",
        "STREAMER_COLLISION",
        "ELECTRODE_ATTACHMENT",
        "BRIDGING",
        "UNRESOLVED",
    ]
    assert contract["single_head_can_be_collision"] is False
    assert "at least two distinct components" in contract["rules"]["STREAMER_COLLISION"]
    assert contract["component_detector"]["does_not_change_C_R3_primary_head"] is True


def test_diagnostics_contract_ek_and_same_cell_semantics():
    contract = load(CONTRACTS / "c_r5_ultrafast_diagnostics_contract.json")
    assert contract["Ek"]["implementation"] == "morrow_lowke_breakdown_field"
    assert contract["Ek"]["semantics"] == "REFERENCE_NORMALIZATION_ONLY"
    assert contract["eta_0"]["E0_SOURCE"] == "STAGE_C_2D_ELECTROSTATIC_REFERENCE"
    assert "same cell and time" in contract["kinetics"]["Pi_RF"]["definition"]
    assert contract["kinetics"]["Pi_RF"]["role"] == "DIAGNOSTIC_CANDIDATE"


def test_scalar_trace_schema_and_accepted_step_timing():
    trace = pd.read_csv(DEV / "ultrafast_event_trace.csv")
    required = {
        "step", "time_s", "dt_s", "topology", "head_count", "head_position_m",
        "head_velocity_m_s", "E_peak_V_m", "eta_peak", "ne_peak_m3",
        "nu_i_at_E_peak_s_1", "tau_i_at_E_peak_s", "sigma_e_at_E_peak_S_m",
        "tau_M_at_E_peak_s", "Pi_RF_at_E_peak", "roi_cell_count",
        "bridge_status", "current_moment_z_A_m",
    }
    assert required <= set(trace.columns)
    assert np.array_equal(trace["step"].to_numpy(), np.arange(1, len(trace) + 1))
    assert np.allclose(
        np.diff(np.r_[0.0, trace["time_s"]]),
        trace["dt_s"],
        rtol=1e-12,
        atol=1e-30,
    )
    assert np.all(np.isfinite(trace["Pi_RF_at_E_peak"]))


def test_trace_preserves_reference_metadata_and_no_interpolation_upgrade():
    trace = pd.read_csv(DEV / "ultrafast_event_trace.csv")
    assert set(trace["Ek_semantics"]) == {"REFERENCE_NORMALIZATION_ONLY"}
    assert set(trace["E0_source"]) == {"STAGE_C_2D_ELECTROSTATIC_REFERENCE"}
    report = load(DEV / "resolution_report.json")
    assert report["raw_temporal_metadata"]["DIAGNOSTIC_TRACE_DT"] == "SOLVER_ACCEPTED_STEP_SPACING"
    assert report["raw_temporal_metadata"]["interpolation_used"] is False
    assert report["kinetic_resolution"]["pulse_width_resolution_claimed"] is False


def test_development_reference_is_conservative_and_not_collision():
    summary = load(DEV / "event_summary.json")
    assert summary["case_classification"] == "DEVELOPMENT_REFERENCE"
    assert summary["head_count_max"] == 1
    assert summary["STREAMER_COLLISION_DEMONSTRATION"] == "NOT_AVAILABLE_IN_REFERENCE_CASE"
    assert summary["scientific_semantics"] == "DEVELOPMENT_REFERENCE_NOT_COLLISION_OR_RF_VALIDATION"


def test_handoff_candidate_does_not_modify_stage_f():
    handoff = load(CONTRACTS / "c_r5_stage_f_handoff_candidate.json")
    assert handoff["name"] == "C_R5_ULTRAFAST_DIAGNOSTIC_HANDOFF"
    assert handoff["consumer"] == "F-R5B"
    assert "new Stage-F source definition" in handoff["excluded"]
    assert "pulse-width trust decision" in handoff["excluded"]


def test_stage_c_frozen_baseline_hashes():
    baseline = load(CONTRACTS / "c_r5_stage_c_frozen_baseline.json")
    for relative, expected in baseline["sha256"].items():
        assert sha256(ROOT / relative) == expected


def test_resource_report_declares_scalar_only_history():
    resource = load(DEV / "resource_report.json")
    assert resource["status"] == "PASS"
    assert resource["steps_accepted"] == 24
    assert resource["history_2d_fields_retained"] is False
    assert resource["overhead_status"] in {"PASS", "OVERHEAD_NOT_RESOLVED_SHORT_RUN"}


def test_scientific_status_preserved():
    status = load(DEV / "status.json")
    assert status["C_R5_DIAGNOSTICS_IMPLEMENTED"] is True
    assert status["MICROGAP_ULTRAFAST_RF_MECHANISM"] == "NOT_VALIDATED"
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert status["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False
    assert status["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"
    assert status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert status["LARGE_SIMULATION_RERUN"] is False
    assert status["NEXT_NODE"] == "F-R5B"
