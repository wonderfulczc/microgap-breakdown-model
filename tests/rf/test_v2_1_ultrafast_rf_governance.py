import hashlib
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
ULTRAFAST = ROOT / "rf/ultrafast"


def load_json(name):
    return json.loads((ULTRAFAST / name).read_text(encoding="utf-8"))


def test_koile_evidence_and_table_values_are_complete():
    card = yaml.safe_load(
        (ROOT / "literature/evidence_cards/KOILE_2021_STREAMER_COLLISION_RF.yaml").read_text(encoding="utf-8")
    )
    benchmark = load_json("koile_2021_benchmark_spec.json")
    assert card["paper_id"] == "KOILE_2021_STREAMER_COLLISION_RF"
    assert card["modification_recommendation"] == "CONTRACT_EXTENSION_BEFORE_IMPLEMENTATION"
    assert len(card["processing_chain"]) == 8
    assert card["project_candidate"] == {
        "name": "Pi_RF",
        "formula": "tau_M / tau_i",
        "status": "DIAGNOSTIC_CANDIDATE",
        "paper_parameter": False,
        "novelty_claim_permitted": False,
    }
    group1 = benchmark["group_1"]["cases"]
    group2 = benchmark["group_2"]["cases"]
    assert [(c["separation_m"], c["pulse_width_s"], c["f3dB_Hz"], c["f10dB_Hz"]) for c in group1] == [
        (0.012, 21.2e-12, 2.98e9, 12.6e9),
        (0.014, 21.7e-12, 2.90e9, 12.0e9),
        (0.020, 20.2e-12, 2.92e9, 12.4e9),
    ]
    assert [(c["Eamb_over_Ek"], c["pulse_width_s"], c["f3dB_Hz"], c["f10dB_Hz"]) for c in group2] == [
        (0.65, 21.7e-12, 2.89e9, 12.0e9),
        (0.85, 16.7e-12, 4.53e9, 18.5e9),
        (1.00, 13.7e-12, 5.98e9, 22.6e9),
        (1.50, 9.6e-12, 6.31e9, 30.0e9),
    ]


def test_diagnostics_contract_fields_units_and_candidate_semantics():
    contract = load_json("ultrafast_rf_diagnostics_contract.json")
    required_units = {
        "eta_0": "1",
        "eta_peak": "1",
        "nu_i": "s^-1",
        "tau_i_s": "s",
        "sigma_e_S_m": "S m^-1",
        "tau_M_s": "s",
        "Pi_RF": "1",
        "dMdt_pulse_width_s": "s",
        "f3dB_Hz": "Hz",
        "f10dB_Hz": "Hz",
        "spectral_slope_dB_decade": "dB decade^-1",
        "event_start_s": "s",
        "event_peak_s": "s",
        "event_end_s": "s",
        "event_time_resolution_s": "s",
        "event_sampling_quality": "enum",
        "trust_status": "enum",
    }
    fields = contract["field_definitions"]
    assert {name: fields[name]["unit"] for name in required_units} == required_units
    assert all({"definition", "source", "computation", "applicability", "not_resolved"} <= set(fields[name]) for name in required_units)
    assert fields["Pi_RF"]["scientific_position"] == "DIAGNOSTIC_CANDIDATE"
    assert fields["Pi_RF"]["koile_original_parameter"] is False
    assert fields["Pi_RF"]["novelty_claim_permitted"] is False
    assert contract["Ek_semantics"] == "NORMALIZATION_REFERENCE_ONLY_NOT_A_SUFFICIENT_MICROGAP_BREAKDOWN_CRITERION"


def test_event_topology_and_temporal_resolution_rule():
    contract = load_json("ultrafast_rf_diagnostics_contract.json")
    assert contract["event_topology"]["enum"] == [
        "PROPAGATION",
        "HEAD_INTERACTION",
        "STREAMER_COLLISION",
        "ELECTRODE_ATTACHMENT",
        "BRIDGING",
        "UNRESOLVED",
    ]
    rule = contract["temporal_resolution_rule"]
    assert rule["minimum_effective_samples"] == 5
    assert rule["preferred_effective_samples"] == 10
    assert rule["below_minimum_status"] == "NOT_RESOLVED_TEMPORAL_RESOLUTION"
    assert rule["interpolation_can_upgrade_status"] is False
    quality = lambda samples: "RESOLVED_PREFERRED" if samples >= 10 else "RESOLVED_MINIMUM" if samples >= 5 else rule["below_minimum_status"]
    assert quality(10) == "RESOLVED_PREFERRED"
    assert quality(5) == "RESOLVED_MINIMUM"
    assert quality(4.999) == "NOT_RESOLVED_TEMPORAL_RESOLUTION"


def test_processing_profile_and_sensitivity_gate():
    profile = load_json("koile_compatible_processing_profile.json")
    assert profile["profile_id"] == "KOILE_COMPATIBLE_PROCESSING_PROFILE"
    assert profile["uniform_resampling"]["method"] == "LINEAR_INTERPOLATION"
    assert profile["smoothing"]["window_samples"] == 10
    assert profile["window"]["type"] == "TUKEY"
    assert profile["window"]["alpha_provenance"] == "PROJECT_BENCHMARK_CANDIDATE_NOT_REPORTED_BY_KOILE_2021"
    assert profile["zero_padding"]["target_duration_s"] == 12e-9
    assert "DOES_NOT_INCREASE_TRUE_FREQUENCY_RESOLUTION" in profile["zero_padding"]["physical_resolution_claim"]
    assert profile["reference_frequency_Hz"] == 1e9
    assert profile["pulse_width"]["definition"] == "FULL_WIDTH_AT_ONE_OVER_E_OF_MAIN_ABS_DMDT_PEAK"
    assert set(profile["sensitivity_profiles"]) == {"RAW", "KOILE_COMPATIBLE", "NO_SMOOTHING", "ALTERNATIVE_WINDOW"}
    assert profile["sensitivity_gate"]["material_change_requires_status"] == "PROCESSING_SENSITIVE"
    assert profile["sensitivity_gate"]["automatic_trusted_status_permitted"] is False


def test_current_moment_two_level_interface_is_preserved():
    contract = load_json("ultrafast_rf_diagnostics_contract.json")
    handoff = contract["stage_f_handoff"]
    assert {"M_t_A_m", "dMdt_A_m_s", "FFT_complex", "ESD", "trust_metadata"} <= set(handoff["required_outputs"])
    assert handoff["level_1_mechanisms_preserved"] == ["HCE", "ACC", "REDIS"]
    assert {"tau_i_s", "tau_M_s", "eta_peak", "ne_m_3", "Pi_RF"} <= set(handoff["level_2_candidate_diagnostics"])


def test_stage_names_and_scientific_status_remain_frozen():
    matrix = load_json("v2_1_stage_status_matrix.json")
    assert matrix["formal_stage_names"] == list("ABCDEFGHIJ")
    assert matrix["new_formal_stage_created"] is False
    assert matrix["internal_mechanism_layer"]["formal_stage"] is False
    assert matrix["signal_path_relation"] == "DISTINCT_NOT_INTERCHANGEABLE"
    assert matrix["scientific_status"] == {
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    }
    release_status = json.loads((ROOT / "packaging/release_preparation_status.json").read_text())
    assert release_status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    stage_i = json.loads((ROOT / "validation/stage_i/final/stage_i_final_status.json").read_text())
    assert stage_i["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert stage_i["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"


def test_frozen_v2_0_contract_hashes_are_unchanged():
    status = load_json("v2_1_r0_status.json")
    assert status["FROZEN_V2_0_BASELINE_PRESERVED"] is True
    for relative, expected in status["frozen_baseline_hashes"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == expected, relative


def test_r0_is_contract_only_and_does_not_authorize_physics_changes():
    status = load_json("v2_1_r0_status.json")
    change = yaml.safe_load((ROOT / "literature/change_requests/CR-0001-v2.1-ultrafast-rf.yaml").read_text())
    assert status["V2_1_R0_STATUS"] == "PASS"
    assert status["SIMULATION_ARCHITECTURE"] == "v2.1"
    assert status["SOFTWARE_VERSION"] == "0.1.0"
    assert status["PHYSICS_SOLVER_CHANGED"] is False
    assert status["LARGE_SIMULATION_RERUN"] is False
    assert status["NEXT_NODE"] == "F-R5A"
    assert change["automatic_code_change_permitted"] is False
    assert change["decision_status"] == "APPROVED_FOR_R0_CONTRACT_ONLY"


def test_source_and_packaged_literature_registries_match():
    source = yaml.safe_load((ROOT / "literature/literature_registry.yaml").read_text())
    packaged = yaml.safe_load((ROOT / "python/streamer_rf/resources/literature/literature_registry.yaml").read_text())
    koile = next(p for p in source["papers"] if p["paper_id"] == "KOILE_2021_STREAMER_COLLISION_RF")
    packaged_koile = next(p for p in packaged["papers"] if p["paper_id"] == "KOILE_2021_STREAMER_COLLISION_RF")
    assert koile == packaged_koile
    assert koile["impact_status"] == "ARCHITECTURE_IMPACT"
    assert koile["frozen_result_risk"] == "NONE_R0_CONTRACT_ONLY"
