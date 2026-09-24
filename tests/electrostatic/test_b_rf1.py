from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from streamer_rf.electrostatic import (
    GeometryCompatibility,
    classify_geometry,
    compare_mesh_levels,
    cross_validate_stage_c,
    extract_descriptor,
    f_r6b_gate,
    load_comsol_export,
    make_stage_c_handoff,
    paper1_descriptor,
    route_geometry,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "rf/b_rf1"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def fixture(level: str = "base") -> tuple[dict, pd.DataFrame]:
    return load_comsol_export(OUT / "fixtures" / level)


def test_comsol_export_schema_units_coordinates_and_hash(tmp_path: Path) -> None:
    metadata, frame = fixture()
    assert metadata["units"] == "SI"
    assert {"x_m", "y_m", "potential_V", "Ex_V_m", "Ey_V_m", "E_mag_V_m"} <= set(frame)

    copied = tmp_path / "export"
    copied.mkdir()
    field_path = copied / "field.csv"
    frame.to_csv(field_path, index=False)
    broken = deepcopy(metadata)
    broken["units"] = "mm-kV"
    broken["source_provenance"]["field_sha256"] = hashlib.sha256(field_path.read_bytes()).hexdigest()
    (copied / "metadata.json").write_text(json.dumps(broken))
    with pytest.raises(ValueError, match="SI exports only"):
        load_comsol_export(copied)

    broken["units"] = "SI"
    broken["coordinate_system"]["handedness"] = "UNDECLARED"
    (copied / "metadata.json").write_text(json.dumps(broken))
    with pytest.raises(ValueError, match="coordinate handedness"):
        load_comsol_export(copied)


def test_hdf5_export_uses_same_contract(tmp_path: Path) -> None:
    metadata, frame = fixture()
    target = tmp_path / "hdf-export"
    target.mkdir()
    field_path = target / "field.h5"
    with h5py.File(field_path, "w") as handle:
        group = handle.create_group("field")
        for column in frame.columns:
            values = frame[column].to_numpy()
            if values.dtype.kind in {"O", "U"}:
                values = values.astype("S")
            group.create_dataset(column, data=values)
    metadata = deepcopy(metadata)
    metadata["field_file"] = field_path.name
    metadata["source_provenance"]["field_sha256"] = hashlib.sha256(field_path.read_bytes()).hexdigest()
    (target / "metadata.json").write_text(json.dumps(metadata))
    loaded_metadata, loaded_frame = load_comsol_export(target)
    assert loaded_metadata["field_file"] == "field.h5"
    assert len(loaded_frame) == len(frame)
    assert extract_descriptor(loaded_metadata, loaded_frame)["roi_statistics"]["sample_count"] > 0


def test_roi_and_electrostatic_descriptor_definitions() -> None:
    metadata, frame = fixture("refined")
    descriptor = extract_descriptor(metadata, frame)
    stats = descriptor["roi_statistics"]
    assert stats["sample_count"] < len(frame)
    roi = metadata["gap_analysis_roi"]
    mask = (
        (frame.x_m >= roi["x_bounds_m"][0]) & (frame.x_m <= roi["x_bounds_m"][1])
        & (frame.y_m >= roi["y_bounds_m"][0]) & (frame.y_m <= roi["y_bounds_m"][1])
    )
    values = frame.loc[mask, "E_mag_V_m"].to_numpy()
    assert stats["E0_max_V_m"] == np.max(values)
    assert stats["E0_mean_V_m"] > 0.0
    nominal = metadata["applied_voltage_V"] / metadata["geometry"]["gap_m"]
    assert np.isclose(descriptor["field_enhancement"]["beta_E"], stats["E0_max_V_m"] / nominal)
    assert np.isclose(
        descriptor["ek_normalization"]["eta_0_p95"],
        stats["E0_p95_V_m"] / metadata["ek_reference"]["Ek_V_m"],
    )
    assert descriptor["ek_normalization"]["semantics"] == "REFERENCE_NORMALIZATION_ONLY"
    assert descriptor["field_enhancement"]["is_breakdown_criterion"] is False


def test_high_field_localization_event_candidate_and_cgap() -> None:
    metadata, frame = fixture("refined")
    descriptor = extract_descriptor(metadata, frame)
    regions = descriptor["high_field_regions"]
    assert regions["q70"]["area_m2"] > regions["q80"]["area_m2"] > regions["q90"]["area_m2"] > 0.0
    assert regions["q70"]["sample_count"] > regions["q80"]["sample_count"] > regions["q90"]["sample_count"]
    assert descriptor["field_localization"]["status"] == "PASS"
    assert descriptor["field_localization"]["L_E_at_Emax_m"] > 0.0
    event = descriptor["static_event_region_candidate"]
    assert event["role"] == "STATIC_EVENT_REGION_CANDIDATE"
    assert event["actual_streamer_initiation_asserted"] is False
    assert descriptor["Cgap"]["role"] == "UPSTREAM_GAP_CAPACITANCE"
    assert descriptor["Cgap"]["unit"] == "F"
    assert descriptor["Cgap"]["reference_terminals"] == ["high_voltage", "ground"]


def test_geometry_compatibility_is_conservative_and_routes_correctly() -> None:
    metadata, _ = fixture()
    assert classify_geometry(metadata) == GeometryCompatibility.AXISYMMETRIC_2D_EXACT
    assert route_geometry(GeometryCompatibility.AXISYMMETRIC_2D_EXACT)["stage_c_allowed"] is True

    foil = deepcopy(metadata)
    foil["geometry"].update({
        "electrode_type": "TRIANGULAR_FOIL_EDGE",
        "axisymmetric": False,
        "alignment": "OFFSET",
        "compatibility_class": "NONAXISYMMETRIC_3D_REQUIRED",
    })
    assert classify_geometry(foil) == GeometryCompatibility.NONAXISYMMETRIC_3D_REQUIRED
    route = route_geometry(classify_geometry(foil))
    assert route["stage_c_allowed"] is False
    assert route["stage_e_r2_required"] is True

    false_exact = deepcopy(foil)
    false_exact["geometry"]["compatibility_class"] = "AXISYMMETRIC_2D_EXACT"
    assert classify_geometry(false_exact) == GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED


def test_mesh_convergence_and_singularity_handling() -> None:
    base = load_json(OUT / "fixture_descriptor_base.json")
    refined = load_json(OUT / "fixture_descriptor_refined.json")
    comparison = compare_mesh_levels(base, refined)
    assert comparison["status"] == "PASS"
    assert comparison["primary_comparison_uses_robust_metrics"] is True

    metadata, frame = fixture()
    singular = deepcopy(metadata)
    singular["geometry"]["singularity"] = {"status": "IDEAL_GEOMETRIC_SINGULARITY", "zero_radius_corner": True}
    descriptor = extract_descriptor(singular, frame)
    assert descriptor["singularity"]["Emax_primary_cross_geometry_metric"] is False
    assert "E0_p95" in descriptor["singularity"]["preferred_metrics"]


def test_stage_c_handoff_and_crossvalidation_boundary() -> None:
    descriptor = load_json(OUT / "fixture_descriptor_refined.json")
    handoff = make_stage_c_handoff(descriptor)
    assert handoff["status"] == "READY"
    assert "self-consistent Stage-C Poisson" in handoff["forbidden_use"]

    stage_c = deepcopy(descriptor)
    stage_c["roi_statistics"]["E0_p95_V_m"] *= 1.01
    stage_c["field_enhancement"]["beta_E_p95"] *= 1.01
    axis = np.linspace(1.0, 2.0, 20)
    result = cross_validate_stage_c(descriptor, stage_c, axis, axis * 1.01)
    assert result["status"] == "PASS"
    assert result["hard_parameter_tuning_to_match_forbidden"] is True


def test_paper1_descriptor_and_f_r6b_gate_reject_fixture() -> None:
    descriptor = load_json(OUT / "fixture_descriptor_refined.json")
    paper = paper1_descriptor(descriptor, material_interface_id="fixture-interface")
    assert paper["descriptor_type"] == "PAPER1_ELECTRODE_DESCRIPTOR"
    assert paper["scientific_result_allowed"] is False
    assert paper["downstream_links"] == {
        "streamer_metrics": None, "rf_metrics": None, "experimental_metrics": None,
    }
    assert f_r6b_gate([descriptor, deepcopy(descriptor)]) is False

    actual_a = deepcopy(descriptor)
    actual_b = deepcopy(descriptor)
    for index, actual in enumerate((actual_a, actual_b), start=1):
        actual["dataset_role"] = "ACTUAL_STAGE_B_EXPORT"
        actual["scientific_result_allowed"] = True
        actual["mesh"]["status"] = "PASS"
        actual["geometry_identity"]["geometry_id"] = f"actual-{index}"
    assert f_r6b_gate([actual_a, actual_b]) is True


def test_generated_status_has_no_scientific_upgrade() -> None:
    status = load_json(OUT / "status.json")
    audit = load_json(OUT / "actual_data_audit.json")
    assert status["B_RF1_ENGINEERING"] == "PASS"
    assert audit["ACTUAL_STAGE_B_DATA"] == "NOT_AVAILABLE"
    assert status["ACTUAL_ELECTRODE_DESCRIPTOR_COUNT"] == 0
    assert status["ACTUAL_ELECTRODE_DESCRIPTOR_READY"] is False
    assert status["F_R6B_ALLOWED"] is False
    assert status["ACTUAL_ELECTRODE_RF_MECHANISM_VALIDATED"] is False
    assert status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"


def test_frozen_stage_c_baseline_and_f_r6a_status_are_preserved() -> None:
    baseline = load_json(ROOT / "rf/c_r5/contracts/c_r5_stage_c_frozen_baseline.json")
    for name, expected in baseline["sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    f_r6a = load_json(ROOT / "rf/f_r6a/status.json")
    assert f_r6a["F_R6A_STATUS"] == "PASS"
    assert f_r6a["F_R6B_ALLOWED"] is False
