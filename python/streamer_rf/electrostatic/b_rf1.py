"""Read-only Stage-B electrostatic import and B-RF1 descriptor extraction."""

from __future__ import annotations

import hashlib
import json
import math
from enum import StrEnum
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np
import pandas as pd


class GeometryCompatibility(StrEnum):
    AXISYMMETRIC_2D_EXACT = "AXISYMMETRIC_2D_EXACT"
    AXISYMMETRIC_2D_APPROXIMATE = "AXISYMMETRIC_2D_APPROXIMATE"
    NONAXISYMMETRIC_3D_REQUIRED = "NONAXISYMMETRIC_3D_REQUIRED"
    GEOMETRY_MAPPING_UNRESOLVED = "GEOMETRY_MAPPING_UNRESOLVED"


_REQUIRED_METADATA = {
    "schema_version", "dataset_role", "geometry", "coordinate_system",
    "field_file", "applied_voltage_V", "gas", "pressure_Pa",
    "temperature_K", "gap_analysis_roi", "mesh", "capacitance",
    "ek_reference", "source_provenance",
}
_REQUIRED_GEOMETRY = {
    "geometry_id", "electrode_type", "geometry_revision", "gap_m",
    "tip_radius_m", "electrode_material", "alignment", "source_file",
    "source_hash", "comsol_model_identifier", "simulation_date",
}
_REQUIRED_FIELD_COLUMNS = {
    "x_m", "y_m", "potential_V", "Ex_V_m", "Ey_V_m", "E_mag_V_m",
    "domain_id",
}
_ALLOWED_ROLES = {"ACTUAL_STAGE_B_EXPORT", "SYNTHETIC_DEVELOPMENT_FIXTURE"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_comsol_export(directory: str | Path) -> tuple[dict[str, Any], pd.DataFrame]:
    """Load and validate a standardized COMSOL export without COMSOL itself."""
    root = Path(directory)
    metadata_path = root / "metadata.json"
    if not metadata_path.is_file():
        raise ValueError("metadata.json is required")
    metadata = json.loads(metadata_path.read_text())
    missing = sorted(_REQUIRED_METADATA - metadata.keys())
    if missing:
        raise ValueError(f"missing metadata fields: {missing}")
    geometry_missing = sorted(_REQUIRED_GEOMETRY - metadata["geometry"].keys())
    if geometry_missing:
        raise ValueError(f"missing geometry identity fields: {geometry_missing}")
    if metadata["dataset_role"] not in _ALLOWED_ROLES:
        raise ValueError("dataset_role is not recognized")
    if metadata.get("units") != "SI":
        raise ValueError("B-RF1 accepts SI exports only")
    if metadata["coordinate_system"].get("handedness") not in {"RIGHT_HANDED", "AXISYMMETRIC_RZ"}:
        raise ValueError("coordinate handedness is missing or unsupported")
    if metadata["ek_reference"].get("implementation") != "morrow_lowke_breakdown_field":
        raise ValueError("Ek must come from the unified project implementation")
    if metadata["ek_reference"].get("semantics") != "REFERENCE_NORMALIZATION_ONLY":
        raise ValueError("Ek semantics must remain reference normalization only")

    field_path = root / metadata["field_file"]
    if field_path.suffix.lower() == ".csv":
        frame = pd.read_csv(field_path)
    elif field_path.suffix.lower() in {".h5", ".hdf5"}:
        with h5py.File(field_path, "r") as handle:
            group = handle[metadata.get("hdf_key", "field")]
            frame = pd.DataFrame({name: group[name][()] for name in group.keys()})
    else:
        raise ValueError("field export must be CSV or HDF5")
    frame = frame.reset_index(drop=True)
    if "domain_id" in frame and frame["domain_id"].dtype.kind == "S":
        frame["domain_id"] = frame["domain_id"].str.decode("utf-8")
    missing_columns = sorted(_REQUIRED_FIELD_COLUMNS - set(frame.columns))
    if missing_columns:
        raise ValueError(f"missing field columns: {missing_columns}")
    if metadata["coordinate_system"].get("dimension") == 3:
        required_3d = {"z_m", "Ez_V_m"}
        if not required_3d <= set(frame.columns):
            raise ValueError("3D export requires z_m and Ez_V_m")
    numeric = [column for column in frame.columns if column != "domain_id"]
    if not np.all(np.isfinite(frame[numeric].to_numpy(dtype=float))):
        raise ValueError("field export contains non-finite numeric values")
    if np.any(frame["E_mag_V_m"].to_numpy() < 0.0):
        raise ValueError("E magnitude cannot be negative")
    components = frame[["Ex_V_m", "Ey_V_m"]].to_numpy(dtype=float)
    if "Ez_V_m" in frame:
        components = np.column_stack((components, frame["Ez_V_m"].to_numpy(dtype=float)))
    component_magnitude = np.linalg.norm(components, axis=1)
    if not np.allclose(component_magnitude, frame["E_mag_V_m"], rtol=5e-6, atol=1e-9):
        raise ValueError("E_mag_V_m is inconsistent with field components")
    declared_hash = metadata["source_provenance"].get("field_sha256")
    if declared_hash and declared_hash != _sha256(field_path):
        raise ValueError("field export SHA-256 mismatch")
    return metadata, frame


def classify_geometry(metadata: Mapping[str, Any]) -> GeometryCompatibility:
    geometry = metadata["geometry"]
    declared = geometry.get("compatibility_class")
    try:
        declared_class = GeometryCompatibility(declared)
    except (ValueError, TypeError):
        return GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED
    axisymmetric = geometry.get("axisymmetric") is True
    alignment = geometry.get("alignment")
    electrode_type = str(geometry.get("electrode_type", "")).upper()
    inherently_3d = any(token in electrode_type for token in ("FOIL_EDGE", "TRIANGULAR_FOIL", "IRREGULAR"))
    if declared_class == GeometryCompatibility.AXISYMMETRIC_2D_EXACT:
        if axisymmetric and alignment in {"COAXIAL", "AXIS_ALIGNED"} and not inherently_3d:
            return declared_class
        return GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED
    if declared_class == GeometryCompatibility.AXISYMMETRIC_2D_APPROXIMATE:
        if geometry.get("mapping_rationale") and not inherently_3d:
            return declared_class
        return GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED
    if declared_class == GeometryCompatibility.NONAXISYMMETRIC_3D_REQUIRED:
        return declared_class
    return GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED


def route_geometry(compatibility: GeometryCompatibility) -> dict[str, Any]:
    stage_c = compatibility in {
        GeometryCompatibility.AXISYMMETRIC_2D_EXACT,
        GeometryCompatibility.AXISYMMETRIC_2D_APPROXIMATE,
    }
    return {
        "stage_c_allowed": stage_c,
        "stage_e_r2_required": compatibility == GeometryCompatibility.NONAXISYMMETRIC_3D_REQUIRED,
        "routing_status": "RESOLVED" if compatibility != GeometryCompatibility.GEOMETRY_MAPPING_UNRESOLVED else "NOT_RESOLVED",
    }


def _roi_mask(frame: pd.DataFrame, roi: Mapping[str, Any]) -> np.ndarray:
    if roi.get("type") != "AXIS_ALIGNED_BOUNDS":
        raise ValueError("only explicit AXIS_ALIGNED_BOUNDS ROI is supported")
    mask = np.ones(len(frame), dtype=bool)
    for axis in ("x", "y", "z"):
        column = f"{axis}_m"
        bounds = roi.get(f"{axis}_bounds_m")
        if column in frame and bounds is not None:
            lower, upper = map(float, bounds)
            if not math.isfinite(lower) or not math.isfinite(upper) or lower >= upper:
                raise ValueError(f"invalid {axis} ROI bounds")
            values = frame[column].to_numpy(dtype=float)
            mask &= (values >= lower) & (values <= upper)
    domains = roi.get("domain_ids")
    if domains:
        mask &= frame["domain_id"].astype(str).isin([str(value) for value in domains]).to_numpy()
    if not np.any(mask):
        raise ValueError("GAP_ANALYSIS_ROI contains no field samples")
    return mask


def _weights(frame: pd.DataFrame, dimension: int) -> tuple[np.ndarray, str]:
    candidate = "cell_volume_m3" if dimension == 3 else "cell_area_m2"
    if candidate in frame:
        values = frame[candidate].to_numpy(dtype=float)
        if np.all(np.isfinite(values)) and np.all(values > 0.0):
            return values, candidate
    return np.ones(len(frame), dtype=float), "UNIFORM_SAMPLE_WEIGHT"


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(values)
    ordered_values = values[order]
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights) - 0.5 * ordered_weights
    cumulative /= np.sum(ordered_weights)
    return float(np.interp(q, cumulative, ordered_values))


def _field_localization(frame: pd.DataFrame, roi_mask: np.ndarray, high_mask: np.ndarray) -> dict[str, Any]:
    selected = roi_mask & high_mask
    l_at_max: float | None = None
    if "grad_E_mag_V_m2" in frame:
        gradient = frame.loc[selected, "grad_E_mag_V_m2"].to_numpy(dtype=float)
        field = frame.loc[selected, "E_mag_V_m"].to_numpy(dtype=float)
        valid = np.isfinite(gradient) & (gradient > 0.0)
        values = field[valid] / gradient[valid]
        roi_indices = np.flatnonzero(roi_mask)
        global_max_index = roi_indices[int(np.argmax(frame.loc[roi_mask, "E_mag_V_m"].to_numpy()))]
        gradient_at_max = float(frame.iloc[global_max_index]["grad_E_mag_V_m2"])
        if gradient_at_max > 0.0:
            l_at_max = float(frame.iloc[global_max_index]["E_mag_V_m"] / gradient_at_max)
    else:
        dimension = 3 if "z_m" in frame else 2
        axes = ["x_m", "y_m"] + (["z_m"] if dimension == 3 else [])
        shape = tuple(frame[axis].nunique() for axis in axes)
        if int(np.prod(shape)) != len(frame):
            return {"status": "NOT_RESOLVED_UNSTRUCTURED_GRID"}
        ordered = frame.sort_values(axes)
        expected = pd.MultiIndex.from_product([np.sort(frame[a].unique()) for a in axes])
        actual = pd.MultiIndex.from_frame(ordered[axes])
        if not actual.equals(expected):
            return {"status": "NOT_RESOLVED_UNSTRUCTURED_GRID"}
        field_grid = ordered["E_mag_V_m"].to_numpy().reshape(shape)
        spacings = [np.diff(np.sort(frame[a].unique())).mean() for a in axes]
        gradients = np.gradient(field_grid, *spacings, edge_order=2)
        grad = np.sqrt(sum(component * component for component in gradients)).ravel()
        ordered_selected = (roi_mask & high_mask)[ordered.index.to_numpy()]
        valid = ordered_selected & np.isfinite(grad) & (grad > 0.0)
        values = field_grid.ravel()[valid] / grad[valid]
        ordered_roi = roi_mask[ordered.index.to_numpy()]
        max_ordered_position = int(np.argmax(np.where(ordered_roi, field_grid.ravel(), -np.inf)))
        if grad[max_ordered_position] > 0.0:
            l_at_max = float(field_grid.ravel()[max_ordered_position] / grad[max_ordered_position])
    if values.size == 0:
        return {"status": "NOT_RESOLVED_ZERO_GRADIENT"}
    return {
        "status": "PASS",
        "definition": "L_E=|E|/|grad|E||",
        "L_E_at_Emax_m": l_at_max,
        "L_E_p50_highfield_m": float(np.quantile(values, 0.50)),
        "L_E_p95_highfield_m": float(np.quantile(values, 0.95)),
    }


def extract_descriptor(metadata: Mapping[str, Any], frame: pd.DataFrame) -> dict[str, Any]:
    """Extract comparable actual-electrode descriptors inside GAP_ANALYSIS_ROI."""
    role = metadata["dataset_role"]
    dimension = int(metadata["coordinate_system"]["dimension"])
    roi_mask = _roi_mask(frame, metadata["gap_analysis_roi"])
    roi = frame.loc[roi_mask].copy()
    weights, weight_source = _weights(roi, dimension)
    field = roi["E_mag_V_m"].to_numpy(dtype=float)
    emax = float(np.max(field))
    emean = float(np.average(field, weights=weights))
    ep95 = _weighted_quantile(field, weights, 0.95)
    voltage = abs(float(metadata["applied_voltage_V"]))
    gap = float(metadata["geometry"]["gap_m"])
    if not math.isfinite(voltage) or voltage <= 0.0 or not math.isfinite(gap) or gap <= 0.0:
        raise ValueError("applied voltage and gap must be positive")
    nominal = voltage / gap
    ek = float(metadata["ek_reference"]["Ek_V_m"])
    if not math.isfinite(ek) or ek <= 0.0:
        raise ValueError("project Ek must be positive")

    axes = ["x_m", "y_m"] + (["z_m"] if dimension == 3 else [])
    high_field: dict[str, Any] = {}
    for q in (0.70, 0.80, 0.90):
        label = f"q{int(q * 100)}"
        selected = field >= q * emax
        selected_weights = weights[selected]
        coordinates = roi.loc[selected, axes].to_numpy(dtype=float)
        centroid = np.average(coordinates, axis=0, weights=selected_weights)
        extent = np.ptp(coordinates, axis=0)
        measure_name = "volume_m3" if dimension == 3 else "area_m2"
        high_field[label] = {
            "threshold_fraction_of_E0_max": q,
            measure_name: float(np.sum(selected_weights)),
            "sample_count": int(np.sum(selected)),
            "centroid_m": centroid.tolist(),
            "extent_m": extent.tolist(),
        }

    q90_global = np.zeros(len(frame), dtype=bool)
    q90_global[np.flatnonzero(roi_mask)] = field >= 0.90 * emax
    localization = _field_localization(frame, roi_mask, q90_global)
    event_region = high_field["q90"] | {
        "role": "STATIC_EVENT_REGION_CANDIDATE",
        "definition": "q90 high-field region containing the E0 maximum candidate",
        "actual_streamer_initiation_asserted": False,
        "distance_to_electrode_m": metadata["gap_analysis_roi"].get("candidate_distance_to_electrode_m"),
        "E0_max_region_V_m": emax,
        "E0_p95_region_V_m": _weighted_quantile(field[field >= 0.9 * emax], weights[field >= 0.9 * emax], 0.95),
    }
    compatibility = classify_geometry(metadata)
    capacitance = metadata["capacitance"]
    if capacitance.get("role") != "UPSTREAM_GAP_CAPACITANCE":
        raise ValueError("Cgap role must remain UPSTREAM_GAP_CAPACITANCE")
    if capacitance.get("unit") != "F" or not math.isfinite(float(capacitance.get("value_F", math.nan))) or float(capacitance["value_F"]) <= 0.0:
        raise ValueError("Cgap must be a positive SI capacitance in F")
    if len(capacitance.get("reference_terminals", [])) != 2 or not capacitance.get("calculation_method"):
        raise ValueError("Cgap reference terminals and calculation method are required")
    result = {
        "schema_version": "1.0",
        "descriptor_type": "B_RF1_ELECTROSTATIC_DESCRIPTOR",
        "dataset_role": role,
        "scientific_result_allowed": role == "ACTUAL_STAGE_B_EXPORT" and metadata["mesh"].get("status") == "PASS",
        "geometry_identity": dict(metadata["geometry"]),
        "compatibility_class": compatibility.value,
        "routing": route_geometry(compatibility),
        "gap_analysis_roi": dict(metadata["gap_analysis_roi"]),
        "roi_statistics": {
            "sample_count": int(np.sum(roi_mask)),
            "measure_weight_source": weight_source,
            "measure_m3" if dimension == 3 else "measure_m2": float(np.sum(weights)),
            "E0_max_V_m": emax,
            "E0_mean_V_m": emean,
            "E0_p95_V_m": ep95,
        },
        "field_enhancement": {
            "E_nominal_V_m": nominal,
            "beta_E": emax / nominal,
            "beta_E_p95": ep95 / nominal,
            "is_breakdown_criterion": False,
        },
        "ek_normalization": {
            **dict(metadata["ek_reference"]),
            "eta_0_max": emax / ek,
            "eta_0_p95": ep95 / ek,
        },
        "high_field_regions": high_field,
        "field_localization": localization,
        "static_event_region_candidate": event_region,
        "Cgap": dict(capacitance),
        "mesh": dict(metadata["mesh"]),
        "singularity": dict(metadata["geometry"].get("singularity", {"status": "NONE_DECLARED"})),
        "provenance": dict(metadata["source_provenance"]),
        "scientific_boundary": "static E-field descriptors are not RF-performance conclusions",
    }
    if result["singularity"].get("status") == "IDEAL_GEOMETRIC_SINGULARITY":
        result["singularity"]["Emax_primary_cross_geometry_metric"] = False
        result["singularity"]["preferred_metrics"] = ["E0_p95", "beta_E_p95", "high_field_extent", "L_E"]
    return result


def compare_mesh_levels(base: Mapping[str, Any], refined: Mapping[str, Any], *, tolerance: float = 0.05) -> dict[str, Any]:
    if base["geometry_identity"]["geometry_id"] != refined["geometry_identity"]["geometry_id"]:
        raise ValueError("mesh levels must describe the same geometry_id")
    def rel(a: float, b: float) -> float:
        return abs(b - a) / max(abs(b), 1e-300)
    measure_key = "volume_m3" if "volume_m3" in base["high_field_regions"]["q90"] else "area_m2"
    changes = {
        "E0_max_relative_difference": rel(base["roi_statistics"]["E0_max_V_m"], refined["roi_statistics"]["E0_max_V_m"]),
        "E0_p95_relative_difference": rel(base["roi_statistics"]["E0_p95_V_m"], refined["roi_statistics"]["E0_p95_V_m"]),
        "beta_E_p95_relative_difference": rel(base["field_enhancement"]["beta_E_p95"], refined["field_enhancement"]["beta_E_p95"]),
        "Cgap_relative_difference": rel(base["Cgap"]["value_F"], refined["Cgap"]["value_F"]),
        "high_field_q90_measure_relative_difference": rel(
            base["high_field_regions"]["q90"][measure_key], refined["high_field_regions"]["q90"][measure_key]
        ),
    }
    robust_pass = all(changes[key] <= tolerance for key in (
        "E0_p95_relative_difference", "beta_E_p95_relative_difference", "Cgap_relative_difference",
        "high_field_q90_measure_relative_difference",
    ))
    emax_sensitive = changes["E0_max_relative_difference"] > tolerance
    return {
        "status": "PASS" if robust_pass else "NOT_RESOLVED",
        "tolerance": tolerance,
        **changes,
        "Emax_status": "E_MAX_MESH_SENSITIVE" if emax_sensitive else "MESH_STABLE_AT_TESTED_LEVELS",
        "primary_comparison_uses_robust_metrics": True,
    }


def make_stage_c_handoff(descriptor: Mapping[str, Any]) -> dict[str, Any]:
    compatibility = GeometryCompatibility(descriptor["compatibility_class"])
    if not route_geometry(compatibility)["stage_c_allowed"]:
        return {"status": "NOT_APPLICABLE", "route_to_stage_e_r2": compatibility == GeometryCompatibility.NONAXISYMMETRIC_3D_REQUIRED}
    geometry = descriptor["geometry_identity"]
    return {
        "status": "READY" if compatibility == GeometryCompatibility.AXISYMMETRIC_2D_EXACT else "APPROXIMATE_REQUIRES_JUSTIFICATION",
        "geometry_id": geometry["geometry_id"],
        "compatibility_class": compatibility.value,
        "gap_m": geometry["gap_m"],
        "tip_radius_m": geometry["tip_radius_m"],
        "tip_profile": geometry.get("tip_profile"),
        "boundary_locations_m": geometry.get("boundary_locations_m"),
        "reference_voltage_V": descriptor["field_enhancement"]["E_nominal_V_m"] * geometry["gap_m"],
        "stage_c_requirement": "reconstruct geometry and solve zero-charge electrostatics",
        "forbidden_use": "do not replace self-consistent Stage-C Poisson with a fixed COMSOL E0 map",
    }


def cross_validate_stage_c(stage_b: Mapping[str, Any], stage_c: Mapping[str, Any],
                           stage_b_axis_field: np.ndarray, stage_c_axis_field: np.ndarray) -> dict[str, Any]:
    if stage_b_axis_field.shape != stage_c_axis_field.shape or stage_b_axis_field.size < 3:
        raise ValueError("axis profiles must have the same shape and at least three samples")
    correlation = float(np.corrcoef(stage_b_axis_field, stage_c_axis_field)[0, 1])
    def rel(left: float, right: float) -> float:
        return abs(left - right) / max(abs(left), 1e-300)
    emax = rel(stage_b["roi_statistics"]["E0_max_V_m"], stage_c["roi_statistics"]["E0_max_V_m"])
    ep95 = rel(stage_b["roi_statistics"]["E0_p95_V_m"], stage_c["roi_statistics"]["E0_p95_V_m"])
    beta = rel(stage_b["field_enhancement"]["beta_E_p95"], stage_c["field_enhancement"]["beta_E_p95"])
    measure_key = "volume_m3" if "volume_m3" in stage_b["high_field_regions"]["q90"] else "area_m2"
    high_field = rel(
        stage_b["high_field_regions"]["q90"][measure_key], stage_c["high_field_regions"]["q90"][measure_key]
    )
    passed = ep95 < 0.05 and beta < 0.05 and correlation >= 0.99
    return {
        "status": "PASS" if passed else "NOT_RESOLVED",
        "E0_max_relative_difference": emax,
        "E0_p95_relative_difference": ep95,
        "beta_E_p95_relative_difference": beta,
        "high_field_q90_measure_relative_difference": high_field,
        "axis_profile_correlation": correlation,
        "initial_acceptance_guidance": {"robust_metric_relative_difference_max": 0.05, "axis_profile_correlation_min": 0.99},
        "hard_parameter_tuning_to_match_forbidden": True,
    }


def paper1_descriptor(descriptor: Mapping[str, Any], *, material_interface_id: str) -> dict[str, Any]:
    identity = descriptor["geometry_identity"]
    stats = descriptor["roi_statistics"]
    return {
        "schema_version": "1.0",
        "descriptor_type": "PAPER1_ELECTRODE_DESCRIPTOR",
        "dataset_role": descriptor["dataset_role"],
        "scientific_result_allowed": descriptor["scientific_result_allowed"],
        "geometry_id": identity["geometry_id"],
        "gap_m": identity["gap_m"],
        "tip_radius_m": identity["tip_radius_m"],
        "electrode_material": identity["electrode_material"],
        "material_interface_id": material_interface_id,
        "E0_max_V_m": stats["E0_max_V_m"],
        "E0_p95_V_m": stats["E0_p95_V_m"],
        "beta_E": descriptor["field_enhancement"]["beta_E"],
        "beta_E_p95": descriptor["field_enhancement"]["beta_E_p95"],
        "eta_0_p95": descriptor["ek_normalization"]["eta_0_p95"],
        "high_field_regions": descriptor["high_field_regions"],
        "field_localization": descriptor["field_localization"],
        "Cgap": descriptor["Cgap"],
        "mesh_status": descriptor["mesh"].get("status"),
        "compatibility_class": descriptor["compatibility_class"],
        "downstream_links": {"streamer_metrics": None, "rf_metrics": None, "experimental_metrics": None},
    }


def f_r6b_gate(descriptors: list[Mapping[str, Any]]) -> bool:
    actual = [item for item in descriptors if item.get("dataset_role") == "ACTUAL_STAGE_B_EXPORT" and item.get("scientific_result_allowed")]
    compatible = any(item.get("compatibility_class") in {
        GeometryCompatibility.AXISYMMETRIC_2D_EXACT.value,
        GeometryCompatibility.AXISYMMETRIC_2D_APPROXIMATE.value,
    } for item in actual)
    return len(actual) >= 2 and compatible
