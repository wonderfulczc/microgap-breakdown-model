"""Strict, read-only FP-1 manifest and COMSOL CSV import."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


class FP1Error(ValueError):
    def __init__(self, code: str, detail: str):
        super().__init__(f"{code}: {detail}")
        self.code = code


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def package_path(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not name or name == "N/A":
        raise FP1Error("MISSING_FILE", str(name))
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise FP1Error("PATH_OUTSIDE_PACKAGE", name)
    if not path.is_file():
        raise FP1Error("MISSING_FILE", name)
    return path


def load_manifest(path: str | Path) -> tuple[dict, Path]:
    path = Path(path)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    schema = json.loads((ROOT / "manifest_schema.json").read_text(encoding="utf-8"))
    missing = set(schema["required"]) - set(manifest)
    if missing:
        raise FP1Error("MISSING_MANIFEST_FIELDS", ", ".join(sorted(missing)))
    unknown = set(manifest) - set(schema["properties"])
    if unknown:
        raise FP1Error("UNKNOWN_MANIFEST_FIELDS", ", ".join(sorted(unknown)))
    for key in ("schema_version", "fp1_manifest_version", "dataset_role", "model_dimension", "coordinate_basis", "coordinate_unit", "normal_convention", "return_reference"):
        rule = schema["properties"][key]
        if "const" in rule and manifest[key] != rule["const"]:
            raise FP1Error("INVALID_MANIFEST_VALUE", key)
        if "enum" in rule and manifest[key] not in rule["enum"]:
            raise FP1Error("INVALID_MANIFEST_VALUE", key)
    if manifest["model_dimension"] == "3D":
        raise FP1Error("UNSUPPORTED_3D_BACKEND", "schema reserved; no validated 3D mapper")
    if manifest["coordinate_basis"] != "r,z":
        raise FP1Error("INVALID_COORDINATES", "2D_axi requires r,z")
    if manifest["unit_voltage_V"] != 1:
        raise FP1Error("INVALID_REFERENCE_VOLTAGE", "FULL static reference must be 1 V")
    if manifest["dataset_role"] == "REAL_COMSOL":
        for key in ("geometry_source", "geometry_source_sha256", "static_mesh_sensitivity_file", "domain_truncation_file"):
            if not manifest.get(key) or manifest[key] == "N/A":
                raise FP1Error("REAL_COMSOL_PROVENANCE_PENDING", key)
    for key in ("geometry_id", "geometry_family", "geometry_model_version", "static_run_id", "boundary_model_id", "port_reference_plane_id", "material_map_id", "outer_boundary_type"):
        if not isinstance(manifest[key], str) or not manifest[key].strip() or manifest[key] == "N/A":
            raise FP1Error("INVALID_MANIFEST_VALUE", key)
    for key in ("gas_region_id", "hv_surface_id", "return_surface_id", "outer_boundary_id", "axis_boundary_id"):
        values = manifest[key]
        if not isinstance(values, list) or (key in ("gas_region_id", "hv_surface_id", "return_surface_id") and not values) or any(type(x) is not int for x in values) or len(set(values)) != len(values):
            raise FP1Error("INVALID_REGION_LABELS", key)
    groups = [set(manifest[key]) for key in ("hv_surface_id", "return_surface_id", "outer_boundary_id", "axis_boundary_id")]
    if any(groups[i] & groups[j] for i in range(4) for j in range(i + 1, 4)):
        raise FP1Error("OVERLAPPING_BOUNDARY_LABELS", "terminal/outer/axis sets must be disjoint")
    for key in ("d_net_m", "C_geom_F"):
        value = manifest[key]
        if type(value) not in (float, int) or not math.isfinite(value) or value <= 0:
            raise FP1Error("INVALID_MANIFEST_VALUE", key)
    if ("epsilon_r_uniform" in manifest) == ("epsilon_field_file" in manifest):
        raise FP1Error("EPSILON_SOURCE_MISMATCH", "specify exactly one uniform epsilon_r or field file")
    if "epsilon_r_uniform" in manifest and (type(manifest["epsilon_r_uniform"]) not in (float, int) or not math.isfinite(manifest["epsilon_r_uniform"]) or manifest["epsilon_r_uniform"] <= 0):
        raise FP1Error("INVALID_EPSILON", "relative permittivity must be positive")
    for key in ("Q_HV_at_1V_C", "U_at_1V_J"):
        if type(manifest[key]) not in (float, int) or not math.isfinite(manifest[key]):
            raise FP1Error("INVALID_MANIFEST_VALUE", key)
    witness = manifest["gap_witness"]
    if not isinstance(witness, dict) or witness.get("gas_connected") is not True or not witness.get("method"):
        raise FP1Error("GAP_WITNESS_PENDING", "facing surfaces and gas connection required")
    try:
        ends = np.asarray([witness["hv_point_m"], witness["return_point_m"]], dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise FP1Error("GAP_WITNESS_PENDING", str(exc)) from exc
    if ends.shape != (2, 2) or not np.isfinite(ends).all() or (ends[:, 0] < 0).any():
        raise FP1Error("INVALID_GAP_WITNESS", "finite r,z points with r>=0 required")
    measured = float(np.linalg.norm(ends[0] - ends[1]))
    if not math.isclose(measured, manifest["d_net_m"], rel_tol=1e-6, abs_tol=1e-12):
        raise FP1Error("GAP_DEFINITION_MISMATCH", f"witness={measured}, declared={manifest['d_net_m']}")
    for key in ("static_solution_provenance", "mesh_mapping", "uncertainty"):
        rule = schema["properties"][key]
        if not isinstance(manifest[key], dict) or set(rule["required"]) - set(manifest[key]):
            raise FP1Error("MISSING_MANIFEST_FIELDS", key)
    if manifest["static_solution_provenance"]["zero_free_charge"] is not True:
        raise FP1Error("NONSTATIC_SOURCE", "weighting reference requires rho=0")
    if manifest["mesh_mapping"]["policy"] not in {"SAME_MESH", "LINEAR_NO_EXTRAPOLATION"}:
        raise FP1Error("INVALID_MAPPING_POLICY", "unsupported policy")
    root = path.parent.resolve()
    for key in ("psi_field_file", "Ew_field_file", "domain_mask_file", "boundary_faces_file", "mesh_metadata_file"):
        package_path(root, manifest[key])
    for key in ("geometry_source", "static_mesh_sensitivity_file", "domain_truncation_file", "epsilon_field_file", "target_poisson_links_file"):
        if manifest.get(key):
            package_path(root, manifest[key])
    hashed = manifest["file_sha256"]
    for key in ("psi_field_file", "Ew_field_file", "domain_mask_file", "boundary_faces_file", "mesh_metadata_file"):
        name = manifest[key]
        if hashed.get(name) != sha256(package_path(root, name)):
            raise FP1Error("FILE_HASH_MISMATCH", name)
    for key in ("epsilon_field_file", "target_poisson_links_file", "static_mesh_sensitivity_file", "domain_truncation_file"):
        if manifest.get(key) and hashed.get(manifest[key]) != sha256(package_path(root, manifest[key])):
            raise FP1Error("FILE_HASH_MISMATCH", manifest[key])
    if manifest.get("geometry_source") and manifest.get("geometry_source_sha256") != sha256(package_path(root, manifest["geometry_source"])):
        raise FP1Error("SOURCE_HASH_MISMATCH", manifest["geometry_source"])
    return manifest, root


def load_field(path: Path, kind: str) -> dict[str, np.ndarray]:
    """Read normalized r,z,psi or r,z,Ewr,Ewz CSV, including COMSOL % headers."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = csv.reader(stream)
        header = None
        metadata = {}
        raw_comsol = False
        data = []
        for row in rows:
            if not row:
                continue
            if header is None:
                if row[0].startswith("%"):
                    if row[0].strip() == "% r" and len(row) > 1 and row[1].strip() == "z":
                        header = [row[0][2:].strip(), *[x.strip() for x in row[1:]]]
                        raw_comsol = True
                    elif len(row) > 1:
                        metadata[row[0][1:].strip()] = row[1].strip()
                    continue
                header = [x.strip() for x in row]
                continue
            data.append(row)
    if header is None or not data:
        raise FP1Error("EMPTY_FIELD", str(path))
    names = ("r_m", "z_m", "psi") if kind == "psi" else ("r_m", "z_m", "Ewr_per_m", "Ewz_per_m")
    if set(names) <= set(header):
        indices = [header.index(x) for x in names]
    elif kind == "psi" and {"r (m)", "z (m)", "V/(1[V]) (1)"} <= set(header):
        indices = [header.index(x) for x in ("r (m)", "z (m)", "V/(1[V]) (1)")]
    elif kind == "Ew" and {"es.Er/(1[V]) (1/m)", "es.Ez/(1[V]) (1/m)"} <= set(header):
        if header[:2] != ["r", "z"]:
            raise FP1Error("INVALID_FIELD_HEADER", "COMSOL coordinate order unclear")
        indices = [0, 1, header.index("es.Er/(1[V]) (1/m)"), header.index("es.Ez/(1[V]) (1/m)")]
    else:
        raise FP1Error("INVALID_FIELD_HEADER", f"{kind}: {header}")
    try:
        values = np.asarray([[float(row[i]) for i in indices] for row in data], dtype=float)
    except (IndexError, ValueError) as exc:
        raise FP1Error("INVALID_FIELD_VALUES", str(path)) from exc
    if not np.isfinite(values).all():
        raise FP1Error("INVALID_FIELD_VALUES", "nonfinite")
    if raw_comsol:
        unit = metadata.get("Length unit")
        if unit not in {"m", "mm", "um"}:
            raise FP1Error("INVALID_COORDINATE_UNIT", f"COMSOL length unit: {unit}")
        factor = {"m": 1, "mm": 1e-3, "um": 1e-6}[unit]
        if kind == "psi":
            native = np.asarray([[float(row[0]), float(row[1])] for row in data]) * factor
            if not np.allclose(native, values[:, :2], rtol=0, atol=1e-12):
                raise FP1Error("COORDINATE_UNIT_MISMATCH", "native and SI coordinate columns differ")
        else:
            values[:, :2] *= factor
            if "es.normE/(1[V]) (1/m)" in header:
                magnitude = np.asarray([float(row[header.index("es.normE/(1[V]) (1/m)")]) for row in data])
                if not np.allclose(np.linalg.norm(values[:, 2:], axis=1), magnitude, rtol=1e-5, atol=1e-8):
                    raise FP1Error("FIELD_COMPONENT_MISMATCH", "COMSOL Ew magnitude")
        if "Nodes" in metadata and int(metadata["Nodes"]) != len(values):
            raise FP1Error("NODE_COUNT_MISMATCH", str(path))
    if (values[:, 0] < -1e-14).any():
        raise FP1Error("INVALID_COORDINATES", "negative radius")
    if kind == "psi" and ((values[:, 2] < -1e-5).any() or (values[:, 2] > 1 + 1e-5).any()):
        raise FP1Error("INVALID_PSI", "1 V/0 V weighting potential outside expected range")
    return {**{name: values[:, i] for i, name in enumerate(names)}, "source_metadata": metadata}


def load_package(path: str | Path) -> dict:
    manifest, root = load_manifest(path)
    psi = load_field(package_path(root, manifest["psi_field_file"]), "psi")
    field = load_field(package_path(root, manifest["Ew_field_file"]), "Ew")
    if len(psi["r_m"]) != len(field["r_m"]) or not np.allclose(np.column_stack((psi["r_m"], psi["z_m"])), np.column_stack((field["r_m"], field["z_m"])), rtol=0, atol=1e-12):
        raise FP1Error("FIELD_COORDINATE_MISMATCH", "psi and Ew must share physical sample coordinates")
    for key in ("Model", "Version"):
        if psi["source_metadata"].get(key) != field["source_metadata"].get(key):
            raise FP1Error("SOURCE_IDENTITY_MISMATCH", key)
    if manifest["dataset_role"] == "REAL_COMSOL":
        provenance = manifest["static_solution_provenance"]
        if psi["source_metadata"].get("Model") != provenance["model_name"] or psi["source_metadata"].get("Version") != provenance["comsol_version"]:
            raise FP1Error("SOURCE_IDENTITY_MISMATCH", "COMSOL model/version disagree with manifest")
    mesh_data = None
    for key in ("domain_mask_file", "mesh_metadata_file"):
        data = json.loads(package_path(root, manifest[key]).read_text(encoding="utf-8"))
        for identity in ("geometry_id", "static_run_id", "boundary_model_id", "port_reference_plane_id"):
            if data.get(identity) != manifest[identity]:
                raise FP1Error("IDENTITY_MISMATCH", f"{key}: {identity}")
        if key == "domain_mask_file":
            codes = data.get("region_codes", {})
            if any(str(region) not in codes for region in manifest["gas_region_id"]):
                raise FP1Error("REGION_LABEL_MISMATCH", "gas region absent from mask")
        else:
            mesh_data = data
    if manifest["mesh_mapping"]["policy"] == "SAME_MESH":
        node_block = mesh_data.get("nodes", {})
        cols = node_block.get("columns", [])
        if not {"r_m", "z_m"} <= set(cols):
            raise FP1Error("INVALID_MESH", "node coordinates missing")
        node_coordinates = np.asarray(node_block["rows"], dtype=float)[:, [cols.index("r_m"), cols.index("z_m")]]
        sample_coordinates = np.column_stack((psi["r_m"], psi["z_m"]))
        if node_coordinates.shape != sample_coordinates.shape or not np.allclose(node_coordinates, sample_coordinates, rtol=0, atol=1e-12):
            raise FP1Error("MESH_IDENTITY_MISMATCH", "SAME_MESH requires matching node coordinates/order")
    return {"manifest": manifest, "root": root, "psi": psi, "field": field}
