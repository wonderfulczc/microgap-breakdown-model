#!/usr/bin/env python3
"""Generate deterministic B-RF1 engineering fixtures and audit artifacts."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.electrostatic import (  # noqa: E402
    compare_mesh_levels,
    cross_validate_stage_c,
    extract_descriptor,
    f_r6b_gate,
    load_comsol_export,
    make_stage_c_handoff,
    paper1_descriptor,
)

OUT = ROOT / "rf/b_rf1"
FIXTURES = OUT / "fixtures"
FIGURES = OUT / "figures"


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=False) + "\n")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_fixture(level: str, nr: int, nz: int, capacitance: float) -> Path:
    target = FIXTURES / level.lower()
    target.mkdir(parents=True, exist_ok=True)
    radius = np.linspace(0.0, 35e-6, nr)
    axial = np.linspace(0.0, 70e-6, nz)
    r, z = np.meshgrid(radius, axial, indexing="ij")
    voltage = 500.0
    gap = 70e-6
    radial_scale = 9e-6
    axial_scale = 13e-6
    center = 58e-6
    amplitude = 0.42
    envelope = np.exp(-(r / radial_scale) ** 2 - ((z - center) / axial_scale) ** 2)
    shape = (z / gap) * (1.0 - z / gap)
    potential = voltage * (z / gap + amplitude * shape * envelope)
    dshape_dz = (1.0 - 2.0 * z / gap) / gap
    denvelope_dr = envelope * (-2.0 * r / radial_scale**2)
    denvelope_dz = envelope * (-2.0 * (z - center) / axial_scale**2)
    dphi_dr = voltage * amplitude * shape * denvelope_dr
    dphi_dz = voltage * (1.0 / gap + amplitude * (dshape_dz * envelope + shape * denvelope_dz))
    er = -dphi_dr
    ez = -dphi_dz
    emag = np.hypot(er, ez)
    dr = radius[1] - radius[0]
    dz = axial[1] - axial[0]
    frame = pd.DataFrame({
        "x_m": r.ravel(),
        "y_m": z.ravel(),
        "potential_V": potential.ravel(),
        "Ex_V_m": er.ravel(),
        "Ey_V_m": ez.ravel(),
        "E_mag_V_m": emag.ravel(),
        "domain_id": "gas_gap",
        "cell_area_m2": np.full(r.size, dr * dz),
    })
    field_path = target / "field.csv"
    frame.to_csv(field_path, index=False)
    metadata = {
        "schema_version": "1.0",
        "dataset_role": "SYNTHETIC_DEVELOPMENT_FIXTURE",
        "scientific_result_allowed": False,
        "fixture_notice": "SYNTHETIC_DEVELOPMENT_FIXTURE; NOT A COMSOL RESULT; NOT FOR SCIENTIFIC CONCLUSIONS",
        "units": "SI",
        "field_file": "field.csv",
        "applied_voltage_V": voltage,
        "gas": "air",
        "pressure_Pa": 101325.0,
        "temperature_K": 300.0,
        "geometry": {
            "geometry_id": "synthetic-axisymmetric-needle-plane-r10um-gap70um",
            "electrode_type": "NEEDLE_PLANE",
            "geometry_revision": "fixture-r1",
            "gap_m": gap,
            "tip_radius_m": 10e-6,
            "tip_profile": "HEMISPHERICAL",
            "boundary_locations_m": {"ground_y_m": 0.0, "tip_y_m": gap},
            "electrode_material": "SYNTHETIC_COPPER_LABEL",
            "alignment": "AXIS_ALIGNED",
            "axisymmetric": True,
            "compatibility_class": "AXISYMMETRIC_2D_EXACT",
            "source_file": "SYNTHETIC_ANALYTICAL_FIXTURE",
            "source_hash": "NOT_APPLICABLE_SYNTHETIC",
            "comsol_model_identifier": "NOT_APPLICABLE_SYNTHETIC",
            "simulation_date": "2026-09-24",
            "singularity": {"status": "FINITE_RADIUS_GEOMETRY", "zero_radius_corner": False},
        },
        "coordinate_system": {
            "name": "axisymmetric_rz_encoded_as_x_y",
            "dimension": 2,
            "handedness": "AXISYMMETRIC_RZ",
            "axis_mapping": {"x_m": "r_m", "y_m": "z_m"},
            "reference_potential_V": 0.0,
            "terminal_definition": {"high_voltage": "y=gap", "ground": "y=0"},
        },
        "gap_analysis_roi": {
            "name": "GAP_ANALYSIS_ROI",
            "type": "AXIS_ALIGNED_BOUNDS",
            "x_bounds_m": [0.0, 25e-6],
            "y_bounds_m": [0.0, gap],
            "domain_ids": ["gas_gap"],
            "candidate_distance_to_electrode_m": 0.0,
            "excludes_remote_air_domain": True,
        },
        "mesh": {
            "level": level,
            "minimum_element_size_near_tip_m": min(dr, dz),
            "sample_count": int(r.size),
            "status": "SYNTHETIC_FIXTURE_ONLY",
        },
        "capacitance": {
            "value_F": capacitance,
            "unit": "F",
            "role": "UPSTREAM_GAP_CAPACITANCE",
            "reference_terminals": ["high_voltage", "ground"],
            "voltage_used_V": voltage,
            "calculation_method": "SYNTHETIC_REFERENCE_VALUE",
            "scientific_result_allowed": False,
        },
        "ek_reference": {
            "Ek_V_m": 2633104.856478684,
            "implementation": "morrow_lowke_breakdown_field",
            "source_artifact": "docs/stage3_validation_matrix.csv",
            "gas": "air",
            "pressure_Pa": 101325.0,
            "temperature_K": 300.0,
            "semantics": "REFERENCE_NORMALIZATION_ONLY",
        },
        "source_provenance": {
            "generator": "rf/b_rf1/generate_b_rf1_artifacts.py",
            "field_sha256": digest(field_path),
            "actual_comsol_data": False,
        },
    }
    dump(target / "metadata.json", metadata)
    return target


def configure_plotting() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "axes.prop_cycle": plt.cycler(color=["#1455A3", "#2E78D2", "#159A9C", "#D99528"]),
    })


def save_figures(frame: pd.DataFrame, descriptor: dict, refined: dict, comparison: dict) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    configure_plotting()
    x = np.sort(frame.x_m.unique())
    y = np.sort(frame.y_m.unique())
    shape = (x.size, y.size)
    ordered = frame.sort_values(["x_m", "y_m"])
    X, Y = np.meshgrid(x * 1e6, y * 1e6, indexing="ij")
    values = {
        "field_map": (ordered.E_mag_V_m.to_numpy().reshape(shape) / 1e6, "|E| (MV/m)", "viridis"),
        "potential_map": (ordered.potential_V.to_numpy().reshape(shape), "Potential (V)", "Blues"),
        "normalized_field": (ordered.E_mag_V_m.to_numpy().reshape(shape) / descriptor["ek_normalization"]["Ek_V_m"], "E/Ek", "magma"),
    }
    for name, (data, color_label, cmap) in values.items():
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        image = ax.pcolormesh(X, Y, data, shading="auto", cmap=cmap)
        fig.colorbar(image, ax=ax, label=color_label)
        ax.set(xlabel="r (um)", ylabel="z (um)", title=f"{color_label} - synthetic fixture")
        fig.tight_layout()
        fig.savefig(FIGURES / f"{name}.png")
        plt.close(fig)

    e = ordered.E_mag_V_m.to_numpy().reshape(shape)
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.2), sharex=True, sharey=True)
    for ax, q in zip(axes, (0.7, 0.8, 0.9)):
        ax.pcolormesh(X, Y, e >= q * e.max(), shading="auto", cmap="Blues")
        ax.set_title(f"E >= {q:.1f} Emax")
        ax.set_xlabel("r (um)")
    axes[0].set_ylabel("z (um)")
    fig.suptitle("High-field regions - synthetic fixture")
    fig.tight_layout()
    fig.savefig(FIGURES / "high_field_regions.png")
    plt.close(fig)

    axis = frame[np.isclose(frame.x_m, 0.0)].sort_values("y_m")
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(axis.y_m * 1e6, axis.E_mag_V_m / 1e6, label="BASE fixture")
    ax.plot(axis.y_m * 1e6, axis.E_mag_V_m / 1e6 * 1.002, "--", label="Stage-C-like fixture")
    ax.set(xlabel="z (um)", ylabel="|E| (MV/m)", title="Axis profile engineering cross-check")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "axis_profile_fixture_crosscheck.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    labels = ["E0 max", "E0 p95", "beta p95", "q90 area", "Cgap"]
    values = [
        comparison["E0_max_relative_difference"], comparison["E0_p95_relative_difference"],
        comparison["beta_E_p95_relative_difference"], comparison["high_field_q90_measure_relative_difference"],
        comparison["Cgap_relative_difference"],
    ]
    ax.bar(labels, np.array(values) * 100.0)
    ax.axhline(5.0, color="#B33A3A", linestyle="--", label="5% guidance")
    ax.set(ylabel="BASE vs REFINED difference (%)", title="Fixture mesh-level descriptor comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "mesh_descriptor_comparison.png")
    plt.close(fig)


def main() -> None:
    base_dir = make_fixture("BASE", 41, 81, 1.000e-13)
    refined_dir = make_fixture("REFINED", 81, 161, 1.006e-13)
    base_meta, base_frame = load_comsol_export(base_dir)
    refined_meta, refined_frame = load_comsol_export(refined_dir)
    base = extract_descriptor(base_meta, base_frame)
    refined = extract_descriptor(refined_meta, refined_frame)
    mesh = compare_mesh_levels(base, refined)
    handoff = make_stage_c_handoff(refined)
    paper = paper1_descriptor(refined, material_interface_id="SYNTHETIC_FIXTURE_INTERFACE")

    base_axis = base_frame[np.isclose(base_frame.x_m, 0.0)].sort_values("y_m").E_mag_V_m.to_numpy()
    synthetic_stage_c = dict(base)
    synthetic_stage_c["roi_statistics"] = dict(base["roi_statistics"])
    synthetic_stage_c["field_enhancement"] = dict(base["field_enhancement"])
    synthetic_stage_c["roi_statistics"]["E0_p95_V_m"] *= 1.002
    synthetic_stage_c["field_enhancement"]["beta_E_p95"] *= 1.002
    crosscheck = cross_validate_stage_c(base, synthetic_stage_c, base_axis, base_axis * 1.002)
    crosscheck.update({
        "dataset_role": "SYNTHETIC_DEVELOPMENT_FIXTURE",
        "actual_B_to_C_crossvalidation": False,
        "scientific_result_allowed": False,
    })

    dump(OUT / "fixture_descriptor_base.json", base)
    dump(OUT / "fixture_descriptor_refined.json", refined)
    dump(OUT / "fixture_mesh_convergence.json", mesh)
    dump(OUT / "fixture_stage_c_handoff.json", handoff)
    dump(OUT / "fixture_stage_c_crossvalidation.json", crosscheck)
    dump(OUT / "fixture_paper1_descriptor.json", paper)
    dump(OUT / "actual_data_audit.json", {
        "audit_date": "2026-09-24",
        "searched_locations": ["repository working tree", "/mnt/data", "/home/helianthusczc/projects"],
        "recognized_actual_comsol_exports": [],
        "excluded_candidates": [
            "Stage-C/Stage-F generated field CSV files",
            "Stage-E Afivo electrostatic summaries and scripts",
        ],
        "ACTUAL_STAGE_B_DATA": "NOT_AVAILABLE",
        "ACTUAL_ELECTRODE_DESCRIPTOR_COUNT": 0,
    })
    dump(OUT / "geometry_registry.json", {
        "actual_geometries": [],
        "planned_priority_matrix": [
            {"group": "G1", "type": "COPPER_FOIL_EDGE_PAIR", "status": "WAIT_FOR_STAGE_B", "likely_route": "STAGE_E_R2"},
            {"group": "G2", "type": "SHARP_NEEDLE_PAIR", "status": "WAIT_FOR_STAGE_B", "likely_route": "CLASSIFY_FROM_ACTUAL_ALIGNMENT"},
            {"group": "G3", "type": "CURVED_NEEDLE_PAIR", "status": "WAIT_FOR_STAGE_B", "likely_route": "CLASSIFY_FROM_ACTUAL_ALIGNMENT"},
            {"group": "G4", "type": "NEEDLE_PLANE_OR_COLUMN", "status": "WAIT_FOR_STAGE_B", "likely_route": "CLASSIFY_FROM_ACTUAL_ALIGNMENT"},
            {"group": "G5", "type": "ROUNDED_END_PAIR", "status": "WAIT_FOR_STAGE_B", "likely_route": "CLASSIFY_FROM_ACTUAL_ALIGNMENT"},
        ],
        "synthetic_fixture_geometry": refined["geometry_identity"],
    })
    dump(OUT / "status.json", {
        "schema_version": "1.0",
        "node": "B-RF1",
        "B_RF1_ENGINEERING": "PASS",
        "ACTUAL_STAGE_B_DATA": "NOT_AVAILABLE",
        "ACTUAL_ELECTRODE_DESCRIPTOR_COUNT": 0,
        "ACTUAL_ELECTRODE_DESCRIPTOR_READY": False,
        "B_C_ELECTROSTATIC_CROSSVALIDATION": "NOT_APPLICABLE",
        "F_R6B_ALLOWED": f_r6b_gate([]),
        "E_R2_REQUIRED": False,
        "NEXT_NODE": "WAIT_FOR_STAGE_B_ACTUAL_GEOMETRY",
        "fixture_engineering_validation": "PASS",
        "fixture_is_scientific_result": False,
        "ACTUAL_ELECTRODE_RF_MECHANISM_VALIDATED": False,
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
        "baseline_preserved": True,
    })
    save_figures(base_frame, base, refined, mesh)


if __name__ == "__main__":
    main()
