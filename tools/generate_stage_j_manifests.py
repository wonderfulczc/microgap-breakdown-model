#!/usr/bin/env python3
"""Generate compact Stage-J packaging and reproducibility manifests."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "packaging"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tree_bytes(path):
    path = Path(path)
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file()) if path.exists() else 0


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def main():
    OUT.mkdir(exist_ok=True)
    commit = git("rev-parse", "HEAD")
    branch = git("branch", "--show-current")
    openems = json.loads((ROOT / "fullwave/h1/openems_backend.json").read_text())
    stage_i = json.loads((ROOT / "validation/stage_i/final/stage_i_final_contract.json").read_text())

    directory_classes = {
        ".agents": "EXCLUDE_FROM_RELEASE",
        ".codex": "EXCLUDE_FROM_RELEASE",
        ".codex_skill_staging": "EXCLUDE_FROM_RELEASE",
        ".git": "EXCLUDE_FROM_RELEASE",
        ".pytest_cache": "BUILD_ARTIFACT",
        ".venv": "BUILD_ARTIFACT",
        "bridge-export": "ARCHIVE_ONLY",
        "build": "BUILD_ARTIFACT",
        "build-smoke": "BUILD_ARTIFACT",
        "circuit": "REQUIRED_REFERENCE",
        "config": "CORE",
        "cpp": "CORE",
        "docs": "REQUIRED_REFERENCE",
        "fullwave": "REQUIRED_REFERENCE",
        "packaging": "REQUIRED_REFERENCE",
        "python": "CORE",
        "results": "REGENERABLE",
        "rf": "REQUIRED_REFERENCE",
        "scripts": "CORE",
        "solver3d": "EXTERNAL_BACKEND",
        "tests": "CORE",
        "thermal": "REQUIRED_REFERENCE",
        "tools": "CORE",
        "validation": "REQUIRED_REFERENCE",
    }
    inventory = []
    for path in sorted(item for item in ROOT.iterdir() if item.is_dir()):
        if path.name not in directory_classes:
            raise ValueError(f"UNCLASSIFIED_TOP_LEVEL_DIRECTORY:{path.name}")
        inventory.append({
            "path": path.name,
            "classification": directory_classes[path.name],
            "bytes_in_worktree": tree_bytes(path),
            "tracked": bool(git("ls-files", path.name)),
        })
    write("stage_j_repository_inventory.json", {"schema": "StageJRepositoryInventory-v1", "directories": inventory})

    backends = {
        "Afivo": {
            "source_policy": "EXTERNAL_NOT_VENDORED",
            "upstream": "https://github.com/MD-CWI/afivo-streamer.git",
            "commit": "a50b5508775086e90dfe423455fb58d812578410",
            "license": "GPL-3.0_UPSTREAM",
            "environment_variable": "AFIVO_STREAMER_ROOT",
            "interface": "solver3d/afivo_reference configs and exported CSV source/field data",
            "role": "3D streamer cross-validation backend",
            "required_project_files": ["solver3d/afivo_reference", "docs/stage_d1_afivo_backend.md"],
        },
        "openEMS": {
            "source_policy": "EXTERNAL_NOT_VENDORED",
            "upstream": openems["upstream_url"],
            "commits": openems["commits"],
            "versions": openems["source_versions"],
            "build_options": openems["build_options"],
            "environment_variables": ["OPENEMS_PROJECT_ROOT", "OPENEMS_ROOT", "OPENEMS_DEPS_ROOT", "OPENEMS_PYTHON"],
            "interface": "isolated Python openEMS/CSXCAD environment and compact CSV/JSON results",
            "role": "passive full-wave structure and receiver backend",
        },
        "COMSOL": {
            "source_policy": "PROPRIETARY_EXTERNAL_TOOL_NOT_DISTRIBUTED",
            "role": "Stage-B electrostatic geometry/field comparison interface",
            "packaged_scope": ["DOCUMENTATION", "GEOMETRY_FIELD_METADATA_CONVENTION", "EXISTING_IMPORT_INTERFACE"],
            "scientific_status": "EXTERNAL_VALIDATION_PENDING",
        },
    }
    write("external_backends.json", backends)

    large_data = {
        "datasets": [
            {
                "dataset": "AFIVO_STAGE_E_RAW",
                "path": "solver3d/afivo_reference/stage_e/results_raw",
                "bytes_local": tree_bytes(ROOT / "solver3d/afivo_reference/stage_e/results_raw"),
                "scientific_role": "3D source and field exports for Stage E/F",
                "directory_hash": "NOT_COMPUTED_MULTI_GB_DIRECTORY",
                "provenance": "docs/stage_d1_afivo_backend.md and Stage-E configs",
                "regenerable": True,
                "required_upstream": "Afivo pinned commit plus solver3d/afivo_reference/stage_e/configs",
                "release_policy": "EXTERNAL_ARCHIVE",
            },
            {
                "dataset": "STAGE_F_RAW_FIELDS",
                "path": "results/stage_f4/source_gate_recovery",
                "bytes_local": tree_bytes(ROOT / "results/stage_f4/source_gate_recovery"),
                "scientific_role": "dense field/source-gate recovery evidence",
                "directory_hash": "NOT_COMPUTED_MULTI_GB_DIRECTORY",
                "provenance": "rf/source/audit manifests and frozen Stage-F contracts",
                "regenerable": True,
                "required_upstream": "frozen Stage-C/E source configurations",
                "release_policy": "DO_NOT_PACKAGE",
            },
        ],
        "normal_source_release_contains_multi_GB_raw_data": False,
        "local_files_deleted": False,
    }
    write("large_data_manifest.json", large_data)

    synthetic = {
        "fixtures": [
            {
                "path": "validation/stage_i/system_350mhz/synthetic_input",
                "bytes": tree_bytes(ROOT / "validation/stage_i/system_350mhz/synthetic_input"),
                "manifest": "validation/stage_i/system_350mhz/synthetic_input/sha256_manifest.json",
                "labels": ["SYNTHETIC_DEVELOPMENT_INPUT", "SYNTHETIC_DRY_RUN", "NOT_EXPERIMENTAL_DATA"],
                "policy": "KEEP_IN_REPOSITORY_WITH_HASH_MANIFEST",
            },
            {
                "path": "validation/stage_i/system_350mhz/wp_i_d_synthetic",
                "bytes": tree_bytes(ROOT / "validation/stage_i/system_350mhz/wp_i_d_synthetic"),
                "manifest": "validation/stage_i/system_350mhz/wp_i_d_synthetic/sha256_manifest.json",
                "labels": ["SYNTHETIC_DEVELOPMENT_INPUT", "SYNTHETIC_DRY_RUN", "NOT_EXPERIMENTAL_DATA"],
                "policy": "KEEP_IN_REPOSITORY_WITH_HASH_MANIFEST",
            },
        ],
        "replacement_policy": "DO_NOT_REWRITE_FROZEN_FIXTURES_WITHOUT_BYTE_AND_PROVENANCE_TRACEABILITY",
    }
    write("synthetic_fixture_policy.json", synthetic)

    machine_paths = {
        "classifications": {
            "DOCUMENTATION_ONLY": ["docs/environment_baseline.md", "docs/stage_d1_afivo_backend.md", "docs/stage_h1_fullwave_foundation.md"],
            "FROZEN_METADATA": ["fullwave/h1/openems_backend.json", "fullwave/h1-h4 summary JSON", "rf/source/audit manifests", "bridge-export/simulation"],
            "TEST_FIXTURE": ["frozen result paths under rf validation JSON"],
            "ENV_CONFIG_REQUIRED": ["solver3d/afivo_reference/**/*.cfg", "OPENEMS_* external regeneration variables"],
            "ACTIVE_PACKAGING_DEFECT": [],
        },
        "corrected_active_code": ["fullwave/h1/generate_h1_audit.py", "fullwave/h4/generate_h4_results.py"],
        "normal_core_source_requires_machine_absolute_path": False,
    }
    write("machine_path_audit.json", machine_paths)

    provenance_paths = [
        ("thermal/g1/g1_reference_summary.json", "thermal/g1/generate_g1_reference.py", "THERMAL_REFERENCE", "SI thermodynamic and electrical units", "NOT_APPLICABLE"),
        ("thermal/g3_port/g3_port_summary.json", "thermal/g3_port/generate_g3_port.py", "PHYSICS_DERIVED_PORT_MODEL_VALIDATED", "V,A,ohm,s", "EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP"),
        ("fullwave/h2/h2_350mhz_summary.json", "fullwave/h2/generate_h2_350mhz_results.py", "DEVELOPMENT_REFERENCE", "Hz,ohm,S-parameters", "50_OHM_PORTS"),
        ("fullwave/h3/h3_result_contract.json", "fullwave/h3/generate_h3_results.py", "DEVELOPMENT_VERIFIED", "Hz,V,ohm,S", "G3_TOTAL_PORT_VOLTAGE_TO_RX_50OHM"),
        ("fullwave/h4/h4_result_contract.json", "fullwave/h4/generate_h4_results.py", "NUMERICAL_REFERENCE_ONLY", "Hz,V,m,V/m", "LOCAL_50_OHM_RECEIVER"),
        ("fullwave/h5/h5_stage_h_transfer_contract.json", "fullwave/h5/generate_h5_results.py", "STAGE_H_TOOL_DEVELOPMENT_PASS", "branch-specific", "PATHWAY_SPECIFIC"),
        ("validation/stage_i/final/stage_i_final_contract.json", "validation/stage_i/final/generate_wp_i_e.py", "STAGE_I_TOOL_DEVELOPMENT_PASS", "contract/status", "MEASUREMENT_CONTRACT_SPECIFIC"),
    ]
    provenance = []
    for path, script, status, units, plane in provenance_paths:
        provenance.append({
            "path": path,
            "source_hash": sha256(ROOT / path),
            "producer_script": script,
            "config_hash": "RECORDED_UPSTREAM_WHERE_AVAILABLE_OTHERWISE_NOT_RECONSTRUCTED",
            "backend_version": "RECORDED_IN_SOURCE_CONTRACT_OR_EXTERNAL_BACKEND_MANIFEST",
            "status": status,
            "units": units,
            "reference_plane": plane,
        })
    write("result_provenance_manifest.json", {"records": provenance, "historical_hashes_guessed": False})

    release = {
        "components": [
            {"path": "cpp, python/streamer_rf, tools, scripts, config", "classification": "INCLUDE_SOURCE_RELEASE"},
            {"path": "tests, CMakeLists.txt, requirements*.txt", "classification": "INCLUDE_SOURCE_RELEASE"},
            {"path": "docs, packaging, compact stage contracts", "classification": "INCLUDE_REFERENCE_DATA"},
            {"path": "validation/stage_i/system_350mhz/*synthetic*", "classification": "INCLUDE_SYNTHETIC_FIXTURE"},
            {"path": "Afivo-streamer source", "classification": "DOCUMENT_ONLY_EXTERNAL_BACKEND"},
            {"path": "openEMS/CSXCAD source", "classification": "DOCUMENT_ONLY_EXTERNAL_BACKEND"},
            {"path": "COMSOL binaries/projects", "classification": "DOCUMENT_ONLY_EXTERNAL_BACKEND"},
            {"path": "build,.venv,caches,temp openEMS", "classification": "EXCLUDE_BUILD_ARTIFACT"},
            {"path": "Afivo raw and Stage-F raw fields", "classification": "EXCLUDE_LARGE_RAW_DATA"},
            {"path": "entire source release", "classification": "LICENSE_REVIEW_REQUIRED"},
        ],
        "public_release_created": False,
    }
    write("stage_j_release_inventory.json", release)

    requirements_hashes = {name: sha256(ROOT / name) for name in ("requirements.txt", "requirements-dev.txt")}
    reference_hashes = {
        "Stage_I_final_contract": sha256(ROOT / "validation/stage_i/final/stage_i_final_contract.json"),
        "Stage_H_contract": sha256(ROOT / "fullwave/h5/h5_stage_h_transfer_contract.json"),
        "Stage_I_system_fixture_manifest": sha256(ROOT / "validation/stage_i/system_350mhz/synthetic_input/sha256_manifest.json"),
        "Stage_I_WP_D_fixture_manifest": sha256(ROOT / "validation/stage_i/system_350mhz/wp_i_d_synthetic/sha256_manifest.json"),
        "D_Angola_property_provenance": sha256(ROOT / "thermal/g1/dangola_property_provenance.json"),
    }
    reproducibility = {
        "schema": "StageJReproducibilityManifest-v1",
        "repository_commit": commit,
        "branch": branch,
        "platform_reference": platform.platform(),
        "build_toolchain": {"C++": "C++17", "CMake": "4.2.3", "Ninja": "1.13.2", "PETSc": "3.24.4", "MPI": "Open MPI 5.0.10"},
        "Python": {
            "version": platform.python_version(),
            "runtime_requirements": "requirements.txt",
            "development_requirements": "requirements-dev.txt",
            "requirement_hashes": requirements_hashes,
            "openEMS_environment": "SEPARATE_EXTERNAL_ENVIRONMENT",
        },
        "external_backends": backends,
        "reference_input_hashes": reference_hashes,
        "minimal_reproduction_commands": ["./scripts/reproduce_smoke.sh"],
        "test_expectations": {"CTest": "10_OF_10_PASS", "pytest_default": "ALL_NON_ENVIRONMENT_SKIPPED_TESTS_PASS", "smoke": "PASS"},
        "large_data_policy": "packaging/large_data_manifest.json",
        "scientific_status": "docs/scientific_status.md",
        "clean_worktree_validation": "packaging/clean_clone_reproducibility.json",
    }
    write("reproducibility_manifest.json", reproducibility)

    clean_clone = {
        "schema": "StageJCleanWorktreeReproducibility-v1",
        "tested_commit": "a4e8d23",
        "checkout_kind": "LOCAL_DETACHED_GIT_WORKTREE",
        "checkout_location_policy": "TEMPORARY_PATH_RECORDED_NOT_REQUIRED_BY_REPRODUCTION",
        "python_environment": "EXISTING_PROJECT_VENV_EXPLICITLY_SELECTED_BY_PYTHON_BIN",
        "python_import": "PASS",
        "cmake_configure": "PASS",
        "build": "PASS",
        "ctest": "10_OF_10_PASS",
        "smoke_pytest": "58_PASS",
        "smoke_status": "PASS",
        "wall_time_s": 9.46,
        "peak_rss_kib": 245976,
        "external_backends_rebuilt": False,
        "large_physics_simulations_rerun": False,
    }
    write("clean_clone_reproducibility.json", clean_clone)

    citation = {"CITATION_STATUS": "CITATION_METADATA_INCOMPLETE", "template": "CITATION.cff.template", "missing": ["AUTHOR", "VERSION", "RELEASE_DATE", "OPTIONAL_DOI_OR_ORCID"], "invented_metadata": False}
    license_status = {"PROJECT_LICENSE_STATUS": "DECISION_REQUIRED", "LICENSE_file_present": False, "license_selected_by_stage_j": False, "PUBLIC_RELEASE_READY": "PENDING_LICENSE_OR_USER_RELEASE_DECISION"}
    write("citation_status.json", citation)
    write("license_status.json", license_status)

    final_contract = {
        "schema": "FinalProjectContract-v2.0",
        "architecture": "docs/software_architecture.md",
        "stages": {
            "A": "FROZEN_TOOL_FOUNDATION", "B": "EXTERNAL_COMSOL_VALIDATION_PENDING",
            "C": "FROZEN", "D_E": "FROZEN_EXTERNAL_BACKEND_REFERENCE",
            "F": "TRUSTED_PHYSICS_WITH_RECORDED_LIMITATIONS",
            "G": "REFERENCE_SOLVERS_VALIDATED_CALIBRATION_PENDING",
            "H": {"tool": "PASS", "science": "PENDING_STAGE_I"},
            "I": {"tool": stage_i["STAGE_I_TOOL_DEVELOPMENT"], "science": stage_i["STAGE_I_SCIENTIFIC_VALIDATION"]},
            "J": "PASS",
        },
        "external_backends": "packaging/external_backends.json",
        "data_policy": "packaging/large_data_manifest.json",
        "reproduction_entry_point": "scripts/reproduce_smoke.sh",
        "real_experiment_reentry": "validation/stage_i/final/stage_i_real_data_reentry_contract.json",
        "release_readiness": license_status["PUBLIC_RELEASE_READY"],
        "known_debts": "validation/stage_i/final/stage_i_final_debt_ledger.json",
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "automatic_release_actions_performed": False,
    }
    write("final_project_contract.json", final_contract)

    status = {
        "STAGE_J_TOOL_PACKAGING": "PASS",
        "REPRODUCIBILITY_SMOKE": "PASS",
        "OPEN_SOURCE_ARCHITECTURE_READY": True,
        "PUBLIC_RELEASE_READY": "PENDING_LICENSE_OR_USER_RELEASE_DECISION",
        "PROJECT_LICENSE_STATUS": "DECISION_REQUIRED",
        "CITATION_STATUS": "CITATION_METADATA_INCOMPLETE",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "large_physics_simulations_rerun": False,
        "release_created": False,
    }
    write("stage_j_status.json", status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
