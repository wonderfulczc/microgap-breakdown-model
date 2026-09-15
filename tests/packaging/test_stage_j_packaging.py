import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGING = ROOT / "packaging"


def load(name):
    return json.loads((PACKAGING / name).read_text())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_top_level_repository_inventory_is_complete():
    inventory = load("stage_j_repository_inventory.json")["directories"]
    recorded = {entry["path"] for entry in inventory}
    actual = {path.name for path in ROOT.iterdir() if path.is_dir()}
    assert recorded == actual
    allowed = {
        "CORE", "REQUIRED_REFERENCE", "REGENERABLE", "EXTERNAL_BACKEND",
        "SYNTHETIC_FIXTURE", "ARCHIVE_ONLY", "BUILD_ARTIFACT", "EXCLUDE_FROM_RELEASE",
    }
    assert {entry["classification"] for entry in inventory} <= allowed


def test_external_backends_are_pinned_and_not_vendored():
    backends = load("external_backends.json")
    assert backends["Afivo"]["commit"] == "a50b5508775086e90dfe423455fb58d812578410"
    assert backends["Afivo"]["source_policy"] == "EXTERNAL_NOT_VENDORED"
    assert backends["openEMS"]["source_policy"] == "EXTERNAL_NOT_VENDORED"
    assert backends["openEMS"]["commits"]["openEMS"] == "8f480d04e0e17a780df28e0ac12e2086be041f9e"
    assert backends["openEMS"]["commits"]["CSXCAD"] == "0458ee11ad711909cb32390b195079c3cff0606e"
    assert backends["openEMS"]["commits"]["fparser"] == "4b9c845b449b520c4b8c5f23c74cd04820084f81"
    assert backends["COMSOL"]["source_policy"] == "PROPRIETARY_EXTERNAL_TOOL_NOT_DISTRIBUTED"


def test_runtime_and_development_requirements_are_pinned_and_separate():
    runtime = (ROOT / "requirements.txt").read_text().splitlines()
    development = (ROOT / "requirements-dev.txt").read_text().splitlines()
    assert runtime
    assert all(re.fullmatch(r"[A-Za-z0-9_-]+==[0-9][A-Za-z0-9.]*", line) for line in runtime)
    assert development[0] == "-r requirements.txt"
    assert development[1].startswith("pytest==")
    assert not any(line.startswith("pytest==") for line in runtime)


def test_smoke_command_is_fail_fast_and_does_not_run_external_backends():
    script = (ROOT / "scripts/reproduce_smoke.sh").read_text()
    assert "set -euo pipefail" in script
    assert "ctest --test-dir" in script
    assert "test_stage_f2_jefimenko.py" in script
    assert "test_stage_g3_port.py" in script
    assert "test_h3_transient.py" in script
    assert "test_stage_i_wp_e.py" in script
    assert "openEMS" not in script
    assert "afivo-streamer" not in script


def test_large_data_is_excluded_without_deletion():
    manifest = load("large_data_manifest.json")
    policies = {entry["dataset"]: entry["release_policy"] for entry in manifest["datasets"]}
    assert policies == {"AFIVO_STAGE_E_RAW": "EXTERNAL_ARCHIVE", "STAGE_F_RAW_FIELDS": "DO_NOT_PACKAGE"}
    assert manifest["normal_source_release_contains_multi_GB_raw_data"] is False
    assert manifest["local_files_deleted"] is False


def test_synthetic_fixtures_cannot_be_mistaken_for_measurements():
    policy = load("synthetic_fixture_policy.json")
    for fixture in policy["fixtures"]:
        assert set(fixture["labels"]) == {
            "SYNTHETIC_DEVELOPMENT_INPUT", "SYNTHETIC_DRY_RUN", "NOT_EXPERIMENTAL_DATA"
        }
        assert len(sha256(ROOT / fixture["manifest"])) == 64


def test_machine_paths_have_no_remaining_active_packaging_defect():
    audit = load("machine_path_audit.json")
    assert audit["classifications"]["ACTIVE_PACKAGING_DEFECT"] == []
    assert audit["normal_core_source_requires_machine_absolute_path"] is False
    h1 = (ROOT / "fullwave/h1/generate_h1_audit.py").read_text()
    assert 'os.environ["OPENEMS_PROJECT_ROOT"]' in h1
    assert 'os.environ["OPENEMS_ROOT"]' in h1


def test_license_and_citation_are_not_invented():
    license_status = load("license_status.json")
    citation = load("citation_status.json")
    assert license_status["PROJECT_LICENSE_STATUS"] == "DECISION_REQUIRED"
    assert license_status["license_selected_by_stage_j"] is False
    assert citation["CITATION_STATUS"] == "CITATION_METADATA_INCOMPLETE"
    assert citation["invented_metadata"] is False


def test_release_inventory_has_all_required_boundaries():
    inventory = load("stage_j_release_inventory.json")
    classes = {item["classification"] for item in inventory["components"]}
    assert classes == {
        "INCLUDE_SOURCE_RELEASE", "INCLUDE_REFERENCE_DATA", "INCLUDE_SYNTHETIC_FIXTURE",
        "DOCUMENT_ONLY_EXTERNAL_BACKEND", "EXCLUDE_BUILD_ARTIFACT", "EXCLUDE_LARGE_RAW_DATA",
        "LICENSE_REVIEW_REQUIRED",
    }
    assert inventory["public_release_created"] is False


def test_result_provenance_hashes_match_files():
    records = load("result_provenance_manifest.json")["records"]
    assert records
    for record in records:
        assert record["source_hash"] == sha256(ROOT / record["path"])
        assert record["producer_script"]
        assert record["units"]


def test_reproducibility_manifest_preserves_stage_i_science_boundary():
    manifest = load("reproducibility_manifest.json")
    assert manifest["branch"] == "stage-j-packaging"
    assert manifest["build_toolchain"]["PETSc"] == "3.24.4"
    assert manifest["minimal_reproduction_commands"] == ["./scripts/reproduce_smoke.sh"]
    status = load("stage_j_status.json")
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert status["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False
    assert status["STAGE_J_TOOL_PACKAGING"] == "PASS"
    assert status["REPRODUCIBILITY_SMOKE"] == "PASS"


def test_clean_worktree_reproduction_evidence_is_recorded():
    evidence = load("clean_clone_reproducibility.json")
    assert evidence["smoke_status"] == "PASS"
    assert evidence["ctest"] == "10_OF_10_PASS"
    assert evidence["smoke_pytest"] == "58_PASS"
    assert evidence["external_backends_rebuilt"] is False
    assert evidence["large_physics_simulations_rerun"] is False


def test_final_project_contract_has_no_release_side_effects():
    contract = load("final_project_contract.json")
    assert contract["schema"] == "FinalProjectContract-v2.0"
    assert contract["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert contract["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False
    assert contract["automatic_release_actions_performed"] is False
    assert contract["stages"]["J"] == "PASS"
    assert contract["real_experiment_reentry"].endswith("stage_i_real_data_reentry_contract.json")


def test_gitignore_excludes_build_cache_and_external_runtime_products():
    ignore = (ROOT / ".gitignore").read_text()
    for pattern in ("build/", "build-smoke/", ".venv/", "__pycache__/", ".pytest_cache/", "*.log", "fullwave/**/raw/"):
        assert pattern in ignore
