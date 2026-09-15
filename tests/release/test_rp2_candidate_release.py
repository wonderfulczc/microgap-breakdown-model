import json
import tarfile
import tomllib
from importlib.resources import files
from pathlib import Path
from zipfile import ZipFile

from streamer_rf.release.audit import (
    CURRENT_VERSION,
    NEXT_CANDIDATE_VERSION,
    audit_archive,
    candidate_version_allowed,
    release_audit,
)


ROOT = Path(__file__).resolve().parents[2]


def test_package_metadata_and_console_entry_point():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert project["name"] == "microgap-rf"
    assert project["version"] == CURRENT_VERSION
    assert project["scripts"]["microgap-rf"] == "streamer_rf.cli:main"
    assert project["requires-python"] == ">=3.11"


def test_runtime_package_data_are_available_by_resource_api():
    package = files("streamer_rf")
    required = [
        "thermal/dangola_coefficients.json",
        "resources/config/unified_case_schema.yaml",
        "resources/literature/literature_registry.yaml",
        "resources/literature/evidence_card_template.yaml",
        "resources/literature/change_request_template.yaml",
        "resources/packaging/final_project_contract.json",
    ]
    assert all(package.joinpath(path).is_file() and package.joinpath(path).read_bytes() for path in required)


def test_wheel_archive_audit_accepts_runtime_and_rejects_external_source(tmp_path):
    good = tmp_path / "good.whl"
    required = [
        "streamer_rf/cli.py", "streamer_rf/thermal/dangola_coefficients.json",
        "streamer_rf/resources/config/unified_case_schema.yaml",
        "streamer_rf/resources/literature/literature_registry.yaml",
        "streamer_rf/resources/literature/evidence_card_template.yaml",
        "streamer_rf/resources/literature/change_request_template.yaml",
        "streamer_rf/resources/packaging/final_project_contract.json",
    ]
    with ZipFile(good, "w") as archive:
        for name in required:
            archive.writestr(name, "x")
    assert audit_archive(good, kind="wheel")["status"] == "PASS"
    bad = tmp_path / "bad.whl"
    with ZipFile(bad, "w") as archive:
        for name in required:
            archive.writestr(name, "x")
        archive.writestr("solver3d/afivo/source.f90", "x")
    result = audit_archive(bad, kind="wheel")
    assert result["status"] == "FAIL" and result["forbidden"]


def test_sdist_archive_exclusion_policy(tmp_path):
    good = tmp_path / "good.tar.gz"
    source = tmp_path / "README.md"
    source.write_text("ok")
    with tarfile.open(good, "w:gz") as archive:
        archive.add(source, arcname="microgap-rf/README.md")
    assert audit_archive(good, kind="sdist")["status"] == "PASS"


def test_license_and_citation_gates_remain_pending_user_input():
    gates = json.loads((ROOT / "packaging/rp2_release_gate.json").read_text())
    citation = json.loads((ROOT / "packaging/citation_metadata_audit.json").read_text())
    assert not (ROOT / "LICENSE").exists()
    assert gates["LICENSE_DECISION"] == "PENDING_USER_DECISION"
    assert gates["CITATION_METADATA"] == "INCOMPLETE_USER_INPUT_REQUIRED"
    assert citation["invented_metadata"] is False
    assert not (ROOT / "CITATION.cff").exists()


def test_scientific_validation_gate_is_immutable():
    gates = json.loads((ROOT / "packaging/rp2_release_gate.json").read_text())
    assert gates["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert gates["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False
    assert gates["SYSTEM_350MHZ_VALIDATION"] == "NOT_MEASURED"
    assert gates["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"


def test_release_audit_has_no_release_side_effect():
    result = release_audit(ROOT)
    assert result["audit_only"] is True
    assert not any(result["side_effects"].values())
    assert result["scientific_status"]["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False


def test_version_state_requires_user_gates_before_rc1():
    gates = json.loads((ROOT / "packaging/rp2_release_gate.json").read_text())
    assert CURRENT_VERSION == "0.1.0.dev0"
    assert NEXT_CANDIDATE_VERSION == "0.1.0rc1"
    assert candidate_version_allowed(gates) is False
    complete = {key: "PASS" for key in ("CLEAN_WHEEL_INSTALL", "SDIST_REBUILD", "LICENSE_DECISION", "CITATION_METADATA", "PRIVACY_AUDIT", "SCIENTIFIC_STATUS_AUDIT", "RELEASE_INVENTORY")}
    assert candidate_version_allowed(complete) is True


def test_external_backends_are_excluded_from_candidate_artifacts():
    inventory = json.loads((ROOT / "packaging/rp2_candidate_release_inventory.json").read_text())
    excluded = " ".join(inventory["exclude"])
    assert all(name in excluded for name in ("Afivo", "openEMS", "COMSOL"))
    assert inventory["release_performed"] is False


def test_synthetic_fixture_labels_and_distribution_advice():
    policy = json.loads((ROOT / "packaging/synthetic_fixture_policy.json").read_text())
    for fixture in policy["fixtures"]:
        assert set(fixture["labels"]) == {"SYNTHETIC_DEVELOPMENT_INPUT", "SYNTHETIC_DRY_RUN", "NOT_EXPERIMENTAL_DATA"}
    inventory = json.loads((ROOT / "packaging/rp2_candidate_release_inventory.json").read_text())
    assert inventory["synthetic_fixture_distribution"] == "KEEP_IN_REPOSITORY_ONLY"


def test_third_party_unknown_is_not_guessed():
    inventory = json.loads((ROOT / "packaging/third_party_license_inventory.json").read_text())
    fparser = next(item for item in inventory["records"] if item["dependency"] == "fparser")
    assert fparser["license"] == "UNKNOWN_REQUIRES_REVIEW"
    assert inventory["status"] == "INCOMPLETE_REVIEW_REQUIRED"
