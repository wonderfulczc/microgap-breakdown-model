import json
import subprocess
import sys
import tarfile
import tomllib
from pathlib import Path

import yaml

from streamer_rf.release.audit import CURRENT_VERSION, release_audit


ROOT = Path(__file__).resolve().parents[2]
ACTIVE_VERSION = "0.1.0"
RC_VERSION = "0.1.0rc1"


def load_json(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_active_release_versions_are_consistent():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    resource = load_json("python/streamer_rf/resources/packaging/release_build_metadata.json")
    assert {project["version"], citation["version"], resource["package_version"], CURRENT_VERSION} == {ACTIVE_VERSION}


def test_source_cli_reports_rc_version():
    result = subprocess.run(
        [sys.executable, "-m", "streamer_rf.cli", "--version"],
        cwd=ROOT,
        env={"PYTHONPATH": str(ROOT / "python")},
        text=True,
        capture_output=True,
        check=True,
    )
    assert result.stdout.strip() == f"microgap-rf {ACTIVE_VERSION}"


def test_changelog_and_release_notes_preserve_scientific_boundary():
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    notes = (ROOT / "release/0.1.0_release_notes.md").read_text(encoding="utf-8")
    assert ACTIVE_VERSION in text and "Initial release" in text
    assert "PENDING_REAL_EXPERIMENT" in text and "PENDING_REAL_EXPERIMENT" in notes
    for debt in ("VNA_MEASUREMENT_PENDING", "FULL_WAVE_LOADING_MISMATCH_HIGH", "FULL_MAXWELL_REFERENCE_PENDING"):
        assert debt in notes


def test_acceptance_matrix_schema_and_no_publish_side_effect():
    record = load_json("release/0.1.0rc1_acceptance.json")
    required = {"VERSION_CONSISTENCY", "LICENSE", "CITATION", "THIRD_PARTY", "WHEEL_BUILD", "SDIST_BUILD", "WHEEL_INSTALL", "SDIST_REBUILD", "CLI_SMOKE", "PYTHON_API", "CTEST", "PYTEST", "REPRODUCIBILITY_SMOKE", "PRIVACY", "SCIENTIFIC_STATUS", "DOCUMENTATION", "EXTERNAL_BACKEND_BOUNDARY"}
    assert set(record["checks"]) == required
    assert set(record["checks"].values()) <= {"PASS", "FAIL", "PENDING"}
    assert record["RC1_PUBLISHED"] is False
    assert record["GITHUB_RELEASE_CREATED"] is False
    assert record["PYPI_PUBLISHED"] is False


def test_rc_artifact_manifest_schema_and_hashes():
    import hashlib

    path = ROOT / "release/0.1.0rc1_artifacts.json"
    if not path.is_file():
        assert "exclude release/0.1.0rc1_artifacts.json" in (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        return
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["version"] == RC_VERSION
    for kind in ("wheel", "sdist"):
        item = record[kind]
        assert len(item["sha256"]) == 64 and item["size_bytes"] > 0
    assert len(record["citation_sha256"]) == 64
    assert record["artifacts_tracked_in_git"] is False
    assert record["RC1_PUBLISHED"] is False


def test_release_audit_reports_rc_without_upgrading_science():
    audit = release_audit(ROOT)
    assert audit["package"]["version"] == ACTIVE_VERSION
    assert audit["gates"]["VERSION"] == ACTIVE_VERSION
    assert audit["gates"]["RC1_ALLOWED"] is True
    assert audit["gates"]["PUBLIC_RELEASE_READY"] == "APPROVED_FOR_GITHUB_RELEASE"
    assert audit["scientific_status"]["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"


def test_manifest_excludes_post_build_audits_and_external_backends():
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "exclude release/0.1.0rc1_artifacts.json" in manifest
    assert "exclude release/0.1.0rc1_acceptance.json" in manifest
    assert "prune solver3d/afivo_reference/stage_e/results_raw" in manifest
