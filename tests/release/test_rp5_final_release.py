import json
import subprocess
import sys
import tomllib
from pathlib import Path

import yaml

from streamer_rf.release.audit import CURRENT_VERSION, FINAL_VERSION, release_audit


ROOT = Path(__file__).resolve().parents[2]
VERSION = "0.1.0"


def test_final_version_metadata_are_consistent():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    resource = json.loads((ROOT / "python/streamer_rf/resources/packaging/release_build_metadata.json").read_text())
    assert {project["version"], citation["version"], resource["package_version"], CURRENT_VERSION, FINAL_VERSION} == {VERSION}
    assert citation["date-released"] == "2026-09-15"


def test_final_release_approval_excludes_pypi():
    approval = json.loads((ROOT / "release/rp5_final_release_approval.json").read_text())
    assert approval["FINAL_RELEASE_APPROVED"] is True
    assert approval["GITHUB_RELEASE_APPROVED"] is True
    assert approval["PYPI_PUBLICATION_APPROVED"] is False


def test_final_release_notes_and_changelog_preserve_science():
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    notes = (ROOT / "release/0.1.0_release_notes.md").read_text(encoding="utf-8")
    assert "## 0.1.0 - Initial release" in changelog
    for status in ("PENDING_REAL_EXPERIMENT", "NOT_MEASURED", "NOT_RESOLVED"):
        assert status in notes
    for debt in ("VNA_MEASUREMENT_PENDING", "FULL_WAVE_LOADING_MISMATCH_HIGH", "FULL_MAXWELL_REFERENCE_PENDING"):
        assert debt in notes


def test_source_cli_and_release_audit_report_final_version_without_side_effects():
    result = subprocess.run(
        [sys.executable, "-m", "streamer_rf.cli", "--version"],
        cwd=ROOT,
        env={"PYTHONPATH": str(ROOT / "python")},
        text=True,
        capture_output=True,
        check=True,
    )
    assert result.stdout.strip() == "microgap-rf 0.1.0"
    audit = release_audit(ROOT)
    assert audit["package"]["version"] == VERSION
    assert audit["gates"]["FINAL_RELEASE_APPROVED"] is True
    assert audit["gates"]["FINAL_VERSION_BUILD_READY"] is True
    assert audit["gates"]["PUBLIC_RELEASE_READY"] == "APPROVED_FOR_GITHUB_RELEASE"
    assert not any(audit["side_effects"].values())


def test_scientific_status_is_not_upgraded_by_release():
    assert release_audit(ROOT)["scientific_status"] == {
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    }
