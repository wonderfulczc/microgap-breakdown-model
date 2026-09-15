import json
import tomllib
from pathlib import Path

import yaml

from streamer_rf.release.audit import CURRENT_VERSION, release_audit, validate_citation_file


ROOT = Path(__file__).resolve().parents[2]


def load_json(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_formal_apache_license_and_pep621_metadata():
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert "Apache License" in license_text and "Version 2.0, January 2004" in license_text
    assert project["license"] == "Apache-2.0"
    assert project["authors"] == [{"name": "Zach"}]
    assert project["version"] == CURRENT_VERSION


def test_formal_citation_preserves_user_metadata_without_invention():
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert validate_citation_file(ROOT / "CITATION.cff")["valid"] is True
    assert citation["authors"] == [{"family-names": "Zach", "affiliation": "College of Integrated Circuits of Southeast University"}]
    assert "orcid" not in citation["authors"][0]
    assert "given-names" not in citation["authors"][0]
    assert citation["date-released"] == "2026-09-15"
    assert "preferred-citation" not in citation


def test_repository_and_future_paper_policy_are_frozen():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    resolved = load_json("release/rp3_user_decision_resolved.json")
    url = "https://github.com/wonderfulczc/microgap-breakdown-model"
    assert project["urls"]["Repository"] == url
    assert resolved["CITATION_POLICY"] == "CURRENT_SOFTWARE_CITATION_WITH_FUTURE_PREFERRED_PAPER"


def test_rc1_gate_is_closed_but_version_and_release_actions_are_unchanged():
    audit = release_audit(ROOT)
    status = load_json("release/rp3_status.json")
    assert audit["gates"]["RC1_ALLOWED"] is True
    assert audit["gates"]["PUBLIC_RELEASE_READY"] == "APPROVED_FOR_GITHUB_RELEASE"
    assert status["CURRENT_VERSION"] == "0.1.0.dev0"
    assert status["release_actions_performed"] is False
    assert not any(audit["side_effects"].values())


def test_scientific_state_is_preserved_exactly():
    science = release_audit(ROOT)["scientific_status"]
    assert science == {
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED",
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    }


def test_stale_candidate_files_are_removed_and_resolution_is_traceable():
    for path in (
        "CITATION.cff.candidate",
        "release/license_candidate_metadata.json",
        "release/release_metadata_user_inputs.json",
        "release/rp3_user_decision_required.json",
    ):
        assert not (ROOT / path).exists()
    assert load_json("release/rp3_user_decision_resolved.json")["status"] == "RESOLVED"
