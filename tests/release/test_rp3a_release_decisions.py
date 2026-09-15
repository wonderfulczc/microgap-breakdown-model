import json
from pathlib import Path

from streamer_rf.release.audit import canonical_repository_url, citation_file_complete, rc1_allowed, release_audit


ROOT = Path(__file__).resolve().parents[2]


def load(path):
    return json.loads((ROOT / path).read_text())


def test_third_party_inventory_schema_and_fparser_evidence():
    inventory = load("packaging/third_party_license_inventory.json")
    required = {"name", "version_or_commit", "license_name", "license_source", "distribution_mode", "linked_or_invoked", "included_in_wheel", "included_in_sdist", "notice_requirement", "redistribution_concern", "confidence"}
    assert inventory["status"] == "PASS"
    assert all(required <= record.keys() for record in inventory["records"])
    fparser = next(record for record in inventory["records"] if record["name"] == "fparser")
    assert fparser["license_status"] == "CONFIRMED"
    assert "docs/fparser.html" in fparser["license_source"]


def test_license_recommendation_is_non_effective_and_requires_user():
    recommendation = load("release/license_recommendation.json")
    assert recommendation["PRIMARY_RECOMMENDATION"]["license"] == "Apache-2.0"
    assert recommendation["SECONDARY_OPTION"]["license"] == "BSD-3-Clause"
    assert recommendation["effective_license_grant"] is False
    assert recommendation["USER_DECISION_REQUIRED"] is True
    assert not (ROOT / "LICENSE").exists()
    assert not (ROOT / "LICENSE.candidate").exists()
    metadata = load("release/license_candidate_metadata.json")
    assert metadata["notice"] == "NOT_EFFECTIVE_UNTIL_USER_APPROVAL"
    assert metadata["effective_license_grant"] is False


def test_citation_candidate_retains_placeholders_and_formal_file_absent():
    text = (ROOT / "CITATION.cff.candidate").read_text()
    assert "<USER_REQUIRED" in text and "<USER_CONFIRM" in text
    assert not (ROOT / "CITATION.cff").exists()
    assert not citation_file_complete(ROOT / "CITATION.cff.candidate")


def test_release_user_input_schema_is_small_and_unresolved():
    record = load("release/rp3_user_decision_required.json")
    assert set(record["questions"]) == {"A_LICENSE", "B_authors_order", "C_repository_URL", "D_optional_ORCID_affiliation", "E_preferred_citation_policy", "F_approve_RC1_later"}
    assert record["automatic_release_action"] is False
    assert not any(record["resolved_gate_inputs"].values())


def test_rc1_gate_rejects_candidate_files_and_missing_approval(tmp_path):
    gates = {"LICENSE_DECISION": "PASS", "THIRD_PARTY_LICENSE_REVIEW": "PASS", "CITATION_METADATA": "PASS", "SCIENTIFIC_STATUS_AUDIT": "PASS", "CLEAN_WHEEL_INSTALL": "PASS"}
    (tmp_path / "LICENSE").write_text("approved text")
    (tmp_path / "CITATION.cff").write_text("cff-version: 1.2.0\nmessage: cite\ntitle: tool\nauthors:\n  - family-names: A\n")
    approved = {"license_approved": True, "repository_url_confirmed": True, "authors_and_order_confirmed": True, "user_approved_rc1": True}
    assert not rc1_allowed(tmp_path, gates, {**approved, "license_approved": False})
    assert not rc1_allowed(tmp_path, gates, {**approved, "user_approved_rc1": False})
    assert rc1_allowed(tmp_path, gates, approved)


def test_repository_url_normalization_does_not_infer_visibility():
    expected = "https://github.com/wonderfulczc/microgap-breakdown-model"
    assert canonical_repository_url("https://github.com/wonderfulczc/microgap-breakdown-model.git") == expected
    assert canonical_repository_url("git@github.com:wonderfulczc/microgap-breakdown-model.git") == expected
    assert load("release/release_metadata_user_inputs.json")["repository"]["REPOSITORY_VISIBILITY"] == "UNKNOWN"


def test_release_audit_preserves_science_and_pending_user_gates():
    result = release_audit(ROOT)
    assert result["gates"]["LICENSE_DECISION"] == "PENDING_USER_DECISION"
    assert result["gates"]["CITATION_METADATA"] == "INCOMPLETE_USER_INPUT_REQUIRED"
    assert result["gates"]["THIRD_PARTY_LICENSE_REVIEW"] == "PASS"
    assert result["gates"]["RC1_ALLOWED"] is False
    assert result["scientific_status"] == {"STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT", "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False, "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED", "NATIVE_RF_350MHZ": "NOT_RESOLVED"}


def test_external_coupling_and_license_risk_schemas():
    coupling = load("release/external_backend_coupling_audit.json")
    assert {item["name"] for item in coupling["backends"]} == {"Afivo-streamer", "openEMS", "CSXCAD", "fparser", "AppCSXCAD"}
    assert all(item["project_artifact_inclusion"] is False for item in coupling["backends"])
    risks = load("release/license_compatibility_audit.json")["candidates"]
    assert "NO_OBVIOUS_CONFLICT" in risks["Apache-2.0"]
    assert "LINKING_REVIEW_REQUIRED" in risks["GPL-3.0-only"]
