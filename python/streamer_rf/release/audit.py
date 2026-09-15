"""Read-only release-candidate gates and version-state policy."""
from __future__ import annotations

import json
import tarfile
from zipfile import ZipFile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from urllib.parse import urlparse

import yaml


CURRENT_VERSION = "0.1.0rc1"
NEXT_CANDIDATE_VERSION = "0.1.0rc1"
FINAL_VERSION = "0.1.0"
FORBIDDEN_ARCHIVE_MARKERS = ("/.git/", "/.venv/", "/build/", "/results/", "results_raw", "__pycache__", "openems_simulation", "solver3d/afivo")
REQUIRED_WHEEL_FILES = (
    "streamer_rf/cli.py",
    "streamer_rf/thermal/dangola_coefficients.json",
    "streamer_rf/resources/config/unified_case_schema.yaml",
    "streamer_rf/resources/literature/literature_registry.yaml",
    "streamer_rf/resources/literature/evidence_card_template.yaml",
    "streamer_rf/resources/literature/change_request_template.yaml",
    "streamer_rf/resources/packaging/final_project_contract.json",
)


def candidate_version_allowed(gates: dict) -> bool:
    required = (
        gates.get("CLEAN_WHEEL_INSTALL") == "PASS",
        gates.get("SDIST_REBUILD") == "PASS",
        gates.get("LICENSE_DECISION") == "PASS",
        gates.get("CITATION_METADATA") == "PASS",
        gates.get("PRIVACY_AUDIT") in {"PASS", "PASS_WITH_SAFE_HISTORICAL_METADATA"},
        gates.get("SCIENTIFIC_STATUS_AUDIT") == "PASS",
        gates.get("RELEASE_INVENTORY") == "PASS",
    )
    return all(required)


def canonical_repository_url(value: str) -> str:
    """Normalize a GitHub remote without inferring repository visibility."""
    value = value.strip()
    if value.startswith("git@github.com:"):
        value = "https://github.com/" + value.removeprefix("git@github.com:")
    parsed = urlparse(value)
    if parsed.scheme != "https" or parsed.netloc.lower() != "github.com":
        raise ValueError("CANONICAL_REPOSITORY_URL_MUST_BE_HTTPS_GITHUB")
    path = parsed.path.removesuffix(".git").strip("/")
    if len(path.split("/")) != 2:
        raise ValueError("INVALID_GITHUB_REPOSITORY_PATH")
    return f"https://github.com/{path}"


def validate_citation_file(path: Path) -> dict:
    """Validate the release-critical subset of CFF 1.2 without inventing metadata."""
    if not path.is_file():
        return {"valid": False, "errors": ["CITATION_FILE_MISSING"]}
    text = path.read_text(encoding="utf-8")
    placeholders = ("<USER_", "NOT_PROVIDED", "REQUIRED_BEFORE_RELEASE")
    errors = []
    try:
        record = yaml.safe_load(text)
    except yaml.YAMLError:
        record = None
        errors.append("INVALID_YAML")
    if not isinstance(record, dict):
        errors.append("CFF_ROOT_MUST_BE_MAPPING")
        return {"valid": False, "errors": errors}
    for field in ("cff-version", "message", "title", "authors", "repository-code", "license", "version"):
        if not record.get(field):
            errors.append(f"MISSING_{field.upper().replace('-', '_')}")
    if str(record.get("cff-version")) != "1.2.0":
        errors.append("UNSUPPORTED_CFF_VERSION")
    authors = record.get("authors")
    if not isinstance(authors, list) or not authors or not all(
        isinstance(author, dict) and (author.get("family-names") or author.get("name"))
        for author in authors
    ):
        errors.append("INVALID_AUTHORS")
    if any(marker in text for marker in placeholders):
        errors.append("UNRESOLVED_PLACEHOLDER")
    return {"valid": not errors, "errors": errors}


def citation_file_complete(path: Path) -> bool:
    return validate_citation_file(path)["valid"]


def rc1_allowed(root: Path, gates: dict, decisions: dict) -> bool:
    formal_license = (root / "LICENSE").is_file() or gates.get("FORMAL_LICENSE_PRESENT") is True
    formal_citation = citation_file_complete(root / "CITATION.cff") or gates.get("FORMAL_CITATION_VALIDATED") is True
    requirements = (
        formal_license and gates.get("LICENSE_DECISION") == "PASS" and decisions.get("license_approved") is True,
        gates.get("THIRD_PARTY_LICENSE_REVIEW") == "PASS" or gates.get("THIRD_PARTY_DOCUMENTED_EXCEPTION") == "ACCEPTED",
        formal_citation and gates.get("CITATION_METADATA") == "PASS",
        decisions.get("repository_url_confirmed") is True,
        decisions.get("authors_and_order_confirmed") is True,
        decisions.get("user_approved_rc1") is True,
        gates.get("SCIENTIFIC_STATUS_AUDIT") == "PASS",
        gates.get("CLEAN_WHEEL_INSTALL") == "PASS",
        gates.get("SDIST_REBUILD") == "PASS",
    )
    return all(requirements)


def audit_archive(path: Path, *, kind: str) -> dict:
    path = Path(path)
    if kind == "wheel":
        with ZipFile(path) as archive:
            names = archive.namelist()
            total = sum(item.file_size for item in archive.infolist())
        missing = [name for name in REQUIRED_WHEEL_FILES if name not in names]
    elif kind == "sdist":
        with tarfile.open(path) as archive:
            members = archive.getmembers()
            names = [item.name for item in members]
            total = sum(item.size for item in members if item.isfile())
        missing = []
    else:
        raise ValueError("UNKNOWN_ARCHIVE_KIND")
    forbidden = [name for name in names if any(marker in f"/{name}" for marker in FORBIDDEN_ARCHIVE_MARKERS)]
    return {"kind": kind, "file_count": len(names), "uncompressed_bytes": total, "required_missing": missing, "forbidden": forbidden, "status": "PASS" if not missing and not forbidden else "FAIL"}


def release_audit(root: Path) -> dict:
    root = Path(root)
    gate_path = root / "packaging/rp2_release_gate.json"
    if not gate_path.is_file():
        gate_path = root / "packaging/release_gate.json"
    gates = json.loads(gate_path.read_text()) if gate_path.is_file() else {}
    gates = dict(gates)
    decisions_path = root / "release/rp3_user_decision_resolved.json"
    if not decisions_path.is_file():
        decisions_path = root / "packaging/rp3_user_decision_resolved.json"
    decisions_record = json.loads(decisions_path.read_text()) if decisions_path.is_file() else {}
    answers = decisions_record.get("resolved_gate_inputs", {})
    formal_license = (root / "LICENSE").is_file() or gates.get("FORMAL_LICENSE_PRESENT") is True
    formal_citation = citation_file_complete(root / "CITATION.cff") or gates.get("FORMAL_CITATION_VALIDATED") is True
    if not formal_license or answers.get("license_approved") is not True:
        gates["LICENSE_DECISION"] = "PENDING_USER_DECISION"
    if not formal_citation:
        gates["CITATION_METADATA"] = "INCOMPLETE_USER_INPUT_REQUIRED"
    inventory_path = root / "packaging/third_party_license_inventory.json"
    if not inventory_path.is_file():
        inventory_path = root / "packaging/third_party_license_status.json"
    if inventory_path.is_file():
        gates["THIRD_PARTY_LICENSE_REVIEW"] = json.loads(inventory_path.read_text()).get("status", "PENDING_EXTERNAL_VERIFICATION")
    gates["RC1_ALLOWED"] = rc1_allowed(root, gates, answers)
    gates["VERSION"] = CURRENT_VERSION
    gates["RC1_BUILD_READY"] = gates["RC1_ALLOWED"] and CURRENT_VERSION == NEXT_CANDIDATE_VERSION
    gates["PUBLIC_RELEASE_READY"] = "AWAITING_FINAL_RELEASE_APPROVAL" if gates["RC1_BUILD_READY"] else "PENDING_LICENSE_OR_USER_RELEASE_DECISION"
    docs = [
        "README.zh-CN.md", "docs/zh/软件架构.md", "docs/zh/CLI架构说明.md",
        "docs/zh/配置与结果合同.md", "docs/zh/科学状态.md", "docs/zh/数据政策.md",
        "docs/zh/实验数据重入.md", "docs/zh/文献维护流程.md",
        "docs/zh/科研模型变更流程.md", "docs/zh/许可证选择决策说明.md",
        "docs/zh/开源许可与知识产权边界.md", "docs/zh/软件引用与学术署名说明.md",
        "docs/zh/RP3用户决策清单.md",
    ]
    if (root / "pyproject.toml").is_file():
        installed_version = CURRENT_VERSION
    else:
        try:
            installed_version = version("microgap-rf")
        except PackageNotFoundError:
            installed_version = CURRENT_VERSION
    return {
        "audit_only": True,
        "side_effects": {"tag_created": False, "release_created": False, "upload_performed": False, "license_selected": False},
        "package": {"project_name": "microgap-rf", "import_name": "streamer_rf", "cli_name": "microgap-rf", "version": installed_version},
        "version_policy": {"CURRENT_VERSION": CURRENT_VERSION, "NEXT_CANDIDATE_VERSION": NEXT_CANDIDATE_VERSION, "FINAL_VERSION": FINAL_VERSION, "candidate_transition_allowed": gates["RC1_ALLOWED"], "RC1_ALLOWED": gates["RC1_ALLOWED"]},
        "required_docs_present": all((root / path).is_file() for path in docs) or gates.get("REQUIRED_DOCS_PRESENT") is True,
        "gates": gates,
        "scientific_status": {"STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT", "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False, "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED", "NATIVE_RF_350MHZ": "NOT_RESOLVED"},
    }
