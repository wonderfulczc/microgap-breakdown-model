"""Read-only release-candidate gates and version-state policy."""
from __future__ import annotations

import json
import tarfile
from zipfile import ZipFile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


CURRENT_VERSION = "0.1.0.dev0"
NEXT_CANDIDATE_VERSION = "0.1.0rc1"
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
        gates.get("PRIVACY_AUDIT") == "PASS",
        gates.get("SCIENTIFIC_STATUS_AUDIT") == "PASS",
        gates.get("RELEASE_INVENTORY") == "PASS",
    )
    return all(required)


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
    gates = json.loads(gate_path.read_text()) if gate_path.is_file() else {}
    docs = [
        "README.zh-CN.md", "docs/zh/软件架构.md", "docs/zh/CLI架构说明.md",
        "docs/zh/配置与结果合同.md", "docs/zh/科学状态.md", "docs/zh/数据政策.md",
        "docs/zh/实验数据重入.md", "docs/zh/文献维护流程.md",
        "docs/zh/科研模型变更流程.md", "docs/zh/许可证选择决策说明.md",
    ]
    try:
        installed_version = version("microgap-rf")
    except PackageNotFoundError:
        installed_version = CURRENT_VERSION
    return {
        "audit_only": True,
        "side_effects": {"tag_created": False, "release_created": False, "upload_performed": False, "license_selected": False},
        "package": {"project_name": "microgap-rf", "import_name": "streamer_rf", "cli_name": "microgap-rf", "version": installed_version},
        "version_policy": {"CURRENT_VERSION": CURRENT_VERSION, "NEXT_CANDIDATE_VERSION": NEXT_CANDIDATE_VERSION, "candidate_transition_allowed": candidate_version_allowed(gates)},
        "required_docs_present": all((root / path).is_file() for path in docs),
        "gates": gates,
        "scientific_status": {"STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT", "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False, "SYSTEM_350MHZ_VALIDATION": "NOT_MEASURED", "NATIVE_RF_350MHZ": "NOT_RESOLVED"},
    }
