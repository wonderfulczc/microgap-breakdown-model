"""Registry-only literature impact analysis; never changes model code."""
from __future__ import annotations

from pathlib import Path

import yaml


STATUSES = {"REFERENCE_ONLY", "PARAMETER_UPDATE_CANDIDATE", "MODEL_UPDATE_CANDIDATE", "ARCHITECTURE_IMPACT", "PENDING_REVIEW"}


def load_registry(path: str | Path) -> dict:
    registry = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if registry.get("schema_version") != "1.0" or not isinstance(registry.get("papers"), list):
        raise ValueError("INVALID_LITERATURE_REGISTRY")
    required = {"paper_id", "title", "authors", "year", "journal", "doi", "url", "topic_tags", "related_stages", "relevance", "review_status", "impact_status", "notes"}
    ids = set()
    for paper in registry["papers"]:
        if not required.issubset(paper) or paper["impact_status"] not in STATUSES:
            raise ValueError("INVALID_LITERATURE_RECORD")
        if paper["paper_id"] in ids:
            raise ValueError("DUPLICATE_PAPER_ID")
        ids.add(paper["paper_id"])
    return registry


def impact_analysis(registry: dict, paper_id: str) -> dict:
    paper = next((item for item in registry["papers"] if item["paper_id"] == paper_id), None)
    if paper is None:
        raise ValueError(f"UNKNOWN_PAPER_ID:{paper_id}")
    recommendation = {
        "REFERENCE_ONLY": "NO_CHANGE",
        "PENDING_REVIEW": "REVIEW",
        "PARAMETER_UPDATE_CANDIDATE": "CHANGE_REQUEST_REQUIRED",
        "MODEL_UPDATE_CANDIDATE": "CHANGE_REQUEST_REQUIRED",
        "ARCHITECTURE_IMPACT": "CHANGE_REQUEST_REQUIRED",
    }[paper["impact_status"]]
    return {
        "paper_id": paper_id,
        "possibly_affected_stages": paper["related_stages"],
        "related_modules": paper.get("related_modules", []),
        "related_tests": paper.get("related_tests", []),
        "frozen_result_risk": paper.get("frozen_result_risk", "REVIEW_REQUIRED"),
        "recommendation": recommendation,
        "automatic_code_change": False,
    }

