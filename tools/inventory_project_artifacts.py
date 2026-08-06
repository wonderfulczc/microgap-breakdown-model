#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    ".md", ".csv", ".json", ".yaml", ".yml", ".html", ".htm", ".txt",
    ".py", ".cpp", ".hpp", ".h", ".cmake", ".svg", ".js", ".css",
}
SCAN_DIRS = ["docs", "results", "config", "cpp", "python", "tests", "tools"]
CANONICAL_DOC_TOKENS = [
    "closure_report", "validation_matrix", "decision_log", "formula_audit",
    "resource", "project_simulation_transfer", "final_project_structure",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_type(path: Path) -> str:
    if path.is_symlink():
        return "symlink"
    if path.suffix:
        return path.suffix.lower().lstrip(".")
    return "no_suffix"


def stage_for(path: str) -> str:
    m = re.search(r"stage([1-5])", path, re.I)
    if m:
        return f"stage{m.group(1)}"
    return ""


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="ignore")
    except Exception:
        return ""


def referenced_paths(files: list[Path]) -> tuple[dict[str, set[str]], list[dict[str, str]]]:
    all_rel = {rel(p) for p in files}
    basename_index: dict[str, set[str]] = defaultdict(set)
    for p in all_rel:
        basename_index[Path(p).name].add(p)
    refs: dict[str, set[str]] = defaultdict(set)
    edges: list[dict[str, str]] = []

    path_pattern = re.compile(r"(?:(?:\./)?(?:docs|results|config|cpp|python|tests|tools|archive|build)/[A-Za-z0-9_./:@%+=,\-]+)")
    for f in files:
        if f.suffix.lower() not in TEXT_SUFFIXES:
            continue
        source = rel(f)
        text = read_text(f)
        found = set(path_pattern.findall(text))
        for token in found:
            target = token.removeprefix("./").rstrip(").,;\"'")
            if target in all_rel:
                refs[target].add(source)
                edges.append({"source": source, "target": target, "reference_type": "path_text"})

        # Registry and manifest CSVs often contain direct path columns.
        if f.name in {"run_registry.csv", "derived_artifact_manifest.csv"} or "validation_matrix" in f.name:
            try:
                with f.open(newline="", errors="ignore") as fp:
                    for row in csv.DictReader(fp):
                        for key, value in row.items():
                            if not value:
                                continue
                            for part in str(value).split(";"):
                                p = part.strip().removeprefix("./")
                                if p in all_rel:
                                    refs[p].add(source)
                                    edges.append({"source": source, "target": p, "reference_type": key or "csv_value"})
            except Exception:
                pass
    return refs, edges


def run_id_for(path: str) -> str:
    try:
        for reg in ROOT.glob("results/stage*/provenance/run_registry.csv"):
            with reg.open(newline="", errors="ignore") as fp:
                for row in csv.DictReader(fp):
                    rf = row.get("result_files", "")
                    if rf == path:
                        return row.get("run_id", "")
    except Exception:
        return ""
    return ""


def classify(path: str, size: int, ref_count: int) -> tuple[str, str, str, str]:
    p = Path(path)
    parts = set(p.parts)
    lower = path.lower()

    if any(part == "__pycache__" for part in p.parts) or p.suffix == ".pyc" or ".pytest_cache" in parts:
        return "cache", "quarantine", "Python/test cache; reproducible and not evidence", "high"
    if p.name.endswith(("~", ".tmp", ".bak", ".swp")) or p.name in {".ninja_log", ".ninja_deps"}:
        return "temporary", "quarantine", "temporary/editor/build log file", "high"
    if path.startswith(("cpp/", "python/", "tools/")):
        return "canonical_source", "keep", "source or workflow tool", "high"
    if path.startswith("config/") or p.name in {"CMakeLists.txt", "README.md"}:
        return "canonical_config", "keep", "configuration/build entrypoint", "high"
    if path.startswith("tests/") or "test_" in p.name:
        return "canonical_test", "keep", "test source", "high"
    if path.startswith("docs/"):
        if any(tok in lower for tok in CANONICAL_DOC_TOKENS):
            return "canonical_evidence", "keep", "report/method/validation documentation", "high"
        if path.startswith("docs/simulation_reconstruction_summary/"):
            return "canonical_result", "keep", "final offline HTML summary artifact", "high"
        return "canonical_evidence", "keep", "documentation retained by default", "medium"
    if "invalidated_results" in lower:
        return "audit_only", "archive", "invalidated Stage 2 results retained only as audit archive", "high"
    if any(tok in lower for tok in ["interrupted", "corrupt", "partial_hung", "paused", "failed_attempt", "aborted"]):
        return "audit_only", "archive", "interrupted/corrupt run retained only as audit archive", "high"
    if path.startswith("results/cleanup/"):
        return "audit_only", "keep", "cleanup audit output", "high"
    if path.startswith("results/") and (
        "run_registry.csv" in path or "derived_artifact_manifest.csv" in path or "validation" in path
    ):
        return "canonical_evidence", "keep", "registry/manifest/validation evidence", "high"
    if path.startswith("results/") and p.suffix.lower() in {".png", ".pdf", ".svg"}:
        return "canonical_figure", "keep", "result figure", "high" if ref_count else "medium"
    if path.startswith("results/"):
        if ref_count > 0 or run_id_for(path):
            return "canonical_result", "keep", "formal result referenced by evidence chain", "high"
        return "uncertain", "keep", "unreferenced result retained by default", "low"
    if path.startswith("build/"):
        return "reproducible_build_artifact", "keep", "current verified build directory retained", "medium"
    if path.startswith(".venv/"):
        return "reproducible_build_artifact", "keep", "project virtual environment retained", "medium"
    return "uncertain", "keep", "not confidently classifiable", "low"


def disk_usage_by_directory(files: list[Path]) -> list[dict[str, str | int]]:
    usage: Counter[str] = Counter()
    for f in files:
        r = rel(f)
        top = r.split("/", 2)[0] if "/" in r else "."
        if top == "results" and len(r.split("/")) > 1:
            top = "/".join(r.split("/")[:2])
        usage[top] += f.stat().st_size
    return [{"directory": k, "size_bytes": v} for k, v in sorted(usage.items())]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="results/cleanup")
    args = ap.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    files = [p for p in ROOT.rglob("*") if p.is_file()]
    refs, edges = referenced_paths(files)
    rows = []
    for f in files:
        r = rel(f)
        stat = f.stat()
        ref_by = sorted(refs.get(r, set()))
        classification, action, reason, confidence = classify(r, stat.st_size, len(ref_by))
        rows.append({
            "path": r,
            "size_bytes": stat.st_size,
            "sha256": sha256(f),
            "modified_time": stat.st_mtime,
            "file_type": file_type(f),
            "stage": stage_for(r),
            "run_id": run_id_for(r),
            "referenced_by": ";".join(ref_by),
            "reference_count": len(ref_by),
            "classification": classification,
            "proposed_action": action,
            "reason": reason,
            "confidence": confidence,
        })

    with (out / "full_inventory.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)
    with (out / "reference_graph.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["source", "target", "reference_type"])
        w.writeheader()
        w.writerows(edges)
    du = disk_usage_by_directory(files)
    with (out / "disk_usage_by_directory.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["directory", "size_bytes"])
        w.writeheader()
        w.writerows(du)
    print(f"Inventory complete: {len(files)} files, {len(edges)} references")


if __name__ == "__main__":
    main()
