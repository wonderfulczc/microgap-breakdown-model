#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_UNAVAILABLE: list[str] = []
DELETED_PATHS: set[str] | None = None


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit(f"Project cleanup integrity validation failed: {message}")


def record_unavailable_evidence(rel: str) -> None:
    EVIDENCE_UNAVAILABLE.append(
        "EVIDENCE_UNAVAILABLE: historical production artifact intentionally "
        f"removed during validated cleanup: {rel}; see "
        "docs/github_archive_assessment.md and docs/github_archive_deletion_manifest.csv"
    )


def deleted_paths() -> set[str]:
    global DELETED_PATHS
    if DELETED_PATHS is None:
        manifest = ROOT / "docs/github_archive_deletion_manifest.csv"
        if manifest.exists():
            DELETED_PATHS = set(pd.read_csv(manifest)["path"].astype(str))
        else:
            DELETED_PATHS = set()
    return DELETED_PATHS


def check_registries() -> None:
    for reg in ROOT.glob("results/stage*/provenance/run_registry.csv"):
        df = pd.read_csv(reg)
        for row in df.itertuples():
            fields = row._asdict()
            result_field = fields.get("result_files") or fields.get("result_file")
            if not result_field or not isinstance(result_field, str):
                continue
            result_paths = [part.strip() for part in re.split(r"[;|]", result_field) if part.strip()]
            sha_field = str(fields.get("result_sha256", "") or "")
            sha_values = [part.strip() for part in re.split(r"[;|]", sha_field) if part.strip()]
            for i, rel in enumerate(result_paths):
                p = ROOT / rel
                require(p.exists(), f"registered result missing: {rel}")
                if len(sha_values) == len(result_paths):
                    require(sha(p) == sha_values[i], f"registry sha mismatch: {rel}")
                elif len(result_paths) == 1 and sha_field:
                    require(sha(p) == sha_field, f"registry sha mismatch: {rel}")


def check_manifests() -> None:
    for mf in ROOT.glob("results/stage*/provenance/derived_artifact_manifest.csv"):
        df = pd.read_csv(mf)
        for row in df.itertuples():
            file_path = getattr(row, "file", "")
            digest = getattr(row, "sha256", "")
            if not file_path or not isinstance(file_path, str):
                continue
            p = ROOT / file_path
            require(p.exists(), f"derived artifact missing: {file_path}")
            require(sha(p) == digest, f"derived artifact sha mismatch: {file_path}")


def check_html_summary() -> None:
    html = ROOT / "docs/simulation_reconstruction_summary/index.html"
    if not html.exists():
        return
    text = html.read_text(errors="ignore")
    require("http://" not in text and "https://" not in text and "//cdn" not in text, "HTML contains external URL dependency")
    for src in re.findall(r'(?:src|href)="([^"]+)"', text):
        if src.startswith("#") or src.startswith("data:"):
            continue
        p = html.parent / src
        require(p.exists(), f"HTML referenced asset missing: {src}")
    source_csv = ROOT / "docs/simulation_reconstruction_summary/data/figure_sources.csv"
    require(source_csv.exists(), "missing figure_sources.csv")
    df = pd.read_csv(source_csv)
    for row in df.itertuples():
        p = ROOT / row.source_file
        if not p.exists() and str(row.source_file).startswith("results/") and (
            not (ROOT / "results").exists() or str(row.source_file) in deleted_paths() or "results" in deleted_paths()
        ):
            record_unavailable_evidence(row.source_file)
            continue
        require(p.exists(), f"figure source missing: {row.source_file}")
        require(sha(p) == row.sha256, f"figure source sha mismatch: {row.source_file}")
    bad = ["invalidated_results", "interrupted_corrupt", "partial_hung", "quarantine"]
    lower = text.lower()
    require(not any(b in lower for b in bad), "HTML references invalidated/corrupt/quarantine content")
    require("1:1复刻完成" not in text and "1:1定量复刻完成" not in text, "HTML contains forbidden reproduction claim")


def check_imports_and_smoke() -> None:
    py = ROOT / ".venv/bin/python"
    code = "import streamer_rf.stage4, streamer_rf.stage5; print('imports ok')"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "python")
    p = subprocess.run([str(py), "-c", code], cwd=ROOT, text=True, capture_output=True, env=env)
    require(p.returncode == 0, f"Python import smoke failed: {p.stderr}")
    exe = ROOT / "build/bin/stage5_cpp_tests"
    if exe.exists():
        p = subprocess.run([str(exe)], cwd=ROOT, text=True, capture_output=True, timeout=120)
        require(p.returncode == 0, f"C++ smoke failed: {p.stdout}\n{p.stderr}")


def main() -> None:
    check_registries()
    check_manifests()
    check_html_summary()
    check_imports_and_smoke()
    for item in EVIDENCE_UNAVAILABLE:
        print(item)
    print("Project cleanup integrity validation passed.")


if __name__ == "__main__":
    main()
