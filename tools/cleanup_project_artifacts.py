#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_inventory() -> list[dict[str, str]]:
    inv = ROOT / "results/cleanup/full_inventory.csv"
    if not inv.exists():
        subprocess.run([sys.executable, "tools/inventory_project_artifacts.py"], cwd=ROOT, check=True)
    with inv.open(newline="") as fp:
        return list(csv.DictReader(fp))


def candidate_dirs() -> dict[str, dict[str, str]]:
    dirs: dict[str, dict[str, str]] = {}
    stage2 = ROOT / "results/stage2/forensic/invalidated_results"
    if stage2.exists():
        dirs[rel(stage2)] = {
            "archive": "stage2_invalidated_results.tar.zst",
            "reason": "Stage 2 invalidated forensic results; not valid final evidence.",
            "replacement": "Stage 2 closure evidence and forensic summary",
        }
    for base, archive_name, replacement in [
        (ROOT / "results/stage4/runs", "stage4_interrupted_runs.tar.zst", "S4-RIGHT-ISOLATED and other formal Stage 4 runs"),
        (ROOT / "results/stage5/runs", "stage5_interrupted_runs.tar.zst", "S5 formal reruns and resource-aware terminations"),
    ]:
        if base.exists():
            matches = []
            for p in base.iterdir():
                low = p.name.lower()
                if p.is_dir() and any(tok in low for tok in ["interrupted", "corrupt", "partial_hung", "paused", "failed", "aborted"]):
                    matches.append(p)
            if matches:
                pseudo = f"{rel(base)}::" + ";".join(rel(p) for p in matches)
                dirs[pseudo] = {
                    "archive": archive_name,
                    "reason": "Interrupted/corrupt/paused run directories replaced by formal evidence.",
                    "replacement": replacement,
                }
    return dirs


def tar_zst(archive: Path, paths: list[Path]) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "w:zst") as tf:
        for p in paths:
            tf.add(p, arcname=p.relative_to(ROOT))


def write_archive_readme(stem: Path, original_dirs: list[str], reason: str, replacement: str) -> None:
    text = f"""# Audit-only archive: {stem.name}

Original directory/directories:

{chr(10).join(f'- `{d}`' for d in original_dirs)}

Reason:

{reason}

Replacement formal evidence:

{replacement}

Archive date: {dt.datetime.now(dt.UTC).isoformat()}

These files must not be used for final scientific results. They are retained only for audit traceability.

To restore:

```bash
tar --zstd -xf {stem.name}
```
"""
    stem.with_suffix(stem.suffix + ".README.md").write_text(text)


def write_dir_manifest(stem: Path, paths: list[Path]) -> None:
    rows = []
    for root in paths:
        for f in root.rglob("*"):
            if f.is_file():
                rows.append({
                    "original_path": rel(f),
                    "size_bytes": f.stat().st_size,
                    "sha256": sha(f),
                })
    with stem.with_suffix(stem.suffix + ".manifest.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["original_path", "size_bytes", "sha256"])
        w.writeheader()
        w.writerows(rows)


def build_candidates(inventory: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates = []
    for row in inventory:
        path = row["path"]
        cls = row["classification"]
        action = row["proposed_action"]
        ref_count = int(row.get("reference_count") or 0)
        if cls in {"cache", "temporary"} and ref_count == 0:
            candidates.append({**row, "deletion_reason": row["reason"], "final_action": "quarantine"})
    return candidates


def archive_audit_only(apply: bool, manifest_rows: list[dict[str, str]]) -> None:
    archive_root = ROOT / "archive/audit_only"
    archive_root.mkdir(parents=True, exist_ok=True)
    archive_manifest = ROOT / "results/cleanup/archive_manifest.csv"
    archive_rows = []
    for key, meta in candidate_dirs().items():
        if "::" in key:
            original_dirs = key.split("::", 1)[1].split(";")
            paths = [ROOT / d for d in original_dirs if (ROOT / d).exists()]
        else:
            original_dirs = [key]
            paths = [ROOT / key] if (ROOT / key).exists() else []
        if not paths:
            continue
        archive = archive_root / meta["archive"]
        if apply:
            tar_zst(archive, paths)
            write_dir_manifest(archive, paths)
            write_archive_readme(archive, original_dirs, meta["reason"], meta["replacement"])
            digest = sha(archive)
            archive.with_suffix(archive.suffix + ".sha256").write_text(f"{digest}  {archive.name}\n")
        else:
            digest = ""
        total_size = sum(f.stat().st_size for p in paths for f in p.rglob("*") if f.is_file())
        archive_rows.append({
            "archive": rel(archive),
            "original_dirs": ";".join(original_dirs),
            "size_bytes": total_size,
            "sha256": digest,
            "reason": meta["reason"],
            "replacement": meta["replacement"],
            "created": str(bool(apply)),
        })
        for p in paths:
            manifest_rows.append({
                "original_path": rel(p),
                "classification": "audit_only",
                "size_bytes": str(total_size),
                "sha256": "",
                "reference_count": "0",
                "replacement_path": rel(archive),
                "deletion_reason": meta["reason"],
                "validation_status": "pending",
                "final_action": "quarantine_directory" if apply else "dry_run_archive_directory",
            })
    with archive_manifest.open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=["archive", "original_dirs", "size_bytes", "sha256", "reason", "replacement", "created"])
        w.writeheader()
        w.writerows(archive_rows)


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["original_path", "classification", "size_bytes", "sha256", "reference_count", "replacement_path", "deletion_reason", "validation_status", "final_action"]
    with path.open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def quarantine_file(path: Path, qroot: Path) -> Path:
    target = qroot / rel(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(target))
    return target


def restore_quarantine(qroot: Path) -> None:
    if not qroot.exists():
        return
    for p in sorted([x for x in qroot.rglob("*") if x.is_file()], key=lambda x: len(x.parts), reverse=True):
        dest = ROOT / p.relative_to(qroot)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(p), str(dest))


def apply_quarantine(candidates: list[dict[str, str]], qroot: Path, archive_only: bool) -> None:
    qroot.mkdir(parents=True, exist_ok=True)
    for row in candidates:
        p = ROOT / row["path"]
        if p.exists() and p.is_file():
            target = quarantine_file(p, qroot)
            row["replacement_path"] = rel(target)
            row["final_action"] = "quarantined"
    if archive_only:
        for key in candidate_dirs():
            if "::" in key:
                dirs = key.split("::", 1)[1].split(";")
            else:
                dirs = [key]
            for d in dirs:
                p = ROOT / d
                if p.exists() and p.is_dir():
                    target = qroot / d
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(p), str(target))


def purge(qroot: Path) -> None:
    if not qroot.exists():
        return
    shutil.rmtree(qroot)


def retained_uncertain(inventory: list[dict[str, str]]) -> None:
    rows = [r for r in inventory if r["classification"] == "uncertain"]
    out = ROOT / "results/cleanup/retained_uncertain.csv"
    if not rows:
        rows = []
    with out.open("w", newline="") as fp:
        fields = ["path", "size_bytes", "sha256", "reason", "confidence"]
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--quarantine-dir", default="results/cleanup/quarantine")
    ap.add_argument("--purge-quarantine", action="store_true")
    ap.add_argument("--manifest", default="results/cleanup/deletion_manifest.csv")
    ap.add_argument("--keep-build", action="store_true")
    ap.add_argument("--archive-audit-only", action="store_true")
    args = ap.parse_args()

    if args.purge_quarantine:
        purge(ROOT / args.quarantine_dir)
        print(f"Purged {args.quarantine_dir}")
        return

    if not args.apply:
        args.dry_run = True
    inventory = read_inventory()
    retained_uncertain(inventory)
    candidates = build_candidates(inventory)
    rows = []
    for c in candidates:
        rows.append({
            "original_path": c["path"],
            "classification": c["classification"],
            "size_bytes": c["size_bytes"],
            "sha256": c["sha256"],
            "reference_count": c["reference_count"],
            "replacement_path": "",
            "deletion_reason": c["deletion_reason"],
            "validation_status": "pending",
            "final_action": "dry_run_quarantine" if args.dry_run else "pending_quarantine",
        })
    if args.archive_audit_only:
        archive_audit_only(args.apply, rows)
    if args.apply:
        apply_quarantine(candidates, ROOT / args.quarantine_dir, args.archive_audit_only)
        # Refresh only cache/temporary manifest rows. Keep audit-only rows and
        # archive manifest produced above intact.
        refreshed = []
        cache_by_path = {c["path"]: c for c in candidates}
        for r in rows:
            c = cache_by_path.get(r["original_path"])
            if c is not None:
                r["replacement_path"] = c.get("replacement_path", "")
                r["final_action"] = c.get("final_action", "quarantined")
            refreshed.append(r)
        rows = refreshed
    write_manifest(ROOT / args.manifest, rows)
    print(f"{'Dry run' if args.dry_run else 'Apply'} complete: {len(rows)} manifest rows")


if __name__ == "__main__":
    main()
