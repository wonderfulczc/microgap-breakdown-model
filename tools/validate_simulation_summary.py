#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / "docs/simulation_reconstruction_summary/index.html"
FIGSRC = ROOT / "docs/simulation_reconstruction_summary/data/figure_sources.csv"
MD = ROOT / "docs/simulation_toolchain_summary_for_project_book.md"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit(f"Simulation summary validation failed: {message}")


def registered_run_ids() -> set[str]:
    ids: set[str] = set()
    for reg in ROOT.glob("results/stage*/provenance/run_registry.csv"):
        df = pd.read_csv(reg)
        if "run_id" in df:
            ids.update(df["run_id"].astype(str))
    ids.add("DOCUMENTATION")
    return ids


def check_html_assets(text: str) -> None:
    require("http://" not in text and "https://" not in text and "//cdn" not in text, "external URL dependency found")
    for src in re.findall(r'(?:src|href)="([^"]+)"', text):
        if src.startswith("#") or src.startswith("data:"):
            continue
        p = HTML.parent / src
        require(p.exists(), f"missing HTML asset: {src}")


def check_figure_sources() -> None:
    require(FIGSRC.exists(), "missing figure_sources.csv")
    df = pd.read_csv(FIGSRC)
    ids = registered_run_ids()
    for row in df.itertuples():
        source = ROOT / row.source_file
        require(source.exists(), f"figure source missing: {row.source_file}")
        require(sha256(source) == row.sha256, f"figure source SHA mismatch: {row.source_file}")
        require(str(row.run_id) in ids, f"figure run_id not registered: {row.run_id}")
    figures = sorted((ROOT / "docs/simulation_reconstruction_summary/assets").glob("figure_*"))
    require(len(figures) >= 13, f"expected at least 13 figures, found {len(figures)}")
    recorded = set(df["figure"].astype(str))
    for fig in figures:
        require(fig.name in recorded, f"figure has no source record: {fig.name}")


def check_numbers(text: str) -> None:
    s4 = pd.read_csv(ROOT / "results/stage4/collision/collision_event.csv").iloc[0]
    s5 = pd.read_csv(ROOT / "results/stage5/trends/case_comparison.csv")
    high = s5[s5["run_id"] == "S5-HIGHFIELD-COLLISION"].iloc[0]
    require(f"{s4.t_collision_s * 1e9:.3f} ns" in text, "Stage 4 collision time mismatch or missing")
    require(f"{s4.field_drop_fraction * 100:.2f}%" in text, "Stage 4 field drop mismatch or missing")
    require(f"{high.t_collision * 1e9:.3f} ns" in text, "Stage 5 highfield collision time mismatch or missing")
    require(f"{high.spectral_centroid / 1e9:.3f} GHz" in text, "Stage 5 spectral centroid mismatch or missing")


def main() -> None:
    require(HTML.exists(), "missing HTML summary")
    require(MD.exists(), "missing Markdown project-book summary")
    text = HTML.read_text(encoding="utf-8")
    md = MD.read_text(encoding="utf-8")
    check_html_assets(text)
    check_figure_sources()
    forbidden_paths = ["invalidated_results", "interrupted_corrupt", "partial_hung", "quarantine"]
    lower = (text + "\n" + md).lower()
    require(not any(s in lower for s in forbidden_paths), "summary references invalidated/corrupt/quarantine content")
    forbidden_claims = ["1:1复刻完成", "1:1定量复刻完成", "shi 2019 figures 1–4 have been quantitatively reproduced"]
    require(not any(s in lower for s in forbidden_claims), "forbidden reproduction claim found")
    check_numbers(text + "\n" + md)
    require("sampling/Nyquist" in text and "strict_delta" in text and "event_local_proxy" in text, "required trust/proxy wording missing")
    print("Simulation summary validation passed.")


if __name__ == "__main__":
    main()
