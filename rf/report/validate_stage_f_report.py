#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "rf/report"


def main() -> None:
    required = [
        REPORT / "stage_f_smoke_paper_report.html",
        REPORT / "figure_manifest.yaml",
        REPORT / "data_manifest.csv",
        REPORT / "summary_matrix.csv",
    ]
    for stem in (
        "Fig_F1_rf_trustworthiness",
        "Fig_F2_jefimenko_validation",
        "Fig_F3_stage_time_frequency",
        "Fig_F4_stage_frequency_fingerprint",
    ):
        required.extend([REPORT / "figures" / f"{stem}.png", REPORT / "figures" / f"{stem}.pdf"])
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise SystemExit(f"missing report artifacts: {missing}")

    html = (REPORT / "stage_f_smoke_paper_report.html").read_text()
    for src in re.findall(r'<img src="([^"]+)"', html):
        target = REPORT / src
        if not target.exists():
            raise SystemExit(f"broken html image link: {src}")
    if "FIG_F5_STATUS = DATA_INSUFFICIENT" not in html:
        raise SystemExit("Fig F5 insufficient-data status is missing")

    stage = pd.read_csv(ROOT / "rf/production/f4_attribution/stage_rf_summary.csv")
    for _, row in stage.iterrows():
        for status_col, energy_col in (
            ("VHF_status", "VHF_energy"),
            ("UHF_status", "UHF_energy"),
            ("1_3GHz_status", "1-3 GHz_energy"),
        ):
            if row[status_col] == "UNTRUSTED" and pd.notna(row[energy_col]):
                raise SystemExit(f"untrusted band has numeric energy: {row['dataset']} {row['stage']} {energy_col}")
    matrix = pd.read_csv(REPORT / "summary_matrix.csv")
    if "ABSENT" in matrix.to_string():
        raise SystemExit("summary matrix must not use ABSENT")
    print("Stage F report validation PASS")


if __name__ == "__main__":
    main()
