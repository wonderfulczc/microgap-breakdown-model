from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

from streamer_rf.circuit import PercolationConfig, analyze_cr4_handoff_run


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate G0 cold-thermal handoff diagnostics from C-R4 output.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", default=Path("circuit/g0_handoff"), type=Path)
    parser.add_argument("--sigma-relative-threshold", default=0.1, type=float)
    args = parser.parse_args()

    t0 = time.perf_counter()
    summary = analyze_cr4_handoff_run(
        args.run_dir,
        args.out_dir,
        PercolationConfig(sigma_relative_threshold=args.sigma_relative_threshold),
    )
    summary["generator_runtime_s"] = time.perf_counter() - t0
    summary["generator_peak_rss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    summary_path = args.out_dir / "g0_handoff_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
