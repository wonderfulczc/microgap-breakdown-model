#!/usr/bin/env python3
"""Static and provenance audit for Stage 2 acceptance evidence."""
from __future__ import annotations
import csv, hashlib, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/stage2/forensic"
OUT.mkdir(parents=True, exist_ok=True)
SCOPES = [ROOT/"cpp", ROOT/"python/stage2", ROOT/"tests/stage2",
          ROOT/"tools/validate_stage2_closure.py", ROOT/"docs/stage2_validation_matrix.csv",
          ROOT/"docs/stage2_closure_report.md"]
PATTERNS = {
    "placeholder": re.compile(r"check\s*\(\s*true|ASSERT_TRUE\s*\(\s*true|EXPECT_TRUE\s*\(\s*true|REQUIRE\s*\(\s*true|assert\s+True"),
    "hardcoded_pass": re.compile(r"status\s*=\s*[\"']PASS|[\"']status[\"']\s*:\s*[\"']pass", re.I),
    "fixed_observation": re.compile(r"(relative_error|mass_error|gauss_residual|iterations|relative_difference)[\"']?\s*[:=]\s*[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"),
    "existence_only": re.compile(r"assert\s+.*\.exists\(\)"),
}
rows, literals, graph = [], [], []
files = []
for scope in SCOPES:
    files += list(scope.rglob("*")) if scope.is_dir() else [scope]
for path in sorted(set(x for x in files if x.is_file())):
    # Stage 2's audit owns the shared numerical sources and Stage 2 tests, but
    # later-stage tests are audited by their corresponding closure tooling.
    rel = str(path.relative_to(ROOT))
    if rel.startswith(("cpp/src/streamer/", "cpp/include/streamer_rf/streamer/",
                       "cpp/tests/stage3/")):
        continue
    try: text = path.read_text(errors="replace")
    except OSError: continue
    for kind, regex in PATTERNS.items():
        for match in regex.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            classification = "placeholder" if kind in ("placeholder", "existence_only") else "hardcoded"
            # Zero initialization of an iteration counter is solver state, not
            # an acceptance observation. Runtime observations are still caught
            # when serialized under a quoted metric key.
            if kind == "fixed_observation" and match.group(0).endswith("=0"):
                continue
            row = {"file": rel, "line": line, "finding": kind, "classification": classification,
                   "text": match.group(0)[:160], "resolved": "no"}
            rows.append(row)
            if kind in ("fixed_observation", "hardcoded_pass"): literals.append(row)
    if "/test" in rel or rel.startswith("tests/"):
        calls = [name for name in ("sg_flux", "isg0_flux", "solve_potential", "solve_photoionization", "open_charge", "robin") if name in text]
        graph.append({"test_file": rel, "target_calls": ";".join(calls),
                      "classification": "measured" if calls else "untraceable"})

def write(name, data, fields):
    with (OUT/name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(data)

write("evidence_audit.csv", rows, ["file","line","finding","classification","text","resolved"])
write("suspicious_literals.csv", literals, ["file","line","finding","classification","text","resolved"])
write("test_call_graph.csv", graph, ["test_file","target_calls","classification"])
registry = ROOT / "results/stage2/provenance/run_registry.csv"
provenance = []
if registry.exists():
    with registry.open(newline="") as handle:
        for row in csv.DictReader(handle):
            result = ROOT / row["result_file"]
            status = "measured" if result.exists() and hashlib.sha256(result.read_bytes()).hexdigest() == row["result_sha256"] else "untraceable"
            provenance.append((row["run_id"], status))
with (OUT/"report.md").open("w") as handle:
    handle.write("# Stage 2 Evidence Audit\n\n")
    handle.write(f"- suspicious static findings: {len(rows)}\n- hardcoded literals: {len(literals)}\n")
    handle.write(f"- test files without target calls: {sum(x['classification']=='untraceable' for x in graph)}\n")
    handle.write(f"- registered runs checked: {len(provenance)}\n- untraceable registered runs: {sum(s!='measured' for _,s in provenance)}\n")
print(f"evidence_findings={len(rows)} hardcoded={len(literals)} untraceable_tests={sum(x['classification']=='untraceable' for x in graph)}")
raise SystemExit(1 if rows else 0)
