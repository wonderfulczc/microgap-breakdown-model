#!/usr/bin/env python3
from __future__ import annotations
import csv,hashlib,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage3/forensic';OUT.mkdir(parents=True,exist_ok=True)
patterns={'placeholder':re.compile(r'(?:check|REQUIRE|ASSERT_TRUE|EXPECT_TRUE)\s*\(\s*true|assert\s+True'),'fixed_pass':re.compile(r'["\']status["\']\s*:\s*["\']PASS',re.I),'hardcoded_observed':re.compile(r'(?:observed|E_max|channel_density|head_velocity|convergence_error)\s*[=:]\s*[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?',re.I)}
files=[p for base in (ROOT/'cpp',ROOT/'python/stage3',ROOT/'tests/stage3',ROOT/'tools/validate_stage3_closure.py') for p in ([base] if base.is_file() else base.rglob('*')) if p.is_file()];rows=[]
for p in files:
 try:text=p.read_text(errors='replace')
 except OSError:continue
 for kind,rx in patterns.items():
  for m in rx.finditer(text):rows.append(dict(file=str(p.relative_to(ROOT)),line=text.count('\n',0,m.start())+1,classification=kind,text=m.group(0)[:120]))
reg=ROOT/'results/stage3/provenance/run_registry.csv';bad_hash=0
if reg.exists():
 for r in csv.DictReader(reg.open()):
  p=ROOT/r['result_files'];bad_hash+=not(p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==r['result_sha256'])
with (OUT/'evidence_audit.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['file','line','classification','text']);w.writeheader();w.writerows(rows)
(OUT/'report.md').write_text(f'# Stage 3 Evidence Audit\n\n- unresolved static findings: {len(rows)}\n- bad registered hashes: {bad_hash}\n')
print(f'evidence_findings={len(rows)} bad_registered_hashes={bad_hash}');raise SystemExit(1 if rows or bad_hash else 0)
