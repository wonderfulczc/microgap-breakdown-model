#!/usr/bin/env python3
"""Automated acceptance gate for the frozen Stage 1 baseline."""
from __future__ import annotations
import csv, hashlib, json, re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
ERRORS = []
def fail(check, reason, fix): ERRORS.append((check, reason, fix))
def run(cmd, check):
    r = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    out = r.stdout + r.stderr
    if r.returncode: fail(check, out.strip()[-1200:], f"Run {' '.join(cmd)} and fix the failures.")
    return out
def counts():
    with (ROOT/'docs/source_registry.csv').open(encoding='utf-8-sig', newline='') as f: sources=sum(1 for _ in csv.DictReader(f))
    with (ROOT/'docs/parameter_registry.csv').open(encoding='utf-8-sig', newline='') as f: params=sum(1 for _ in csv.DictReader(f))
    equations=len(re.findall(r'^## EQ-\d{4}\b',(ROOT/'docs/equation_registry.md').read_text(encoding='utf-8'),re.M))
    unknowns=len(re.findall(r'^\| UNK-\d{4} \|',(ROOT/'docs/unknowns.md').read_text(encoding='utf-8'),re.M))
    for name,actual,minimum in [('sources',sources,10),('parameters',params,145),('equations',equations,35),('unknowns',unknowns,19)]:
        if actual<minimum: fail('registry counts',f'{name}: {actual} < {minimum}','Restore the frozen registry baseline.')
def manifest():
    p=ROOT/'docs/stage1_artifact_manifest.csv'
    if not p.exists(): fail('manifest','manifest missing','Generate it.'); return
    with p.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    if not rows: fail('manifest','manifest empty','Register all Stage 1 artifacts.')
    for row in rows:
        rel=row.get('path',''); target=ROOT/Path(rel)
        if not target.is_file(): fail('manifest',f'missing: {rel}','Restore it or re-freeze.'); continue
        if hashlib.sha256(target.read_bytes()).hexdigest()!=row.get('sha256'): fail('manifest',f'SHA256 mismatch: {rel}','Review change and deliberately regenerate manifest.')
def metrics():
    d=json.loads((ROOT/'results/stage_A2/summary.json').read_text(encoding='utf-8'))
    expected={'t_peak_s':7.758250938877745e-09,'I_peak_A_m':0.3693459484952068,'quad_max_relative_error':1.2375722305025981e-11,'fft_trusted_max_relative_error':0.0002226720457802179,'energy_consistency_relative_error':0.0}
    for k,v in expected.items():
        a=d.get(k)
        if a is None or abs(a-v)>max(1e-15,abs(v)*1e-10): fail('summary metrics',f'{k}={a!r}, expected {v!r}','Restore validated Stage 1.2 baseline.')
    m={x['metric']:x['value'] for x in d.get('I0_metrics',[])}
    if abs(m.get('spectrum_amplitude_ratio',0)-1.1)>1e-10 or abs(m.get('ESD_ratio',0)-1.21)>1e-10: fail('amplitude sensitivity','ratios changed','Restore 1.10 amplitude and 1.21 ESD ratios.')
def scope():
    prohibited={'isg0flux.cpp','poissonsolver.cpp','photoionizationsp3.cpp','streamersolver.cpp','isg0_solver.cpp','poisson_solver.cpp','sp3_solver.cpp','streamer_solver.cpp','case_i_solver.cpp','collision_solver.cpp'}
    # Later-stage implementations are allowed after Stage 1 closure. Scope leakage
    # means a prohibited solver is placed in a Stage 1-owned tree or result set.
    roots=[ROOT/'python/streamer_rf',ROOT/'python/reproduce_stage_A2.py',ROOT/'results/stage_A2']
    found=[p.relative_to(ROOT).as_posix() for base in roots for p in ([base] if base.is_file() else base.rglob('*')) if p.is_file() and p.name.lower() in prohibited]
    if found: fail('scope exclusions',f'prohibited solver files: {found}','Move them out of Stage 1.')
    report=(ROOT/'docs/stage1_closure_report.md').read_text(encoding='utf-8').lower()
    for term in ['isg-0','poisson','sp3','petsc/mpi','case i','figure 2','figure 3 inset','green collision fft','purple isolated-streamer','figure 4b']:
        if term not in report: fail('closure report exclusions',f'missing: {term}','Add exclusion and re-freeze.')
def main():
    out=run([sys.executable,'tools/validate_registry.py'],'registry validator')
    if '0 error(s), 0 warning(s)' not in out: fail('registry validator','not zero errors/warnings','Resolve registry findings.')
    out=run([sys.executable,'-m','pytest','tests/test_lifecycle.py','tests/test_fourier.py','tests/test_esd.py','tests/test_stage_A2_regression.py','-v'],'Stage 1.2 tests')
    passed=max([int(x) for x in re.findall(r'(\d+) passed',out)] or [0])
    if passed<24 or re.search(r'\b(?:failed|skipped|error)s?\b',out,re.I): fail('Stage 1.2 tests',f'required 24 clean passes; parsed {passed}','Fix failures/skips.')
    counts(); metrics(); scope(); manifest()
    if ERRORS:
        for c,r,f in ERRORS: print(f'FAIL [{c}]\n  reason: {r}\n  suggested fix: {f}')
        return 1
    print('Stage 1 closure validation passed.'); return 0
if __name__=='__main__': raise SystemExit(main())
