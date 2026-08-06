#!/usr/bin/env python3
from __future__ import annotations
import csv,hashlib,json,subprocess
from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];failures=[]
def check(x,m):
 if not x:failures.append(m)
def run(cmd):
 p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True);check(p.returncode==0,f"command failed: {' '.join(cmd)}\n{p.stdout}\n{p.stderr}");return p
run(['.venv/bin/python','tools/validate_stage1_closure.py']);s2=run(['.venv/bin/python','tools/validate_stage2_closure.py']);check('Stage 2 closure validation passed with traceable measured evidence.' in s2.stdout,'Stage 2 anti-fraud gate text missing');run(['cmake','-S','.','-B','build','-G','Ninja','-DCMAKE_BUILD_TYPE=Release']);run(['cmake','--build','build','--parallel']);run(['ctest','--test-dir','build','--output-on-failure']);run(['.venv/bin/python','-m','pytest','tests/stage3','-q']);run(['.venv/bin/python','tools/audit_stage3_evidence.py'])
formula=(ROOT/'docs/stage3_formula_audit.md').read_text();check('Morrow–Lowke 1997 PDF p.14' in formula and 'A1–A11' in formula,'Morrow-Lowke source/formulas not audited')
required={'coarse_ml','baseline_ml','fine_ml','coarse_paper','coarse_no_sp3',*[f'nref_{x}' for x in ('1e-14','1e-12','1e-10','1e-8','1e-6')],*[f'dt_{x}' for x in ('0.5','1.0','1.5')],*[f'eta_{x}' for x in ('0','1e-4','1e-3')],'domain_r1.5',*[f'mpi_{x}' for x in ('1','2','4')]}
regp=ROOT/'results/stage3/provenance/run_registry.csv';check(regp.exists(),'run registry missing');run_ids=set()
if regp.exists():
 rows=list(csv.DictReader(regp.open()));run_ids={r['run_id'].removeprefix('S3-').lower() for r in rows}
 for r in rows:
  for k in ('run_id','timestamp_utc','executable','executable_sha256','config','config_sha256','mpi_ranks','command','stdout','stderr','result_files','result_sha256'):check(bool(r.get(k)),f"empty provenance field {k}")
  p=ROOT/r['result_files'];check(p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==r['result_sha256'],f"result hash mismatch {p}")
  e=ROOT/r['executable'];check(e.exists() and hashlib.sha256(e.read_bytes()).hexdigest()==r['executable_sha256'],f"executable hash mismatch {e}")
check(all(any(name.replace('.','_') in rid for rid in run_ids) for name in required),'one or more required independent runs missing')
summaryp=ROOT/'results/stage3/closure/run_summary.csv';check(summaryp.exists(),'measured run summary missing')
if summaryp.exists():
 d=pd.read_csv(summaryp);names=set(d.run);check(required<=names,'required analyzed runs missing');core=d[d.run.isin(['coarse_ml','baseline_ml','fine_ml'])];check(len(core)==3 and core.formed.all(),'double-headed streamer criteria failed');check(core.head_definitions_agree.all(),'head definitions differ by more than two cells');check((core.photoionization_ahead>0).all(),'positive-head photoionization missing');check(core.max_conservation_residual.max()<1e-6,'conservation residual threshold failed');check(np.isfinite(core.select_dtypes('number')).all().all(),'nonfinite core metric')
 nr=d[d.run.str.startswith('nref_')].sort_values('run');stable=False
 for i in range(max(0,len(nr)-2)):
  q=nr.iloc[i:i+3];cols=['E_max_V_m','channel_density_m_3','lower_head_z_m','upper_head_z_m','total_electrons','minimum_dt_s'];stable|=all((q[c].max()-q[c].min())/max(abs(q[c]).max(),1e-300)<.01 for c in cols)
 check(stable,'no three-point n_ref stable interval');photo=d[d.run.isin(['coarse_ml','coarse_no_sp3'])];check(len(photo)==2 and photo.iloc[photo.run.tolist().index('coarse_no_sp3')].upper_head_z_m<photo.iloc[photo.run.tolist().index('coarse_ml')].upper_head_z_m,'SP3-off positive-head response not weaker')
 mpi=d[d.run.str.startswith('mpi_')];
 if len(mpi)==3:
  base=mpi.iloc[0];check(max(abs(mpi.E_max_V_m-base.E_max_V_m)/abs(base.E_max_V_m))<1e-9,'MPI Emax mismatch');check(max(abs(mpi.total_electrons-base.total_electrons)/abs(base.total_electrons))<1e-9,'MPI electron mismatch')
 else:check(False,'MPI 1/2/4 results missing')
for pattern in ('*case_i*','*collision*','*current_moment*','*fft*','*esd*','*radiation*'):check(not any((ROOT/'results/stage3').rglob(pattern)),f'forbidden output {pattern}')
for p in ('docs/stage3_validation_matrix.csv','docs/stage3_closure_report_v2.md','docs/stage4_entry_plan.md'):check((ROOT/p).exists(),f'missing {p}')
if failures:print('\n'.join('FAIL: '+x for x in failures));raise SystemExit(1)
print('Stage 3 closure validation passed with traceable measured evidence.')
