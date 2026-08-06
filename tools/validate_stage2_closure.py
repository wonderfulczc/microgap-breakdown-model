#!/usr/bin/env python3
"""Close Stage 2 only when measured evidence and its provenance are intact."""
from __future__ import annotations
import csv,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; failures=[]
def check(ok,msg):
 if not ok: failures.append(msg)
def run(cmd):
 p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True);check(p.returncode==0,f"command failed: {' '.join(cmd)}\n{p.stdout}\n{p.stderr}");return p
run(['.venv/bin/python','tools/validate_registry.py'])
run(['.venv/bin/python','tools/validate_stage1_closure.py'])
run(['cmake','-S','.','-B','build','-G','Ninja','-DCMAKE_BUILD_TYPE=Release'])
run(['cmake','--build','build','--parallel'])
run(['ctest','--test-dir','build','--output-on-failure'])
run(['.venv/bin/python','-m','pytest','tests/stage2','-q'])
audit=run(['.venv/bin/python','tools/audit_stage2_evidence.py'])
summary_path=ROOT/'results/stage2/integration/stage2_summary.json';check(summary_path.exists(),'missing measured summary')
if summary_path.exists():
 s=json.loads(summary_path.read_text());w1=s['wp2_1'];w2=s['wp2_2'];w3=s['wp2_3']
 for ok,msg in [(w1['cpp_python']<1e-11,'WP2.1 C++/Python'),(w1['mass']<1e-10,'WP2.1 mass'),(w1['strong_isg']<.05 and w1['strong_isg']<w1['strong_sg'],'WP2.1 strong gradient'),(w2['minimum_order']>=1.8,'Poisson order'),(w2['axis_ring_error']<.005,'axis ring'),(w2['off_axis_ring_error']<.01,'off-axis ring'),(w2['gaussian_error']<.01,'Gaussian boundary'),(w2['gauss_residual']<1e-8,'Gauss law'),(w2['boundary_update']>0,'OpenCharge update'),(w3['minimum_robin_order']>=1.8,'Robin order'),(w3['max_final_residual']<1e-8,'fixed point'),(w3['max_sp3_error']<.1,'SP3 integral'),(s['mpi_max_difference']<1e-11,'MPI')]:check(ok,msg)
registry_path=ROOT/'results/stage2/provenance/run_registry.csv';check(registry_path.exists(),'missing run registry');run_ids=set()
if registry_path.exists():
 rows=list(csv.DictReader(registry_path.open()));run_ids={r['run_id'] for r in rows}
 for r in rows:
  for key in ('timestamp_utc','git_commit_or_worktree_hash','executable_path','executable_sha256','config_path','config_sha256','mpi_ranks','petsc_version','command_line','stdout_log','stderr_log','exit_code','result_file','result_sha256'):check(bool(r[key]),f"empty registry field {key} in {r['run_id']}")
  result=ROOT/r['result_file'];check(result.exists(),f"missing registered result {result}")
  if result.exists():check(hashlib.sha256(result.read_bytes()).hexdigest()==r['result_sha256'],f"result hash mismatch {result}")
  exe=ROOT/r['executable_path'];check(exe.exists() and hashlib.sha256(exe.read_bytes()).hexdigest()==r['executable_sha256'],f"executable hash mismatch {exe}")
  config=ROOT/r['config_path'];check(config.exists() and hashlib.sha256(config.read_bytes()).hexdigest()==r['config_sha256'],f"config hash mismatch {config}")
  check((ROOT/r['stdout_log']).exists() and (ROOT/r['stderr_log']).exists(),f"missing logs {r['run_id']}")
 expected={f'S2-{solver}-MPI-{n:02d}' for solver in ('POISSON','SP3') for n in (1,2,4)};check(expected<=run_ids,'independent 1/2/4-rank runs missing')
matrix_path=ROOT/'docs/stage2_validation_matrix.csv';check(matrix_path.exists(),'missing validation matrix')
if matrix_path.exists():
 for r in csv.DictReader(matrix_path.open()):
  check(not(r['status']=='PASS' and r['evidence_status']=='invalid'),f"invalid PASS {r['validation_id']}")
  check(bool(r['run_id']) and bool(r['result_sha256']),f"untraceable matrix row {r['validation_id']}")
source=(ROOT/'cpp/src/poisson.cpp').read_text()+(ROOT/'cpp/src/sp3.cpp').read_text();check('MPI_Allreduce' in source and 'RingQuadratureMode::CellCenterRing' in source and 'open_charge_boundary' in source,'OpenCharge implementation absent');check('beta2' in source and 'final_residual' in source and 'solve_axisymmetric_robin' in source,'coupled Robin implementation absent')
check('SUPERSEDED AND INVALIDATED' in (ROOT/'docs/stage2_closure_report.md').read_text(),'old report not invalidated')
check('STAGE 3 NOT STARTED' in (ROOT/'docs/stage3_closure_report.md').read_text(),'Stage 3 status not rolled back')
for pattern in ('*case_i*','*collision*','*current_moment*','*radiation*'):
 check(not any((ROOT/'results/stage2').rglob(pattern)),f'forbidden Stage 2 output: {pattern}')
if failures:
 print('\n'.join('FAIL: '+x for x in failures));raise SystemExit(1)
print('Stage 2 closure validation passed with traceable measured evidence.')
