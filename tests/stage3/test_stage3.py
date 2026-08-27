from pathlib import Path
import sys
import hashlib
import numpy as np,pandas as pd,pytest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from streamer_rf.streamer.morrow_lowke_reference import evaluate
from streamer_rf.streamer.diagnostics import axis_heads,trajectory_metrics
EVIDENCE_REASON="EVIDENCE_UNAVAILABLE: historical production artifact intentionally removed during validated cleanup; see docs/github_archive_assessment.md and docs/github_archive_deletion_manifest.csv"
def require_evidence(path):
 if not path.exists(): pytest.skip(EVIDENCE_REASON)
def test_reference_nonnegative():
 N=101325/(1.380649e-23*300)
 for E in np.geomspace(1e4,4e7,100):
  q=evaluate(E,N);assert all(np.isfinite(list(q.values())));assert all(v>=0 for v in q.values())
def test_reference_junctions_bounded():
 N=101325/(1.380649e-23*300)
 for x in [2.6e-17,1e-16,1.05e-15,1.5e-15,2e-15]:
  E=x*N/1e4;l=evaluate(E*(1-1e-9),N);r=evaluate(E*(1+1e-9),N);assert abs(l['mobility']-r['mobility'])/max(l['mobility'],r['mobility'])<1e-7
@pytest.mark.evidence
def test_measured_cpp_python_agreement():
 p=ROOT/'results/stage3/transport/cpp_python_comparison.csv';require_evidence(p);d=pd.read_csv(p);assert d.relative_difference.max()<1e-11
@pytest.mark.evidence
def test_provenance_hashes():
 p=ROOT/'results/stage3/provenance/run_registry.csv';require_evidence(p);d=pd.read_csv(p);assert len(d)>0
 for r in d.itertuples():assert hashlib.sha256((ROOT/r.result_files).read_bytes()).hexdigest()==r.result_sha256
def test_no_forbidden_results():
 names=' '.join(str(p).lower() for p in (ROOT/'results/stage3').rglob('*'));assert 'case_i' not in names and 'collision' not in names and 'current_moment' not in names
def test_independent_head_identification():
 z=np.arange(20.)*1e-5;rows=[]
 for j,x in enumerate(z):rows.append(dict(i=0,j=j,r_m=5e-6,z_m=x,E_V_m=1+10*np.exp(-((j-4)/1.2)**2)+9*np.exp(-((j-15)/1.2)**2),rho_C_m_3=-np.exp(-((j-4)/1.2)**2)+np.exp(-((j-15)/1.2)**2)))
 h=axis_heads(pd.DataFrame(rows));assert h['lower']['charge']<0<h['upper']['charge'];assert h['lower']['charge_z']==pytest.approx(4e-5);assert h['upper']['charge_z']==pytest.approx(15e-5)
def test_sustained_double_head_criterion():
 t=np.arange(6.)*1e-10;d=pd.DataFrame(dict(time_s=t,lower_head_z_m=5e-3-np.arange(6)*2e-5,upper_head_z_m=5e-3+np.arange(6)*2e-5,lower_head_rho_C_m_3=-np.ones(6),upper_head_rho_C_m_3=np.ones(6)));q=trajectory_metrics(d,2e-5);assert q['formed'] and q['lower_velocity_m_s']<0<q['upper_velocity_m_s']
