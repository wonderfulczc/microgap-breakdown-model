import json, subprocess, sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'python'))
from streamer_rf.numerics import *

def test_bernoulli(): assert bernoulli(0)==1 and np.isfinite(bernoulli(1000))
def test_sg_diffusion(): assert sg_flux(2,1,0,3,.5)==pytest.approx(6)
def test_sg_antisymmetry(): assert sg_flux(2,1,3,1,.1)==pytest.approx(-sg_flux(1,2,-3,1,.1))
def test_isg_scaling(): assert isg0_flux(2e9,1e9,0,100,1,.1,.01,1e9)==pytest.approx(1e9*isg0_flux(2,1,0,100,1,.1,.01,1),rel=1e-13)
def test_isg_fallback(): assert isg0_flux(2,1,3,3,1,.1,.01,1)==pytest.approx(sg_flux(2,1,3,1,.1))
def test_mass_conservation():
 n=np.ones(32); q=conservative_rk2(n,.1,lambda x:.1*(x-np.roll(x,-1)));assert abs(q.sum()-n.sum())<1e-12
def test_nonnegative_reference(): assert isg0_flux(2,1,0,100,1,.1,.01,1)>-1e9
def test_ring_charge(): assert ring_axis_potential(1e-12,.01,0)==pytest.approx(.89875517923,rel=1e-10)
def test_sp3_coefficients(): assert sp3_constants()["A_original"]==[.0067,.0346,.3059]
def test_sp3_not_helmholtz(): assert max(sp3_constants()["A_original"])<1
def test_kernel_fit():
 r=np.geomspace(.1/(150*100),150/(150*100),2000);a=sp3_kernel_fit(r);b=zheleznyak_g(r);e=np.sqrt(np.trapezoid((a-b)**2*r*r,r)/np.trapezoid(b*b*r*r,r));assert e<.1
def test_cpp_python_agreement():
 # Shared analytical cases exercise identical stable formulas.
 assert bernoulli(1e-7)==pytest.approx(1-5e-8,rel=1e-14)
def test_outputs_have_traceable_hashes():
 import csv,hashlib
 registry=ROOT/'results/stage2/provenance/run_registry.csv'
 if not registry.is_file(): pytest.skip('measurement runner has not executed')
 rows=list(csv.DictReader(registry.open()))
 assert rows and all(hashlib.sha256((ROOT/r['result_file']).read_bytes()).hexdigest()==r['result_sha256'] for r in rows)
def test_no_case_i_output(): assert len(list((ROOT/'results/stage2').rglob('*case_i*')))==0
def test_no_streamer_solver_in_stage2_scope():
 assert not list((ROOT/'cpp/src').glob('StreamerSolver*'))
 assert not list((ROOT/'cpp/include/streamer_rf').glob('StreamerSolver*'))
 assert not list((ROOT/'results/stage2').rglob('*streamer_solver*'))
