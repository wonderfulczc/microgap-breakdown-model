from pathlib import Path
import sys
import json
from copy import deepcopy
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from streamer_rf.fullwave.foundation import *


def geometry(): return json.loads((ROOT/'fullwave/h1/h1_reference_geometry.json').read_text())


def test_backend_pin_and_runtime_status_are_independent():
    r=json.loads((ROOT/'fullwave/h1/openems_backend.json').read_text())
    assert validate_backend(r)
    # Valid source hashes alone must never promote an untested installation.
    r['installation_status']='BUILD_BLOCKED_MISSING_DEPENDENCIES'
    assert not validate_backend(r)
    r['commits']['openEMS']='missing'
    with pytest.raises(ValueError): validate_backend(r)


def test_reference_geometry_and_port_polarity():
    g=geometry(); assert validate_geometry(g)
    g['positive_current_direction']='GAP_TO_EXTERNAL'
    with pytest.raises(ValueError): validate_geometry(g)


def test_production_pending():
    g=geometry();g['kind']='PRODUCTION_GEOMETRY'
    with pytest.raises(ValueError,match='PENDING_STAGE_B'):validate_geometry(g)


@pytest.mark.parametrize('partition',['GEOMETRY_RESOLVED_CGAP','BOTH','UNKNOWN'])
def test_partition_exclusive(partition):
    g=geometry();g['cgap_partition']=partition
    with pytest.raises(ValueError):validate_geometry(g)


def test_no_duplicate_capacitance():
    g=geometry();g['add_duplicate_lumped_cgap']=True
    with pytest.raises(ValueError):validate_geometry(g)


@pytest.mark.parametrize('extent',[[np.nan,1,1],[0,1,1],[-1,1,1],[1,2]])
def test_invalid_geometry(extent):
    g=geometry();g['domain_extent_m']=extent
    with pytest.raises(ValueError):validate_geometry(g)


def test_wavelength_SI_conversion():
    assert wavelength_cell(1e9)==pytest.approx(.299792458/20)
    assert wavelength_cell(1e9,20,4)==pytest.approx(wavelength_cell(1e9)/2)


def test_graded_mesh():
    a=graded_axis(.06,[.02,.03,.04],.002,.0002)
    d=np.diff(a)
    assert d.max()<=.002
    assert np.max(np.maximum(d[1:]/d[:-1],d[:-1]/d[1:]))<=1.4+1e-12
    assert .02 in a


def test_microscopic_edge_rejected():
    with pytest.raises(ValueError,match='MICROSCOPIC'):graded_axis(.06,[70e-6],.002,.0002)


def test_CFL_cost():
    assert cfl_dt([1e-3]*3)==pytest.approx(1e-3/C0/np.sqrt(3))
    c=cost_estimate([10,20,30],[1e-3]*3,{'window':1e-9})
    assert c['cells']==6000 and c['estimated_memory_bytes']==6000*256
    assert c['steps']['window']>500


def test_multiband():
    assert multiband_decision(20e-9,1e-6)['decision']=='MULTIBAND_FULLWAVE_RUNS_REQUIRED'


def test_G3_contract_readability():
    s=json.loads((ROOT/'thermal/g3_port/g3_port_summary.json').read_text())
    assert s['contract']['reference_plane_id']==REFERENCE_PLANE
    assert s['method_status']=='PHYSICS_DERIVED_PORT_MODEL_VALIDATED'


def test_runtime_mesh_reporting():
    axes=[np.linspace(0,.02,21)]*3
    r=mesh_report(axes,2e9)
    assert r['cells']==8000
    assert r['estimated_CFL_dt_s']==pytest.approx(1e-3/C0/np.sqrt(3))


@pytest.mark.parametrize('axes',[[[0,0,1]]*3,[[0,np.nan,1]]*3,[[0,1e-6,.02]]*3])
def test_runtime_bad_mesh_rejected(axes):
    with pytest.raises(ValueError):mesh_report(axes,2e9)
