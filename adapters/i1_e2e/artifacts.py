import hashlib
import json
import csv
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation/i1_e2e/v1_20261020'
P1=ROOT/'validation/i1_p1/implementation_20261020'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def dump(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def table(p,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('x',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def protect():
    OUT.mkdir(parents=True,exist_ok=False)
    old=json.loads((P1/'protection_after.json').read_text())['after_sha256']
    paths=set(old)
    for folder in ('adapters/i1_p1','validation/i1_p1'):
        paths.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).rglob('*') if p.is_file())
    dump(OUT/'protection_before.json',dict(sha256={p:sha(ROOT/p) for p in sorted(paths)},
        head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
        status=subprocess.check_output(['git','status','--porcelain'],text=True)))
    (OUT/'GENERIC_DEV_before.diff').write_bytes(subprocess.check_output(['git','diff','--binary']))
    for p in ('docs/i1_research_status_and_roadmap.md','docs/i1_open_issues_and_gates.md'):
        target=OUT/'baseline'/p;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/p).read_bytes())
    dump(OUT/'execution_plan.json',dict(wall_limit_s=2400,max_scenarios=4,
        order=['synthetic_e2e','rom_R3','rom_R300','rom_gap100um'],
        source='P1 synthetic x=0.1mm, sigma=1e-6, epsPET=3 epsFEP=2; unchanged in all cases',
        synthetic=dict(end_s=100e-9,dt_s=20e-12,R_ohm=5.,Cp_F=.2e-12),
        rom=dict(end_s=200e-12,dt_s=2e-12,cells=12,ne_np0_m3=1e16,nn0_m3=0.,
            radius_m=50e-6,gap_m=70e-6,mu_p=2e-4,D_p=5e-6,gamma_i=0.,
            ion_parameter_class='C_SYNTHETIC_NOT_CALIBRATED',electron_coefficients='EXISTING_NATIVE_PARITY_TRANSCRIPTION',
            timestep_safety=.15,max_substeps_per_trial=20000,negative_density='REJECT_NO_FLOOR'),
        checks=dict(graph_charge_relative=1e-10,cell_charge_relative=1e-8,
            source_quadrature_relative=1e-3,receiver_algebra_relative=1e-12,
            circuit_KVL_V=1e-7,energy_relative_exploratory=1e-3),
        no_autoretry_failed_scenario=True,no_full_or_g1_steps=True))
    print('Protection established',len(paths))

if __name__=='__main__':protect()
