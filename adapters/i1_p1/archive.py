"""Exclusive P1 artifacts and preservation of the uncommitted worktree."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'validation/i1_p1/implementation_20261020'
OLD = ROOT/'validation/i1_cold_physics_v3_5/implementation_20261020'
AUTHORIZED = ('adapters/i1_cold_physics_v3_5/scenario.py',
              'docs/i1_research_status_and_roadmap.md','docs/i1_open_issues_and_gates.md')


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def dump(p, data):
    with Path(p).open('x') as f:
        json.dump(data,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')


def ref(p):
    p=Path(p).resolve()
    return dict(file=str(p.relative_to(ROOT)),sha256=sha(p))


def load(r):
    p=(ROOT/r['file']).resolve()
    if not p.is_relative_to(ROOT) or not p.is_file() or sha(p)!=r['sha256']:
        raise ValueError('SOURCE_IDENTITY_MISMATCH')
    return json.loads(p.read_text())


def table(p, rows):
    with Path(p).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def protect():
    OUT.mkdir(parents=True,exist_ok=False)
    paths=set(json.loads((OLD/'protection_after.json').read_text())['after_sha256'])
    for folder in ('adapters/i1_cold_physics_v3_5','validation/i1_cold_physics_v3_5','docs'):
        paths.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).rglob('*')
                     if p.is_file() and '__pycache__' not in str(p))
    dump(OUT/'protection_before.json',dict(sha256={p:sha(ROOT/p) for p in sorted(paths)},
        head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
        status=subprocess.check_output(['git','status','--porcelain'],text=True)))
    (OUT/'GENERIC_DEV_before.diff').write_bytes(subprocess.check_output(['git','diff','--binary']))
    for name in AUTHORIZED:
        dest=OUT/'baseline'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dest)
    dump(OUT/'execution_plan.json',dict(wall_limit_s=900,compile_separate=True,
        dynamic_cases_max=2,end_s=1e-8,dt_s=1e-12,
        dielectric_reference=dict(PET=3.,FEP=2.,adhesive=3.,identity='SYNTHETIC_ENGINEERING_NOT_LITERATURE'),
        sigma_eff_C_m2=1e-6,sigma_identity='SYNTHETIC_ENGINEERING_NOT_DEVICE',
        x_references_m=[1e-4,1e-3],dynamic_TENG_x_m=1e-4,
        mathematical_R_ohm=200.,synthetic_Cp_F=.2e-12,
        amplifier_reference_V=1.,synthetic_Rbleed_ohm=1e6,
        dynamic_amplifier_impedance='IDEAL_AMPLIFIER_LIMIT',
        state_reference_relative_target=.01,charge_KVL_energy_scaled_target=1e-8,
        ownership_single_trial_target=1e-10,
        normalization='fixed initial voltage, Cext*V, Cext*V^2 and V/sqrt(L/Cext); retain absolute errors',
        independent_reference='constant augmented ODE matrix exponential, not interface history',
        cold_steps=0,thermal_steps=0,propagation_calls=0,
        missing_assessment='I1_parameter_assessment_20261020.md NOT_FOUND_IN_WORKSPACE',
        no_retry_failed_dynamic_case=True))
    print('Protected',len(paths),'files')


if __name__=='__main__': protect()
