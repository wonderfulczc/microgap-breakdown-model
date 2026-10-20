"""Exclusive outputs and source identities; never publish into historical runs."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'validation/i1_cold_physics_v3_5/implementation_20261020'
PREVIOUS = ROOT / 'validation/i1_cold_physics_v3_4/implementation_20261019'
DOCS = ('docs/i1_research_status_and_roadmap.md', 'docs/i1_open_issues_and_gates.md')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def table(path, rows):
    with Path(path).open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def file_ref(path):
    p = Path(path).resolve()
    return {'file': str(p.relative_to(ROOT)), 'sha256': sha(p)}


def protect():
    OUT.mkdir(parents=True, exist_ok=False)
    paths = set(json.loads((PREVIOUS/'protection_after.json').read_text())['after_sha256'])
    for folder in ('adapters', 'validation/i1_cold_physics_v3_4', 'python/streamer_rf',
                   'cpp/include', 'cpp/src', 'docs'):
        paths.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).rglob('*')
                     if p.is_file() and '__pycache__' not in str(p)
                     and 'i1_cold_physics_v3_5' not in str(p))
    paths.update(DOCS)
    for folder in ROOT.glob('build_i1_cold_physics*'):
        paths.update(str(p.relative_to(ROOT)) for p in folder.glob('*') if p.is_file()
                     and (p.suffix == '.a' or p.name.startswith('i1_')))
    dump(OUT/'protection_before.json', dict(sha256={p: sha(ROOT/p) for p in sorted(paths)},
         head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
         branch=subprocess.check_output(['git', 'branch', '--show-current'], text=True).strip(),
         status=subprocess.check_output(['git', 'status', '--porcelain'], text=True)))
    (OUT/'GENERIC_DEV_before.diff').write_bytes(subprocess.check_output(['git', 'diff', '--binary']))
    for name in DOCS:
        dest = OUT/'baseline'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, dest)
    dump(OUT/'execution_plan.json', dict(wall_limit_s=900, compile_separate=True,
         cold_accepted_max=3, cold_step_calls_max=12, planned_cold_calls=3,
         thermal_accepted_max=4, thermal_trial_call_cap=4, thermal_dt_cap_s=5e-11,
         thermal_energy_exchange_relative_limit=1e-3, thermal_mass_relative_limit=1e-10,
         thermal_current_area_relative_limit=1e-12, cell_charge_relative_limit=1e-6,
         global_cell_charge_relative_limit=1e-8, time_identity='exact anchor plus local time',
         no_new_propagation=True, no_automatic_retry=True, thermal_fixture='ARCHIVED_G1_1MM',
         schema_validation='Explicit semantic checks; no claim of generic Draft meta-schema certification'))
    print('Protected files:', len(paths))


if __name__ == '__main__':
    protect()
