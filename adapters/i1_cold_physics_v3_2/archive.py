from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation/i1_cold_physics_v3_2/implementation_20261019'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def protect():
    OUT.mkdir(parents=True,exist_ok=False)
    prior=json.loads((ROOT/'validation/i1_cold_physics_v3_1/implementation_20261019/protection_after.json').read_text())
    paths=set(prior['after_sha256'])
    for folder in ('adapters/i1_cold_physics_v3_1','validation/i1_cold_physics_v3_1','build_i1_cold_physics_v3_1','python/streamer_rf/rf/source'):
        paths.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p))
    for name in paths:
        if name.startswith(('cpp/','python/','adapters/')) or name=='CMakeLists.txt' or name.startswith('docs/'):
            p=OUT/'baseline'/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,p)
    dump('protection_before.json',dict(sha256={p:sha(ROOT/p) for p in sorted(paths)},
      head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
      status=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)))
    (OUT/'GENERIC_DEV_before.diff').write_bytes(subprocess.check_output(['git','diff','--binary'],cwd=ROOT))
    print(len(paths))
if __name__=='__main__':protect()
