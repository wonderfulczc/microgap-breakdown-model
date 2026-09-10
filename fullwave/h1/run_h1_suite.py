"""Run or reuse exactly the four authorized H1 numerical configurations."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

CASES={'baseline':(.5,8),'coarse':(.75,8),'fine':(.25,8),'pml10':(.5,10)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--raw-root',type=Path,required=True)
    p.add_argument('--case',choices=list(CASES)+['all'],default='all')
    a=p.parse_args();a.raw_root.mkdir(parents=True,exist_ok=True)
    script=Path(__file__).with_name('run_h1_reference.py')
    for case in CASES if a.case=='all' else [a.case]:
        result=a.raw_root/case/'result.json'
        if result.exists():
            stored=json.loads(result.read_text())
            if stored['script_sha256']!=hashlib.sha256(script.read_bytes()).hexdigest():
                raise ValueError('EXISTING_RESULT_CODE_HASH_DIFFERS')
            print('Reusing',case,flush=True)
            continue
        step,pml=CASES[case]
        with (a.raw_root/f'{case}.log').open('x') as log:
            subprocess.run([sys.executable,str(script),'--case',case,'--output',str(a.raw_root/case),
                '--spacing-mm',str(step),'--pml',str(pml)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
        print(case,result,flush=True)


if __name__=='__main__':main()
