"""Preserve three scope CSVs and expose units/quality before waveform interpretation."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation/stage_i/i1_waveform_review_20261016'
FILES=('signal.csv','difference-high-voltage.csv','current.csv')


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def write(name,value):
    with (OUT/name).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False)


def table(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    with (OUT/name).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT).decode()


def header(path):
    with path.open('rb') as f:
        prefix=f.read(4096)
    encoding='utf-8-sig' if prefix.startswith(b'\xef\xbb\xbf') else 'ascii'
    prefix.decode(encoding,errors='strict')
    metadata={}
    raw=[]
    with path.open(encoding=encoding,newline='') as f:
        for number,row in enumerate(csv.reader(f,delimiter=','),1):
            raw.append(row)
            if len(row)>=2 and row[0].strip().lower()=='second' and row[1].strip().lower()=='value':
                return dict(encoding=encoding,delimiter=',',header_lines=number,
                            raw_header_rows=raw,scope_metadata=metadata,
                            newline='CRLF' if b'\r\n' in prefix else 'LF')
            if row:
                metadata[row[0].strip()]=','.join(row[1:]).strip().strip(',')
    raise ValueError('DATA_HEADER_NOT_IDENTIFIED')


def prepare():
    OUT.mkdir(exist_ok=False,parents=True)
    prior=ROOT/'validation/i1_part_vi/spark_formation_20261016/protection_before.json'
    paths={ROOT/p for p in json.loads(prior.read_text())['files']}
    for base in ('validation/i1_part_vi/spark_formation_20261016','adapters/i1_spark_formation'):
        paths.update(p for p in (ROOT/base).rglob('*') if p.is_file() and '__pycache__' not in str(p))
    paths.update(ROOT/p for p in FILES)
    paths.update(ROOT/p for p in ('validation/stage_i/stage_i_measurement_contract.json',
          'validation/stage_i/stage_i_processing_policy.json','validation/stage_i/stage_i_scope_contract.json',
          'docs/stage_h5_signal_evidence_chain.md','docs/i1_part_vi_port_circuit_energy_contract.md'))
    paths.update(ROOT/p for p in git('diff','--name-only').splitlines())
    write('protection_before.json',dict(branch=git('branch','--show-current').strip(),HEAD=git('rev-parse','HEAD').strip(),
          files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths) if p.is_file()}))
    for name,value in [('generic_dev_before.diff',git('diff','--binary')),('git_status_before.txt',git('status','--short','--branch'))]:
        with (OUT/name).open('x') as f:
            f.write(value)
    write('analysis_plan.json',dict(files=FILES,raw_bytes_read_only=True,time_limit_s=600,
          amplitudes='instrument-export declared units, independent probe/calibration unknown',
          baseline='first1% up to100000 uniformly sampled finite observations, median/MAD',
          primary_window='strongest absolute departure from baseline: -2us to +10us, clipped to record',
          anchor_audit='also record largest raw increment and trigger-nearest time; no cross-record peak alignment',
          thresholds_sigma=[4,6,8],nominal_sigma=6,
          peak_rule='signed peaks height and prominence >= threshold; no expected-frequency distance',
          spectrum='event window; demean documented; periodic Hann and rectangular comparison; no padding',
          synchronization='compare original trigger timestamps/serial/time grids; no arbitrary shift',
          resistance='no pointwise V/I, no Toepler/RW identification without branch/calibration/same-event evidence',
          forbidden_calls=['PDE','G1','FULL','PETSc','COMSOL','fullwave','RLC evolution','formation replay']))


def run():
    started=time.monotonic()
    records={}
    for name in FILES:
        if time.monotonic()-started>550:
            raise TimeoutError('OFFLINE_BUDGET')
        path=ROOT/name
        info=header(path)
        meta=info['scope_metadata']
        if meta.get('Horizontal Units')!='s':
            raise ValueError('TIME_UNIT_NOT_CONFIRMED')
        values=np.loadtxt(path,delimiter=',',skiprows=info['header_lines'],encoding=info['encoding'],ndmin=2)
        if values.shape[1]!=2:
            raise ValueError('UNEXPECTED_DATA_COLUMNS')
        t,y=values[:,0],values[:,1]
        finite=np.isfinite(t)&np.isfinite(y)
        n=len(t)
        dt=np.diff(t)
        declared=float(meta['Sample Interval'])
        baseline_sample=y[:min(n//100,100000)]
        baseline_sample=baseline_sample[np.isfinite(baseline_sample)]
        baseline=float(np.median(baseline_sample))
        noise=float(1.482602218505602*np.median(np.abs(baseline_sample-baseline)))
        k=int(np.nanargmax(np.abs(y-baseline)))
        dy=np.diff(y)
        jump=int(np.nanargmax(np.abs(dy)))+1
        dt_median=float(np.median(dt))
        start=max(0,int(np.searchsorted(t,t[k]-2e-6)))
        end=min(n,int(np.searchsorted(t,t[k]+10e-6,side='right')))
        with (OUT/(name+'.window.npz')).open('xb') as f:
            np.savez(f,t=t[start:end],y=y[start:end],original_start_index=start,peak_original_index=k)
        # Preserve actual min/max sample times in each block, not an averaging trace.
        block=max(1,n//5000)
        m=(n//block)*block
        chunks=y[:m].reshape(-1,block)
        base=np.arange(len(chunks))*block
        ids=np.unique(np.r_[base+np.nanargmin(chunks,axis=1),base+np.nanargmax(chunks,axis=1),n-1])
        with (OUT/(name+'.overview.npz')).open('xb') as f:
            np.savez(f,t=t[ids],y=y[ids],original_indices=ids)
        quant_sample=np.unique(y[::max(1,n//200000)])
        qdiff=np.diff(quant_sample)
        positive=qdiff[qdiff>0]
        qstep=float(np.median(positive)) if len(positive) else None
        ymin,ymax=float(np.nanmin(y)),float(np.nanmax(y))
        at_min,at_max=y==ymin,y==ymax
        def run_length(mask):
            changes=np.diff(np.r_[False,mask,False].astype(np.int8))
            a,b=np.flatnonzero(changes==1),np.flatnonzero(changes==-1)
            return int(np.max(b-a)) if len(a) else 0
        trigger_index=int(meta['Trigger index'].split(':')[-1])-1
        info.update(path=name,sha256=sha(path),bytes=path.stat().st_size,samples=n,
             finite_pair_count=int(np.count_nonzero(finite)),nonfinite_pair_count=int(np.count_nonzero(~finite)),
             declared_record_length=int(meta['Record Length'].split(':')[-1]),
             t_start_s=float(t[0]),t_end_s=float(t[-1]),sample_interval_median_s=dt_median,
             sample_interval_min_s=float(np.min(dt)),sample_interval_max_s=float(np.max(dt)),
             max_interval_deviation_from_metadata_s=float(np.max(np.abs(dt-declared))),
             duplicate_times=int(np.count_nonzero(dt==0)),decreasing_times=int(np.count_nonzero(dt<0)),
             interval_gap_candidates=int(np.count_nonzero(dt>1.5*declared)),
             sample_rate_Hz=1/declared,nyquist_Hz=.5/declared,
             exported_amplitude_unit=meta.get('Vertical Units'),unit_status='EXPORT_DECLARED_NOT_INDEPENDENTLY_CALIBRATED',
             calibrated_engineering_amplitude=False,probe_factor=None,probe_model=None,probe_bandwidth_Hz=None,
             input_loading_ohm=None,polarity_mapping=None,
             baseline_median_exported=baseline,baseline_noise_MAD_sigma_exported=noise,
             baseline_samples=len(baseline_sample),baseline_window_end_s=float(t[len(baseline_sample)-1]),
             y_min_exported=ymin,y_max_exported=ymax,global_departure_peak_index=k,
             global_departure_peak_time_s=float(t[k]),largest_increment_time_s=float(t[jump]),
             largest_increment_exported=float(dy[jump-1]),
             trigger_sample_original_index=trigger_index,
             trigger_index_sample_time_s=float(t[trigger_index]) if 0<=trigger_index<n else None,
             median_positive_unique_step_sampled=qstep,
             minimum_value_count=int(np.count_nonzero(at_min)),maximum_value_count=int(np.count_nonzero(at_max)),
             longest_minimum_plateau_samples=run_length(at_min),longest_maximum_plateau_samples=run_length(at_max),
             saturation_status='RAILS_AND_PROBE_RANGE_UNKNOWN; extrema repetitions are flags, not saturation proof',
             truncation_status='REQUIRES_SELECTED_EVENT_BOUNDARY_AUDIT',
             selected_window=dict(original_start_index=start,original_stop_index_exclusive=end,
                 start_s=float(t[start]),end_s=float(t[end-1]),anchor='strongest raw absolute departure',
                 clipped_left=start==0,clipped_right=end==n),
             measurement_role='TX_VOLTAGE_USER_REPORTED' if name.startswith('difference') else
                              'TX_CURRENT_USER_REPORTED' if name=='current.csv' else 'OSCILLATORY_SIGNAL_MEASUREMENT_LOCATION_UNKNOWN')
        records[name]=info
        print(json.dumps(dict(file=name,samples=n,dt_s=dt_median,range_s=[float(t[0]),float(t[-1])],
              baseline=baseline,noise=noise,minimum=ymin,maximum=ymax,peak_time_s=float(t[k]),jump_time_s=float(t[jump]),
              window_samples=end-start,elapsed_s=time.monotonic()-started)),flush=True)
        del values,t,y,dt,dy,chunks,at_min,at_max,finite
    write('ingested_metadata.json',records)
    write('ingest_cost.json',dict(wall_s=time.monotonic()-started,physical_solver_calls=0))


if __name__=='__main__':
    import sys
    if sys.argv[1:]==['prepare']:
        prepare()
    elif sys.argv[1:]==['run']:
        run()
    else:
        raise SystemExit('prepare then run once')
