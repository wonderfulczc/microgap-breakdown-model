"""Noise-qualified burst catalog; preserves same-frame coordinates without shifts."""
import hashlib
import json
import time

import numpy as np

from .ingest import ROOT,OUT,FILES,write,table


def burst_groups(indices,maximum_gap):
    if not len(indices):
        return []
    cuts=np.flatnonzero(np.diff(indices)>maximum_gap)+1
    return [(int(a[0]),int(a[-1])) for a in np.split(indices,cuts)]


def main():
    started=time.monotonic()
    metadata=json.loads((OUT/'ingested_metadata.json').read_text())
    catalog=[]
    details={}
    signal=metadata['signal.csv']
    # Registered independent raw time range, not a peak-alignment operation.
    common_center=signal['global_departure_peak_time_s']
    for name in FILES:
        m=metadata[name]
        data=np.loadtxt(ROOT/name,delimiter=',',skiprows=m['header_lines'],encoding=m['encoding'])
        t,y=data[:,0],data[:,1]
        baseline=m['baseline_median_exported']
        sigma=m['baseline_noise_MAD_sigma_exported']
        if sigma<=0:
            details[name]=dict(status='NO_ROBUST_NOISE_SCALE',reason='No arbitrary epsilon or default threshold')
            continue
        maximum_gap=round(200e-9/float(m['scope_metadata']['Sample Interval']))
        groups=burst_groups(np.flatnonzero(np.abs(y-baseline)>=6*sigma),maximum_gap)
        dominant=None
        peak=m['global_departure_peak_index']
        for number,(left,right) in enumerate(groups):
            strongest=int(left+np.argmax(np.abs(y[left:right+1]-baseline)))
            is_dominant=left<=peak<=right
            row=dict(file=name,burst_id=number,onset_index=left,last_above_threshold_index=right,
                     onset_time_s=float(t[left]),last_above_threshold_time_s=float(t[right]),
                     strongest_index=strongest,strongest_time_s=float(t[strongest]),
                     strongest_departure_exported=float(y[strongest]-baseline),
                     amplitude_threshold_exported=6*sigma,merge_gap_s=200e-9,
                     dominant=is_dominant,clipped_left=left==0,clipped_right=right==len(y)-1)
            catalog.append(row)
            if is_dominant:
                dominant=row
        details[name]=dict(noise_scale=sigma,burst_count_nominal=len(groups),dominant=dominant,
                          time_axis_float64_sha256=hashlib.sha256(np.ascontiguousarray(t).view(np.uint8)).hexdigest(),
                          algorithm='abs(raw-baseline)>=6*MAD_sigma; merge gaps<=200ns; no frequency targeting')
        if name!='difference-high-voltage.csv' and dominant:
            # One burst only; include baseline and a tail, never absorb a later burst.
            l=max(0,dominant['onset_index']-round(100e-9/float(m['scope_metadata']['Sample Interval'])))
            r=min(len(y),dominant['last_above_threshold_index']+round(100e-9/float(m['scope_metadata']['Sample Interval']))+1)
            subsequent=[g[0] for g in groups if g[0]>dominant['last_above_threshold_index']]
            if subsequent:
                r=min(r,subsequent[0])
            with (OUT/(name+'.burst.npz')).open('xb') as f:
                np.savez(f,t=t[l:r],y=y[l:r],original_start_index=l,
                         nominal_onset_s=dominant['onset_time_s'],nominal_last_above_s=dominant['last_above_threshold_time_s'])
        if name in ('signal.csv','difference-high-voltage.csv'):
            l=int(np.searchsorted(t,common_center-300e-9))
            r=int(np.searchsorted(t,common_center+1200e-9,side='right'))
            with (OUT/(name+'.common_frame.npz')).open('xb') as f:
                np.savez(f,t=t[l:r],y=y[l:r],original_start_index=l,
                         identity='same metadata capture; original time, no deskew or peak alignment')
        print(json.dumps(dict(file=name,groups=len(groups),dominant=dominant)),flush=True)
        del data,t,y
    table('raw_burst_catalog.csv',catalog)
    write('event_catalog.json',details)
    write('event_extraction_cost.json',dict(wall_s=time.monotonic()-started,new_physical_solver_calls=0))


if __name__=='__main__':
    main()
