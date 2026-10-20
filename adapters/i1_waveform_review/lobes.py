"""Hysteretic excursion lobes distinguish cycles from within-lobe noise extrema."""
import json
import time
import numpy as np
from .ingest import OUT,write,table


def excursions(t,centered,threshold):
    ids=np.flatnonzero(np.abs(centered)>=threshold)
    if len(ids)==0: return []
    signs=np.sign(centered[ids])
    cuts=np.r_[0,np.flatnonzero(np.diff(signs)!=0)+1,len(ids)]
    result=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        chosen=ids[a:b]; k=chosen[np.argmax(np.abs(centered[chosen]))]
        result.append(dict(time_s=float(t[k]),departure_exported=float(centered[k]),
            sign=int(np.sign(centered[k])),first_above_threshold_s=float(t[chosen[0]]),
            last_above_threshold_s=float(t[chosen[-1]])))
    return result


def main():
    started=time.monotonic()
    meta=json.loads((OUT/'ingested_metadata.json').read_text())
    choices=json.loads((OUT/'window_selection_and_noise.json').read_text())
    rows=[]; extrema=[]; reductions=[]
    for name in ['signal.csv','current.csv']:
        with np.load(OUT/(name+'.burst.npz')) as z: t,y=z['t'],z['y']; onset=float(z['nominal_onset_s'])
        v=y-meta[name]['baseline_median_exported']; sigma=choices[name]['conservative_sigma']
        for a,b in [(0,100),(100,500),(0,500)]:
            mask=(t>=onset+a*1e-9)&(t<onset+b*1e-9)
            for factor in [4,6,8]:
                lobes=excursions(t[mask],v[mask],factor*sigma)
                row=dict(file=name,window_ns=f'{a}_{b}',noise_multiplier=factor,sign_reversals=max(0,len(lobes)-1),
                    cycle_definition='consecutive same-polarity dominant excursion maxima with intervening opposite excursion; no frequency filter')
                for sign,label in [(1,'positive'),(-1,'negative')]:
                    selected=[x for x in lobes if x['sign']==sign]
                    intervals=np.diff([x['time_s'] for x in selected])
                    median=float(np.median(intervals)) if len(intervals) else None
                    iqr=float(np.quantile(intervals,.75)-np.quantile(intervals,.25)) if len(intervals) else None
                    row[label+'_complete_candidate_cycles']=len(intervals)
                    row[label+'_median_period_ns']=median*1e9 if median else None
                    row[label+'_period_IQR_ns']=iqr*1e9 if iqr is not None else None
                    stable=len(intervals)>=3 and iqr/median<=.2
                    row[label+'_single_mode_status']='REGULAR_LOCAL_EXCURSIONS' if stable else 'MULTIMODE_OR_INSUFFICIENT_CYCLES'
                    decrement=None
                    if stable:
                        amplitude=np.abs([x['departure_exported'] for x in selected])
                        ratios=np.log(amplitude[:-1]/amplitude[1:])
                        if np.count_nonzero(ratios>0)>=.8*len(ratios): decrement=float(np.median(ratios))
                    row[label+'_log_decrement']=decrement
                rows.append(row)
                extrema.extend(dict(file=name,window_ns=f'{a}_{b}',noise_multiplier=factor,**x) for x in lobes)
        early=np.abs(v[(t>=onset)&(t<onset+100e-9)])
        late=np.abs(v[(t>=onset+100e-9)&(t<onset+500e-9)])
        reductions.append(dict(file=name,early_max_exported=float(early.max()),late_max_exported=float(late.max()),
            late_over_early_peak=float(late.max()/early.max()) if early.max()>0 else None,
            early_RMS_exported=float(np.sqrt(np.mean(early**2))),late_RMS_exported=float(np.sqrt(np.mean(late**2))),
            status='WINDOW_DESCRIPTIVE_NOT_A_FITTED_DECAY_CONSTANT'))
    table('excursion_cycle_sensitivity.csv',rows)
    table('excursion_extrema.csv',extrema)
    write('early_late_amplitude_comparison.json',reductions)
    write('excursion_analysis_cost.json',dict(wall_s=time.monotonic()-started))
    print(json.dumps([r for r in rows if r['noise_multiplier']==6],indent=2))

if __name__=='__main__': main()
