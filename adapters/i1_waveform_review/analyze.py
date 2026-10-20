"""Offline, raw-scale waveform diagnostics; no resistance or circuit fitting."""
import json
import time
import numpy as np
from scipy.signal import find_peaks, windows
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ingest import ROOT,OUT,FILES,write,table


def spectrum(t,y,kind='hann'):
    n=len(y)
    dt=float(np.median(np.diff(t)))
    w=windows.hann(n,sym=False) if kind=='hann' else np.ones(n)
    a=2*np.abs(np.fft.rfft((y-np.mean(y))*w))/np.sum(w)
    a[0]/=2
    if n%2==0: a[-1]/=2
    f=np.fft.rfftfreq(n,dt)
    k=1+int(np.argmax(a[1:]))
    return f,a,dict(samples=n,first_sample_s=float(t[0]),last_sample_s=float(t[-1]),
        observation_Ndt_s=n*dt,bin_spacing_Hz=1/(n*dt),nyquist_Hz=.5/dt,
        window=kind,demean=True,zero_padding=False,
        ENBW_Hz=float(np.sum(w*w)/np.sum(w)**2/dt),
        dominant_nonDC_bin_Hz=float(f[k]),dominant_bin_amplitude_exported=float(a[k]),
        amplitude_normalization='single-sided 2*abs(rfft)/sum(window); DC/Nyquist not doubled',
        frequency_status='BIN_DIAGNOSTIC_NOT_SINGLE_MODE_IDENTIFICATION')


def peaks(t,y,baseline,sigma,multiplier):
    if sigma<=0: return dict(status='NO_ROBUST_NOISE_SCALE'),[]
    out=[]
    values={}
    periods=[]
    for sign,label in [(1,'positive'),(-1,'negative')]:
        ids,prop=find_peaks(sign*(y-baseline),height=multiplier*sigma,prominence=multiplier*sigma)
        gaps=np.diff(t[ids])
        median=float(np.median(gaps)) if len(gaps) else None
        iqr=float(np.quantile(gaps,.75)-np.quantile(gaps,.25)) if len(gaps) else None
        stable=len(gaps)>=3 and median>0 and iqr/median<=.2
        values.update({label+'_peaks':len(ids),label+'_adjacent_same_sign_intervals':max(0,len(ids)-1),
            label+'_median_interval_s':median,label+'_period_IQR_s':iqr,
            label+'_median_interval_inverse_Hz':1/median if median else None,
            label+'_period_status':'APPROXIMATELY_REGULAR_IN_SELECTED_WINDOW' if stable else 'MULTICOMPONENT_OR_INSUFFICIENT_CYCLES'})
        decrement=None
        if stable:
            amp=prop['peak_heights']
            logs=np.log(amp[:-1]/amp[1:])
            # Do not infer one decaying mode from alternating growing/shrinking peaks.
            if np.count_nonzero(logs>0)>=.8*len(logs):
                decrement=float(np.median(logs))
        values[label+'_median_log_decrement']=decrement
        values[label+'_decay_status']='FINITE_WINDOW_DESCRIPTIVE_ONLY' if decrement is not None else 'NO_RELIABLE_SINGLE_EXPONENTIAL_DECAY'
        for idx,h,p in zip(ids,prop['peak_heights'],prop['prominences']):
            out.append(dict(time_s=float(t[idx]),polarity=label,departure_exported=float(sign*h),prominence_exported=float(p)))
        periods.extend(gaps.tolist())
    labels=np.zeros(len(y),dtype=int)
    labels[y-baseline>=multiplier*sigma]=1
    labels[y-baseline<=-multiplier*sigma]=-1
    qualified=labels[labels!=0]
    crossings=int(np.count_nonzero(np.diff(qualified)!=0)) if len(qualified)>1 else 0
    values.update(noise_multiplier=multiplier,threshold_exported=multiplier*sigma,
                  hysteretic_sign_reversals=crossings,
                  count_warning='Intervals are candidate cycles, not guaranteed complete single-mode cycles')
    return values,out


def main():
    started=time.monotonic()
    meta=json.loads((OUT/'ingested_metadata.json').read_text())
    events=json.loads((OUT/'event_catalog.json').read_text())
    write('diagnostic_window_plan.json',dict(
        reason='Raw inspection shows multiscale bursts and a slow CH3 voltage background',
        fast_windows_ns=[[0,100],[100,500],[0,500]],
        anchor='first nominal 6sigma crossing of dominant burst, each original record separately',
        voltage_windows='Original CH1/CH3 common frame; CH3 original 0..60us slow response',
        voltage_AC='Optional least-squares affine trend removal only; not bandpass or physical branch separation',
        noise='max(first1% MAD_sigma, raw 100ns pre-onset MAD_sigma), multiplier4/6/8',
        regularity='>=3 same-polarity intervals and IQR/median<=0.2 diagnostic, not a physics gate',
        decay='only if regular and >=80% of consecutive amplitude ratios decay; otherwise unavailable',
        new_physical_model_calls=0))
    rows=[]; peakrows=[]; specrows=[]; specdetails=[]; selections={}
    # Read only one existing scalar file for a missing slow-time window.
    name='difference-high-voltage.csv'; m=meta[name]
    raw=np.loadtxt(ROOT/name,delimiter=',',skiprows=m['header_lines'],encoding=m['encoding'])
    mask=(raw[:,0]>=0)&(raw[:,0]<=60e-6)
    slow=raw[mask].copy(); del raw
    with (OUT/'difference-high-voltage.csv.slow.npz').open('xb') as f:
        np.savez(f,t=slow[:,0],y=slow[:,1])
    for name in FILES:
        m=meta[name]; baseline=m['baseline_median_exported']; global_sigma=m['baseline_noise_MAD_sigma_exported']
        if name!='difference-high-voltage.csv':
            with np.load(OUT/(name+'.burst.npz')) as z: t,y=z['t'],z['y']; onset=float(z['nominal_onset_s'])
            pre=y[t<onset]
            local=float(1.482602218505602*np.median(np.abs(pre-np.median(pre))))
            sigma=max(global_sigma,local)
            windows_to_use=[]
            for a,b in [(0,100),(100,500),(0,500)]:
                mask=(t>=onset+a*1e-9)&(t<onset+b*1e-9)
                windows_to_use.append((f'{a}_{b}ns',t[mask],y[mask],baseline,sigma))
            selections[name]=dict(original_onset_s=onset,global_noise_sigma=global_sigma,local_pre_noise_sigma=local,
                conservative_sigma=sigma,burst_count_merge200ns=events[name]['burst_count_nominal'],
                warning='200ns grouping can merge renewed excitations; fixed first500ns used, not all clusters as one decay')
        else:
            with np.load(OUT/(name+'.common_frame.npz')) as z: t,y=z['t'],z['y']
            # Affine subtraction is disclosed; large nonlinear slow background remains.
            coord=(t-t.mean())/(t[-1]-t[0])
            trend=np.polyval(np.polyfit(coord,y,1),coord)
            windows_to_use=[('common_frame_raw',t,y,baseline,global_sigma),
                            ('common_frame_affine_residual',t,y-trend,0.,global_sigma),
                            ('slow_0_60us_raw',slow[:,0],slow[:,1],baseline,global_sigma)]
            selections[name]=dict(warning='Nonzero slow background; noise extrema are not complete RF cycles; affine residual is not gas voltage',
                response_right_censored=events[name]['dominant']['clipped_right'])
        for label,t,y,base,sigma in windows_to_use:
            for k in [4,6,8]:
                result,p=peaks(t,y,base,sigma,k)
                rows.append(dict(file=name,window_id=label,unit_exported=m['exported_amplitude_unit'],
                    first_sample_s=float(t[0]),last_sample_s=float(t[-1]),
                    baseline_exported=base,y_min_exported=float(y.min()),y_max_exported=float(y.max()),**result))
                peakrows.extend(dict(file=name,window_id=label,noise_multiplier=k,**x) for x in p)
            for win in ['hann','rectangular']:
                f,a,s=spectrum(t,y,win)
                specdetails.append(dict(file=name,window_id=label,**s))
                use=(f>0)&(f<=200e6)
                specrows.extend(dict(file=name,window_id=label,window=win,frequency_Hz=float(ff),amplitude_exported=float(aa)) for ff,aa in zip(f[use],a[use]))
    table('waveform_characteristics.csv',rows)
    table('peak_candidates_threshold_sensitivity.csv',peakrows)
    table('spectral_bins.csv',specrows)
    write('spectral_methods_and_results.json',specdetails)
    write('window_selection_and_noise.json',selections)
    # Raw first500ns plus actual signed peak candidates; distinct record clocks stay distinct.
    fig,axes=plt.subplots(3,2,figsize=(15,11),layout='constrained')
    for r,name in enumerate(FILES):
        if name!='difference-high-voltage.csv':
            with np.load(OUT/(name+'.burst.npz')) as z: t,y=z['t'],z['y']; onset=float(z['nominal_onset_s'])
            mask=(t>=onset-30e-9)&(t<=onset+500e-9); t,y=t[mask],y[mask]
            label='0_500ns'
        else:
            with np.load(OUT/(name+'.common_frame.npz')) as z: t,y=z['t'],z['y']
            label='common_frame_raw'
        axes[r,0].plot(t*1e6,y,lw=.7)
        candidates=[p for p in peakrows if p['file']==name and p['window_id']==label and p['noise_multiplier']==6]
        axes[r,0].scatter([p['time_s']*1e6 for p in candidates],
            [p['departure_exported']+meta[name]['baseline_median_exported'] for p in candidates],s=9,color='red')
        axes[r,0].set_title(name+' / original record time, no deskew')
        axes[r,0].set_xlabel('us'); axes[r,0].set_ylabel('exported '+meta[name]['exported_amplitude_unit'])
        for win in ['hann','rectangular']:
            s=[p for p in specrows if p['file']==name and p['window_id']==label and p['window']==win]
            axes[r,1].plot([p['frequency_Hz']/1e6 for p in s],[p['amplitude_exported'] for p in s],label=win,lw=.8)
        axes[r,1].set_title('Demeaned finite-window amplitude, no zero padding')
        axes[r,1].set_xlabel('MHz');axes[r,1].set_ylabel('exported unit amplitude');axes[r,1].legend()
        for ax in axes[r]: ax.grid(alpha=.2)
    fig.savefig(OUT/'oscillation_windows_peaks_and_spectra.png',dpi=160);plt.close(fig)
    write('analysis_cost.json',dict(wall_s=time.monotonic()-started,physical_solver_calls=0,parameter_identification_runs=0))
    print(json.dumps(dict(rows=len(rows),spectra=len(specdetails),wall_s=time.monotonic()-started)))

if __name__=='__main__': main()
