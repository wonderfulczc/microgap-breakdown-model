from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .lifecycle import current_moment,current_moment_derivative,peak_time
from .fourier import derivative_spectrum_hz
from .esd import esd_per_hz

def save(fig,path):
    path=Path(path); fig.tight_layout(); fig.savefig(path.with_suffix('.png'),dpi=220); fig.savefig(path.with_suffix('.pdf')); plt.close(fig)
def make_all(cfg,out):
    out=Path(out); cases={k:{'alpha':float(cfg[k]['alpha']),'beta':float(cfg[k]['beta']),'I0':float(cfg[k]['I_CM0']),'T0':float(cfg[k]['T0'])} for k in ('primary','figure3_caption_variant','slow_decay','later_transition')}; p=cases['primary']
    t=np.linspace(-2e-9,70e-9,3000); tn=t*1e9
    fig,ax=plt.subplots(figsize=(7,4.4)); ax.plot(tn,current_moment(t,**p),label='primary'); ax.axvline(p['T0']*1e9,ls='--',label='T0'); ax.axvline(peak_time(**p)*1e9,ls=':',label='actual t_peak'); ax.set(xlabel='Time (ns)',ylabel='I_CM (A·m)'); ax.legend(); ax.grid(alpha=.25); save(fig,out/'figure3a_lifecycle_primary')
    fig,ax=plt.subplots(figsize=(7,4.4)); ax.plot(tn,current_moment_derivative(t,**p)); ax.axhline(0,color='k',lw=.7); ax.set(xlabel='Time (ns)',ylabel='dI_CM/dt (A·m/s)'); ax.grid(alpha=.25); save(fig,out/'figure3b_derivative_primary')
    fig,ax=plt.subplots(figsize=(7,4.4));
    for k in ('primary','figure3_caption_variant'): ax.plot(tn,current_moment(t,**cases[k]),label=f"I0 = {cases[k]['I0']:.2f} A·m")
    ax.set(xlabel='Time (ns)',ylabel='I_CM (A·m)'); ax.legend(); ax.grid(alpha=.25); save(fig,out/'figure3_I0_conflict_comparison')
    f=np.logspace(6,11,1800); colors=cfg.get('plot_curve_identity_colors',{})
    fig,ax=plt.subplots(figsize=(7,4.8));
    for k,label in [('primary','primary (blue-equivalent)'),('slow_decay','slow decay (red-equivalent)'),('later_transition','later transition (yellow-equivalent)')]:
        y=derivative_spectrum_hz(f,**cases[k]); mask=y>=1e-18; ax.loglog(f[mask]/1e9,y[mask],label=label,color=colors.get(k))
    ax.set_ylim(1e-18,2); ax.set(xlabel='Frequency (GHz)',ylabel='ω |Ĩ_CM| (A·m)',title='Partial reconstruction of Shi 2019 Figure 4a:\nanalytical lifecycle-model curves only.'); ax.legend(); ax.grid(which='both',alpha=.22); save(fig,out/'figure4a_analytic_model_curves')
    f_esd=np.logspace(6,np.log10(3e10),1200); fig,ax=plt.subplots(figsize=(7,4.8));
    for k in ('primary','slow_decay','later_transition'): ax.loglog(f_esd,np.maximum(esd_per_hz(f_esd,**cases[k]),1e-60),label=k,color=colors.get(k))
    ax.set_ylim(1e-60,None)
    ax.set(xlabel='Frequency (Hz)',ylabel='One-sided ESD (J/Hz)'); ax.legend(); ax.grid(which='both',alpha=.22); save(fig,out/'esd_analytic_cases')




