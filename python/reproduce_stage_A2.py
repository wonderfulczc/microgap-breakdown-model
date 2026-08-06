#!/usr/bin/env python
from pathlib import Path
import argparse,csv,json,platform,sys
import numpy as np
from scipy.integrate import quad,simpson
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve(); ROOT=HERE.parents[1]; sys.path.insert(0,str(HERE.parent))
from streamer_rf.config import load_config,params,cross_check_registry
from streamer_rf.lifecycle import *
from streamer_rf.fourier import *
from streamer_rf.esd import *
from streamer_rf.plotting import make_all

def write_csv(path,fields,rows):
 path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf-8') as h:w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
def quad_transform(w,p):
 q=p['alpha']+p['beta']; a=p['alpha']/q; k=w/q; L=-max(80,38/a); R=max(80,38/(1-a))
 if abs(k)<5:
  def base(y): return np.exp(a*y-np.logaddexp(0,y))
  re=quad(base,L,R,weight='cos',wvar=k,epsabs=2e-12,epsrel=2e-12,limit=1000)[0]
  im=-quad(base,L,R,weight='sin',wvar=k,epsabs=2e-12,epsrel=2e-12,limit=1000)[0]
 else:
  # Shift below the real-time contour without crossing the nearest pole at -i*pi.
  delta=min(0.1,2/abs(k)); cshift=np.pi-delta; phase=np.exp(-1j*a*cshift)
  def g(u):
   if u<=0: return phase*np.exp(a*u)/(1+np.exp(u-1j*cshift))
   return phase*np.exp((a-1)*u+1j*cshift)/(1+np.exp(-u+1j*cshift))
  rc=quad(lambda u:g(u).real,L,R,weight='cos',wvar=k,epsabs=1e-11,limit=1500)[0]
  rs=quad(lambda u:g(u).real,L,R,weight='sin',wvar=k,epsabs=1e-11,limit=1500)[0]
  ic=quad(lambda u:g(u).imag,L,R,weight='cos',wvar=k,epsabs=1e-11,limit=1500)[0]
  iss=quad(lambda u:g(u).imag,L,R,weight='sin',wvar=k,epsabs=1e-11,limit=1500)[0]
  re=(rc+iss)*np.exp(-abs(k)*cshift); im=(ic-rs)*np.exp(-abs(k)*cshift)
  if k<0: im=-im
 return p['I0']/q*(re+1j*im)*np.exp(-1j*w*p['T0'])
def extrema(p):
 tp=peak_time(**p); span=max(25/p['alpha'],25/p['beta']); f=lambda t: current_moment_derivative(t,**p)
 pos=minimize_scalar(lambda t:-f(t),bounds=(tp-span,tp),method='bounded',options={'xatol':1e-18}); neg=minimize_scalar(f,bounds=(tp,tp+span),method='bounded',options={'xatol':1e-18}); return f(pos.x),f(neg.x)
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); ap.add_argument('--output',required=True); args=ap.parse_args(); cfg=load_config(args.config); cross_check_registry(cfg,ROOT/'docs/parameter_registry.csv'); out=Path(args.output); data=out/'data'; val=out/'validation'; figs=out/'figures'; [x.mkdir(parents=True,exist_ok=True) for x in (data,val,figs)]
 names=('primary','figure3_caption_variant','slow_decay','later_transition'); C={n:params(cfg[n]) for n in names}; p=C['primary']
 tp=peak_time(**p); ip=peak_value(**p); pp,nn=extrema(p)
 write_csv(val/'peak_parameter_check.csv',['case','T0_s','t_peak_s','t_peak_minus_T0_s','I0_A_m','I_peak_A_m','I_peak_over_I0'],[{'case':n,'T0_s':C[n]['T0'],'t_peak_s':peak_time(**C[n]),'t_peak_minus_T0_s':peak_time(**C[n])-C[n]['T0'],'I0_A_m':C[n]['I0'],'I_peak_A_m':peak_value(**C[n]),'I_peak_over_I0':peak_value(**C[n])/C[n]['I0']} for n in names])
 t=np.linspace(-2e-9,100e-9,5101); rows=[]
 for ti,ii,di in zip(t,current_moment(t,**p),current_moment_derivative(t,**p)): rows.append({'time_s':ti,'I_CM_A_m':ii,'dI_CM_dt_A_m_per_s':di})
 write_csv(data/'lifecycle_primary.csv',rows[0].keys(),rows)
 rows=[]
 for ti in t:
  z={'time_s':ti}
  for n in names:z[n+'_I_CM_A_m']=current_moment(ti,**C[n]);z[n+'_dI_dt_A_m_per_s']=current_moment_derivative(ti,**C[n])
  rows.append(z)
 write_csv(data/'lifecycle_all_cases.csv',rows[0].keys(),rows)
 f=np.logspace(6,11,2001); sr=[]; er=[]
 for fi in f:
  s={'frequency_Hz':fi};e={'frequency_Hz':fi}
  for n in names:s[n+'_derivative_spectrum_A_m']=derivative_spectrum_hz(fi,**C[n]);e[n+'_ESD_J_per_Hz']=esd_per_hz(fi,**C[n])
  sr.append(s);er.append(e)
 write_csv(data/'spectrum_analytic_all_cases.csv',sr[0].keys(),sr);write_csv(data/'esd_analytic_all_cases.csv',er[0].keys(),er)
 bands=cfg['frequency_bands_hz']; br=[]
 dense=np.logspace(np.log10(3e7),np.log10(3e10),30001)
 for n in names:
  ee=esd_per_hz(dense,**C[n])
  for b,(lo,hi) in bands.items(): br.append({'case':n,'band':b,'f_low_Hz':lo,'f_high_Hz':hi,'energy_J':integrate_band_energy(dense,ee,float(lo),float(hi))})
 write_csv(data/'band_energy_summary.csv',br[0].keys(),br)
 vf=np.array([1e6,1e7,3e7,1e8,3e8,1e9,3e9,1e10,3e10]); qr=[]
 for fi in vf:
  w=2*np.pi*fi; az=current_moment_transform_omega(w,**p); qz=quad_transform(w,p); qr.append({'frequency_Hz':fi,'analytic_real_A_m_s':az.real,'analytic_imag_A_m_s':az.imag,'quad_real_A_m_s':qz.real,'quad_imag_A_m_s':qz.imag,'relative_error':abs(az-qz)/abs(az)})
 write_csv(val/'fourier_validation.csv',qr[0].keys(),qr)
 tg=build_fft_time_grid(**p,max_frequency_hz=30e9); fg,zd=discrete_fft_derivative(tg,**p,method='analytic'); fm=np.interp(vf,fg,np.abs(zd)); ana=derivative_spectrum_hz(vf,**p); res=fg[1]-fg[0]; ny=fg[-1]; trusted=(vf>=5*res)&(vf<=ny/10)&(ana>=np.max(ana)*1e-6); fr=[{'frequency_Hz':x,'fft_A_m':y,'analytic_A_m':z,'relative_error':abs(y-z)/z,'trusted':bool(k),'frequency_resolution_Hz':res,'nyquist_Hz':ny} for x,y,z,k in zip(vf,fm,ana,trusted)]
 write_csv(val/'fft_validation.csv',fr[0].keys(),fr)
 # FFT plot: mark untrusted points explicitly
 import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
 ff=np.logspace(np.log10(max(res,1e6)),np.log10(3e10),600); interp=np.interp(ff,fg,np.abs(zd)); fig,ax=plt.subplots(figsize=(7,4.8));ax.loglog(ff,derivative_spectrum_hz(ff,**p),label='analytic');ax.loglog(ff,interp,'--',label='FFT, analytic derivative sampled');ax.scatter(vf[~trusted],ana[~trusted],marker='x',label='validation points marked untrusted');ax.set_ylim(1e-16,2);ax.set(xlabel='Frequency (Hz)',ylabel='Derivative spectrum (A·m)');ax.legend();ax.grid(which='both',alpha=.22);fig.savefig(figs/'fft_vs_analytic_validation.png',dpi=220,bbox_inches='tight');fig.savefig(figs/'fft_vs_analytic_validation.pdf',bbox_inches='tight');plt.close(fig)
 # rad/s vs Hz energy consistency over full requested analytical range
 wf=2*np.pi*dense; Ehz=simpson(esd_per_hz(dense,**p),x=dense); Eom=2*simpson(esd_per_angular_frequency(wf,**p),x=wf); consistency=abs(Ehz-Eom)/Ehz
 p40=C['figure3_caption_variant']; pp40,nn40=extrema(p40); spec_ratio=derivative_spectrum_hz(1e9,**p)/derivative_spectrum_hz(1e9,**p40); esd_ratio=esd_per_hz(1e9,**p)/esd_per_hz(1e9,**p40)
 im=[{'metric':'I_peak_ratio_0.44_over_0.40','value':ip/peak_value(**p40),'unit':'1'},{'metric':'positive_derivative_peak_ratio','value':pp/pp40,'unit':'1'},{'metric':'negative_derivative_peak_magnitude_ratio','value':abs(nn)/abs(nn40),'unit':'1'},{'metric':'spectrum_amplitude_ratio','value':spec_ratio,'unit':'1'},{'metric':'ESD_ratio','value':esd_ratio,'unit':'1'}]
 write_csv(val/'I0_conflict_metrics.csv',im[0].keys(),im)
 make_all(cfg,figs)
 summary={'python':sys.version.split()[0],'platform':platform.platform(),'t_peak_s':tp,'I_peak_A_m':ip,'t_peak_minus_T0_s':tp-p['T0'],'positive_derivative_peak':pp,'negative_derivative_peak':nn,'quad_max_relative_error':max(x['relative_error'] for x in qr),'fft_trusted_max_relative_error':max(x['relative_error'] for x in fr if x['trusted']),'energy_consistency_relative_error':consistency,'fft_points':len(tg),'fft_resolution_Hz':res,'I0_metrics':im,'band_energies':br}
 (val/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()






