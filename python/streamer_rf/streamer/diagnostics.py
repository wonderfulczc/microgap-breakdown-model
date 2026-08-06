from __future__ import annotations
import re
import numpy as np,pandas as pd
from scipy.signal import find_peaks
E_CHARGE=1.602176634e-19

def _parabolic_position(z:np.ndarray,y:np.ndarray,idx:int)->float:
 if idx<=0 or idx>=len(y)-1:return float(z[idx])
 ym,y0,yp=float(y[idx-1]),float(y[idx]),float(y[idx+1]);den=ym-2*y0+yp
 if abs(den)<1e-300:return float(z[idx])
 delta=.5*(ym-yp)/den
 delta=float(np.clip(delta,-1,1))
 return float(z[idx]+delta*(z[1]-z[0]))

def _contiguous_runs(indices:np.ndarray):
 if len(indices)==0:return []
 cuts=np.flatnonzero(np.diff(indices)>1)+1
 return np.split(indices,cuts)

def _front_for_side(z,E,rho,ne,side,z0,sigma,E0,rho_fraction,previous_z=None):
 if side=='lower':
  side_mask=z<z0-2*sigma
 else:
  side_mask=z>z0+2*sigma
 ix=np.flatnonzero(side_mask)
 if len(ix)==0:return None
 grad=np.abs(np.gradient(ne,z))
 side_peak=max(float(np.max(np.abs(rho[ix]))),1e-300)
 grad_peak=max(float(np.max(grad[ix])),1e-300)
 mask=side_mask&(E>E0)&(np.abs(rho)>=rho_fraction*side_peak)&(grad>=0.02*grad_peak)
 runs=_contiguous_runs(np.flatnonzero(mask))
 if not runs:
  mask=side_mask&(np.abs(rho)>=rho_fraction*side_peak)
  runs=_contiguous_runs(np.flatnonzero(mask))
 if not runs:return None
 if previous_z is not None:
  region=min(runs,key=lambda r:abs(float(np.mean(z[r]))-previous_z))
 else:
  region=max(runs,key=lambda r:float(np.max(np.abs(rho[r]))))
 er=region[np.argmax(E[region])]
 rr=region[np.argmax(np.abs(rho[region]))]
 gr=region[np.argmax(grad[region])]
 zE=_parabolic_position(z,E,er);zr=_parabolic_position(z,np.abs(rho),rr);zg=_parabolic_position(z,grad,gr)
 sep=max(abs(zE-zr),abs(zE-zg),abs(zr-zg))
 return dict(region_min_z=float(z[region[0]]),region_max_z=float(z[region[-1]]),
             z_E=float(zE),z_rho=float(zr),z_gradient=float(zg),
             z_E_grid=float(z[er]),z_rho_grid=float(z[rr]),z_gradient_grid=float(z[gr]),
             charge=float(rho[rr]),field=float(E[er]),gradient=float(grad[gr]),
             max_pairwise_separation=float(sep),region_cell_count=int(len(region)))

def axis_heads(field:pd.DataFrame,z0:float|None=None,sigma:float=1e-4,E0:float|None=None,rho_fraction:float=.2,previous:dict|None=None):
 a=field[field.i==field.i.min()].sort_values('z_m');z=a.z_m.to_numpy();E=a.E_V_m.to_numpy();rho=a.rho_C_m_3.to_numpy();ne=a.ne_m_3.to_numpy() if 'ne_m_3' in a else np.ones_like(z)
 if z0 is None:z0=.5*(float(z[0])+float(z[-1]))
 if E0 is None:E0=float(np.nanpercentile(E,10))
 out={}
 for side in ('lower','upper'):
  prev=None if previous is None or side not in previous else previous[side].get('z_rho')
  h=_front_for_side(z,E,rho,ne,side,z0,sigma,E0,rho_fraction,prev)
  if h is None:
   mask=z<z0 if side=='lower' else z>z0;ix=np.flatnonzero(mask);qr=ix[np.argmax(np.abs(rho[ix]))];h=dict(region_min_z=float(z[qr]),region_max_z=float(z[qr]),z_E=float(z[qr]),z_rho=float(z[qr]),z_gradient=float(z[qr]),z_E_grid=float(z[qr]),z_rho_grid=float(z[qr]),z_gradient_grid=float(z[qr]),charge=float(rho[qr]),field=float(E[qr]),gradient=0.0,max_pairwise_separation=0.0,region_cell_count=1)
  # Compatibility keys used by older tests and reports.
  h['field_z']=h['z_E'];h['charge_z']=h['z_rho'];h['definition_difference']=abs(h['z_E']-h['z_rho'])
  out[side]=h
 return out

def head_radius(field:pd.DataFrame,z_head:float,quantity='ne_m_3'):
 zs=np.sort(field.z_m.unique());z=zs[np.argmin(abs(zs-z_head))];d=field[field.z_m==z].sort_values('r_m');y=np.abs(d[quantity].to_numpy());r=d.r_m.to_numpy();half=.5*y.max();below=np.flatnonzero(y<=half);return float(r[below[0]]) if len(below) else float(r[-1])

def trajectory_metrics(history:pd.DataFrame,dr:float):
 if len(history)<6:return dict(formed=False,reason='insufficient history')
 lo=np.asarray(history.lower_head_z_m.to_numpy(),dtype=float);hi=np.asarray(history.upper_head_z_m.to_numpy(),dtype=float);t=np.asarray(history.time_s.to_numpy(),dtype=float)
 good=np.isfinite(lo)&np.isfinite(hi)&np.isfinite(t)
 lo=lo[good];hi=hi[good];t=t[good]
 if len(t):
  _, rev_keep=np.unique(t[::-1],return_index=True)
  keep=np.sort(len(t)-1-rev_keep)
  lo=lo[keep];hi=hi[keep];t=t[keep]
 if len(t)<6 or len(np.unique(t))<2:return dict(formed=False,reason='insufficient valid trajectory')
 vlo=np.gradient(lo,t);vhi=np.gradient(hi,t);disp_lo=abs(lo[-1]-lo[0]);disp_hi=abs(hi[-1]-hi[0])
 opposite=np.sign(history.lower_head_rho_C_m_3.iloc[-1])==-np.sign(history.upper_head_rho_C_m_3.iloc[-1])
 tail=max(5,len(t)//3)
 lower_mono=float(np.mean(np.diff(lo[-tail:])<=0)) if tail>1 else 0.0
 upper_mono=float(np.mean(np.diff(hi[-tail:])>=0)) if tail>1 else 0.0
 lower_v=float(np.median(vlo[-tail:]));upper_v=float(np.median(vhi[-tail:]))
 threshold=max(5*dr,1e-4)
 sustained=opposite and disp_lo>=threshold and disp_hi>=threshold and lower_mono>=.8 and upper_mono>=.8 and abs(lower_v)>1e4 and abs(upper_v)>1e4
 return dict(formed=bool(sustained),lower_displacement_m=float(disp_lo),upper_displacement_m=float(disp_hi),lower_velocity_m_s=lower_v,upper_velocity_m_s=upper_v,opposite_charge=bool(opposite),lower_monotonic_ratio=lower_mono,upper_monotonic_ratio=upper_mono)

def field_criteria(field:pd.DataFrame,dr:float,z0:float|None=None,sigma:float=1e-4,E0:float|None=None):
 h=axis_heads(field,z0=z0,sigma=sigma,E0=E0);mid=.5*(field.z_m.min()+field.z_m.max());channel=field[(abs(field.z_m-mid)<=2*dr)&(field.r_m<=2*dr)];quasi=float(np.max(np.abs(channel.rho_C_m_3))/(E_CHARGE*np.max(channel.ne_m_3)+1e-300));positive=max((h[k] for k in h),key=lambda q:q['charge']);ahead=field[(field.z_m>positive['z_rho']+dr)&(field.r_m<=2*dr)];ok=[]
 for side,v in h.items():
  radius=max(head_radius(field,v['z_rho']),dr)
  ok.append(v['max_pairwise_separation']<=max(5*dr,.5*radius) and v['region_cell_count']>=1)
 return dict(heads=h,quasineutral_ratio=quasi,photoionization_ahead=float(ahead.Sph_m_3_s_1.max()) if len(ahead) else 0.,definitions_agree=all(ok),max_head_separation=max(v['max_pairwise_separation'] for v in h.values()))

def field_step(path):
 m=re.search(r'fields_(\d+)\.csv$',str(path))
 return int(m.group(1)) if m else None
