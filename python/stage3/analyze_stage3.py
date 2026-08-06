#!/usr/bin/env python3
"""Independent postprocessing of measured Stage 3 executable outputs."""
from __future__ import annotations
import json,sys,hashlib
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'python'))
from streamer_rf.streamer.diagnostics import axis_heads,field_criteria,field_step,head_radius,trajectory_metrics

def _sha(p:Path):return hashlib.sha256(p.read_bytes()).hexdigest()

def _field_snapshots(path:Path,hist:pd.DataFrame):
 rows=[];prev=None
 for p in sorted(path.glob('fields_*.csv'),key=lambda q:(field_step(q) is None,field_step(q) or 10**12,str(q))):
  if p.name=='fields_final.csv':continue
  step=0 if p.name=='fields_000000.csv' else field_step(p)
  if step is None or step==0:
   time=0.0
  else:
   m=hist[hist.step==step]
   if m.empty:continue
   time=float(m.time_s.iloc[-1])
  f=pd.read_csv(p);heads=axis_heads(f,previous=prev);prev=heads
  rows.append(dict(source_file=str(p.relative_to(ROOT)),source_sha256=_sha(p),step=step,time_s=time,
                   lower_head_z_m=heads['lower']['z_rho'],upper_head_z_m=heads['upper']['z_rho'],
                   lower_z_E_m=heads['lower']['z_E'],upper_z_E_m=heads['upper']['z_E'],
                   lower_z_rho_m=heads['lower']['z_rho'],upper_z_rho_m=heads['upper']['z_rho'],
                   lower_z_gradient_m=heads['lower']['z_gradient'],upper_z_gradient_m=heads['upper']['z_gradient'],
                   lower_head_rho_C_m_3=heads['lower']['charge'],upper_head_rho_C_m_3=heads['upper']['charge'],
                   lower_head_separation_m=heads['lower']['max_pairwise_separation'],
                   upper_head_separation_m=heads['upper']['max_pairwise_separation'],
                   lower_region_min_z=heads['lower']['region_min_z'],lower_region_max_z=heads['lower']['region_max_z'],
                   upper_region_min_z=heads['upper']['region_min_z'],upper_region_max_z=heads['upper']['region_max_z']))
 if (path/'fields_final.csv').exists():
  p=path/'fields_final.csv';f=pd.read_csv(p);heads=axis_heads(f,previous=prev);last=hist.iloc[-1]
  rows.append(dict(source_file=str(p.relative_to(ROOT)),source_sha256=_sha(p),step=int(last.step),time_s=float(last.time_s),
                   lower_head_z_m=heads['lower']['z_rho'],upper_head_z_m=heads['upper']['z_rho'],
                   lower_z_E_m=heads['lower']['z_E'],upper_z_E_m=heads['upper']['z_E'],
                   lower_z_rho_m=heads['lower']['z_rho'],upper_z_rho_m=heads['upper']['z_rho'],
                   lower_z_gradient_m=heads['lower']['z_gradient'],upper_z_gradient_m=heads['upper']['z_gradient'],
                   lower_head_rho_C_m_3=heads['lower']['charge'],upper_head_rho_C_m_3=heads['upper']['charge'],
                   lower_head_separation_m=heads['lower']['max_pairwise_separation'],
                   upper_head_separation_m=heads['upper']['max_pairwise_separation'],
                   lower_region_min_z=heads['lower']['region_min_z'],lower_region_max_z=heads['lower']['region_max_z'],
                   upper_region_min_z=heads['upper']['region_min_z'],upper_region_max_z=heads['upper']['region_max_z']))
 d=pd.DataFrame(rows).drop_duplicates(subset=['time_s','source_file']).sort_values('time_s')
 if len(d):d.to_csv(path/'head_trajectory.csv',index=False)
 return d

def analyze(path):
 h=pd.read_csv(path/'scalar_history.csv');f=pd.read_csv(path/'fields_final.csv');r=np.sort(f.r_m.unique());z=np.sort(f.z_m.unique());dr=float(np.median(np.diff(r)));traj=_field_snapshots(path,h)
 if len(traj)>=6:
  q=trajectory_metrics(traj,dr)
  lower_z=float(traj.lower_z_rho_m.iloc[-1]);upper_z=float(traj.upper_z_rho_m.iloc[-1])
  lower_sep=float(traj.lower_head_separation_m.iloc[-1]);upper_sep=float(traj.upper_head_separation_m.iloc[-1])
 else:
  q=trajectory_metrics(h,dr);lower_z=float(h.lower_head_z_m.iloc[-1]);upper_z=float(h.upper_head_z_m.iloc[-1]);lower_sep=np.nan;upper_sep=np.nan
 fc=field_criteria(f,dr);mid=.5*(z.min()+z.max());channel=f[(abs(f.z_m-mid)<=2*dr)&(f.r_m<=2*dr)];last=h.iloc[-1]
 return dict(run=path.name,dr_m=dr,time_s=float(last.time_s),E_max_V_m=float(last.E_max_V_m),ne_max_m_3=float(last.ne_max_m_3),channel_density_m_3=float(channel.ne_m_3.mean()),lower_head_z_m=lower_z,upper_head_z_m=upper_z,lower_head_velocity_m_s=q.get('lower_velocity_m_s',np.nan),upper_head_velocity_m_s=q.get('upper_velocity_m_s',np.nan),lower_displacement_m=q.get('lower_displacement_m',np.nan),upper_displacement_m=q.get('upper_displacement_m',np.nan),lower_monotonic_ratio=q.get('lower_monotonic_ratio',np.nan),upper_monotonic_ratio=q.get('upper_monotonic_ratio',np.nan),formed=q.get('formed',False),opposite_charge=q.get('opposite_charge',False),head_definitions_agree=fc['definitions_agree'],lower_head_separation_m=lower_sep,upper_head_separation_m=upper_sep,photoionization_ahead=fc['photoionization_ahead'],quasineutral_ratio=fc['quasineutral_ratio'],lower_radius_m=head_radius(f,lower_z),upper_radius_m=head_radius(f,upper_z),total_electrons=float(last.total_electrons),total_charge_C=float(last.total_charge_C),minimum_dt_s=float(h.dt_s.min()),maximum_dt_s=float(h.dt_s.max()),median_dt_s=float(h.dt_s.median()),max_conservation_residual=float(h.charge_source_residual.max()),steps=len(h))

def main():
 root=ROOT/'results/stage3';[(root/x).mkdir(parents=True,exist_ok=True) for x in ('closure','convergence','sensitivity','mpi')];dirs=[p for p in (root/'single_seed').iterdir() if p.is_dir() and (p/'fields_final.csv').exists() and (p/'scalar_history.csv').exists()];rows=[analyze(p) for p in dirs];summary=pd.DataFrame(rows);summary.to_csv(root/'closure/run_summary.csv',index=False)
 def select(prefix):return summary[summary.run.str.startswith(prefix)] if len(summary) else summary
 mesh=summary[summary.run.isin(['stage3_verification_baseline_20um','stage3_verification_baseline_10um','stage3_verification_baseline_5um','coarse_ml','baseline_ml','fine_ml'])].copy();mesh.to_csv(root/'convergence/mesh_convergence.csv',index=False) if len(mesh) else None
 select('nref_').to_csv(root/'sensitivity/n_ref_sensitivity.csv',index=False);select('eta_').to_csv(root/'sensitivity/open_boundary_sensitivity.csv',index=False);summary[summary.run.str.contains('sp3',case=False,regex=False)].to_csv(root/'sensitivity/photoionization_sensitivity.csv',index=False);summary[summary.run.isin(['coarse_ml','coarse_paper','coarse_ml_recovery'])].to_csv(root/'sensitivity/background_field_sensitivity.csv',index=False);select('dt_').to_csv(root/'convergence/timestep_sensitivity.csv',index=False);select('domain_').to_csv(root/'convergence/domain_sensitivity.csv',index=False);select('mpi_').to_csv(root/'mpi/rank_consistency.csv',index=False)
 metrics={'runs':rows,'maximum_conservation_residual':float(summary.max_conservation_residual.max()) if len(summary) else None};(root/'closure/stage3_summary.json').write_text(json.dumps(metrics,indent=2,default=lambda x:bool(x) if isinstance(x,np.bool_) else float(x))+'\n');print(json.dumps(metrics,indent=2,default=str))
if __name__=='__main__':main()
