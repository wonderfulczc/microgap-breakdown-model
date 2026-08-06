#!/usr/bin/env python3
"""Stage 3 executable orchestrator; Python derives metrics but never invents observations."""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,shlex,subprocess
from pathlib import Path
import numpy as np,pandas as pd,yaml,h5py
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
import sys;sys.path.insert(0,str(ROOT/'python'))
from streamer_rf.streamer.morrow_lowke_reference import evaluate
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(run_id,cmd,exe,config,ranks,output,allow_nonzero=False):
 logs=output/'provenance/logs';logs.mkdir(parents=True,exist_ok=True);so=logs/f'{run_id}.stdout.log';se=logs/f'{run_id}.stderr.log';p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True);so.write_text(p.stdout);se.write_text(p.stderr)
 if p.returncode and not allow_nonzero:raise RuntimeError(f'{run_id} failed: {se}')
 return dict(run_id=run_id,timestamp_utc=datetime.datetime.now(datetime.UTC).isoformat(),executable=str(Path(exe).relative_to(ROOT)),executable_sha256=sha(exe),config=str(Path(config).relative_to(ROOT)),config_sha256=sha(config),mpi_ranks=ranks,command=shlex.join(map(str,cmd)),exit_code=p.returncode,stdout=str(so.relative_to(ROOT)),stderr=str(se.relative_to(ROOT)))
def csv_to_h5(path):
 d=pd.read_csv(path);target=path.with_suffix('.h5')
 with h5py.File(target,'w') as h:
  for c in d.columns:h.create_dataset(c,data=d[c].to_numpy(),compression='gzip',shuffle=True)
  h.attrs['source_csv_sha256']=sha(path)
 return target
def main():
 a=argparse.ArgumentParser();a.add_argument('--config-root',default='config/stage3');a.add_argument('--build-dir',default='build');a.add_argument('--output',default='results/stage3');a.add_argument('--transport-only',action='store_true');a.add_argument('--run-name');x=a.parse_args();out=ROOT/x.output;out.mkdir(parents=True,exist_ok=True);cfg=ROOT/x.config_root/'baseline.yaml';registry=[];generated={};exe=ROOT/x.build_dir/'bin/stage3_transport';raw=out/'transport/morrow_lowke_coefficients.csv';raw.parent.mkdir(parents=True,exist_ok=True);registry.append(run('S3-TRANSPORT-001',[str(exe),str(raw)],exe,cfg,1,out));d=pd.read_csv(raw);N=101325/(1.380649e-23*300);rows=[]
 for r in d.itertuples():
  q=evaluate(r.electric_field_V_m,N)
  for key,col in [('mobility','mobility_m2_V_s'),('diffusion','diffusion_m2_s'),('ionization_frequency','nu_i_s_1'),('attachment_two_body_frequency','nu_a2_s_1'),('attachment_three_body_frequency','nu_a3_s_1')]:
   cv=getattr(r,col);rows.append(dict(electric_field_V_m=r.electric_field_V_m,quantity=key,cpp=cv,python=q[key],relative_difference=abs(cv-q[key])/max(abs(cv),abs(q[key]),1.)))
 comp=pd.DataFrame(rows);comp.to_csv(out/'transport/cpp_python_comparison.csv',index=False);text=(ROOT/registry[0]['stdout']).read_text();ek=float(text.split('breakdown_field_V_m=')[1].split()[0]);pd.DataFrame([dict(breakdown_field_V_m=ek,paper_field_V_m=3.2e6,relative_difference=abs(ek-3.2e6)/3.2e6)]).to_csv(out/'transport/breakdown_field.csv',index=False)
 fig,ax=plt.subplots();ax.loglog(d.E_over_N_Td,d.nu_i_s_1,label='ionization');ax.loglog(d.E_over_N_Td,d.nu_a2_s_1+d.nu_a3_s_1,label='attachment');ax.set_xlabel('E/N (Td)');ax.set_ylabel('frequency (s$^{-1}$)');ax.legend();fig.tight_layout();fig.savefig(out/'transport/coefficient_curves.png');fig.savefig(out/'transport/coefficient_curves.pdf');plt.close(fig);generated['S3-TRANSPORT-001']=list((out/'transport').glob('*'))
 if not x.transport_only:
  specs=yaml.safe_load((ROOT/x.config_root/'runs.yaml').read_text())['runs'];specs=[s for s in specs if x.run_name is None or s['name']==x.run_name]
  if x.run_name and not specs:raise ValueError(f'unknown run {x.run_name}')
  solver=ROOT/x.build_dir/'bin/stage3_run'
  for s in specs:
   nr=round(float(s['r_max_m'])/float(s['dr_m']));nz=round(float(s['z_max_m'])/float(s['dr_m']));n0=float(s.get('n0',1e20));field=float(s.get('field_factor',1.5))*ek if s['field']=='ml' else float(s.get('field_value',4.8e6));dest=out/'single_seed'/s['name'];dest.mkdir(parents=True,exist_ok=True);cmd=['mpirun','-np',str(s['mpi_ranks']),str(solver),str(dest),str(nr),str(nz),str(float(s['t_end_s'])),str(field),str(n0*float(s['n_ref_ratio'])),str(int(s['sp3'])),str(float(s['eta'])),str(float(s['r_max_m'])),str(float(s['z_max_m'])),str(float(s.get('dt_scale',1.0))),str(n0),str(float(s.get('sigma_m',1e-4))),str(int(s.get('max_steps',50000))),str(float(s.get('physical_stop_after_s',0.0)))]
   has_resume=bool(s.get('resume_checkpoint'))
   if has_resume:
    cmd += [str(ROOT/s['resume_checkpoint']),str(int(s.get('initial_step',0)))]
   if 'minimum_dt_s' in s:
    if not has_resume:
     cmd += ['', '0']
    cmd += [str(float(s['minimum_dt_s']))]
   if 'output_interval_s' in s:
    if 'minimum_dt_s' not in s:
     if not has_resume:
      cmd += ['', '0']
     cmd += [str(float(s.get('minimum_dt_s',1e-18)))]
    cmd += [str(float(s['output_interval_s']))]
   entry=run('S3-'+s['name'].upper().replace('-','_').replace('.','_'),cmd,solver,ROOT/x.config_root/'runs.yaml',s['mpi_ranks'],out,allow_nonzero=True);registry.append(entry)
   for p in dest.glob('fields_*.csv'):csv_to_h5(p)
   hist=pd.read_csv(dest/'scalar_history.csv');hist.assign(lower_head_velocity_m_s=np.gradient(hist.lower_head_z_m,hist.time_s),upper_head_velocity_m_s=np.gradient(hist.upper_head_z_m,hist.time_s)).to_csv(dest/'head_trajectory.csv',index=False);hist[['time_s','charge_source_residual','total_electrons','total_charge_C']].to_csv(dest/'conservation.csv',index=False);hist[['time_s','dt_s','controller']].to_csv(dest/'timestep_controller.csv',index=False);hist[['time_s','poisson_iterations','sp3_ksp_iterations','sp3_boundary_iterations']].to_csv(dest/'solver_iterations.csv',index=False);generated[entry['run_id']]=list(dest.glob('*'))
 reg=[]
 for base in registry:
  for p in generated[base['run_id']]:
   if p.is_file():reg.append({**base,'result_files':str(p.relative_to(ROOT)),'result_sha256':sha(p)})
 path=out/'provenance/run_registry.csv';path.parent.mkdir(parents=True,exist_ok=True);old=pd.read_csv(path).to_dict('records') if path.exists() else [];replaced={x['run_id'] for x in registry};old=[r for r in old if r['run_id'] not in replaced];pd.DataFrame(old+reg).to_csv(path,index=False);print(json.dumps(dict(breakdown_field_V_m=ek,max_cpp_python_difference=comp.relative_difference.max()),indent=2))
if __name__=='__main__':main()
