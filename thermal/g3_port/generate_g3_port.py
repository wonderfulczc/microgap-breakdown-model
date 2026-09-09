"""Recover circuit-resolution data from frozen accepted G2 intervals, not thermal replay."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"python"))
from streamer_rf.circuit.model import CircuitParameters, SeriesRLCGapCircuit
from streamer_rf.circuit.signals import ConductanceProfile, PiecewiseLinearSignal
from streamer_rf.thermal.port import PortTransientContract, port_frame, energy_metrics, dynamic_impedance


def main():
    start=time.perf_counter()
    out=Path(__file__).parent
    source=ROOT/"thermal/g2/g2_primary_timeseries.csv"
    summary_path=ROOT/"thermal/g2/g2_reference_summary.json"
    g2=json.loads(summary_path.read_text())
    if g2["framework_status"]!="BIDIRECTIONAL_REFERENCE_COUPLING_VALIDATED":
        raise ValueError("G2_REFERENCE_NOT_VALIDATED")
    p=g2["parameters"]
    saved=pd.read_csv(source)
    native=[]
    pieces=[]
    errors=[]
    for n in range(1,len(saved)):
        prev,row=saved.iloc[n-1],saved.iloc[n]
        t0,t1=float(prev.time_s),float(row.time_s)
        params=CircuitParameters(p["L_H"],p["Cext_F"],p["Cgap_F"],p["Rs_ohm"],
            ConductanceProfile.constant(1/row.interval_Rsp_ohm),PiecewiseLinearSignal.constant(p["Vs_V"]))
        kernel=SeriesRLCGapCircuit(params)
        sol=solve_ivp(kernel.rhs,(t0,t1),[prev.I_L_A,prev.V_channel_V],method="DOP853",
            rtol=2e-10,atol=[1e-13,1e-11],max_step=(t1-t0)/4,dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        errors.append(np.abs(sol.y[:,-1]-[row.I_L_A,row.V_channel_V]))
        f=port_frame(kernel.evaluate_solution(sol.t,sol.y[0],sol.y[1]))
        f["interval_index"]=n
        f["sample_side"]="INTERIOR"
        f.loc[f.index[0],"sample_side"]="RIGHT"
        f.loc[f.index[-1],"sample_side"]="LEFT"
        native.append(f)
        pieces.append((t0,t1,kernel,sol))
    original=pd.concat(native,ignore_index=True)
    spacing=np.diff(original.time_s.to_numpy())
    # Ignore zero-width jump pairs and the roundoff-sized last G2 remainder.
    meaningful=spacing[spacing>100*np.spacing(float(saved.time_s.iloc[-1]))]
    typical=float(np.median(meaningful))
    duration=float(saved.time_s.iloc[-1]-saved.time_s.iloc[0])
    count=int(np.ceil(duration/typical))
    times=np.linspace(saved.time_s.iloc[0],saved.time_s.iloc[-1],count+1)
    frames=[]
    for n,(t0,t1,kernel,sol) in enumerate(pieces):
        select=(times>=t0)&((times<t1) if n<len(pieces)-1 else (times<=t1))
        t=times[select]
        if not len(t):
            continue
        y=sol.sol(t)
        frames.append(port_frame(kernel.evaluate_solution(t,y[0],y[1])))
    uniform=pd.concat(frames,ignore_index=True)
    # Use one common gate per exported waveform, not a different peak per interval.
    for frame in (original,uniform):
        z,valid,gate=dynamic_impedance(frame.V_port_V,frame.I_port_A)
        frame["Z_dynamic_ohm"],frame["Z_dynamic_valid"],frame["impedance_current_gate_A"]=z,valid,gate
    contract=PortTransientContract("G2_primary_150_LAGGED",uniform,duration/count,
        {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (source,summary_path)}).validate()
    a,b=energy_metrics(original),energy_metrics(uniform)
    comparison={k:abs(b[k]-a[k])/abs(a[k]) for k in ("port_energy_J","conduction_energy_J")}
    for key in ("V_port_V","I_port_A"):
        comparison["peak_"+key]=abs(uniform[key].abs().max()-original[key].abs().max())/original[key].abs().max()
    q=g2["results"]["primary"]["Q_sp_circuit_J"]
    errors=np.asarray(errors)
    endpoint=dict(max_I_L_error_A=float(errors[:,0].max()),max_V_error_V=float(errors[:,1].max()))
    passed=(max(comparison.values())<=.01 and b["energy_residual_relative"]<=.01 and
        abs(b["conduction_energy_J"]-q)/q<=.01 and endpoint["max_V_error_V"]<1e-7 and endpoint["max_I_L_error_A"]<1e-9)
    contract["G3_status"]="PHYSICS_DERIVED_PORT_MODEL_VALIDATED" if passed else "PORT_CHECKS_FAILED"
    contract["energy_closure_status"]="PASS" if passed else "MINIMAL_FIX_REQUIRED"
    contract["valid_time_interval_s"]=[float(times[0]),float(times[-1])]
    contract["uniform_filename"]="g3_port_uniform.csv"
    result=dict(contract=contract,endpoint_recovery=endpoint,native_energy=a,uniform_energy=b,
        resampling_relative_errors=comparison,conduction_vs_G2_relative=abs(b["conduction_energy_J"]-q)/q,
        native_median_internal_dt_s=typical,native_rows=len(original),
        uniform_ranges={k:[float(uniform[k].min()),float(uniform[k].max())] for k in
            ("V_port_V","I_port_A","I_conduction_A","I_gap_displacement_A","Rsp_ohm","P_port_W")},
        impedance_valid_fraction=float(uniform.Z_dynamic_valid.mean()),
        impedance_p10_median_p90_ohm=uniform.loc[uniform.Z_dynamic_valid,"Z_dynamic_ohm"].quantile([.1,.5,.9]).tolist(),
        KCL_max_A=float(uniform.KCL_residual_A.abs().max()),
        method_status="PHYSICS_DERIVED_PORT_MODEL_VALIDATED" if passed else "PORT_CHECKS_FAILED",
        energy_closure_status="PASS" if passed else "MINIMAL_FIX_REQUIRED",
        runtime_s=time.perf_counter()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    original.to_csv(out/"g3_port_timeseries.csv",index=False)
    uniform.to_csv(out/"g3_port_uniform.csv",index=False)
    (out/"g3_port_summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
