"""Small passive R-port fixture; run with the isolated openEMS Python.

No G3 transient is injected. Raw output stays under a caller-owned external root.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'python'))
from streamer_rf.fullwave.foundation import mesh_report


def run(output, case, spacing_mm, pml):
    start_time=time.perf_counter()
    if output.exists():
        raise ValueError("EXISTING_RAW_CASE_REQUIRES_EXPLICIT_REUSE")
    output.mkdir(parents=True)
    fdtd=openEMS(NrTS=20000,EndCriteria=1e-6)
    fdtd.SetGaussExcite(1e9,1e9)
    fdtd.SetBoundaryCond([f'PML_{pml}']*6)
    csx=ContinuousStructure();fdtd.SetCSX(csx)
    grid=csx.GetGrid();grid.SetDeltaUnit(1e-3)
    # PEC terminal strips connect a source port to an independent 100-ohm load.
    metal=csx.AddMetal('terminals')
    for z in [-.5,.5]:
        metal.AddBox([-3,-.5,z],[3,.5,z],priority=10)
    load=csx.AddLumpedElement('reference_load',ny=2,caps=True,R=100.)
    load.AddBox([2,-.5,-.5],[2,.5,.5],priority=20)
    port=fdtd.AddLumpedPort(1,50.,[-2,-.5,-.5],[-2,.5,.5],'z',1.,priority=20)
    for axis,lines in [('x',[-3,-2,0,2,3]),('y',[-.5,0,.5]),('z',[-.5,0,.5])]:
        grid.AddLine(axis,lines)
    grid.SmoothMeshLines('all',spacing_mm,1.4)
    for axis in 'xyz':grid.AddLine(axis,[-30,30])
    grid.SmoothMeshLines('all',2.,1.4)
    # Keep inner-domain dimensions fixed; sensitivity changes absorber thickness.
    axes=[]
    for axis in 'xyz':
        lines=np.asarray(grid.GetLines(axis),float)
        extra_left=lines[0]-np.arange(pml,0,-1)*(lines[1]-lines[0])
        extra_right=lines[-1]+np.arange(1,pml+1)*(lines[-1]-lines[-2])
        grid.AddLine(axis,np.r_[extra_left,extra_right])
        axes.append(np.asarray(grid.GetLines(axis),float))
    mesh=mesh_report([a*1e-3 for a in axes],2e9)
    probe=csx.AddDump('E_frequency',dump_type=10,file_type=1,frequency=[2e8])
    probe.AddBox([-10,0,-10],[10,0,10])
    fdtd.Run(str(output),cleanup=False,numThreads=2)
    freq=np.linspace(1e8,5e8,81)
    port.CalcPort(str(output),freq)
    zin=port.uf_tot/port.if_tot
    s11=port.uf_ref/port.uf_inc
    # Yee voltage/current probes are time-staggered; align only on shared support.
    keep=(port.u_time>=port.i_time[0])&(port.u_time<=port.i_time[-1])
    t=port.u_time[keep]
    u=port.ut_tot[keep]
    i=np.interp(t,port.i_time,port.it_tot)
    uinc=.5*(u+50*i);uref=.5*(u-50*i)
    def integrate(p):return float(np.sum(.5*(p[1:]+p[:-1])*np.diff(t)))
    energies=dict(incident_J=integrate(uinc**2/50),reflected_J=integrate(uref**2/50),accepted_J=integrate(u*i))
    idx=int(np.argmin(abs(freq-2e8)))
    columns=np.column_stack([freq,zin.real,zin.imag,s11.real,s11.imag,port.P_inc,port.P_ref,port.P_acc])
    np.savetxt(output/'port_response.csv',columns,delimiter=',',header='frequency_Hz,Z_real_ohm,Z_imag_ohm,S11_real,S11_imag,P_inc,P_ref,P_acc',comments='')
    fields={}
    for path in output.glob('E_frequency*.h5'):
        with h5py.File(path) as h:
            def inspect(name,obj):
                if isinstance(obj,h5py.Dataset) and name.startswith('FieldData'):
                    data=np.asarray(obj)
                    if data.dtype.kind in 'fci':
                        fields[name]=dict(shape=list(data.shape),maximum_abs=float(np.max(abs(data))),finite=bool(np.all(np.isfinite(data))))
            h.visititems(inspect)
    result=dict(case=case,spacing_mm=spacing_mm,pml=pml,raw_directory=str(output),
        **mesh,
        reference_frequency_Hz=float(freq[idx]),reference_load_ohm=100.,Z_real_ohm=float(zin[idx].real),Z_imag_ohm=float(zin[idx].imag),
        S11_abs=float(abs(s11[idx])),Z_relative_error=float(abs(zin[idx]-100)/100),
        response_finite=bool(np.all(np.isfinite(columns))),
        max_power_balance_relative=float(np.max(abs(port.P_inc-port.P_ref-port.P_acc))/np.max(port.P_inc)),
        max_S11_abs=float(np.max(abs(s11))),min_accepted_power=float(np.min(port.P_acc)),
        time_domain_energies=energies,voltage_probe_samples=len(port.u_time),
        current_probe_samples=len(port.i_time),simulated_time_s=float(port.u_time[-1]),
        spectral_power_semantics="PULSE_FOURIER_AMPLITUDE_SQUARED_DIAGNOSTIC_NOT_CW_WATTS",
        field_datasets=fields,runtime_s=time.perf_counter()-start_time,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        raw_bytes=sum(f.stat().st_size for f in output.rglob('*') if f.is_file()),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',required=True)
    p.add_argument('--spacing-mm',type=float,required=True)
    p.add_argument('--pml',type=int,default=8)
    a=p.parse_args();run(a.output,a.case,a.spacing_mm,a.pml)
