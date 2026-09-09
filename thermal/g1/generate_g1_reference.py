"""One synthetic thermal case at three grids; compact evidence only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"python"))
from streamer_rf.thermal.properties import ATM, EquilibriumAirProperties
from streamer_rf.thermal.radial import RadialGrid, RadialThermalSolver, thermal_initialization_from_g0


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,default=ROOT/"thermal/g1")
    args=parser.parse_args()
    out=args.out_dir
    out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    air=EquilibriumAirProperties()
    temperatures=np.array([300,1000,3000,5000,10000,15000,30000])
    values=air.evaluate(temperatures,ATM)
    sanity=pd.DataFrame(dict(T_K=temperatures,**values))
    sanity.to_csv(out/"g1_property_sanity.csv",index=False)
    dh=(air.evaluate(temperatures+.001,ATM)["h_J_kg"]-air.evaluate(temperatures-.001,ATM)["h_J_kg"])/.002
    inverse_T,inverse_p=air.invert(values["rho_kg_m3"],values["e_J_kg"],temperatures*.9,ATM*1.1)
    property_check=dict(all_selected_values_valid=bool(np.all(values["valid"])),
        max_cp_vs_dh_relative=float(np.max(np.abs(values["cp_J_kgK"]-dh)/dh)),
        max_inverse_T_relative=float(np.max(np.abs(inverse_T-temperatures)/temperatures)),
        max_inverse_p_relative=float(np.max(np.abs(inverse_p-ATM)/ATM)),
        published_figure_check="Scale-only check against indexed Figures6/7 axes; no digitization. Image endpoint unavailable; pointwise curve comparison not claimed.")
    t_end=1e-8
    grid_rows=[]
    runs=[]
    for n in (100,150,225):
        t0=time.perf_counter()
        grid=RadialGrid(n,5e-4)
        f=np.exp(-(grid.centers/5e-5)**2)
        drive=lambda t: .05*(1+.5*np.sin(np.pi*t/t_end)**2)
        solver=RadialThermalSolver(grid,1e-3,300+9700*f,ATM*(1+f),properties=air,current=drive)
        records=[solver.diagnostics()]
        initial_p=solver.p.copy()
        if n == 150:
            pd.DataFrame(solver.profile()).to_csv(out/"g1_profile_initial.csv",index=False)
        for label,target in (("mid",t_end/2),("final",t_end)):
            while solver.time < target:
                solver.step(target-solver.time)
                records.append(solver.diagnostics())
                if solver.steps > 20000:
                    raise RuntimeError("reference exceeded bounded step budget")
            d=solver.diagnostics()
            grid_rows.append(dict(radial_cells=n,sample=label,**{k:d[k] for k in
                ("time_s","T_axis_K","channel_radius_m","Rsp_ohm","total_energy_J")}))
            if n==150:
                pd.DataFrame(solver.profile()).to_csv(out/f"g1_profile_{label}.csv",index=False)
        frame=pd.DataFrame(records)
        disturbance=np.abs(solver.p-initial_p)
        outer=grid.centers > .8*grid.radius_m
        run=dict(radial_cells=n,accepted_steps=solver.steps,rejected_steps=solver.rejected_steps,
                 runtime_s=time.perf_counter()-t0,final=solver.diagnostics(),
                 max_joule_closure=float(frame.joule_closure_relative.max()),
                 max_energy_residual_relative=float(np.max(np.abs(frame.energy_residual_J))/solver.initial_energy),
                 max_energy_residual_over_Joule=float(np.max(np.abs(frame.energy_residual_J))/solver.QJ),
                 outer_shell_pressure_change_fraction=float(np.max(disturbance[outer])/ATM),
                 radius_initial_m=float(frame.channel_radius_m.iloc[0]),
                 radius_final_m=float(frame.channel_radius_m.iloc[-1]),
                 channel_radius_max_step_change_m=float(np.max(np.abs(np.diff(frame.channel_radius_m)))))
        runs.append(run)
        if n==150:
            frame.to_csv(out/"g1_reference_timeseries.csv",index=False)
    pd.DataFrame(grid_rows).to_csv(out/"g1_grid_sensitivity.csv",index=False)
    sensitivity={}
    for key in ("T_axis_K","channel_radius_m","Rsp_ohm","total_energy_J"):
        v=[r["final"][key] for r in runs]
        sensitivity[key]=dict(values=v,relative_medium_fine=abs(v[2]-v[1])/abs(v[2]),
                              coarse_medium_change=abs(v[1]-v[0]),medium_fine_change=abs(v[2]-v[1]))
    g0=json.loads((ROOT/"circuit/g0_handoff/g0_handoff_summary.json").read_text())
    source_files=[ROOT/"python/streamer_rf/thermal"/name for name in ("properties.py","radial.py","dangola_coefficients.json")]
    summary=dict(model_status="LTE_REFERENCE_MODEL",LTE_applicability="LTE_APPLICABILITY_PENDING_CALIBRATION",
        radiation_status="RADIATION_MODEL_NOT_ENABLED",initial_condition_status="NUMERICAL_REFERENCE_INITIAL_CONDITION",
        property_source="PUBLISHED_EQUILIBRIUM_AIR_FIT",coefficient_sha256=air.coefficient_sha256,
        reference=dict(ambient_pressure_Pa=ATM,ambient_temperature_K=300.,Rmax_m=5e-4,L_channel_m=1e-3,
                       initial_T="300 + 9700 exp(-(r/50um)^2) K",initial_p="101325*(1+exp(-(r/50um)^2)) Pa",
                       current="0.05*(1+0.5*sin(pi*t/10ns)^2) A",end_time_s=t_end,primary_cells=150),
        property_check=property_check,
        g0_initialization=thermal_initialization_from_g0(g0),Q_required_status="NOT_AVAILABLE",
        runs=runs,grid_sensitivity=sensitivity,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        runtime_s=time.perf_counter()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    usable=(all(r["max_energy_residual_relative"]<1e-8 and r["max_joule_closure"]<1e-12 and
                r["final"]["outward_velocity_max_m_s"]>0 and r["radius_final_m"]>r["radius_initial_m"] for r in runs)
            and all(v["relative_medium_fine"]<.1 or v["medium_fine_change"]<v["coarse_medium_change"] for v in sensitivity.values()))
    summary["reference_verdict"]="PASS" if usable else "MINIMAL_FIX_REQUIRED"
    summary["output_size_bytes"]=sum(p.stat().st_size for p in out.glob("g1_*.csv"))
    (out/"g1_reference_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps(summary,indent=2,allow_nan=False))
    return 0 if usable else 1


if __name__=="__main__":
    raise SystemExit(main())
