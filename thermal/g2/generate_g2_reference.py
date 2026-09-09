"""One canonical/G1 reference, plus only the requested bounded comparisons."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.circuit.model import CircuitParameters
from streamer_rf.circuit.signals import ConductanceProfile, PiecewiseLinearSignal
from streamer_rf.thermal.radial import RadialGrid, RadialThermalSolver
from streamer_rf.thermal.properties import ATM
from streamer_rf.thermal.coupling import ThermalCircuitCoupling


def reference(cells=150, interval=5e-11, feedback=True):
    grid = RadialGrid(cells, 5e-4)
    f = np.exp(-(grid.centers/5e-5)**2)
    thermal = RadialThermalSolver(grid, 1e-3, 300+9700*f, ATM*(1+f))
    params = CircuitParameters(1e-6, .8e-12, .2e-12, 2.,
        ConductanceProfile.constant(2e-4), PiecewiseLinearSignal.constant(0.))
    return ThermalCircuitCoupling(thermal, params, (0., 100.), interval, feedback=feedback)


def metrics(s):
    rows = pd.DataFrame(s.rows)
    final = s.rows[-1]
    denom = max(s.initial_circuit_energy, abs(s.Wsource))
    return dict(peak_I_A=s.peak_current, final_T_axis_K=final["T_axis_K"],
        final_Rsp_ohm=final["Rsp_ohm"], Q_sp_circuit_J=s.Qsp, Q_sp_thermal_J=s.thermal.QJ,
        circuit_internal_steps=s.internal_steps, circuit_rhs_evaluations=s.nfev,
        thermal_accepted_steps=s.thermal.steps, coupling_steps=len(s.rows)-1,
        discarded_intervals=s.discarded_intervals, final_V_channel_V=final["V_channel_V"],
        circuit_energy_residual_max_J=float(rows.circuit_energy_residual_J.abs().max()),
        circuit_energy_residual_normalized=float(rows.circuit_energy_residual_J.abs().max()/denom),
        cross_domain_relative=abs(s.Qsp-s.thermal.QJ)/s.Qsp if s.Qsp else 0.,
        thermal_energy_residual_max_J=float(rows.energy_residual_J.abs().max()),
        final_radius_m=final["channel_radius_m"], final_sigma_axis_S_m=final["sigma_axis_S_m"],
        final_total_energy_J=final["total_energy_J"],
        minimum_T_K=float(rows.T_min_K.min()), minimum_p_Pa=float(rows.p_min_Pa.min()),
        minimum_rho_kg_m3=float(rows.rho_min_kg_m3.min()),
        max_interval_power_relative=float(((rows.P_sp_circuit_W-rows.P_sp_thermal_W).abs()/
            rows.P_sp_circuit_W.where(rows.P_sp_circuit_W>0)).max()))


def main():
    start = time.perf_counter()
    out = Path(__file__).parent
    results, frames = {}, {}
    for name, cells, dt, feedback in [("primary",150,5e-11,True),
            ("half_interval",150,2.5e-11,True), ("finer_grid",225,5e-11,True),
            ("one_way",150,5e-11,False)]:
        begin = time.perf_counter()
        s = reference(cells, dt, feedback)
        if name == "primary":
            pd.DataFrame(s.thermal.profile()).to_csv(out/"g2_profile_initial.csv",index=False)
        s.run(1e-8)
        frames[name] = pd.DataFrame(s.rows)
        frames[name].to_csv(out/f"g2_{name}_timeseries.csv",index=False)
        results[name] = metrics(s)
        results[name].update(runtime_s=time.perf_counter()-begin, cells=cells, coupling_cap_s=dt)
        results[name]["process_peak_RSS_KiB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if name == "primary":
            pd.DataFrame(s.thermal.profile()).to_csv(out/"g2_profile_final.csv",index=False)
        print(name, results[name], flush=True)
    comparisons = {}
    for name in ("half_interval","finer_grid","one_way"):
        comparisons[name] = {key: abs(results[name][key]-results["primary"][key])/abs(results["primary"][key])
            for key in ("peak_I_A","Q_sp_circuit_J","final_T_axis_K","final_Rsp_ohm")}
        a, b = frames["primary"], frames[name]
        interp = np.interp(a.time_s,b.time_s,b.I_A)
        comparisons[name]["current_waveform_normalized_L2"] = float(np.linalg.norm(interp-a.I_A)/np.linalg.norm(a.I_A))
    frozen = ["python/streamer_rf/thermal/radial.py", "python/streamer_rf/thermal/properties.py",
              "python/streamer_rf/thermal/dangola_coefficients.json", "python/streamer_rf/circuit/model.py",
              "python/streamer_rf/circuit/signals.py", "thermal/g1/g1_reference_summary.json"]
    summary = dict(reference_status="NUMERICAL_REFERENCE_INITIAL_CONDITION", mode="LAGGED",
        circuit_provenance="circuit/validation/generate_stage_g1_validation.py: conductive_gap",
        thermal_provenance="thermal/g1/generate_g1_reference.py",
        parameters=dict(L_H=1e-6,Cext_F=.8e-12,Cgap_F=.2e-12,Rs_ohm=2.,Vs_V=0.,
            initial_Vgap_V=100.,initial_I_L_A=0.,Rmax_m=5e-4,L_channel_m=1e-3,end_time_s=1e-8),
        results=results,comparisons=comparisons,
        status_flags=["HANDOFF_CALIBRATION_PENDING","G0_THERMAL_INITIALIZATION_UNRESOLVED","THERMAL_REFERENCE_SOLVER_VALIDATED",
            "PRODUCTION_THERMAL_RLC_COUPLING_NOT_CALIBRATED","LTE_APPLICABILITY_PENDING_CALIBRATION",
            "RADIATION_MODEL_NOT_ENABLED","THERMAL_RSP_GRID_SENSITIVITY_PRESENT"],
        frozen_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in frozen},
        runtime_s=time.perf_counter()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    a = results["primary"]
    passed = (a["cross_domain_relative"] <= .01 and a["max_interval_power_relative"] < 1e-12
        and a["circuit_energy_residual_normalized"] < 1e-7
        and a["thermal_energy_residual_max_J"]/a["Q_sp_thermal_J"] < 1e-6
        and all(a[key] > 0 for key in ("minimum_T_K","minimum_p_Pa","minimum_rho_kg_m3")) and all(
        comparisons["half_interval"][key] <= .05 for key in
        ("peak_I_A","Q_sp_circuit_J","final_T_axis_K","final_Rsp_ohm")))
    summary["coupling_acceptance"] = "PASS" if passed else "MINIMAL_FIX_REQUIRED"
    summary["framework_status"] = "BIDIRECTIONAL_REFERENCE_COUPLING_VALIDATED" if passed else "COUPLING_NOT_RESOLVED"
    (out/"g2_reference_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
