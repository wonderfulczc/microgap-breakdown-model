"""G3 gap-side port diagnostics, using the canonical circuit branch currents."""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


def dynamic_impedance(voltage, current, absolute_floor_A=1e-12, relative_floor=1e-8):
    v, i = np.asarray(voltage, float), np.asarray(current, float)
    if (v.shape != i.shape or v.ndim != 1 or not v.size or
            not np.all(np.isfinite(v)) or not np.all(np.isfinite(i)) or
            not np.isfinite(absolute_floor_A) or absolute_floor_A < 0 or
            not np.isfinite(relative_floor) or not 0 <= relative_floor < 1):
        raise ValueError("INVALID_PORT_WAVEFORM_OR_GATE")
    gate = max(absolute_floor_A, relative_floor*float(np.max(np.abs(i))))
    valid = np.abs(i) > gate
    with np.errstate(over="ignore", invalid="ignore"):
        z = np.divide(v, i, out=np.full_like(v, np.nan), where=valid)
    valid &= np.isfinite(z)
    z[~valid] = np.nan
    return z, valid, gate


def port_frame(solution):
    c = solution.currents
    total = solution.I_L_A-c.I_Cext_A
    z, valid, gate = dynamic_impedance(solution.Vgap_V, total)
    g = solution.Gb_S
    return pd.DataFrame(dict(time_s=solution.time_s, V_port_V=solution.Vgap_V,
        I_port_A=total, I_conduction_A=c.I_gap_cond_A,
        I_gap_displacement_A=c.I_gap_disp_A, I_external_capacitor_A=c.I_Cext_A,
        I_series_A=solution.I_L_A,
        Rsp_ohm=np.divide(1.,g,out=np.full_like(g,np.nan),where=g>0),
        Z_dynamic_ohm=z,Z_dynamic_valid=valid,impedance_current_gate_A=gate,
        waveform_valid=np.ones(len(total),bool),
        KCL_residual_A=total-c.I_gap_cond_A-c.I_gap_disp_A,
        P_port_W=solution.Vgap_V*total,
        P_conduction_W=solution.Vgap_V*c.I_gap_cond_A,
        E_Cgap_J=solution.cgap_energy_J))


def energy_metrics(frame):
    t=frame.time_s.to_numpy()
    if len(t)<2 or np.any(~np.isfinite(t)) or np.any(np.diff(t)<0) or t[-1]<=t[0]:
        raise ValueError("INVALID_PORT_TIME")
    # Repeated native timestamps represent opposite sides of a held-R jump.
    # Their zero-width trapezoid does not blend a jump across an interval.
    def integral(values):
        v=np.asarray(values,float)
        if not np.all(np.isfinite(v)):
            raise ValueError("NONFINITE_PORT_POWER")
        return float(np.sum(.5*(v[1:]+v[:-1])*np.diff(t)))
    port=integral(frame.P_port_W)
    cond=integral(frame.P_conduction_W)
    delta=float(frame.E_Cgap_J.iloc[-1]-frame.E_Cgap_J.iloc[0])
    residual=port-cond-delta
    scale=max(abs(port),abs(cond),abs(delta))
    return dict(port_energy_J=port,conduction_energy_J=cond,Cgap_energy_change_J=delta,
                energy_residual_J=residual,energy_residual_relative=abs(residual)/scale if scale else 0.)


@dataclass
class PortTransientContract:
    source_case_id: str
    waveform: pd.DataFrame
    dt_port_s: float
    metadata: dict = field(default_factory=dict)

    def validate(self):
        f=self.waveform
        core=["time_s","V_port_V","I_port_A","I_conduction_A","I_gap_displacement_A",
              "I_external_capacitor_A","I_series_A"]
        required=core+["Rsp_ohm","Z_dynamic_ohm","Z_dynamic_valid"]
        if not self.source_case_id or any(k not in f for k in required) or len(f)<2:
            raise ValueError("INVALID_PORT_CONTRACT")
        if not np.all(np.isfinite(f[core].to_numpy())) or not np.isfinite(self.dt_port_s) or self.dt_port_s<=0:
            raise ValueError("INVALID_PORT_CONTRACT")
        if not np.allclose(np.diff(f.time_s),self.dt_port_s,rtol=1e-9,atol=0):
            raise ValueError("PORT_NOT_UNIFORMLY_SAMPLED")
        scale=max(float(np.max(np.abs(f.I_port_A))),1e-300)
        if np.max(np.abs(f.I_port_A-f.I_conduction_A-f.I_gap_displacement_A))>1e-12*scale:
            raise ValueError("PORT_KCL_FAILED")
        if np.max(np.abs(f.I_port_A-f.I_series_A+f.I_external_capacitor_A))>1e-12*scale:
            raise ValueError("PORT_REFERENCE_PLANE_FAILED")
        r=f.Rsp_ohm.to_numpy()
        if np.any(np.isinf(r)) or np.any(r[np.isfinite(r)]<=0) or np.any(np.isnan(r)&(f.I_conduction_A!=0)):
            raise ValueError("INVALID_CHANNEL_RESISTANCE")
        z, valid, _=dynamic_impedance(f.V_port_V,f.I_port_A)
        if not np.array_equal(f.Z_dynamic_valid,valid) or not np.allclose(f.Z_dynamic_ohm,z,equal_nan=True):
            raise ValueError("INVALID_IMPEDANCE_CONTRACT")
        return dict(source_case_id=self.source_case_id,reference_plane_id="EXTERNAL_CEXT_TO_GAP_CGAP_PARALLEL_GSP",
            voltage_reference="GAP_NODE_MINUS_GROUND",positive_current_direction="EXTERNAL_TO_GAP",
            dt_port_s=self.dt_port_s,sample_count=len(f),Nyquist_Hz=.5/self.dt_port_s,
            sampling_status="BOUNDED_CIRCUIT_DENSE_OUTPUT",
            G0_status="HANDOFF_CALIBRATION_PENDING",
            G1_status="THERMAL_REFERENCE_SOLVER_VALIDATED",
            G2_status="BIDIRECTIONAL_REFERENCE_COUPLING_VALIDATED",
            G2_production_status="PRODUCTION_THERMAL_RLC_COUPLING_NOT_CALIBRATED",
            production_calibration_status="PRODUCTION_PORT_TRANSIENT_NOT_CALIBRATED",
            interface_constraints=["H_GAP_CAPACITANCE_PARTITION_REQUIRED",
                "FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED","LTE_APPLICABILITY_PENDING_CALIBRATION",
                "THERMAL_RSP_GRID_SENSITIVITY_PRESENT","RADIATION_MODEL_NOT_ENABLED"],
            provenance=self.metadata)
