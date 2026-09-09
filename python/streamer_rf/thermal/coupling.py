"""G2 lagged, transactional coupling to the unchanged canonical circuit RHS."""
from __future__ import annotations

from copy import copy
from dataclasses import replace

import numpy as np
from scipy.integrate import solve_ivp

from streamer_rf.circuit.model import SeriesRLCGapCircuit
from streamer_rf.circuit.signals import ConductanceProfile
from .radial import joule_heating, thermal_initialization_from_g0


def resistance_conductance(resistance_ohm):
    if not np.isfinite(resistance_ohm) or resistance_ohm <= 0:
        raise ValueError("INVALID_THERMAL_RESISTANCE")
    g = 1. / resistance_ohm
    if not np.isfinite(g) or g <= 0:
        raise ValueError("INVALID_THERMAL_CONDUCTANCE")
    return ConductanceProfile.constant(g, source_kind="G1_THERMAL_RESISTANCE")


def require_g0_initialization(summary):
    result = thermal_initialization_from_g0(summary)
    if not result["valid"]:
        raise ValueError(result["status"])


def circuit_interval(params, state, t0, t1):
    """Reuse canonical RHS; augment only its three power quadratures.

    DOP853 internal nodes are transient. Integral states avoid output-cadence
    trapezoidal loss error; the canonical diagnostic evaluator checks KCL/KVL.
    """
    state = np.asarray(state, float)
    if state.shape != (2,) or not np.all(np.isfinite(state)) or not np.isfinite(t0+t1) or t1 <= t0:
        raise ValueError("INVALID_CIRCUIT_INTERVAL")
    params.conductance.require_supports(t0, t1)
    params.source_voltage.require_supports(t0, t1)
    kernel = SeriesRLCGapCircuit(params)

    def rhs(t, y):
        i, v = y[:2]
        return np.r_[kernel.rhs(t, y[:2]), params.conductance(t)*v*v,
                     params.Rs_ohm*i*i, params.source_voltage(t)*i]

    sol = solve_ivp(rhs, (t0, t1), np.r_[state, 0., 0., 0.], method="DOP853",
                    rtol=2e-10, atol=[1e-13, 1e-11, 1e-22, 1e-22, 1e-22],
                    max_step=(t1-t0)/4)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise ValueError("CIRCUIT_INTEGRATION_FAILED")
    d = kernel.evaluate_solution(sol.t, sol.y[0], sol.y[1])
    return dict(state=sol.y[:2, -1], Qsp=float(sol.y[2, -1]), Qfixed=float(sol.y[3, -1]),
                Wsource=float(sol.y[4, -1]), stored=float(d.stored_energy_J[-1]),
                internal_steps=len(sol.t)-1, nfev=sol.nfev,
                peak_current=float(np.max(np.abs(d.currents.I_gap_cond_A))),
                kcl_max=float(np.max(np.abs(d.KCL_residual_A))),
                kvl_max=float(np.max(np.abs(d.KVL_residual_V))))


class ThermalCircuitCoupling:
    """One circuit interval per accepted G1 step; no mutation on rejected trials.

    FIXED_RSP_REPLAY is deliberately one-way and may have a cross-domain energy
    mismatch once the actual thermal resistance differs from the fixed circuit R.
    """

    def __init__(self, thermal, params, initial_state, coupling_dt_s, *, feedback=True):
        if not np.isfinite(coupling_dt_s) or coupling_dt_s <= 0:
            raise ValueError("INVALID_COUPLING_INTERVAL")
        state = np.asarray(initial_state, float)
        if state.shape != (2,) or not np.all(np.isfinite(state)) or thermal.time != 0 or thermal.steps != 0:
            raise ValueError("INVALID_INITIAL_COUPLING_STATE")
        self.thermal, self.params, self.state = copy(thermal), params, state.copy()
        self.thermal.current = lambda t: 0.
        self.dt_cap, self.feedback = coupling_dt_s, bool(feedback)
        self.fixed_R = self.thermal.diagnostics()["Rsp_ohm"]
        resistance_conductance(self.fixed_R)
        self.initial_circuit_energy = .5*params.L_H*state[0]**2 + .5*params.C_total_F*state[1]**2
        self.Qsp = self.Qfixed = self.Wsource = 0.
        self.internal_steps = self.nfev = self.discarded_intervals = 0
        self.peak_current = abs(state[1]/self.fixed_R)
        self.rows = [self._row(self.fixed_R, 0., 0., 0.)]

    def _row(self, used_R, dt, pc, pt):
        d = self.thermal.diagnostics()
        # G1's instantaneous drive is zero outside trial steps. Do not export
        # those placeholder powers as if they described the coupled interval.
        for key in ("PJ_volume_W", "PJ_circuit_W", "Ez_V_m", "joule_closure_relative"):
            d.pop(key)
        r = d["Rsp_ohm"] if self.feedback else self.fixed_R
        i, v = self.state
        circuit = SeriesRLCGapCircuit(replace(self.params, conductance=resistance_conductance(r)))
        canonical = circuit.evaluate_solution(np.array([self.thermal.time]), np.array([i]), np.array([v]))
        stored = float(canonical.stored_energy_J[0])
        balance = stored-self.initial_circuit_energy+self.Qsp+self.Qfixed-self.Wsource
        d.update(I_A=v/r, I_L_A=i, V_channel_V=v, Rsp_circuit_ohm=r, Gsp_S=1/r,
                 interval_Rsp_ohm=used_R, interval_dt_s=dt,
                 P_sp_circuit_W=pc, P_sp_thermal_W=pt,
                 Q_sp_circuit_J=self.Qsp, Q_sp_thermal_J=self.thermal.QJ,
                 Q_fixed_R_J=self.Qfixed, source_work_J=self.Wsource,
                 inductor_energy_J=float(canonical.inductor_energy_J[0]),
                 cext_energy_J=float(canonical.cext_energy_J[0]),
                 cgap_energy_J=float(canonical.cgap_energy_J[0]),
                 circuit_energy_J=stored, thermal_energy_J=self.thermal.total_energy(),
                 circuit_energy_residual_J=balance,
                 cross_domain_energy_residual_J=self.Qsp-self.thermal.QJ,
                 coupling_status="LAGGED" if self.feedback else "FIXED_RSP_REPLAY")
        return d

    def step(self, maximum_dt=None):
        if maximum_dt is not None and (not np.isfinite(maximum_dt) or maximum_dt <= 0):
            raise ValueError("INVALID_COUPLING_INTERVAL")
        limit = self.dt_cap if maximum_dt is None else min(self.dt_cap, maximum_dt)
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError("INVALID_COUPLING_INTERVAL")
        base = self.thermal
        used_R = base.diagnostics()["Rsp_ohm"] if self.feedback else self.fixed_R
        params = replace(self.params, conductance=resistance_conductance(used_R))
        # Probe the unchanged G1 hydro/conduction limit without a guessed current.
        probe = copy(base)
        probe.current = lambda t: 0.
        dt = min(limit, probe._rhs()[1])
        for _ in range(24):
            c = circuit_interval(params, self.state, base.time, base.time+dt)
            if c["Qsp"] < 0:
                raise ValueError("NEGATIVE_CONDUCTIVE_ENERGY")
            rms = np.sqrt(c["Qsp"]/(used_R*dt))
            trial = copy(base)
            trial.current = lambda t, value=rms: value
            props = base.primitive()[3]
            independent = joule_heating(props["sigma_S_m"], base.grid.volumes, base.length, rms)
            accepted_dt = trial.step(dt)
            if accepted_dt < dt:
                self.discarded_intervals += 1
                dt = accepted_dt
                continue
            break
        else:
            raise ValueError("COUPLING_ACCEPTANCE_FAILED")
        # Both subsystem states and all budgets are committed together.
        trial.current = lambda t: 0.
        self.thermal, self.state = trial, c["state"]
        self.Qsp += c["Qsp"]
        self.Qfixed += c["Qfixed"]
        self.Wsource += c["Wsource"]
        self.internal_steps += c["internal_steps"]
        self.nfev += c["nfev"]
        self.peak_current = max(self.peak_current, c["peak_current"], abs(self.state[1]/(
            trial.diagnostics()["Rsp_ohm"] if self.feedback else self.fixed_R)))
        row = self._row(used_R, dt, c["Qsp"]/dt, independent["PJ_volume_W"])
        row.update(KCL_max_A=c["kcl_max"], KVL_max_V=c["kvl_max"], I_rms_interval_A=rms)
        self.rows.append(row)
        return dt

    def run(self, end_time_s):
        if not np.isfinite(end_time_s) or end_time_s <= self.thermal.time:
            raise ValueError("INVALID_END_TIME")
        while self.thermal.time < end_time_s:
            self.step(end_time_s-self.thermal.time)
        return self.rows
