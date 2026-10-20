"""Only the authorized synthetic series L-Rs-Cext-port circuit, midpoint ODE."""
from dataclasses import dataclass, replace
import math
import time
from uuid import uuid4

import numpy as np

from .ports import (CapacitanceOwner, PortError, TrialResult, finite, positive,
                    require_capacitance_assembly)


@dataclass(frozen=True)
class CircuitParameters:
    L_H: float = 100e-9
    Cext_F: float = 101e-12
    Rs_ohm: float = 3.

    def __post_init__(self):
        positive(L=self.L_H, Cext=self.Cext_F)
        finite(self.Rs_ohm)
        if self.Rs_ohm < 0:
            raise PortError("NEGATIVE_SERIES_RESISTANCE")


@dataclass(frozen=True)
class SyntheticSource:
    voltage_V: float
    source_kind: str = "SYNTHETIC_IDEAL_PRESCRIBED_CONSTANT"

    def __post_init__(self):
        if self.source_kind != "SYNTHETIC_IDEAL_PRESCRIBED_CONSTANT":
            raise PortError("SOURCE_MODEL_UNAVAILABLE")
        finite(self.voltage_V)


@dataclass(frozen=True)
class Scales:
    V: float
    T: float
    I: float
    Q: float
    W: float

    @classmethod
    def from_system(cls, params, port, i0, vc0):
        # Fixed before stepping. Never normalize by a cancelling instantaneous value.
        v = 1.
        t = 2*math.pi*math.sqrt(params.L_H*params.Cext_F)
        i = v/math.sqrt(params.L_H/params.Cext_F)
        q = max(params.Cext_F, port.Cp, port.Cd or 0.)*v
        initial = (.5*params.L_H*i0*i0 + .5*params.Cext_F*vc0*vc0
                   + .5*port.Cp*port.accepted.voltage_V**2
                   + .5*(port.Cd or 0.)*port.accepted.x_V**2)
        w = max(initial, max(params.Cext_F, port.Cp, port.Cd or 0.)*v*v)
        return cls(v, t, i, q, w)


@dataclass(frozen=True)
class CircuitAccepted:
    token: str
    version: int
    port_token: str
    time_s: float
    i_A: float
    Vc_V: float
    source_work_J: float = 0.
    series_loss_J: float = 0.


@dataclass(frozen=True)
class CircuitTrial:
    base_token: str
    dt_s: float
    i1_A: float
    Vc1_V: float
    port_trial: TrialResult
    iterations: int
    scaled_residual: tuple


class SeriesPortCircuit:
    def __init__(self, params, port, source, *, i0_A, Vc0_V,
                 reference_plane_id="DIFFERENTIAL_H_TO_R",
                 circuit_assembles_Cp=None):
        finite(i0_A, Vc0_V)
        if port.accepted.time_s != 0 or port.accepted.version != 0:
            raise PortError("FRESH_PORT_REQUIRED")
        if reference_plane_id != port.descriptor.reference_plane_id:
            raise PortError("INCOMPATIBLE_REFERENCE_PLANE")
        if port.kind == "R" and port.accepted.voltage_V != port.R*i0_A:
            raise PortError("INCONSISTENT_RESISTOR_INITIAL_VOLTAGE")
        expected = port.descriptor.capacitance_owner is CapacitanceOwner.SPLIT
        if circuit_assembles_Cp is None:
            circuit_assembles_Cp = expected
        require_capacitance_assembly(port.descriptor,
                                    circuit_assembles_Cp=circuit_assembles_Cp)
        self.params, self.port, self.source = params, port, source
        self.scales = Scales.from_system(params, port, i0_A, Vc0_V)
        self._instance = uuid4().hex
        self._accepted = CircuitAccepted(self._token(0), 0, port.accepted.token,
                                        0., i0_A, Vc0_V)
        self.initial_energy = self.total_energy()
        self.records = []
        self.trial_evaluations = 0

    def _token(self, version):
        return f"{self._instance}:{version}"

    @property
    def accepted(self):
        return self._accepted

    def total_energy(self, i=None, vc=None, vp=None, x=None):
        s, p, c = self.accepted, self.port.accepted, self.params
        i = s.i_A if i is None else i
        vc = s.Vc_V if vc is None else vc
        vp = p.voltage_V if vp is None else vp
        x = p.x_V if x is None else x
        return (.5*c.L_H*i*i+.5*c.Cext_F*vc*vc
                +.5*self.port.Cp*vp*vp+.5*(self.port.Cd or 0.)*x*x)

    def _residual(self, dt, i1, vp1):
        s, c = self.accepted, self.params
        if s.port_token != self.port.accepted.token:
            raise PortError("COMPONENT_CIRCUIT_DESYNCHRONIZED")
        self.trial_evaluations += 1
        trial = self.port.evaluate_trial(self.port.request(dt, vp1))
        selected_q = trial.delta_Q_selected_C+trial.delta_q_transfer_C
        if self.port.descriptor.capacitance_owner is CapacitanceOwner.SPLIT:
            # The circuit, not a second component charge, owns this Cp contribution.
            assembled_q = selected_q+self.port.Cp*(vp1-self.port.accepted.voltage_V)
        else:
            assembled_q = selected_q
        im = .5*(s.i_A+i1)
        vc1 = s.Vc_V+dt*im/c.Cext_F
        vm = .5*(self.port.accepted.voltage_V+vp1)
        r = np.array([(c.L_H*(i1-s.i_A)/dt+c.Rs_ohm*im
                       +.5*(s.Vc_V+vc1)+vm-self.source.voltage_V)/self.scales.V,
                      (im-assembled_q/dt)/self.scales.I])
        jac = np.array([[(c.L_H/dt+.5*c.Rs_ohm+dt/(4*c.Cext_F))/self.scales.V,
                         .5/self.scales.V],
                        [.5/self.scales.I,
                         -trial.d_deltaq_port_d_voltage_n1_F/dt/self.scales.I]])
        return r, jac, trial, vc1

    def propose(self, dt_s, *, max_iterations=8):
        positive(dt=dt_s)
        if not 1 <= max_iterations <= 8:
            raise PortError("ITERATION_BUDGET_INVALID")
        y = np.array([self.accepted.i_A, self.port.accepted.voltage_V])
        for n in range(1, max_iterations+1):
            r, jac, trial, vc1 = self._residual(dt_s, *y)
            if np.max(np.abs(r)) <= 1e-10:
                return CircuitTrial(self.accepted.token, dt_s, float(y[0]), vc1,
                                    trial, n, tuple(float(v) for v in r))
            y -= np.linalg.solve(jac, r)
        raise PortError("SYNTHETIC_LINEAR_SOLVE_DID_NOT_CONVERGE")

    def commit(self, trial):
        if trial.base_token != self.accepted.token:
            raise PortError("STALE_OR_FOREIGN_CIRCUIT_TRIAL")
        finite(trial.i1_A, trial.Vc1_V)
        r, _, verified, vc1 = self._residual(trial.dt_s, trial.i1_A,
                                            trial.port_trial.candidate.voltage_V)
        if (np.max(np.abs(r)) > 1e-10 or verified != trial.port_trial
                or trial.Vc1_V != vc1):
            raise PortError("UNCONVERGED_OR_ALTERED_COUPLED_TRIAL")
        self.port.validate_commit(verified, converged=True)
        s, p, c = self.accepted, self.port.accepted, self.params
        im = .5*(s.i_A+trial.i1_A)
        source_work = trial.dt_s*self.source.voltage_V*im
        series_loss = trial.dt_s*c.Rs_ohm*im*im
        e0 = self.total_energy()
        e1 = self.total_energy(trial.i1_A, trial.Vc1_V,
                               verified.candidate.voltage_V, verified.candidate.x_V)
        new_state = CircuitAccepted(self._token(s.version+1), s.version+1,
                                    verified.candidate.token, verified.candidate.time_s,
                                    trial.i1_A, trial.Vc1_V,
                                    s.source_work_J+source_work,
                                    s.series_loss_J+series_loss)
        ledger = dict(time_s=new_state.time_s, dt_s=trial.dt_s,
                      i_A=trial.i1_A, Vc_V=trial.Vc1_V,
                      Vp_V=verified.candidate.voltage_V, x_V=verified.candidate.x_V,
                      total_energy_J=e1, delta_energy_J=e1-e0,
                      source_work_J=source_work, series_loss_J=series_loss,
                      port_loss_J=verified.resistor_dissipation_J,
                      delta_q_port_C=verified.delta_q_port_C,
                      charge_residual_C=trial.dt_s*im-verified.delta_q_port_C,
                      port_work_J=verified.port_work_J,
                      component_work_J=verified.component_work_J,
                      component_U_n_J=verified.U_n_J, component_U_n1_J=verified.U_n1_J,
                      circuit_Cp_U_n_J=verified.external_Cp_U_n_J,
                      circuit_Cp_U_n1_J=verified.external_Cp_U_n1_J,
                      port_energy_residual_J=verified.total_energy_residual_J,
                      component_energy_residual_J=verified.component_energy_residual_J,
                      energy_residual_J=e1-e0-source_work+series_loss+verified.resistor_dissipation_J,
                      kvl_residual_V=float(r[0]*self.scales.V),
                      current_residual_A=float(r[1]*self.scales.I),
                      source_work_cumulative_J=new_state.source_work_J,
                      loss_cumulative_J=new_state.series_loss_J+verified.candidate.cumulative_dissipation_J,
                      port_charge_cumulative_C=verified.candidate.cumulative_charge_C,
                      port_work_cumulative_J=verified.candidate.cumulative_port_work_J,
                      iterations=trial.iterations,
                      accepted_version=new_state.version,
                      quadrature_id=verified.quadrature_id)
        # All validation and ledger construction precede either accepted-state assignment.
        self.port.commit(verified, converged=True)
        self._accepted = new_state
        self.records.append(ledger)
        return ledger

    def advance(self, dt_s, *, fail_after_trial=False):
        trial = self.propose(dt_s)
        if fail_after_trial:
            self.port.reject(trial.port_trial)
            raise PortError("CONTROLLED_REJECTION_AFTER_INTERNAL_X")
        return self.commit(trial)

    def run(self, end_s, dt_s, *, deadline=None):
        positive(end=end_s, dt=dt_s)
        start = self.accepted.time_s
        if end_s <= start:
            raise PortError("INVALID_END_TIME")
        # Endpoint-based scheduling avoids a roundoff-only extra final step.
        span = end_s-start
        ratio = span/dt_s
        nearest = round(ratio)
        count = nearest if nearest > 0 and abs(span-nearest*dt_s) <= 4*math.ulp(span) else math.ceil(ratio)
        for k in range(1, count+1):
            if deadline is not None and time.monotonic() >= deadline:
                raise TimeoutError("FP3_TOTAL_EXECUTION_RESOURCE_LIMIT")
            target = end_s if k == count else min(start+k*dt_s, end_s)
            if target > self.accepted.time_s:
                self.advance(target-self.accepted.time_s)
        return self.records
