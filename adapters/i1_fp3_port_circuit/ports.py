"""I1-PORT-1.0 synthetic, differential two-terminal transaction protocol."""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from enum import Enum
import math
from typing import Protocol
from uuid import uuid4


class PortError(ValueError):
    pass


class CapacitanceOwner(str, Enum):
    TOTAL = "TOTAL_CHARGE"
    SPLIT = "SPLIT_GEOMETRIC_CAPACITANCE"


class ChargeScope(str, Enum):
    TOTAL = "COMPLETE_Q_H"
    REMAINING = "Q_H_MINUS_CP_V"


class EnergyScope(str, Enum):
    TOTAL = "COMPONENT_INCLUDES_CP"
    REMAINING = "COMPONENT_EXCLUDES_CP"


UNITS = (("time", "s"), ("voltage", "V"), ("current", "A"),
         ("charge", "C"), ("energy", "J"), ("capacitance", "F"))


def finite(*values):
    if not all(math.isfinite(v) for v in values):
        raise PortError("NONFINITE_INPUT")


def positive(**values):
    finite(*values.values())
    if any(v <= 0 for v in values.values()):
        raise PortError("NONPOSITIVE_PARAMETER")


@dataclass(frozen=True)
class PortDescriptor:
    provider_id: str
    geometry_id: str
    operator_id: str
    reference_plane_id: str
    capacitance_owner: CapacitanceOwner
    charge_scope: ChargeScope
    energy_scope: EnergyScope
    model_identity: tuple
    schema_version: str = "I1-FP3-ENGINEERING-1.0"
    units: tuple = UNITS
    capabilities: tuple = ("synthetic_terminal_charge", "synthetic_transfer_current",
                           "synthetic_component_energy", "synthetic_transactions")
    validity_domain: str = "SYNTHETIC_LINEAR_RC_ONLY_NO_PLASMA"
    source_identity: str = "SYNTHETIC_TEST_PARAMETERS_NOT_COMSOL"

    def __post_init__(self):
        if not all((self.provider_id, self.geometry_id, self.operator_id,
                    self.reference_plane_id)):
            raise PortError("MISSING_IDENTITY")
        if self.units != UNITS:
            raise PortError("INCOMPATIBLE_UNITS")
        if not isinstance(self.capacitance_owner, CapacitanceOwner):
            raise PortError("AMBIGUOUS_CAPACITANCE_OWNER")
        expected = (ChargeScope.TOTAL, EnergyScope.TOTAL) if self.capacitance_owner is CapacitanceOwner.TOTAL else (ChargeScope.REMAINING, EnergyScope.REMAINING)
        if (self.charge_scope, self.energy_scope) != expected:
            raise PortError("DUPLICATE_OR_INCONSISTENT_CAPACITANCE_OWNERSHIP")


@dataclass(frozen=True)
class AcceptedState:
    token: str
    version: int
    time_s: float
    voltage_V: float
    x_V: float
    cumulative_charge_C: float = 0.
    cumulative_port_work_J: float = 0.
    cumulative_component_work_J: float = 0.
    cumulative_dissipation_J: float = 0.


@dataclass(frozen=True)
class TrialInput:
    accepted_token: str
    t_n_s: float
    dt_s: float
    voltage_n_V: float
    voltage_n1_V: float
    geometry_id: str
    operator_id: str
    reference_plane_id: str
    voltage_representation: str = "LINEAR_ENDPOINTS_IMPLICIT_MIDPOINT"


@dataclass(frozen=True)
class TrialResult:
    input: TrialInput
    candidate: AcceptedState
    delta_Q_selected_C: float
    delta_q_external_Cp_C: float
    delta_q_transfer_C: float
    delta_q_port_C: float
    U_n_J: float
    U_n1_J: float
    external_Cp_U_n_J: float
    external_Cp_U_n1_J: float
    port_work_J: float
    component_work_J: float
    resistor_dissipation_J: float
    transfer_current_mid_A: float
    d_deltaq_port_d_voltage_n1_F: float
    component_energy_residual_J: float
    total_energy_residual_J: float
    quadrature_id: str = "IMPLICIT_MIDPOINT_DISCRETE_POWER_NOT_CONTINUOUS_EXACT_WORK"
    provenance: str = "SYNTHETIC_LINEAR_RC_NOT_PLASMA"


class PortProvider(Protocol):
    descriptor: PortDescriptor

    @property
    def accepted(self) -> AcceptedState: ...
    def evaluate_trial(self, request: TrialInput) -> TrialResult: ...
    def commit(self, result: TrialResult, *, converged: bool) -> AcceptedState: ...
    def reject(self, result: TrialResult) -> None: ...


class SyntheticPort:
    """Pure trials; only commit mutates the immutable accepted-state reference."""

    def __init__(self, *, kind, Cp_F=0., R_ohm=None, Cd_F=None,
                 owner=CapacitanceOwner.TOTAL, voltage0_V=0., x0_V=0.,
                 geometry_id="SYNTHETIC_TWO_TERMINAL", operator_id="LINEAR_RC_V1",
                 reference_plane_id="DIFFERENTIAL_H_TO_R"):
        finite(Cp_F, voltage0_V, x0_V)
        if Cp_F < 0:
            raise PortError("INVALID_CP")
        if kind == "R":
            positive(R=R_ohm)
            if Cp_F != 0 or x0_V != 0:
                raise PortError("INVALID_RESISTOR_STATE")
        elif kind == "C":
            positive(Cp=Cp_F)
            if x0_V != 0:
                raise PortError("INVALID_CAPACITOR_STATE")
        elif kind == "MEMORY_RC":
            positive(Cp=Cp_F, Cd=Cd_F, Rd=R_ohm)
        else:
            raise PortError("UNSUPPORTED_SYNTHETIC_PORT")
        self.kind, self.Cp, self.R, self.Cd = kind, Cp_F, R_ohm, Cd_F
        self._instance = uuid4().hex
        total = owner is CapacitanceOwner.TOTAL
        self.descriptor = PortDescriptor(
            kind, geometry_id, operator_id, reference_plane_id, owner,
            ChargeScope.TOTAL if total else ChargeScope.REMAINING,
            EnergyScope.TOTAL if total else EnergyScope.REMAINING,
            (kind, Cp_F, R_ohm, Cd_F))
        self._accepted = AcceptedState(self._token(0), 0, 0., voltage0_V, x0_V)

    def _token(self, version):
        return f"{self._instance}:{version}"

    @property
    def accepted(self):
        return self._accepted

    def request(self, dt_s, voltage_n1_V):
        s, d = self.accepted, self.descriptor
        return TrialInput(s.token, s.time_s, dt_s, s.voltage_V, voltage_n1_V,
                          d.geometry_id, d.operator_id, d.reference_plane_id)

    def evaluate_trial(self, request):
        s, d = self.accepted, self.descriptor
        if request.accepted_token != s.token:
            raise PortError("STALE_OR_FOREIGN_TOKEN")
        if (request.geometry_id, request.operator_id, request.reference_plane_id) != (d.geometry_id, d.operator_id, d.reference_plane_id):
            raise PortError("INSTANCE_OR_REFERENCE_MISMATCH")
        finite(request.t_n_s, request.dt_s, request.voltage_n_V, request.voltage_n1_V)
        if request.dt_s <= 0 or request.t_n_s != s.time_s or request.voltage_n_V != s.voltage_V:
            raise PortError("INVALID_TIME_LAYER_OR_DT")
        if request.voltage_representation != "LINEAR_ENDPOINTS_IMPLICIT_MIDPOINT":
            raise PortError("UNSUPPORTED_QUADRATURE")
        dt, v1 = request.dt_s, request.voltage_n1_V
        vm = .5 * (s.voltage_V + v1)
        x1, j, dj_dv1 = s.x_V, 0., 0.
        if self.kind == "R":
            j, dj_dv1 = vm / self.R, .5 / self.R
        elif self.kind == "MEMORY_RC":
            # Cd*(x1-x0)/dt=(Vbar-(x0+x1)/2)/Rd.
            a = dt / (2 * self.R * self.Cd)
            x1 = ((1-a)*s.x_V + 2*a*vm)/(1+a)
            j = (vm - .5*(s.x_V+x1))/self.R
            dj_dv1 = .5/(self.R*(1+a))
        cp_q = self.Cp * (v1 - s.voltage_V)
        total = d.capacitance_owner is CapacitanceOwner.TOTAL
        selected_q, external_q = (cp_q, 0.) if total else (0., cp_q)
        transfer_q = dt*j
        port_q = selected_q + external_q + transfer_q
        cp0, cp1 = .5*self.Cp*s.voltage_V**2, .5*self.Cp*v1**2
        cd0 = .5*self.Cd*s.x_V**2 if self.Cd is not None else 0.
        cd1 = .5*self.Cd*x1**2 if self.Cd is not None else 0.
        u0, u1 = (cp0+cd0, cp1+cd1) if total else (cd0, cd1)
        ec0, ec1 = (0., 0.) if total else (cp0, cp1)
        work, component_work = vm*port_q, vm*(selected_q+transfer_q)
        loss = dt*self.R*j*j if self.R is not None else 0.
        candidate = AcceptedState(self._token(s.version+1), s.version+1,
                                 s.time_s+dt, v1, x1,
                                 s.cumulative_charge_C+port_q,
                                 s.cumulative_port_work_J+work,
                                 s.cumulative_component_work_J+component_work,
                                 s.cumulative_dissipation_J+loss)
        result = TrialResult(request, candidate, selected_q, external_q, transfer_q,
                             port_q, u0, u1, ec0, ec1, work, component_work, loss,
                             j, self.Cp+dt*dj_dv1,
                             component_work-(u1-u0)-loss,
                             work-((u1+ec1)-(u0+ec0))-loss)
        finite(*[v for v in asdict(candidate).values() if isinstance(v, (float, int))],
               work, loss, port_q, result.total_energy_residual_J)
        return result

    def validate_commit(self, result, *, converged):
        if not converged:
            raise PortError("UNCONVERGED_TRIAL")
        if result != self.evaluate_trial(result.input):
            raise PortError("ALTERED_OR_FOREIGN_TRIAL")

    def commit(self, result, *, converged):
        self.validate_commit(result, converged=converged)
        self._accepted = result.candidate
        return self._accepted

    def reject(self, result):
        self.validate_commit(result, converged=True)
        # Pure evaluate_trial has no accepted histories or budgets to undo.


class IdealResistorPort(SyntheticPort):
    def __init__(self, R_test, **kwargs):
        super().__init__(kind="R", R_ohm=R_test, **kwargs)


class IdealCapacitorPort(SyntheticPort):
    def __init__(self, Cp, **kwargs):
        super().__init__(kind="C", Cp_F=Cp, **kwargs)


class SyntheticMemoryRCPort(SyntheticPort):
    def __init__(self, Cp, Rd, Cd, **kwargs):
        super().__init__(kind="MEMORY_RC", Cp_F=Cp, R_ohm=Rd, Cd_F=Cd, **kwargs)


def differential_voltage(V_H, V_R):
    finite(V_H, V_R)
    return V_H-V_R


def require_capacitance_assembly(descriptor, *, circuit_assembles_Cp):
    expected = descriptor.capacitance_owner is CapacitanceOwner.SPLIT
    if circuit_assembles_Cp is not expected:
        raise PortError("DUPLICATE_OR_MISSING_CP_ASSEMBLY")
