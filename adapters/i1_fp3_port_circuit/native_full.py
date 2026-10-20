"""Serial native FULL protocol adapter; local numerical use must be explicit."""
from dataclasses import dataclass
import json
import os
import select
import subprocess
import time

from .ports import (AcceptedState, TrialInput, TrialResult, PortDescriptor, PortError,
                    CapacitanceOwner as Owner, ChargeScope, EnergyScope, finite)


@dataclass(frozen=True)
class NativeAcceptedState(AcceptedState):
    x_V: float | None = None
    cumulative_dissipation_J: float | None = None
    physical_hash: str = ""
    history_hash: str = ""


@dataclass(frozen=True)
class NativeTrialResult(TrialResult):
    # A cold electrical-work ledger is not a passive lumped resistor or gas heat.
    resistor_dissipation_J: float | None
    transfer_current_mid_A: float | None
    d_deltaq_port_d_voltage_n1_F: float | None
    native_trial_id: int = 0
    native_metrics: tuple = ()


class NativeFullProvider:
    stable_evaluation = False

    def __init__(self, binary, output_dir, *, owner=Owner.TOTAL,
                 operation_scope="PHYSICAL_TERMINAL", wall_budget_s=1800,
                 stable_evaluation=False):
        if operation_scope != "LOCAL_NUMERICAL_TEST":
            raise PortError("UNAVAILABLE_CAPABILITY_GLOBAL_PHYSICAL_PORT")
        self.owner = owner
        self.stable_evaluation = stable_evaluation
        self._begun = time.monotonic()
        self._deadline = self._begun+wall_budget_s
        self._stderr = output_dir.with_suffix(".stderr.log").open("x")
        self._protocol_log = output_dir.with_suffix(".stdout.log").open("x")
        self._process = subprocess.Popen(["timeout", str(wall_budget_s)+"s", str(binary), str(output_dir)],
                                          stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                          stderr=self._stderr, text=True, bufsize=1,
                                          env=dict(os.environ, FP3_STABLE_PORT="1" if stable_evaluation else "0"))
        self.rpc_records = []
        self.native_step_calls = 0
        self._pending = {}
        info = self._receive()
        self.rpc_records.append(dict(command="INITIALIZE", reply=info))
        self.info = info
        self._accepted = self._make_state(info)
        self.Cp = info["metrics"]["C_local_F"]
        self.descriptor = self._descriptor(output_dir)

    def _descriptor(self, output_dir):
        import hashlib
        identity = (output_dir/"operator_identity.txt").read_text().strip()
        operator = "ACTIVE_FV_NETWORK_V1:"+hashlib.sha256(identity.encode()).hexdigest()
        return PortDescriptor(
            "NATIVE_FULL_COLD", "B0_R02_A30_g70_Nominal", operator, "LOCAL_HV_RETURN_ACTIVE_FV",
            self.owner, ChargeScope.TOTAL if self.owner is Owner.TOTAL else ChargeScope.REMAINING,
            EnergyScope.TOTAL if self.owner is Owner.TOTAL else EnergyScope.REMAINING,
            ("STAGE_C_NATIVE", 200, 400, "q_B_active_operator", "ZERO_NEUMANN", self.Cp),
            capabilities=("provide_local_numerical_port_charge", "serial_native_transactions", "local_electrical_energy_ledger"),
            validity_domain="LOCAL_DISCRETE_DOMAIN; 499--501V; SP3_OFF; NATIVE_COLD; NO_GLOBAL_PORT_OR_HEAT",
            source_identity="NATIVE_FULL_COLD_NOT_COMPLETE_MICROPHYSICS")

    def _receive(self):
        remaining = self._deadline-time.monotonic()
        if remaining <= 0 or not select.select([self._process.stdout], [], [], remaining)[0]:
            self._process.kill()
            raise PortError("NATIVE_WALL_BUDGET_EXHAUSTED")
        line = self._process.stdout.readline()
        if not line:
            raise PortError("NATIVE_PROCESS_EOF_CHECK_STDERR")
        self._protocol_log.write(line);self._protocol_log.flush()
        # This machine's MPI/UCX initialization can print warnings on stdout.
        # Preserve them verbatim; only the explicitly tagged reply is protocol data.
        while not line.startswith('{"ok":'):
            if time.monotonic() >= self._deadline:
                raise PortError("NATIVE_WALL_BUDGET_EXHAUSTED")
            line = self._process.stdout.readline()
            if not line:
                raise PortError("NATIVE_PROCESS_EOF_CHECK_STDERR")
            self._protocol_log.write(line);self._protocol_log.flush()
        reply = json.loads(line)
        self.native_step_calls = reply.get("step_calls", self.native_step_calls)
        if self.native_step_calls > 120:
            raise PortError("NATIVE_STEP_BUDGET_EXHAUSTED")
        return reply

    def rpc(self, command):
        self._process.stdin.write(command+"\n")
        self._process.stdin.flush()
        reply = self._receive()
        self.rpc_records.append(dict(command=command, reply=reply))
        if not reply["ok"]:
            raise PortError(reply["error"])
        return reply

    @staticmethod
    def _make_state(info, prior=None, result=None):
        return NativeAcceptedState(token=info["token"], version=info["version"],
            time_s=info["time_s"], voltage_V=info["voltage_V"],
            physical_hash=info["physical_hash"], history_hash=info["history_hash"],
            cumulative_charge_C=prior.cumulative_charge_C+result.delta_q_port_C if prior else 0.,
            cumulative_port_work_J=prior.cumulative_port_work_J+result.port_work_J if prior else 0.,
            cumulative_component_work_J=prior.cumulative_component_work_J+result.component_work_J if prior else 0.)

    @property
    def accepted(self):
        return self._accepted

    def request(self, dt_s, voltage_n1_V):
        s, d = self.accepted, self.descriptor
        return TrialInput(s.token, s.time_s, dt_s, s.voltage_V, voltage_n1_V,
                          d.geometry_id, d.operator_id, d.reference_plane_id,
                          "LINEAR_ENDPOINTS_NATIVE_EULER_SUBCYCLED")

    def evaluate_trial(self, request, *, inject_failure=False):
        s, d = self.accepted, self.descriptor
        finite(request.t_n_s, request.dt_s, request.voltage_n_V, request.voltage_n1_V)
        if request.accepted_token != s.token:
            raise PortError("STALE_OR_FOREIGN_TOKEN")
        if (request.geometry_id, request.operator_id, request.reference_plane_id) != (d.geometry_id, d.operator_id, d.reference_plane_id):
            raise PortError("INSTANCE_OR_REFERENCE_MISMATCH")
        if (request.t_n_s != s.time_s or request.voltage_n_V != s.voltage_V
                or not 0 < request.dt_s <= 20e-15 or not 499 <= request.voltage_n1_V <= 501
                or request.voltage_representation != "LINEAR_ENDPOINTS_NATIVE_EULER_SUBCYCLED"):
            raise PortError("OUT_OF_LOCAL_NUMERICAL_SCOPE")
        info = self.rpc(f"TRIAL {s.token} {s.time_s:.17g} {request.dt_s:.17g} {s.voltage_V:.17g} {request.voltage_n1_V:.17g} {int(inject_failure)}")
        m = info["metrics"]
        if info["accepted_hash"] != s.physical_hash:
            raise PortError("ACCEPTED_NATIVE_STATE_CONTAMINATED")
        total = self.owner is Owner.TOTAL
        cpq = self.Cp*(request.voltage_n1_V-s.voltage_V)
        cp0, cp1 = .5*self.Cp*s.voltage_V**2, .5*self.Cp*request.voltage_n1_V**2
        dq = m["stable_Deltaq_port_C"] if self.stable_evaluation else m["Deltaq_port_C"]
        work = m["stable_Wport_J"] if self.stable_evaluation else m["Wport_native_J"]
        cap_work = .5*(s.voltage_V+request.voltage_n1_V)*cpq if self.stable_evaluation else cp1-cp0
        un0 = m["stable_U_remaining_start_J"] if self.stable_evaluation else m["U_E_start_J"]-cp0
        un1 = m["stable_U_remaining_end_J"] if self.stable_evaluation else m["U_E_end_J"]-cp1
        qselected = dq-m["Deltaq_Hg_C"] if self.stable_evaluation else m["q_B_active_end_C"]-m["q_B_active_start_C"]
        if not total and self.stable_evaluation:
            qselected = m["stable_DeltaQ_rho_C"]
        elif not total:
            qselected -= cpq
        residual = m["stable_energy_residual_J"] if self.stable_evaluation else m["energy_residual_J"]
        candidate = NativeAcceptedState(
            token="CANDIDATE:"+str(info["trial_id"]), version=s.version+1,
            time_s=m["time_n1_s"], voltage_V=request.voltage_n1_V,
            cumulative_charge_C=s.cumulative_charge_C+dq,
            cumulative_port_work_J=s.cumulative_port_work_J+work,
            cumulative_component_work_J=s.cumulative_component_work_J+work-(0. if total else cap_work),
            physical_hash=info["physical_hash"], history_hash=info["history_hash"])
        result = NativeTrialResult(request, candidate,
            qselected,
            0. if total else cpq, m["Deltaq_Hg_C"], dq,
            (un0+(cp0 if total else 0.) if self.stable_evaluation else m["U_E_start_J"]-(0. if total else cp0)),
            (un1+(cp1 if total else 0.) if self.stable_evaluation else m["U_E_end_J"]-(0. if total else cp1)),
            0. if total else cp0, 0. if total else cp1,
            work, work-(0. if total else cap_work),
            None, None, None, residual, residual,
            quadrature_id="SUM_NATIVE_SUBSTEP_ENDPOINT_MEAN_V_TIMES_ACTUAL_CHARGE",
            provenance="LOCAL_NUMERICAL_TEST_NATIVE_FULL_GAS_HEAT_UNAVAILABLE",
            native_trial_id=info["trial_id"], native_metrics=tuple(sorted(m.items())))
        self._pending[result.native_trial_id] = result
        return result

    def validate_commit(self, trial, *, converged):
        if not converged or trial.input.accepted_token != self.accepted.token:
            raise PortError("STALE_FOREIGN_OR_UNCONVERGED_TRIAL")
        if self._pending.get(trial.native_trial_id) is not trial:
            raise PortError("ALTERED_OR_FOREIGN_TRIAL")

    def commit(self, trial, *, converged):
        self.validate_commit(trial, converged=converged)
        info = self.rpc(f"COMMIT {self.accepted.token} {trial.native_trial_id}")
        prior = self.accepted
        self._accepted = self._make_state(info, prior, trial)
        self._pending.clear()
        return self.accepted

    def reject(self, trial):
        self.validate_commit(trial, converged=True)
        info = self.rpc(f"REJECT {self.accepted.token} {trial.native_trial_id}")
        if (info["physical_hash"], info["history_hash"]) != (self.accepted.physical_hash, self.accepted.history_hash):
            raise PortError("ROLLBACK_HISTORY_CONTAMINATED")
        del self._pending[trial.native_trial_id]

    def reset_initial_branch(self, owner):
        info = self.rpc("RESET_INITIAL")
        self.owner = owner
        self._accepted = self._make_state(info)
        from dataclasses import replace
        self.descriptor = replace(self.descriptor, capacitance_owner=owner,
            charge_scope=ChargeScope.TOTAL if owner is Owner.TOTAL else ChargeScope.REMAINING,
            energy_scope=EnergyScope.TOTAL if owner is Owner.TOTAL else EnergyScope.REMAINING)
        self._pending.clear()

    def require(self, capability):
        if capability not in self.descriptor.capabilities:
            raise PortError("UNAVAILABLE_CAPABILITY:"+capability)

    def close(self):
        if self._process.poll() is None:
            self._process.stdin.write("QUIT\n");self._process.stdin.flush()
            self._process.wait(timeout=30)
        self._stderr.close()
        self._protocol_log.close()
        if self._process.returncode != 0:
            raise PortError("NATIVE_PROCESS_NONZERO_EXIT")
