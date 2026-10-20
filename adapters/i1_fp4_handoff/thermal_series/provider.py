"""Lagged LTE column / floating series circuit; no changes to native G1.

Trial workspaces deep-copy the Python G1 state, not a reconstructed initial
condition. All physical histories survive; rejected work still costs resources.
"""
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import pickle
import time
from uuid import uuid4

import numpy as np
from scipy.integrate import solve_ivp

from adapters.i1_fp3_port_circuit.ports import (
    CapacitanceOwner, ChargeScope, EnergyScope, PortDescriptor,
    require_capacitance_assembly,
)
from adapters.i1_fp3_port_circuit.circuit import CircuitParameters, SyntheticSource


def zero_current(t):
    return 0.


class Budget:
    def __init__(self):
        self.start = time.monotonic()
        self.steps = self.trials = self.discarded = self.native_rejections = 0

    def check(self):
        if time.monotonic()-self.start >= 600 or self.steps >= 2000:
            raise RuntimeError('RESOURCE_LIMIT')


@dataclass
class Trial:
    token: str
    thermal: object
    state: np.ndarray
    row: dict
    sums: dict
    signature: str = ''

def digest(*objects):
    return hashlib.sha256(pickle.dumps(objects, protocol=5)).hexdigest()


def interval(params, source, state, resistance, cp, dt):
    """New series topology, with independent charge and power quadratures.

    Quadratures: source work, Rs loss, column loss, port/column charge,
    integral Vc, integral Vp. Integral KVL and capacitor balances are checked
    against endpoint states rather than merely evaluating RHS identities.
    """
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(resistance) or resistance <= 0:
        raise ValueError('INVALID_INTERVAL_OR_COLUMN')
    L, C, Rs = params.L_H, params.Cext_F, params.Rs_ohm
    n = 3 if cp else 2
    initial = state[:n]

    def rhs(t, y):
        i, vc = y[:2]
        vp = y[2] if cp else resistance*i
        col = vp/resistance
        core = [(source.voltage_V-Rs*i-vc-vp)/L, i/C]
        if cp:
            core.append((i-col)/cp)
        return np.r_[core, source.voltage_V*i, Rs*i*i,
                     resistance*col*col, i, col, vc, vp]

    atol = [1e-13, 1e-11]+([1e-11] if cp else [])+[1e-23]*3+[1e-23]*2+[1e-21]*2
    sol = solve_ivp(rhs, (0., dt), np.r_[initial, np.zeros(7)],
                    method='DOP853', rtol=2e-10, atol=atol, max_step=dt/4)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError('CIRCUIT_TRIAL_FAILED')
    end = sol.y[:n, -1]
    end = np.r_[end, resistance*end[0]] if not cp else end
    ws, qr, qj, qi, qc, ivc, ivp = sol.y[n:, -1]
    vp0 = state[2] if cp else resistance*state[0]
    def energy(s):
        return .5*L*s[0]**2+.5*C*s[1]**2+.5*cp*s[2]**2
    denergy = energy(end)-energy(state)
    voltages = sol.y[2] if cp else resistance*sol.y[0]
    return dict(state=end, W_source_J=ws, Q_Rs_J=qr, QJ_circuit_J=qj,
                q_port_C=qi, q_column_C=qc, delta_q_cp_C=cp*(end[2]-vp0),
                E_electrical_J=energy(end), delta_E_electrical_J=denergy,
                R_circuit_J=denergy+qr+qj-ws,
                KCL_integrated_C=qi-qc-cp*(end[2]-vp0),
                Cext_charge_residual_C=C*(end[1]-state[1])-qi,
                KVL_integrated_Vs=L*(end[0]-state[0])+Rs*qi+ivc+ivp-source.voltage_V*dt,
                peak_i_A=float(np.max(np.abs(sol.y[0]))),
                peak_Icol_A=float(np.max(np.abs(voltages/resistance))),
                nfev=sol.nfev)


class ThermalSeriesProvider:
    def __init__(self, solver, cp, owner, budget):
        if not np.isfinite(cp) or cp < 0:
            raise ValueError('INVALID_SYNTHETIC_CP')
        self.params = CircuitParameters(L_H=100e-9, Cext_F=101e-12, Rs_ohm=3.)
        self.source = SyntheticSource(10.)
        self.cp, self.budget = cp, budget
        self.descriptor = PortDescriptor(
            provider_id='G1_GIVEN_LTE_SERIES', geometry_id='GIVEN_RADIAL_1MM',
            operator_id='NATIVE_G1_UNMODIFIED_LAGGED_RMS', reference_plane_id='COLUMN_ENDS_NO_ELECTRODE_DROP',
            capacitance_owner=owner,
            charge_scope=ChargeScope.TOTAL if owner is CapacitanceOwner.TOTAL else ChargeScope.REMAINING,
            energy_scope=EnergyScope.TOTAL if owner is CapacitanceOwner.TOTAL else EnergyScope.REMAINING,
            model_identity=(('Cp_F', cp),), capabilities=('given_lte_column', 'thermal_transactions'),
            validity_domain='NUMERICAL_ENGINEERING_REFERENCE_NOT_B0',
            source_identity='SYNTHETIC_IDEAL_10V_STEP')
        require_capacitance_assembly(self.descriptor, circuit_assembles_Cp=owner is CapacitanceOwner.SPLIT)
        self.thermal = deepcopy(solver)
        self.thermal.current = zero_current
        self.state = np.zeros(3)
        self.token = uuid4().hex
        self.sums = {k: 0. for k in ('W_source_J', 'Q_Rs_J', 'QJ_circuit_J', 'q_port_C', 'q_column_C')}
        self.pending = {}
        self.initial_energy = self.thermal.total_energy()

    def semantic_hash(self):
        return digest(self.thermal.__dict__, self.state, self.sums, self.token)

    def evaluate_trial(self, token, dt, *, fail_after_thermal=False):
        if token != self.token:
            raise ValueError('STALE_OR_FOREIGN_TOKEN')
        self.budget.check()
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('INVALID_DT')
        self.budget.trials += 1
        base = self.thermal
        used = base.diagnostics()['Rsp_ohm']
        dt = min(dt, base._rhs()[1])  # zero drive outside trials
        for _ in range(16):
            self.budget.check()
            c = interval(self.params, self.source, self.state, used, self.cp, dt)
            rms = np.sqrt(c['QJ_circuit_J']/(used*dt))
            trial = deepcopy(base)
            trial.current = lambda t, value=rms: value
            actual = trial.step(dt)
            self.budget.steps += 1
            self.budget.native_rejections += trial.rejected_steps-base.rejected_steps
            if fail_after_thermal:
                self.budget.discarded += 1
                raise RuntimeError('CONTROLLED_REJECTION_AFTER_NATIVE_STEP')
            if actual < dt:
                self.budget.discarded += 1
                dt = actual
                continue
            break
        else:
            raise RuntimeError('SHORTENING_RETRY_LIMIT')
        trial.current = zero_current
        d = trial.diagnostics()
        e0, e1 = base.total_energy(), trial.total_energy()
        qjt = trial.QJ-base.QJ
        cond, hydro = trial.Qconduction-base.Qconduction, trial.Qhydro-base.Qhydro
        rt = e1-e0-qjt+cond+hydro
        rx = c['QJ_circuit_J']-qjt
        total = c['delta_E_electrical_J']+(e1-e0)+c['Q_Rs_J']+cond+hydro-c['W_source_J']
        exchange = max(abs(c['W_source_J']), abs(c['Q_Rs_J'])+abs(c['QJ_circuit_J'])+abs(cond)+abs(hydro))
        # Partition is applied before assembly; it never changes the physical Cp.
        split = self.descriptor.capacitance_owner is CapacitanceOwner.SPLIT
        selected_q = 0. if split else c['delta_q_cp_C']
        external_q = c['delta_q_cp_C'] if split else 0.
        cp_energy = .5*self.cp*c['state'][2]**2
        row = {k: float(v) for k, v in c.items() if k != 'state'}
        row.update(time_s=trial.time, dt_s=dt, i_port_A=c['state'][0], I_L_A=c['state'][0],
                   Vc_V=c['state'][1], Vp_interval_end_V=c['state'][2],
                   Icol_interval_end_A=c['state'][2]/used, Icol_RMS_A=rms,
                   I_Cp_interval_end_A=(c['state'][0]-c['state'][2]/used) if self.cp else 0.,
                   Icol_updated_state_A=(c['state'][2]/d['Rsp_ohm']) if self.cp else c['state'][0],
                   Vp_updated_state_V=c['state'][2] if self.cp else c['state'][0]*d['Rsp_ohm'],
                   R_used_ohm=used, R_end_ohm=d['Rsp_ohm'], T_axis_K=trial.T[0],
                   p_axis_Pa=trial.p[0], u_max_m_s=float(np.max(np.abs(trial.U[:, 1]/trial.U[:, 0]))),
                   mass_kg=float(np.dot(trial.U[:, 0], trial.grid.volumes)*trial.length),
                   mass_residual_kg=d['mass_residual_kg'], E_G1_J=e1, delta_E_G1_J=e1-e0,
                   QJ_thermal_J=qjt, Q_conduction_out_J=cond, Q_hydro_out_J=hydro,
                   R_thermal_J=rt, R_cross_J=rx, R_total_J=total, exchange_scale_J=exchange,
                   identity_residual_J=total-(c['R_circuit_J']+rt-rx),
                   delta_Q_selected_C=selected_q, delta_q_external_Cp_C=external_q,
                   assembled_charge_residual_C=selected_q+external_q+c['q_column_C']-c['q_port_C'],
                   component_energy_J=e1+(0. if split else cp_energy),
                   external_Cp_energy_J=cp_energy if split else 0.,
                   total_stored_energy_J=e1+c['E_electrical_J'],
                   native_steps=trial.steps, native_rejected_steps=trial.rejected_steps)
        sums = {k: self.sums[k]+row[k] for k in self.sums}
        row.update({'cumulative_'+k: v for k, v in sums.items()})
        row.update(cumulative_QJ_thermal_J=trial.QJ, cumulative_Q_conduction_out_J=trial.Qconduction,
                   cumulative_Q_hydro_out_J=trial.Qhydro, cumulative_mass_out_kg=trial.mass_out)
        row['native_semantic_sha256'] = digest(trial.__dict__)
        row['cumulative_R_total_J'] = (e1-self.initial_energy+c['E_electrical_J']+
                                       sums['Q_Rs_J']+trial.Qconduction+trial.Qhydro-sums['W_source_J'])
        if exchange <= 0 or max(abs(row[k])/exchange for k in
                                ('R_circuit_J', 'R_thermal_J', 'R_cross_J', 'R_total_J')) > 1e-3:
            raise RuntimeError('INTERVAL_ENERGY_GATE_FAILED')
        result = Trial(token, trial, c['state'], row, sums)
        result.signature = digest(trial.__dict__, result.state, row, sums)
        self.pending[id(result)] = result
        return result

    def reject(self, result):
        if self.pending.pop(id(result), None) is not result:
            raise ValueError('UNKNOWN_TRIAL')
        self.budget.discarded += 1

    def commit(self, result):
        if result.token != self.token or self.pending.get(id(result)) is not result:
            raise ValueError('STALE_FOREIGN_OR_DUPLICATE_COMMIT')
        if result.signature != digest(result.thermal.__dict__, result.state, result.row, result.sums):
            raise ValueError('MUTATED_TRIAL')
        self.thermal, self.state, self.sums = result.thermal, result.state, result.sums
        self.pending.clear()
        self.token = uuid4().hex
        return dict(result.row)
