from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, solve_ivp

from .signals import ConductanceProfile, PiecewiseLinearSignal


@dataclass(frozen=True)
class CircuitParameters:
    L_H: float
    Cext_F: float
    Cgap_F: float
    Rs_ohm: float
    conductance: ConductanceProfile
    source_voltage: PiecewiseLinearSignal
    cgap_source: str = "synthetic"

    def __post_init__(self) -> None:
        for name in ("L_H", "Cext_F", "Cgap_F"):
            value = float(getattr(self, name))
            if value < 0.0 or not np.isfinite(value):
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.L_H <= 0.0:
            raise ValueError("L_H must be positive")
        if self.Cext_F + self.Cgap_F <= 0.0:
            raise ValueError("Cext_F + Cgap_F must be positive")
        if self.Rs_ohm < 0.0 or not np.isfinite(self.Rs_ohm):
            raise ValueError("Rs_ohm must be finite and nonnegative")

    @property
    def C_total_F(self) -> float:
        return self.Cext_F + self.Cgap_F


@dataclass(frozen=True)
class GapCurrents:
    I_gap_cond_A: np.ndarray
    I_gap_disp_A: np.ndarray
    I_gap_total_A: np.ndarray
    I_Cext_A: np.ndarray


@dataclass(frozen=True)
class CircuitSolution:
    time_s: np.ndarray
    I_L_A: np.ndarray
    Vgap_V: np.ndarray
    dI_L_dt_Aps: np.ndarray
    dVgap_dt_Vps: np.ndarray
    source_voltage_V: np.ndarray
    Gb_S: np.ndarray
    currents: GapCurrents
    KCL_residual_A: np.ndarray
    KVL_residual_V: np.ndarray
    inductor_energy_J: np.ndarray
    cext_energy_J: np.ndarray
    cgap_energy_J: np.ndarray
    gap_conductive_loss_J: np.ndarray
    resistor_loss_J: np.ndarray
    source_energy_J: np.ndarray
    energy_balance_residual_J: np.ndarray

    @property
    def stored_energy_J(self) -> np.ndarray:
        return self.inductor_energy_J + self.cext_energy_J + self.cgap_energy_J

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "time_s": self.time_s,
                "I_L_A": self.I_L_A,
                "I_source_A": self.I_L_A,
                "Vgap_V": self.Vgap_V,
                "dI_L_dt_Aps": self.dI_L_dt_Aps,
                "dVgap_dt_Vps": self.dVgap_dt_Vps,
                "source_voltage_V": self.source_voltage_V,
                "Gb_S": self.Gb_S,
                "I_gap_cond_A": self.currents.I_gap_cond_A,
                "I_gap_disp_A": self.currents.I_gap_disp_A,
                "I_gap_total_A": self.currents.I_gap_total_A,
                "I_Cext_A": self.currents.I_Cext_A,
                "KCL_residual_A": self.KCL_residual_A,
                "KVL_residual_V": self.KVL_residual_V,
                "inductor_energy_J": self.inductor_energy_J,
                "cext_energy_J": self.cext_energy_J,
                "cgap_energy_J": self.cgap_energy_J,
                "gap_conductive_loss_J": self.gap_conductive_loss_J,
                "resistor_loss_J": self.resistor_loss_J,
                "source_energy_J": self.source_energy_J,
                "energy_balance_residual_J": self.energy_balance_residual_J,
            }
        )


class SeriesRLCGapCircuit:
    """Reference ODE for Vs -> Rs -> L -> node, with Cext || Cgap || Gb(t) to ground."""

    def __init__(self, params: CircuitParameters):
        self.params = params

    def rhs(self, time_s: float, state: np.ndarray) -> np.ndarray:
        i_l, v_gap = float(state[0]), float(state[1])
        vs = float(self.params.source_voltage(time_s))
        g = float(self.params.conductance(time_s))
        di_dt = (vs - self.params.Rs_ohm * i_l - v_gap) / self.params.L_H
        dv_dt = (i_l - g * v_gap) / self.params.C_total_F
        return np.array([di_dt, dv_dt], dtype=float)

    def solve(
        self,
        *,
        t_span: tuple[float, float],
        initial_state: tuple[float, float],
        output_dt_s: float,
        rtol: float = 1e-9,
        atol: float = 1e-12,
        max_step_s: float | None = None,
    ) -> CircuitSolution:
        t0, t1 = float(t_span[0]), float(t_span[1])
        if t1 <= t0:
            raise ValueError("t_span end must be greater than start")
        if output_dt_s <= 0.0 or not np.isfinite(output_dt_s):
            raise ValueError("output_dt_s must be finite and positive")
        self.params.source_voltage.require_supports(t0, t1)
        self.params.conductance.require_supports(t0, t1)
        n = int(np.floor((t1 - t0) / output_dt_s)) + 1
        t_eval = t0 + np.arange(n) * output_dt_s
        if t_eval[-1] < t1:
            t_eval = np.r_[t_eval, t1]
        max_step = max_step_s if max_step_s is not None else output_dt_s
        sol = solve_ivp(
            self.rhs,
            (t0, t1),
            np.asarray(initial_state, dtype=float),
            t_eval=t_eval,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            method="DOP853",
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        return self.evaluate_solution(sol.t, sol.y[0], sol.y[1])

    def evaluate_solution(self, time_s: np.ndarray, I_L_A: np.ndarray, Vgap_V: np.ndarray) -> CircuitSolution:
        t = np.asarray(time_s, dtype=float)
        i = np.asarray(I_L_A, dtype=float)
        v = np.asarray(Vgap_V, dtype=float)
        vs = np.asarray(self.params.source_voltage(t), dtype=float)
        g = np.asarray(self.params.conductance(t), dtype=float)
        di_dt = (vs - self.params.Rs_ohm * i - v) / self.params.L_H
        dv_dt = (i - g * v) / self.params.C_total_F

        i_gap_cond = g * v
        i_gap_disp = self.params.Cgap_F * dv_dt
        i_cext = self.params.Cext_F * dv_dt
        i_gap_total = i_gap_cond + i_gap_disp
        kcl = i - i_cext - i_gap_total
        kvl = vs - self.params.Rs_ohm * i - self.params.L_H * di_dt - v

        e_l = 0.5 * self.params.L_H * i * i
        e_cext = 0.5 * self.params.Cext_F * v * v
        e_cgap = 0.5 * self.params.Cgap_F * v * v
        p_gap = g * v * v
        p_rs = self.params.Rs_ohm * i * i
        p_source = vs * i
        gap_loss = np.r_[0.0, cumulative_trapezoid(p_gap, t)]
        rs_loss = np.r_[0.0, cumulative_trapezoid(p_rs, t)]
        source_energy = np.r_[0.0, cumulative_trapezoid(p_source, t)]
        stored = e_l + e_cext + e_cgap
        balance = (stored - stored[0]) + gap_loss + rs_loss - source_energy

        return CircuitSolution(
            time_s=t,
            I_L_A=i,
            Vgap_V=v,
            dI_L_dt_Aps=di_dt,
            dVgap_dt_Vps=dv_dt,
            source_voltage_V=vs,
            Gb_S=g,
            currents=GapCurrents(i_gap_cond, i_gap_disp, i_gap_total, i_cext),
            KCL_residual_A=kcl,
            KVL_residual_V=kvl,
            inductor_energy_J=e_l,
            cext_energy_J=e_cext,
            cgap_energy_J=e_cgap,
            gap_conductive_loss_J=gap_loss,
            resistor_loss_J=rs_loss,
            source_energy_J=source_energy,
            energy_balance_residual_J=balance,
        )


def coupling_error(time_s: np.ndarray, Vgap_circuit: np.ndarray, Vgap_reference: np.ndarray) -> float:
    t = np.asarray(time_s, dtype=float)
    vc = np.asarray(Vgap_circuit, dtype=float)
    vr = np.asarray(Vgap_reference, dtype=float)
    if t.ndim != 1 or vc.shape != t.shape or vr.shape != t.shape:
        raise ValueError("time and voltage arrays must be matching 1D arrays")
    denom = np.max(np.abs(vr))
    if denom <= 0.0:
        return float(np.max(np.abs(vc - vr)))
    return float(np.max(np.abs(vc - vr)) / denom)
