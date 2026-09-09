"""First-order cylindrical Euler/conduction reference solver with prescribed I.

No thermal calibration, external circuit dynamics or cold PDE coupling lives here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .properties import ATM, EquilibriumAirProperties


@dataclass(frozen=True)
class RadialGrid:
    cells: int
    radius_m: float

    def __post_init__(self):
        if not isinstance(self.cells, int) or self.cells < 3 or not np.isfinite(self.radius_m) or self.radius_m <= 0:
            raise ValueError("INVALID_RADIAL_GRID")

    @property
    def faces(self):
        return np.linspace(0, self.radius_m, self.cells + 1)

    @property
    def centers(self):
        return .5 * (self.faces[:-1] + self.faces[1:])

    @property
    def volumes(self):
        return np.pi * np.diff(self.faces**2)

    @property
    def areas(self):
        return 2*np.pi*self.faces

    @property
    def dr(self):
        return self.radius_m/self.cells


def prescribed_current_from_csv(path, extrapolation="invalid"):
    from streamer_rf.circuit.signals import PiecewiseLinearSignal
    signal = PiecewiseLinearSignal.from_csv(path, "current_A", extrapolation=extrapolation)

    def current(time_s):
        if not np.isfinite(time_s):
            raise ValueError("INVALID_CURRENT_SAMPLE_TIME")
        signal.require_supports(time_s, time_s)
        return signal(time_s)

    return current


def joule_heating(sigma, areas, length_m, current_A):
    sigma, areas = np.asarray(sigma, float), np.asarray(areas, float)
    if (sigma.shape != areas.shape or sigma.ndim != 1 or
            np.any(~np.isfinite(sigma) | (sigma < 0) | ~np.isfinite(areas) | (areas <= 0)) or
            not np.isfinite(length_m) or length_m <= 0 or not np.isfinite(current_A)):
        raise ValueError("INVALID_CONDUCTIVE_STATE")
    g = float(np.dot(sigma, areas))
    if g == 0:
        if current_A != 0:
            raise ValueError("PRESCRIBED_CURRENT_WITHOUT_CONDUCTANCE")
        return dict(G_length_S_m=0., Rsp_ohm=np.inf, Ez_V_m=0., qJ_W_m3=np.zeros_like(sigma),
                    PJ_volume_W=0., PJ_circuit_W=0., joule_closure_relative=0.)
    resistance, field = length_m/g, current_A/g
    q = sigma*field**2
    volume, circuit = float(length_m*np.dot(q, areas)), float(current_A**2*resistance)
    if not np.all(np.isfinite(q)) or not np.isfinite(circuit):
        raise ValueError("JOULE_POWER_NONFINITE")
    return dict(G_length_S_m=g, Rsp_ohm=resistance, Ez_V_m=field, qJ_W_m3=q,
                PJ_volume_W=volume, PJ_circuit_W=circuit,
                joule_closure_relative=abs(volume-circuit)/max(volume, circuit) if circuit else 0.)


def conduction_power(grid, T, k, outer_temperature=None):
    """Outward power/length at faces; harmonic k, no axis heat flux."""
    T, k = np.asarray(T, float), np.asarray(k, float)
    if T.shape != (grid.cells,) or k.shape != T.shape or np.any(~np.isfinite(T) | ~np.isfinite(k) | (k < 0)):
        raise ValueError("INVALID_CONDUCTION_STATE")
    faces = np.zeros(grid.cells+1)
    denom = k[:-1]+k[1:]
    kh = np.divide(2*k[:-1]*k[1:], denom, out=np.zeros_like(denom), where=denom > 0)
    faces[1:-1] = -grid.areas[1:-1]*kh*np.diff(T)/grid.dr
    if outer_temperature is not None:
        if not np.isfinite(outer_temperature) or outer_temperature <= 0:
            raise ValueError("INVALID_OUTER_TEMPERATURE")
        faces[-1] = -grid.areas[-1]*k[-1]*(outer_temperature-T[-1])/(.5*grid.dr)
    return -np.diff(faces)/grid.volumes, faces


class RadialThermalSolver:
    radiation_status = "RADIATION_MODEL_NOT_ENABLED"

    def __init__(self, grid: RadialGrid, length_m: float, T, p, *, properties=None,
                 current: Callable[[float], float] = lambda t: 0., boundary="ambient",
                 ambient_T=300., ambient_p=ATM, cfl=.3, conduction=True):
        if boundary not in ("ambient", "closed") or not 0 < cfl <= .4 or not np.isfinite(length_m) or length_m <= 0:
            raise ValueError("INVALID_SOLVER_CONFIGURATION")
        self.grid, self.length, self.current = grid, length_m, current
        self.properties = properties if properties is not None else EquilibriumAirProperties()
        self.T = np.broadcast_to(np.asarray(T, float), (grid.cells,)).copy()
        self.p = np.broadcast_to(np.asarray(p, float), (grid.cells,)).copy()
        self.boundary, self.ambient_T, self.ambient_p = boundary, ambient_T, ambient_p
        self.cfl, self.conduction_enabled = cfl, conduction
        initial = self.properties.evaluate(self.T, self.p)
        ambient = self.properties.evaluate(ambient_T, ambient_p)
        if not np.all(initial["valid"]) or not np.all(ambient["valid"]):
            raise ValueError("INITIAL_PROPERTY_INVALID")
        rho, e = initial["rho_kg_m3"], initial["e_J_kg"]
        if np.any(e <= 0):
            raise ValueError("INITIAL_ENERGY_NONPOSITIVE")
        self.U = np.stack((rho, np.zeros_like(rho), rho*e), axis=-1)
        self.ambient_U = np.array([float(ambient["rho_kg_m3"]), 0., float(ambient["rho_kg_m3"]*ambient["e_J_kg"])])
        self.time = 0.
        self.steps = self.rejected_steps = 0
        self.QJ = self.Qconduction = self.Qhydro = self.mass_out = 0.
        self.initial_energy = self.total_energy()
        self.initial_mass = float(np.dot(rho, grid.volumes)*length_m)

    def total_energy(self):
        return float(np.dot(self.U[:, 2], self.grid.volumes)*self.length)

    def primitive(self):
        rho = self.U[:, 0]
        if np.any(~np.isfinite(self.U)) or np.any(rho <= 0):
            raise ValueError("INVALID_CONSERVED_STATE")
        u = self.U[:, 1]/rho
        e = self.U[:, 2]/rho - .5*u**2
        if np.any(e <= 0):
            raise ValueError("NONPOSITIVE_INTERNAL_ENERGY")
        props = self.properties.evaluate(self.T, self.p)
        if not np.all(props["valid"]):
            raise ValueError("PROPERTY_STATE_INVALID")
        return rho, u, e, props

    def _rhs(self):
        grid = self.grid
        rho, u, e, props = self.primitive()
        sound, cv = self.properties.acoustic_properties(self.T, self.p)
        left_ghost = self.U[0].copy()
        left_ghost[1] *= -1
        right_ghost = self.ambient_U.copy() if self.boundary == "ambient" else self.U[-1].copy()
        if self.boundary == "closed":
            right_ghost[1] *= -1
        U = np.vstack((left_ghost, self.U, right_ghost))
        outer_p = self.ambient_p if self.boundary == "ambient" else self.p[-1]
        p = np.r_[self.p[0], self.p, outer_p]
        outer_sound = self.properties.acoustic_properties(self.ambient_T, self.ambient_p)[0] if self.boundary == "ambient" else sound[-1]
        speed = np.r_[sound[0], sound, outer_sound]
        vel = U[:, 1]/U[:, 0]
        F = np.stack((U[:, 1], U[:, 1]*vel+p, (U[:, 2]+p)*vel), axis=-1)
        wave = np.maximum(np.abs(vel[:-1])+speed[:-1], np.abs(vel[1:])+speed[1:])
        flux = .5*(F[:-1]+F[1:])-.5*wave[:, None]*np.diff(U, axis=0)
        rhs = -np.diff(grid.areas[:, None]*flux, axis=0)/grid.volumes[:, None]
        rhs[:, 1] += self.p*np.diff(grid.areas)/grid.volumes
        drive = joule_heating(props["sigma_S_m"], grid.volumes, self.length, float(self.current(self.time)))
        rhs[:, 2] += drive["qJ_W_m3"]
        heat_faces = np.zeros(grid.cells+1)
        if self.conduction_enabled:
            heat, heat_faces = conduction_power(grid, self.T, props["k_W_mK"],
                                               self.ambient_T if self.boundary == "ambient" else None)
            rhs[:, 2] += heat
        # Geometric CFL bounds are based on both face areas, including the axis cell.
        dt_hydro = self.cfl*np.min(grid.volumes/(grid.areas[:-1]+grid.areas[1:])/(np.abs(u)+sound))
        k = props["k_W_mK"]
        dt_heat = .1*np.min(rho*cv*grid.dr**2/np.maximum(k, np.finfo(float).tiny)) if self.conduction_enabled else np.inf
        positive_q = drive["qJ_W_m3"] > 0
        dt_joule = .05*np.min((rho*e)[positive_q]/drive["qJ_W_m3"][positive_q]) if np.any(positive_q) else np.inf
        return rhs, min(dt_hydro, dt_heat, dt_joule), drive, heat_faces[-1]*self.length, flux[-1, 2]*grid.areas[-1]*self.length, flux[-1, 0]*grid.areas[-1]*self.length

    def step(self, maximum_dt):
        if not np.isfinite(maximum_dt) or maximum_dt <= 0:
            raise ValueError("INVALID_TIMESTEP")
        rhs, dt_limit, drive, heat_out, hydro_out, mass_out = self._rhs()
        dt = min(maximum_dt, dt_limit)
        for _ in range(16):
            trial = self.U + dt*rhs
            try:
                if not np.all(np.isfinite(trial)) or np.any(trial[:, 0] <= 0):
                    raise ValueError("NONPOSITIVE_DENSITY")
                u = trial[:, 1]/trial[:, 0]
                e = trial[:, 2]/trial[:, 0] - .5*u**2
                if np.any(e <= 0):
                    raise ValueError("NONPOSITIVE_INTERNAL_ENERGY")
                T, p = self.properties.invert(trial[:, 0], e, self.T, self.p)
                break
            except ValueError:
                self.rejected_steps += 1
                dt *= .5
        else:
            raise ValueError("TIMESTEP_REJECTED_PROPERTY_OR_POSITIVITY")
        self.U, self.T, self.p = trial, T, p
        self.QJ += dt*drive["PJ_volume_W"]
        self.Qconduction += dt*heat_out
        self.Qhydro += dt*hydro_out
        self.mass_out += dt*mass_out
        self.time += dt
        self.steps += 1
        return dt

    def diagnostics(self):
        rho, u, e, props = self.primitive()
        d = joule_heating(props["sigma_S_m"], self.grid.volumes, self.length, float(self.current(self.time)))
        d.pop("qJ_W_m3")
        kinetic = float(np.dot(.5*rho*u**2, self.grid.volumes)*self.length)
        excess = self.T-self.ambient_T
        # Connected-from-axis 10%-of-peak temperature-excess is a numerical radius.
        radius = np.nan
        if excess[0] > 0 and np.max(excess) > 0:
            threshold = .1*np.max(excess)
            outside = np.flatnonzero(excess < threshold)
            if outside.size and outside[0] > 0:
                i = outside[0]
                radius = float(np.interp(threshold, [excess[i], excess[i-1]], [self.grid.centers[i], self.grid.centers[i-1]]))
        mass = float(np.dot(rho, self.grid.volumes)*self.length)
        d.update(time_s=self.time, accepted_steps=self.steps, rejected_steps=self.rejected_steps,
                 T_axis_K=float(self.T[0]), T_max_K=float(np.max(self.T)), T_min_K=float(np.min(self.T)),
                 p_axis_Pa=float(self.p[0]), p_min_Pa=float(np.min(self.p)), rho_axis_kg_m3=float(rho[0]),
                 rho_min_kg_m3=float(np.min(rho)), sigma_axis_S_m=float(props["sigma_S_m"][0]),
                 sigma_max_S_m=float(np.max(props["sigma_S_m"])), channel_radius_m=radius,
                 outward_velocity_max_m_s=float(np.max(u)), internal_energy_J=self.total_energy()-kinetic,
                 kinetic_energy_J=kinetic, total_energy_J=self.total_energy(), cumulative_Joule_J=self.QJ,
                 conduction_boundary_loss_J=self.Qconduction, hydro_boundary_loss_J=self.Qhydro,
                 radiation_loss_J=0., energy_residual_J=self.total_energy()-self.initial_energy-self.QJ+self.Qconduction+self.Qhydro,
                 mass_residual_kg=mass-self.initial_mass+self.mass_out)
        return d

    def profile(self):
        rho, u, _, props = self.primitive()
        return dict(r_m=self.grid.centers, T_K=self.T.copy(), p_Pa=self.p.copy(), rho_kg_m3=rho.copy(),
                    u_m_s=u.copy(), sigma_S_m=props["sigma_S_m"].copy())


def thermal_initialization_from_g0(summary):
    """Cold scalar observables cannot specify a calibrated radial thermal state."""
    candidate = summary.get("G1_handoff_candidate_contract", {})
    if not candidate.get("state_valid") or candidate.get("thermal_profile_status") == "NOT_AVAILABLE":
        return {"valid": False, "status": "G0_THERMAL_INITIALIZATION_UNRESOLVED"}
    return {"valid": False, "status": "CALIBRATED_RADIAL_PROFILE_REQUIRED"}


def validate_calibrated_radial_state(grid, T, p, source_id, source_reference, properties):
    if not source_id or not source_reference or np.shape(T) != (grid.cells,) or np.shape(p) != (grid.cells,):
        raise ValueError("INVALID_CALIBRATED_RADIAL_CONTRACT")
    if not np.all(properties.evaluate(T, p)["valid"]):
        raise ValueError("CALIBRATED_PROFILE_OUT_OF_RANGE")
    return {"T_K": np.array(T, copy=True), "p_Pa": np.array(p, copy=True),
            "source_id": source_id, "source_reference": source_reference}


def required_thermal_energy(mass_kg, cold_e_J_kg, target_e_J_kg):
    """Internal-energy difference for explicitly matched material masses.

    The caller must supply the target. Expansion work/losses are not included.
    """
    m, cold, target = np.broadcast_arrays(*map(lambda x: np.asarray(x, float), (mass_kg, cold_e_J_kg, target_e_J_kg)))
    if np.any(~np.isfinite(m) | ~np.isfinite(cold) | ~np.isfinite(target) | (m <= 0) | (target < cold)):
        raise ValueError("INVALID_THERMAL_TARGET")
    return float(np.sum(m*(target-cold)))
