"""D'Angola et al. (2008), equations 48-59 and tables 21-23,25-26.

Coefficient provenance is recorded in thermal/g1/dangola_property_provenance.json.
Array operations allow conservative-state inversion without a tabulated EOS.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit

ATM = 101325.0
R_GAS = 8.31446261815324
CAL_G_TO_J_KG = 4184.0
COEFFICIENT_FILE = Path(__file__).with_name("dangola_coefficients.json")


def gaussian(T, c, delta):
    return np.exp(-np.square((T - c) / delta))


def sigmoid(T, c, delta):
    return expit(2 * (T - c) / delta)


def xi(T, a, c, delta, w):
    return a - c * np.exp(-np.power(T / delta, w))


def phi(T, a, w):
    return a * np.power(T, w)


class EquilibriumAirProperties:
    model_status = "LTE_REFERENCE_MODEL"
    applicability_status = "LTE_APPLICABILITY_PENDING_CALIBRATION"

    def __init__(self):
        raw = COEFFICIENT_FILE.read_bytes()
        self.coefficient_sha256 = hashlib.sha256(raw).hexdigest()
        self.tables = json.loads(raw)
        for key, count in (("21", 6), ("22", 7), ("23", 11), ("25", 7), ("26", 10)):
            table = self.tables[key]
            if len(table["terms"]) != count:
                raise ValueError("incomplete published coefficient table")
            for term in table["terms"]:
                if len(term) != 3 or any(not np.all(np.isfinite(row)) for row in term):
                    raise ValueError("invalid published coefficient row")

    @staticmethod
    def _poly(row, x):
        return np.polynomial.polynomial.polyval(x, row)

    def _sum(self, key, T, x):
        table = self.tables[key]
        values = []
        result = np.zeros_like(T, dtype=float)
        for i, (ar, cr, dr) in enumerate(table["terms"]):
            mode = table.get("amplitude_modes", ["exp"] * len(table["terms"]))[i]
            if mode == "dependent":
                a = values[2] + values[3] + values[4] - values[0] - values[1] - values[5]
            else:
                a = self._poly(ar, x)
                if mode.endswith("exp"):
                    a = np.exp(a)
                if mode.startswith("negative"):
                    a = -a
            values.append(a)
            c, d = self._poly(cr, x), self._poly(dr, x)
            if not (key == "26" and i == 0):
                c, d = np.exp(c), np.exp(d)
            fn = gaussian if (key == "23" and i >= 5) or (key == "26" and i >= 6) else sigmoid
            result = result + a * fn(T, c, d)
        return result

    def _thermo(self, T, p):
        x = np.log(p / ATM)
        M = self.tables["21"]["base"][0][0] - self._sum("21", T, x)
        c1, c2 = (self._poly(row, x) for row in self.tables["22"]["base"])
        h = CAL_G_TO_J_KG * (c1 * T + c2 * T**2 + self._sum("22", T, x))
        rho = p * M / (R_GAS * T)
        return M, rho, h, h - p / rho

    def evaluate(self, T_K, p_Pa):
        T, p = np.broadcast_arrays(np.asarray(T_K, float), np.asarray(p_Pa, float))
        inside = np.isfinite(T) & np.isfinite(p) & (T >= 50) & (T <= 60000) & (p >= .01 * ATM) & (p <= 100 * ATM)
        # Substitute only for evaluation safety; invalid cells are returned as NaN.
        Ts, ps = np.where(inside, T, 300), np.where(inside, p, ATM)
        x, logT = np.log(ps / ATM), np.log(Ts)
        with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
            M, rho, h, e = self._thermo(Ts, ps)
            c0, c1 = (self._poly(row, x) for row in self.tables["23"]["base"])
            cp = CAL_G_TO_J_KG * (c0 + c1 * Ts + self._sum("23", Ts, x))
            params = [np.exp(self._poly(row, x)) for row in self.tables["25"]["xi_log"]]
            sigma = np.exp(xi(logT, *params) + self._sum("25", Ts, x))
            k = np.exp(self._poly(self.tables["26"]["base"], x) + self._sum("26", logT, x))
        output = dict(molar_mass_kg_mol=M, rho_kg_m3=rho, h_J_kg=h, e_J_kg=e,
                      cp_J_kgK=cp, k_W_mK=k, sigma_S_m=sigma)
        physical = (M > 0) & (rho > 0) & (cp > 0) & (k >= 0) & (sigma >= 0)
        for value in output.values():
            physical &= np.isfinite(value)
        valid = inside & physical
        for key, value in output.items():
            output[key] = np.where(valid, value, np.nan)
        output["valid"] = valid
        output["status"] = np.where(~inside, "PROPERTY_OUT_OF_RANGE", np.where(physical, "VALID", "PROPERTY_FIT_INVALID"))
        return output

    def thermo_jacobian(self, T, p):
        """Bounded finite differences of the SAME rho,e EOS, not the cp fit."""
        T, p = np.broadcast_arrays(np.asarray(T, float), np.asarray(p, float))
        if np.any((T < 50) | (T > 60000) | (p < .01*ATM) | (p > 100*ATM)):
            raise ValueError("PROPERTY_OUT_OF_RANGE")
        loT, hiT = np.maximum(50, T*(1-1e-5)), np.minimum(60000, T*(1+1e-5))
        lop, hip = np.maximum(.01*ATM, p*(1-1e-5)), np.minimum(100*ATM, p*(1+1e-5))
        _, r0, _, e0 = self._thermo(loT, p)
        _, r1, _, e1 = self._thermo(hiT, p)
        _, r2, _, e2 = self._thermo(T, lop)
        _, r3, _, e3 = self._thermo(T, hip)
        return (r1-r0)/(hiT-loT), (r3-r2)/(hip-lop), (e1-e0)/(hiT-loT), (e3-e2)/(hip-lop)

    def acoustic_properties(self, T, p):
        _, rho, _, _ = self._thermo(T, p)
        rt, rp, et, ep = self.thermo_jacobian(T, p)
        det = rt*ep - rp*et
        # dp/drho|e + (p/rho^2) dp/de|rho, from de = p/rho^2 drho at ds=0.
        a2 = (-et + rt*p/rho**2) / det
        cv = et - ep*rt/rp
        if np.any(~np.isfinite(a2) | ~np.isfinite(cv) | (a2 <= 0) | (cv <= 0)):
            raise ValueError("EOS_DERIVATIVE_INVALID")
        return np.sqrt(a2), cv

    def invert(self, rho, e, T_guess=300.0, p_guess=ATM):
        """Bounded damped Newton; per-cell bounded least-squares fallback.

        Success requires a forward EOS residual below 2e-8, not just optimizer
        termination. No extrapolation or modification of conserved targets.
        """
        rho, e, T, p = np.broadcast_arrays(*map(lambda a: np.asarray(a, float), (rho, e, T_guess, p_guess)))
        T, p = T.copy(), p.copy()
        if np.any(~np.isfinite(rho) | ~np.isfinite(e) | ~np.isfinite(T) | ~np.isfinite(p) | (rho <= 0)):
            raise ValueError("EOS_INVERSION_INVALID_INPUT")
        T, p = np.clip(T, 50, 60000), np.clip(p, .01*ATM, 100*ATM)
        scale = np.maximum(np.abs(e), 1.0)
        for _ in range(18):
            _, rc, _, ec = self._thermo(T, p)
            err = np.maximum(np.abs(rc-rho)/rho, np.abs(ec-e)/scale)
            if np.all(err < 2e-10):
                break
            rt, rp, et, ep = self.thermo_jacobian(T, p)
            det = rt*ep-rp*et
            dT = ((rho-rc)*ep-rp*(e-ec))/det
            dp = (rt*(e-ec)-(rho-rc)*et)/det
            damping = np.ones_like(T)
            for _ in range(20):
                tn, pn = T+damping*dT, p+damping*dp
                inside = (tn >= 50) & (tn <= 60000) & (pn >= .01*ATM) & (pn <= 100*ATM)
                _, rn, _, en = self._thermo(np.where(inside, tn, T), np.where(inside, pn, p))
                improved = inside & (np.maximum(np.abs(rn-rho)/rho, np.abs(en-e)/scale) <= err)
                if np.all(improved | (err < 2e-10)):
                    break
                damping = np.where(improved | (err < 2e-10), damping, damping*.5)
            T = np.where(err < 2e-10, T, np.clip(T+damping*dT, 50, 60000))
            p = np.where(err < 2e-10, p, np.clip(p+damping*dp, .01*ATM, 100*ATM))
        _, rc, _, ec = self._thermo(T, p)
        bad = np.maximum(np.abs(rc-rho)/rho, np.abs(ec-e)/scale) > 2e-8
        for index in zip(*np.where(np.atleast_1d(bad))):
            idx = index if T.ndim else ()
            def residual(y):
                _, r, _, energy = self._thermo(*np.exp(y))
                return [(r-rho[idx])/rho[idx], (energy-e[idx])/scale[idx]]
            sol = least_squares(residual, np.log([T[idx], p[idx]]),
                                bounds=(np.log([50, .01*ATM]), np.log([60000, 100*ATM])),
                                xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200)
            if np.max(np.abs(residual(sol.x))) > 2e-8:
                raise ValueError("EOS_INVERSION_FAILED")
            T[idx], p[idx] = np.exp(sol.x)
        if not np.all(self.evaluate(T, p)["valid"]):
            raise ValueError("EOS_INVERSION_PROPERTY_INVALID")
        return T, p
