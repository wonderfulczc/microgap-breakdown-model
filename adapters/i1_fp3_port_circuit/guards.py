"""Capability refusal, not an invented plasma port or FULL state lifting."""
from dataclasses import dataclass
import math

from .ports import PortError, finite


@dataclass(frozen=True)
class C1Context:
    geometry_id: str = "B0_R02_A30_g70_Nominal"
    domain_id: str = "LOCAL_A_FULL_GAS"
    voltage_V: float = 500.
    excitation: str = "CONSTANT_FROM_START"
    pressure_Pa: float = 101325.
    temperature_K: float = 300.
    sp3: bool = False
    width_m: float = 12e-6
    density_m3: float = 1e16
    z0_m: float = 0.


class C1Guard:
    k_star = 9134636280.4
    capabilities = dict(predict_internal_GR=True,
                        provide_physical_terminal_current=False,
                        provide_gas_heating=False,
                        accept_validated_circuit_feedback=False,
                        provide_full_state_lifting=False)

    def __init__(self, G0_S, context=C1Context(), *, elapsed_from_origin_s=0.):
        self._check_scope(context, elapsed_from_origin_s)
        finite(G0_S)
        if G0_S < 0 or elapsed_from_origin_s != 0:
            raise PortError("OUT_OF_VALIDATED_SCOPE")
        self.G0_S, self.context = G0_S, context
        self._last_absolute_time_s = 0.

    @staticmethod
    def _check_scope(c, t):
        finite(t, c.voltage_V, c.width_m, c.density_m3, c.z0_m,
               c.pressure_Pa, c.temperature_K)
        supported_points = {(12e-6, 1e16), (10.8e-6, 1e16), (13.2e-6, 1e16),
                            (12e-6, .5e16), (12e-6, 2e16)}
        if (c.geometry_id != "B0_R02_A30_g70_Nominal" or c.domain_id != "LOCAL_A_FULL_GAS"
                or c.voltage_V != 500. or c.excitation != "CONSTANT_FROM_START"
                or c.pressure_Pa != 101325. or c.temperature_K != 300. or c.sp3
                or c.z0_m != 0 or (c.width_m, c.density_m3) not in supported_points
                or not 0 <= t <= 10e-12):
            raise PortError("OUT_OF_VALIDATED_SCOPE")

    def predict_internal_GR(self, absolute_time_s):
        self._check_scope(self.context, absolute_time_s)
        if absolute_time_s < self._last_absolute_time_s:
            raise PortError("CLOCK_RESET_NOT_ALLOWED")
        self._last_absolute_time_s = absolute_time_s
        return self.G0_S*math.exp(self.k_star*absolute_time_s)

    def request(self, capability):
        if capability == "provide_full_state_lifting":
            raise PortError("STATE_LIFTING_UNAVAILABLE")
        raise PortError("UNAVAILABLE_CAPABILITY")

    def evaluate_trial(self, *args, **kwargs):
        raise PortError("UNAVAILABLE_CAPABILITY")


class FullProviderUnavailable:
    """Deliberately has no legal native state, charge, energy, or rollback."""
    def evaluate_trial(self, *args, **kwargs):
        raise PortError("FULL_PROVIDER_NOT_IMPLEMENTED")

    def lift_from_GR(self, *args, **kwargs):
        raise PortError("STATE_LIFTING_UNAVAILABLE")
