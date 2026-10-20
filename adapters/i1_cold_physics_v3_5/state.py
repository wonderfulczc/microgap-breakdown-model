"""Explicit state identities across independent segments; no lifting or resume."""
from dataclasses import dataclass, asdict
import math
from .assembly import require


def missing(reason):
    return dict(value=None, status='UNAVAILABLE', reason=reason)


@dataclass
class ScenarioStateEnvelope:
    global_time: float
    segment_local_time: float
    time_anchor_s: float
    segment_id: str
    physics_mode: str
    identity: dict
    electrical_circuit_state: dict
    source_state: dict
    port_charge_and_current: dict
    field_energy: dict
    material_energy_if_available: dict
    accumulated_work_and_losses: dict
    stage_capabilities: dict
    checkpoint_provenance: dict
    native_state: dict
    prior_cold_history: dict
    capacitance_owner: str
    circuit_assembles_Cp: bool
    schema_version: str = 'I1_SCENARIO_STATE_ENVELOPE_1.0'

    def to_dict(self):
        require(self.physics_mode in ('COLD_FULL', 'GIVEN_THERMAL_SERIES'), 'UNSUPPORTED_PHYSICS_MODE')
        require(all(math.isfinite(t) and t >= 0 for t in
                    (self.global_time, self.segment_local_time, self.time_anchor_s)), 'INVALID_TIME')
        require(self.global_time == self.time_anchor_s+self.segment_local_time, 'TIME_ANCHOR_MISMATCH')
        require(bool(self.segment_id), 'SEGMENT_ID_REQUIRED')
        require(all(self.identity.get(k) for k in ('geometry_id', 'reference_plane_id', 'operator_id', 'mesh_id')),
                'STATE_IDENTITY_MISSING')
        require(self.capacitance_owner in ('TOTAL_CHARGE', 'SPLIT_GEOMETRIC_CAPACITANCE') and
                self.circuit_assembles_Cp == (self.capacitance_owner == 'SPLIT_GEOMETRIC_CAPACITANCE'),
                'DOUBLE_CAPACITANCE_OWNERSHIP')
        if self.physics_mode == 'GIVEN_THERMAL_SERIES':
            require(self.prior_cold_history == dict(value=None, status='NOT_APPLICABLE_INDEPENDENT_REFERENCE'),
                    'COLD_HISTORY_MAPPING_UNAVAILABLE')
            require(self.checkpoint_provenance['kind'] == 'DIAGNOSTIC_SNAPSHOT', 'THERMAL_RESUME_NOT_IMPLEMENTED')
        value = asdict(self)
        value['units'] = {'time': 's', 'energy': 'J', 'charge': 'C', 'current': 'A', 'voltage': 'V'}
        return value


def refuse_automatic_handoff(_state=None):
    require(False, 'MISSING_PHYSICAL_CLOSURE',
            'No source-bound material energy, mass/momentum/profile projection, LTE criterion, '
            'common energy zero or matched field/capacitor/source inventories; GR/Wcond cannot lift T/p/u')
