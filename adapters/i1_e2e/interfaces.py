"""Common result contract; differing time layers and unavailable energy stay explicit."""
from dataclasses import dataclass,asdict
import csv
import json
from pathlib import Path
from .artifacts import sha


@dataclass(frozen=True)
class ChannelStateEnvelope:
    physics_mode: str
    time_s: float
    voltage_V: float
    currents: dict
    charge: dict
    energy: dict
    state_source: dict
    scope: str
    version: str='I1_CHANNEL_PROVIDER_RESULT_1.0'

    def to_dict(self):return asdict(self)


def legacy_envelope(path):
    path=Path(path);s=json.loads(path.read_text());mode=s['physics_mode']
    if mode not in ('COLD_FULL','GIVEN_THERMAL_SERIES'):raise ValueError('LEGACY_MODE_NOT_SUPPORTED')
    v=s['source_state']['Vp_V'] if mode=='COLD_FULL' else s['electrical_circuit_state']['Vp_interval_end_V']
    return ChannelStateEnvelope(mode,s['global_time'],v,
        currents=s['port_charge_and_current'],charge=s['port_charge_and_current'],
        energy=dict(field=s['field_energy'],material=s['material_energy_if_available'],
                    interval_or_cumulative_ledger=s['accumulated_work_and_losses']),
        state_source=dict(path=str(path),sha256=sha(path),identity=s['identity'],checkpoint=s['checkpoint_provenance']),
        scope='LEGACY_LOCAL_PROVIDER_NO_AUTOMATIC_TRANSITION').to_dict()


def predictive_envelope(path):
    path=Path(path);s=json.loads((path/'scenario_inputs.json').read_text())
    with (path/'time_series.csv').open() as f:rows=list(csv.DictReader(f))
    r=rows[-1];mode=s['channel']['kind']
    num=lambda k:float(r[k]) if r.get(k) not in (None,'') else None
    return ChannelStateEnvelope(mode,num('time_s'),num('Vp_V'),
        currents=dict(i_port_endpoint_A=num('i_A'),i_port_mid_A=num('I_port_mid_A'),
            electrode_conduction_interval_A=num('Icond_H_mid_A') if 'Icond_H_mid_A' in r else num('Icol_mid_A'),
            displacement_interval_A=num('Idisplacement_H_mid_A') if 'Idisplacement_H_mid_A' in r else num('I_Cp_mid_A'),
            time_layer='Endpoint circuit i distinct from accepted-interval boundary averages'),
        charge=dict(interval_port_C=num('delta_q_port_C'),interval_electrode_storage_C=num('delta_Q_selected_C'),
                    source_transferred_charge_C=num('q_T_C'),Cext_charge_C=num('Q_Cext_C')),
        energy=dict(field_J=num('field_energy_J'),electrical_total_J=num('E_total_J'),
            interval_port_work_J=num('W_port_J'),interval_cond_work_J=num('Wcond_face_J') if 'Wcond_face_J' in r else num('W_column_J'),
            material_J=None,material_reason='No material energy equation in this provider'),
        state_source=dict(path=str(path/'scenario_inputs.json'),sha256=sha(path/'scenario_inputs.json')),
        scope='CONDITIONAL_CHANNEL_MODEL_NOT_FULL_GLOBAL_PORT').to_dict()
