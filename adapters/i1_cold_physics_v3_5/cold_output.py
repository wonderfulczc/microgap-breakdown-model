"""Read actual accepted native outputs and use the existing Stage F consumer."""
import csv
import json
import numpy as np
from streamer_rf.rf.source.accepted_species import load_accepted_species_packet
from streamer_rf.rf.source.cell_charge import load_cell_ledger, audit_cell_charge
from streamer_rf.rf.source.conductor_geometry import I1ConductorGeometryAdapter
from streamer_rf.rf.source.global_packet import GlobalSourcePacket, UNITS
from streamer_rf.rf.source.integrals import current_moment
from streamer_rf.rf.jefimenko.axisymmetric import rotate_axisymmetric_record
from .artifacts import ROOT, dump, sha, file_ref, table
from .assembly import require
from .state import ScenarioStateEnvelope, missing


def read_values(path):
    with path.open() as f:
        rows = list(csv.reader(f))[1:]
    return {k: float(v) for k, v in rows}


def consume_cold(runner, spec, data, out):
    native = out/'native'
    execution = json.loads((native/'execution_native.json').read_text())
    require(execution['normal_exit'] and execution['native_step_calls_this_process'] == 3
            and execution['accepted_step_id'] == 3, 'NATIVE_CALL_OR_EXIT_MISMATCH')
    geometry = I1ConductorGeometryAdapter().from_native(native,
        ROOT/'adapters/i1_cold_physics_v3_1/recovery_native.cpp',
        ROOT/'adapters/comsol_fp2/b0_geometry_boundary_config.json')
    dump(out/'actual_b0_geometry.json', geometry)
    weak = [read_values(native/f'weak_step_{step}.csv') for step in (1, 2, 3)]
    rows = []
    for step, w in enumerate(weak, 1):
        index = json.loads((native/f'cell_step_{step}_index.json').read_text())
        index['files'] = {key: dict(file=p.name, sha256=sha(p)) for key, p in
                          [('cells', native/f'cell_step_{step}.csv'), ('faces', native/f'cell_step_{step}_faces.csv')]}
        sealed = native/f'cell_step_{step}_sealed.json'
        dump(sealed, index)
        _, cells, faces = load_cell_ledger(sealed)
        _, audit = audit_cell_charge(cells, faces, geometry['mesh_or_surface_mapping']['nr'])
        require(audit['max_cell_exchange_relative'] <= 1e-6 and audit['L1_relative'] <= 1e-8,
                'NATIVE_CELL_CHARGE_CONTINUITY_FAILED')
        rows.append(dict(accepted_step=step, segment_id=spec['segment_id'], local_time_s=w['time_n1_s'],
            global_time_s=spec['time_anchor_s']+w['time_n1_s'], dt_s=w['dt_s'],
            delta_q_port_C=w['stable_Deltaq_port_C'], I_port_interval_average_A=w['stable_Deltaq_port_C']/w['dt_s'],
            Wport_stable_J=w['stable_Wport_J'], Wcond_stable_J=w['stable_Wcond_J'],
            DeltaU_field_stable_J=w['stable_DeltaU_E_J'], work_residual_J=w['stable_energy_residual_J'],
            cell_charge_max_exchange_relative=audit['max_cell_exchange_relative'],
            global_charge_L1_relative=audit['L1_relative'], physics='CONDITIONAL_SYNTHETIC_PHYSICS'))
    template = ROOT/'validation/i1_cold_physics_v3_2/implementation_20261019/stage_f_source_packet.json'
    manifest = json.loads(template.read_text())
    manifest.update(source_hash=sha(ROOT/'cpp/src/streamer/StreamerSolver.cpp'), native_mesh_id=index['mesh_id'],
        native_operator_id=index['operator_id'], pre_time_s=index['pre_time_s'], post_time_s=index['post_time_s'])
    manifest['files'] = {name: sha(native/name) for name in ('species_center_source.csv', 'species_resolved_current_faces.csv')}
    dump(native/'stage_f_source_packet.json', manifest)
    records = {s: load_accepted_species_packet(native/'stage_f_source_packet.json', s)
               for s in ('electron', 'positive_ion', 'total')}
    gm = dict(version='I1-GLOBAL-SOURCE-1.1', geometry_id=manifest['geometry_id'], mesh_id=manifest['mesh_id'],
        operator_id=manifest['operator_id'], boundary_id='B0_LOCAL_ZERO_NEUMANN', reference_plane_id=manifest['reference_plane_id'],
        accepted_state_id=spec['segment_id']+':3', species_model_id=manifest['species_model_id'], parameter_source_class='C',
        source_hash=manifest['source_hash'], build_hash=data['build_sha256'], config_hash=sha(out/'executed_input.json'),
        units=UNITS, time_interval=[index['pre_time_s'], index['post_time_s']], time_layer='ACCEPTED_INTERVAL_PRE_FLUX',
        includes_displacement_current=False, secondary_emission_in_net_electron_flux=True,
        secondary_emission_added_as_source=False, negative_ion_transport='IMMOBILE_BY_MODEL',
        capacitance_owner=data['capacitance_mode'], capacitance_copies=1,
        cell_charge_ledger=dict(index_path=str(sealed), index_sha256=sha(sealed), accepted_step_id=3, nr=200),
        external_return=None, surface_charge=None, external_branches=None,
        unavailable={'external_return': 'UNKNOWN_EXPERIMENTAL_GEOMETRY', 'surface_charge': 'NO_GLOBAL_SURFACE_SOLUTION',
                     'external_operator': 'GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE'})
    consumed = runner.source('B0_LOCAL_CONDUCTOR_SOURCE', GlobalSourcePacket(gm, local_record=records['total']), geometry)
    dump(out/'global_source_packet.json', gm)
    moments = {s: np.asarray(current_moment(rotate_axisymmetric_record(r, 8))) for s, r in records.items()}
    require(np.allclose(moments['total'], moments['electron']+moments['positive_ion'], rtol=1e-12, atol=1e-28),
            'SPECIES_CURRENT_MOMENT_MISMATCH')
    table(out/'stage_f_consumer_audit.csv', [dict(species=s, Mx_A_m=float(v[0]), My_A_m=float(v[1]), Mz_A_m=float(v[2]),
          consumer='accepted_species -> global assembler -> axisymmetric rotation -> current_moment', propagation=False)
          for s, v in moments.items()])
    w = weak[-1]
    cumulative = read_values(native/'step_3_cumulative.csv')
    cumulative['stable_Wport_J'] = sum(a['stable_Wport_J'] for a in weak)
    cumulative['stable_Wcond_J'] = sum(a['stable_Wcond_J'] for a in weak)
    cp = w['C_local_F']
    energy = .5*cp*w['Vp_n1_V']**2+w['stable_U_remaining_end_J']
    identity = {k: index[k] for k in ('geometry_id', 'mesh_id', 'operator_id', 'reference_plane_id')}
    envelope = ScenarioStateEnvelope(
        global_time=spec['time_anchor_s']+w['time_n1_s'], segment_local_time=w['time_n1_s'],
        time_anchor_s=spec['time_anchor_s'], segment_id=spec['segment_id'], physics_mode='COLD_FULL', identity=identity,
        electrical_circuit_state=missing('Prescribed-voltage cold segment; no RLC dynamic state'),
        source_state=dict(id=data['electrical_source_id'], Vp_V=500., qualification='CONDITIONAL_SYNTHETIC_PHYSICS',
                          q_T=missing('Not a TENG source')),
        port_charge_and_current=dict(method_id='q_B_active_operator_stable_V2', scope='LOCAL_DISCRETE_DOMAIN',
            Q_H_C=w['end_Cp_V_C']+w['end_Q_rho_projection_C'], delta_q_port_C=w['stable_Deltaq_port_C'],
            I_port_interval_average_A=w['stable_Deltaq_port_C']/w['dt_s'], I_Hg_interval_average_A=w['Deltaq_Hg_C']/w['dt_s'],
            cumulative_q_port_C=cumulative['stable_port_charge_C'], I_col=missing('No uniform thermal column'),
            time_layer='Accepted endpoint charge; PRE_FLUX interval average currents'),
        field_energy=dict(U_field_local_stable_J=energy, Cp_F=cp,
                          U_remaining_J=w['stable_U_remaining_end_J'], scope='SAME_LOCAL_ACTIVE_NETWORK'),
        material_energy_if_available=missing('Cold model has no material energy inventory'),
        accumulated_work_and_losses=cumulative,
        stage_capabilities=dict(local_species_source=True, local_port=True, global_port=False,
                                gas_heat='UNAVAILABLE', AUTO='MISSING_PHYSICAL_CLOSURE', RF='NOT_VALIDATED'),
        checkpoint_provenance=dict(kind='ACCEPTED_CHECKPOINT', payload=file_ref(native/'step_3_accepted.bin'),
              config=file_ref(out/'executed_input.json'), status='NATIVE_ACCEPTED_PAYLOAD_NOT_NEW_RESTART_TEST',
              resume_entry='Existing V3.1 seal/verify required; envelope alone is not a restore token'),
        native_state=dict(fields=file_ref(native/'step_3_state.bin'),
              field_order=['ne', 'np', 'nn', 'rho', 'phi', 'Er', 'Ez', 'Emag', 'Sph'],
              field_shape=[400, 200], mechanism_cells=file_ref(native/'cell_step_3.csv'),
              mechanism_faces=file_ref(native/'cell_step_3_faces.csv'), cumulative=file_ref(native/'step_3_cumulative.csv')),
        prior_cold_history=dict(value=None, status='NEW_INDEPENDENT_COLD_SEGMENT'),
        capacitance_owner=data['capacitance_mode'], circuit_assembles_Cp=data['circuit_assembles_Cp'])
    return envelope.to_dict(), rows
