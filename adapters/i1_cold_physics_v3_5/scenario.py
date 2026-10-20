"""Unified dispatcher over the existing cold and thermal providers, no new RHS."""
import argparse
import copy
from dataclasses import replace
import json
from pathlib import Path
import time
import numpy as np

from adapters.i1_cold_physics_v3_4.scenario import I1ScenarioRunner
from adapters.i1_fp3_port_circuit.ports import CapacitanceOwner
from adapters.i1_fp4_handoff.given_thermal import GivenThermalInitialStateAdapter
from adapters.i1_fp4_handoff.thermal_series.provider import ThermalSeriesProvider, Budget
from .artifacts import ROOT, dump, table, file_ref
from .assembly import I1PhysicalAssemblyPackage, bound_json, require
from .state import ScenarioStateEnvelope, missing, refuse_automatic_handoff

VERSION = 'I1_MULTIPHASE_SCENARIO_1.0'
SERIES = {'L_H': 100e-9, 'Rs_ohm': 3., 'Cext_F': 101e-12}
INITIAL = {'i_A': 0., 'Vc_V': 0., 'Vp_V': 0.}
THERMAL_SOURCE = {'id': 'SYNTHETIC_IDEAL_10V_STEP', 'voltage_V': 10.,
                  'time_basis': 'SEGMENT_LOCAL', 'before_V': 0., 'at_and_after_V': 10.}


class SmokeBudget(Budget):
    def __init__(self, deadline):
        super().__init__()
        self.deadline = deadline

    def check(self):
        require(time.monotonic() < self.deadline and self.steps < 4, 'RESOURCE_LIMIT')


class I1MultiphaseScenarioRunner(I1ScenarioRunner):
    """Two independent IVPs; metadata never grants a cold-to-hot transition."""
    def preflight_mode(self, spec):
        if spec.get('mode') == 'AUTO_COLD_TO_THERMAL':
            refuse_automatic_handoff()
        if spec.get('mode') == 'SOURCE_NETWORK_REFERENCE':
            from adapters.i1_p1.network import preflight
            return preflight(spec)
        if 'device_case' in spec or 'requested_electrical_source' in spec:
            from adapters.i1_p1.device import compatibility
            compatibility(spec)
        require(spec.get('version') == VERSION, 'SCENARIO_VERSION_MISMATCH')
        require(spec.get('mode') in ('COLD_FULL', 'GIVEN_THERMAL_SERIES'), 'UNSUPPORTED_PHYSICS_MODE')
        require(spec.get('qualification') == 'SYNTHETIC_ENGINEERING', 'REAL_DEVICE_PARAMETERS_MISSING')
        require(isinstance(spec.get('segment_id'), str) and bool(spec['segment_id']), 'SEGMENT_ID_REQUIRED')
        anchor = spec.get('time_anchor_s')
        require(isinstance(anchor, (float, int)) and np.isfinite(anchor) and anchor >= 0, 'INVALID_TIME_ANCHOR')
        require(spec.get('anchor_identity') == 'SYNTHETIC_CLOCK_NOT_HANDOFF_TIME', 'TIME_ANCHOR_IDENTITY')
        require(spec.get('prior_cold_history') is None, 'COLD_HISTORY_MAPPING_UNAVAILABLE')
        assembly = I1PhysicalAssemblyPackage(bound_json(spec.get('assembly')))
        assembly.validate()
        if spec.get('requested_source_scope') == 'FULL_GLOBAL_PHYSICAL_SOURCE':
            assembly.require_physical_source()
        limit = spec.get('wall_budget_s')
        require(isinstance(limit, (float, int)) and np.isfinite(limit) and 0 < limit <= 900, 'RESOURCE_LIMIT')
        if spec['mode'] == 'COLD_FULL':
            require(spec.get('requested_source_scope') == 'B0_LOCAL_CONDUCTOR_SOURCE', 'UNSUPPORTED_SOURCE_SCOPE')
            require(spec.get('accepted_steps') == 3, 'FIXED_NATIVE_DRIVER_THREE_STEPS_ONLY')
            data = bound_json(spec.get('native_scenario'))
            require(data['end_time_s'] == 6e-14 and data['max_native_steps'] == 3
                    and data['max_step_s'] == 2e-14, 'NATIVE_RESOURCE_IDENTITY')
            geometry = bound_json(assembly.data['local_geometry'])
            require(data['geometry_id'] == geometry['geometry_source_id'], 'GEOMETRY_IDENTITY_MISMATCH')
            binary = ROOT/spec['native_binary']['file']
            require(file_ref(binary) == spec['native_binary'], 'BUILD_IDENTITY_MISMATCH')
            super().preflight(data, binary)
            return data
        require(spec.get('requested_source_scope') == 'THERMAL_COLUMN_REFERENCE', 'UNSUPPORTED_SOURCE_SCOPE')
        require(spec.get('current_branch') == 'COLUMN_CONDUCTION_CURRENT', 'CURRENT_BRANCH_MISMATCH')
        require(spec.get('circuit') == SERIES and spec.get('circuit_initial') == INITIAL
                and spec.get('source') == THERMAL_SOURCE, 'UNSUPPORTED_SERIES_FIXTURE',
                'Existing provider implements this explicit ideal-source fixture, not arbitrary TENG/source initial states')
        require(spec.get('thermal_geometry_relation') == 'INDEPENDENT_REFERENCE_NOT_DEVICE_GEOMETRY', 'THERMAL_GEOMETRY_SCOPE')
        steps = spec.get('accepted_steps')
        require(type(steps) is int and 1 <= steps <= 4, 'THERMAL_STEP_LIMIT')
        dt = spec.get('maximum_dt_s')
        require(isinstance(dt, (float, int)) and np.isfinite(dt) and 0 < dt <= 5e-11, 'THERMAL_TIME_DOMAIN')
        cp = spec.get('Cp_F')
        require(isinstance(cp, (float, int)) and np.isfinite(cp) and cp >= 0
                and spec.get('Cp_identity') == 'SYNTHETIC_NOT_B0_OR_FULL_LOCAL', 'CAPACITANCE_IDENTITY')
        owner = spec.get('capacitance_owner')
        require(owner in [v.value for v in CapacitanceOwner] and
                spec.get('circuit_assembles_Cp') == (owner == CapacitanceOwner.SPLIT.value), 'DOUBLE_CAPACITANCE_OWNERSHIP')
        package = bound_json(spec.get('thermal_initial_package'))
        # Construction package supplies provenance and initial EOS state only. The
        # series provider owns every later Joule callback, using its accepted RMS.
        require(package['time']['thermal_tau0_s'] == 0, 'RESUME_NOT_IMPLEMENTED')
        GivenThermalInitialStateAdapter().validate(package)
        runtime = json.loads((ROOT/package['identity']['source_manifest_uri']).read_text())
        require(steps*dt <= runtime['run_interval_local_s'][1], 'THERMAL_TIME_SUPPORT')
        return package

    def run_mode(self, spec, output):
        if spec.get('mode') == 'SOURCE_NETWORK_REFERENCE':
            from adapters.i1_p1.network import run
            return run(spec, output)
        started = time.monotonic()
        data = self.preflight_mode(spec)
        out = Path(output)
        out.mkdir(parents=True, exist_ok=False)
        dump(out/'executed_input.json', spec)
        deadline = started + spec['wall_budget_s']
        if spec['mode'] == 'COLD_FULL':
            record = super().run_segment(data, ROOT/spec['native_binary']['file'], out/'native',
                                         timeout_s=max(.01, deadline-time.monotonic()))
            require(record['returncode'] == 0, 'NATIVE_EXECUTION_STOPPED_NO_RETRY', record['reason'])
            from .cold_output import consume_cold
            envelope, rows = consume_cold(self, spec, data, out)
            dump(out/'state_envelope.json', envelope)
            table(out/'cold_mode_integration_audit.csv', rows)
            result = dict(mode=spec['mode'], accepted_steps=len(rows), native_process=record,
                          science_wall_s=time.monotonic()-started)
        else:
            result = self._thermal(spec, data, out, deadline)
            result['science_wall_s'] = time.monotonic()-started
        dump(out/'execution_record.json', result)
        return result

    def _thermal(self, spec, package, out, deadline):
        package = copy.deepcopy(package)
        package['time']['global_t_h_s'] = spec['time_anchor_s']
        package['time']['new_budget_segment_id'] = spec['segment_id']
        dump(out/'derived_independent_initial_package.json', package)
        state = GivenThermalInitialStateAdapter().construct(package)
        dump(out/'initial_inventory.json', state.initial_inventory)
        dump(out/'initial_semantic_audit.json', state.semantic_audit)
        budget = SmokeBudget(deadline)
        provider = ThermalSeriesProvider(state.solver, spec['Cp_F'], CapacitanceOwner(spec['capacitance_owner']), budget)
        provider.descriptor = replace(provider.descriptor, geometry_id=package['identity']['geometry_id'],
             model_identity=(('Cp_F', spec['Cp_F']), ('ell_ch_m', state.solver.length)))
        rows = []
        for step in range(spec['accepted_steps']):
            trial = provider.evaluate_trial(provider.token, spec['maximum_dt_s'])
            row = provider.commit(trial)
            row.update(segment_id=spec['segment_id'], accepted_step=step+1,
                       global_time_s=spec['time_anchor_s']+row['time_s'])
            rows.append(row)
            table(out/f'accepted_{step+1}.csv', [row])
        table(out/'given_thermal_mode_integration_audit.csv', rows)
        s = provider.thermal
        diag = s.diagnostics()
        profile = s.profile()
        sigma = profile['sigma_S_m']
        gl = float(np.dot(sigma, s.grid.volumes))
        icol = rows[-1]['Icol_updated_state_A']
        jz = sigma*icol/gl
        error = float(np.dot(jz, s.grid.volumes)-icol)
        require(icol != 0 and abs(error/icol) <= 1e-12, 'COLUMN_CURRENT_RECONSTRUCTION')
        require(abs(diag['mass_residual_kg'])/s.initial_mass <= 1e-10, 'THERMAL_MASS_BUDGET')
        np.savez(out/'thermal_accepted_state.npz', U=s.U, T_K=s.T, p_Pa=s.p,
                 u_r_m_s=profile['u_m_s'], sigma_S_m=sigma, r_m=s.grid.centers,
                 radial_faces_m=s.grid.faces, Jz_A_m2=jz, area_m2=s.grid.volumes)
        table(out/'thermal_column_current_profile.csv', [dict(r_m=float(r), sigma_S_m=float(v),
              Jz_A_m2=float(j), area_m2=float(a)) for r, v, j, a in zip(s.grid.centers, sigma, jz, s.grid.volumes)])
        dump(out/'thermal_column_source.json', dict(identity='THERMAL_COLUMN_REFERENCE', time_layer='ACCEPTED_POST_UPDATED_R',
             time_s=s.time, ell_ch_m=s.length, G_length_S_m=gl, Rcol_ohm=s.length/gl,
             Icol_A=icol, current_area_residual_A=error, profile=file_ref(out/'thermal_accepted_state.npz'),
             rho_free_C_m3=None, reason='G1 does not evolve charge; no fabricated zero-charge SourceRecord',
             global_source='NOT_READY', propagation='NOT_RUN', assumption='AXIALLY_UNIFORM_LTE_LOCAL_COLUMN',
             reference_plane_id=provider.descriptor.reference_plane_id))
        row = rows[-1]
        identity = {k: package['identity'][k] for k in ('geometry_id', 'mesh_id', 'operator_id')}
        identity['reference_plane_id'] = provider.descriptor.reference_plane_id
        identity['initial_package_reference_plane'] = package['identity']['reference_plane_id']
        identity['reference_mapping'] = 'Same column endpoints; prescribed-initializer drive replaced by lagged-series RMS, no global map'
        envelope = ScenarioStateEnvelope(
            global_time=spec['time_anchor_s']+s.time, segment_local_time=s.time, time_anchor_s=spec['time_anchor_s'],
            segment_id=spec['segment_id'], physics_mode=spec['mode'], identity=identity,
            electrical_circuit_state=dict(i_A=row['i_port_A'], Vc_V=row['Vc_V'],
                Vp_interval_end_V=row['Vp_interval_end_V'], Vp_updated_state_V=row['Vp_updated_state_V'],
                Cext_charge_C=SERIES['Cext_F']*row['Vc_V'], Cp_charge_C=spec['Cp_F']*row['Vp_updated_state_V'],
                E_electrical_J=row['E_electrical_J']),
            source_state=dict(**spec['source'], q_T=missing('Ideal prescribed source, not TENG')),
            port_charge_and_current=dict(i_port_A=row['i_port_A'], Icol_A=icol, Icol_interval_end_A=row['Icol_interval_end_A'],
                Icol_RMS_A=row['Icol_RMS_A'], I_Cp_interval_end_A=row['I_Cp_interval_end_A'],
                cumulative_q_port_C=provider.sums['q_port_C'], current_time_layer='Endpoint signed currents; interval RMS for Joule only'),
            field_energy=missing('No resolved electrostatic field in G1; synthetic Cp energy in electrical inventory once'),
            material_energy_if_available=dict(U_material_J=diag['internal_energy_J'], K_gas_J=diag['kinetic_energy_J'],
                E_G1_J=s.total_energy(), mass_kg=row['mass_kg'], energy_zero_id=package['energy_convention']['energy_zero_id']),
            accumulated_work_and_losses=dict(**provider.sums, QJ_thermal_J=s.QJ, Qconduction_out_J=s.Qconduction,
                Qhydro_out_J=s.Qhydro, mass_out_kg=s.mass_out, radiation='MODEL_DISABLED', electrode_axial='NOT_MODELED'),
            stage_capabilities=dict(local_thermal_column=True, full_global_source=False, RF='NOT_VALIDATED',
                                    AUTO='MISSING_PHYSICAL_CLOSURE', gas_charge='UNAVAILABLE'),
            checkpoint_provenance=dict(kind='DIAGNOSTIC_SNAPSHOT', resume='NOT_IMPLEMENTED_FOR_THERMAL',
                                       source_input=file_ref(out/'derived_independent_initial_package.json')),
            native_state=file_ref(out/'thermal_accepted_state.npz'),
            prior_cold_history=dict(value=None, status='NOT_APPLICABLE_INDEPENDENT_REFERENCE'),
            capacitance_owner=spec['capacitance_owner'], circuit_assembles_Cp=spec['circuit_assembles_Cp'])
        dump(out/'state_envelope.json', envelope.to_dict())
        return dict(mode=spec['mode'], accepted_steps=s.steps, thermal_trial_step_calls=budget.steps,
                    native_rejections=budget.native_rejections, discarded_trials=budget.discarded,
                    time_s=s.time, source_current_area_residual_A=error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scenario', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(I1MultiphaseScenarioRunner().run_mode(json.loads(args.scenario.read_text()), args.output), indent=2))


if __name__ == '__main__':
    main()
