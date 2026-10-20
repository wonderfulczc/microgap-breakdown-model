"""Provenance-bound standalone T/p + exactly zero initial velocity adapter.

No PDE, EOS or time integrator is duplicated. No constructor-state overwrite.
The frozen FP4 schema is retained; extra runtime settings live in its hashed
source_manifest_uri and current.source_uri sidecars.
"""
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal
import hashlib
import json
from pathlib import Path

import numpy as np

from streamer_rf.thermal.properties import EquilibriumAirProperties
from streamer_rf.thermal.radial import RadialGrid, RadialThermalSolver

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
SCHEMA = BASE/'cold_thermal_handoff_schema.json'


class HandoffError(ValueError):
    def __init__(self, code, reason, field=''):
        self.code, self.reason, self.field = code, reason, field
        super().__init__(f'{code}: {field}: {reason}')

    def as_dict(self):
        return dict(code=self.code, reason=self.reason, field=self.field)


def reject(code, reason, field=''):
    raise HandoffError(code, reason, field)


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_json(path):
    def real(token):
        value=float(token)
        if not np.isfinite(value) or (value==0 and Decimal(token)!=0):
            reject('NONREPRESENTABLE_NUMERIC_INPUT','nonzero underflow or overflow is not an explicit zero')
        return value
    try:
        return json.loads(Path(path).read_text(), parse_float=real, parse_constant=lambda x: reject('NONFINITE_INPUT', x))
    except (OSError, json.JSONDecodeError) as e:
        reject('SOURCE_UNAVAILABLE', str(e), str(path))


def number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        reject('INVALID_NUMERIC_INPUT', 'finite number required', name)
    if positive and value <= 0:
        reject('INVALID_NUMERIC_INPUT', 'strictly positive required', name)
    return float(value)


def required(obj, keys, name):
    if not isinstance(obj, dict):
        reject('INPUT_CONTRACT_INVALID', 'object required', name)
    absent = set(keys)-obj.keys()
    if absent:
        reject('REQUIRED_FIELD_MISSING', ','.join(sorted(absent)), name)


def array(value, name):
    def bad_element(x):
        if isinstance(x,(list,tuple)):
            return any(bad_element(y) for y in x)
        if isinstance(x,np.ndarray):
            return x.dtype.kind not in 'fiu'
        return isinstance(x,(bool,np.bool_)) or not isinstance(x,(int,float,np.number))
    if bad_element(value):
        reject('INVALID_NUMERIC_INPUT','numeric arrays required without boolean/string coercion',name)
    try:
        raw=np.asarray(value)
        if raw.dtype.kind not in 'fiu':
            reject('INVALID_NUMERIC_INPUT','numeric arrays required, not booleans or numeric strings',name)
        a = raw.astype(float)
    except (TypeError, ValueError) as e: reject('INVALID_NUMERIC_INPUT', str(e), name)
    if not np.all(np.isfinite(a)):
        reject('NONFINITE_INPUT', 'finite numeric array required', name)
    return a


@dataclass
class GivenThermalState:
    solver: RadialThermalSolver
    initial_inventory: dict
    semantic_audit: dict
    t_anchor_s: float
    segment_id: str

    def global_time(self):
        return self.t_anchor_s+self.solver.time


class GivenThermalInitialStateAdapter:
    """Only standalone radial LTE reference initialization, never a checkpoint."""
    def __init__(self, root=ROOT):
        self.root = Path(root).resolve()

    def path(self, uri):
        if not isinstance(uri, str) or not uri.strip():
            reject('SOURCE_UNAVAILABLE', 'nonempty local source URI required')
        p = (self.root/uri).resolve()
        if not p.is_relative_to(self.root) or not p.is_file():
            reject('SOURCE_UNAVAILABLE', 'source must exist within declared source root', uri)
        return p

    def validate(self, package):
        try:
            return self._validate(package)
        except HandoffError:
            raise
        except (KeyError,TypeError,ValueError,IndexError,AttributeError,OverflowError) as e:
            reject('INPUT_CONTRACT_INVALID',str(e))

    def _validate(self, package):
        if not isinstance(package, dict):
            reject('INPUT_CONTRACT_INVALID', 'object required')
        if package.get('entry_mode') == 'AUTO_COLD_TO_THERMAL' or any(
                k in package for k in ('G_R', 'W_cond', 'W_port', 'ne', 'np', 'nn')):
            reject('MISSING_PHYSICAL_CLOSURE', 'cold electrical/species fields cannot determine a thermal state')
        if package.get('entry_mode') != 'GIVEN_THERMAL_INITIAL_STATE':
            reject('INPUT_CONTRACT_INVALID', 'standalone given state required', 'entry_mode')
        for key in ('history','thermal_initial','current','port_partition','control_volumes','inventories'):
            if key in package and not isinstance(package[key],dict):
                reject('INPUT_CONTRACT_INVALID','object required',key)
        h = package.get('history', {})
        if h.get('restart_policy') != 'NEW_STANDALONE_REFERENCE_NOT_COLD_RESTART' or any(
                h.get(k) is not None for k in ('accepted_state_token', 'cold_snapshot_uri', 'circuit_state_uri', 'previous_ledger_uri', 'numerical_corrections_uri')):
            reject('RESUME_NOT_IMPLEMENTED', 'no checkpoint/history continuation')
        if package.get('control_volumes', {}).get('cold_control_volume_id') is not None:
            reject('COLD_HISTORY_MAPPING_UNAVAILABLE', 'no inherited cold material, mass or charge')
        initial = package.get('thermal_initial', {})
        if package.get('target_backend') != 'CURRENT_G1_TP_ZERO_VELOCITY' or initial.get('representation') != 'T_P_U' or any(k in initial for k in ('U_h', 'rho_m_kg_m3', 'momentum_kg_m2s', 'total_material_J_m3')):
            reject('UNSUPPORTED_INITIAL_STATE_MODE', 'only T/p + explicit u_r=0 is implemented')
        if 'u_r_m_s' not in initial:
            reject('UNKNOWN_INITIAL_VELOCITY', 'missing is unknown, not zero', 'thermal_initial.u_r_m_s')
        current = package.get('current', {})
        if current.get('role') != 'COLUMN_CONDUCTION_NOT_TOTAL_PORT' or current.get('I_col', {}).get('unit') != 'A':
            reject('CURRENT_BRANCH_MISMATCH', 'only independently declared column conduction in A')
        partition = package.get('port_partition', {})
        if partition.get('topology_id') != 'STANDALONE_PRESCRIBED_ICOL' or partition.get('capacitance_owner') != 'NOT_APPLICABLE_STANDALONE_PRESCRIBED_ICOL':
            reject('THERMAL_SERIES_COUPLING_NOT_IMPLEMENTED', 'old parallel G2 and new series circuit are not connected by this initializer')
        if package.get('inventories', {}).get('cold_material_energy', {}).get('value') is not None or package.get('inventories', {}).get('net_charge', {}).get('value') is not None:
            reject('COLD_HISTORY_MAPPING_UNAVAILABLE', 'cold material/charge inheritance is not provided')
        schema = read_json(SCHEMA)
        required(package, schema['required'], 'package')
        if set(package)-schema['properties'].keys():
            reject('INPUT_CONTRACT_INVALID', 'unknown top-level fields')
        if package['contract_version'] != 'I1-FP4-HANDOFF-1.0':
            reject('INPUT_CONTRACT_INVALID', 'unsupported contract version')
        for key in ('identity', 'control_volumes', 'time', 'energy_convention', 'port_partition', 'boundary', 'current', 'inventories', 'history'):
            spec = schema['properties'][key]
            required(package[key], spec['required'], key)
            if set(package[key])-spec['properties'].keys():
                reject('INPUT_CONTRACT_INVALID', 'unknown section fields', key)
        required(initial, ('T_K', 'p_Pa', 'u_r_m_s', 'source_id', 'velocity_policy'), 'thermal_initial')
        if set(initial)-schema['$defs']['primitive_profile']['properties'].keys():
            reject('UNSUPPORTED_INITIAL_STATE_MODE', 'only primitive profile fields accepted')
        # All values used below come from hash-checked files, including the runtime
        # sidecars. Merely supplying a nonempty provenance string is not sufficient.
        sources = {}
        if not isinstance(package['provenance'], list) or not package['provenance']:
            reject('SOURCE_UNAVAILABLE', 'provenance list required')
        for item in package['provenance']:
            required(item, ('source_id', 'uri', 'source_type', 'sha256'), 'provenance')
            if item['source_type'] not in schema['$defs']['provenance']['properties']['source_type']['enum']:
                reject('SOURCE_UNAVAILABLE', 'unknown source type')
            p = self.path(item['uri'])
            if sha(p) != item['sha256']:
                reject('SOURCE_HASH_MISMATCH', 'actual SHA256 differs', item['uri'])
            if item['source_id'] in sources:
                reject('SOURCE_UNAVAILABLE', 'duplicate source ID')
            sources[item['source_id']] = (p, item['source_type'])
        def source_by_uri(uri):
            path = self.path(uri)
            if path not in [p for p, _ in sources.values()]:
                reject('SOURCE_UNAVAILABLE', 'URI lacks hash-bound provenance', uri)
            return read_json(path)
        config = source_by_uri(package['identity']['source_manifest_uri'])
        required(config, ('grid', 'native_options', 'run_interval_local_s', 'current_time_basis', 'material', 'claims', 'energy_zero_id', 'source_type'), 'runtime_config')
        if config['source_type'] not in ('SYNTHETIC_REFERENCE', 'EXTERNAL_GIVEN_STATE'):
            reject('SOURCE_UNAVAILABLE', 'independent thermal input source type required')
        allowed_claims = {'NUMERICAL_REFERENCE_INITIAL_CONDITION', 'LTE_ASSUMED_FOR_REFERENCE', 'NOT_B0_THERMAL_CALIBRATION'}
        if set(config['claims']) != allowed_claims or not allowed_claims.issubset(package['limitations']):
            reject('UNSUPPORTED_PHYSICAL_CLAIM', 'not real B0, automatic breakdown or LTE experiment validation')
        if config['material'] != 'EQUILIBRIUM_AIR':
            reject('EOS_VERSION_MISMATCH', 'current backend supports its existing equilibrium-air EOS only')
        grid_config = config['grid']
        required(grid_config, ('cells', 'radius_m'), 'grid')
        if type(grid_config['cells']) is not int:
            reject('INVALID_GRID', 'integer cell count required')
        try:
            grid = RadialGrid(grid_config['cells'], number(grid_config['radius_m'], 'radius_m', positive=True))
        except ValueError as e:
            reject('INVALID_GRID', str(e))
        c = package['control_volumes']
        ell = number(c['ell_ch_m'], 'ell_ch_m', positive=True)
        if c['profile_sampling'] != 'CELL_CENTER':
            reject('GRID_MAPPING_UNSUPPORTED', 'no averaging, smoothing, extrapolation or conservative remapping')
        faces = array(c['radial_faces_m'], 'radial_faces_m')
        if faces.shape != (grid.cells+1,) or not np.all(np.isfinite(faces)) or not np.array_equal(faces, grid.faces):
            reject('GRID_MAPPING_UNSUPPORTED', 'faces must exactly match the native uniform radial grid')
        if c['mapping_scope'] != 'STANDALONE_GIVEN_REFERENCE_NO_COLD_MAP':
            reject('COLD_HISTORY_MAPPING_UNAVAILABLE', 'no cold-domain mapping')
        arrays = {}
        for name in ('T_K', 'p_Pa', 'u_r_m_s'):
            try: a = array(initial[name],name)
            except HandoffError as e: reject('INVALID_PROFILE',e.reason,name)
            if a.shape != (grid.cells,) or not np.all(np.isfinite(a)):
                reject('INVALID_PROFILE', 'finite one-dimensional native-cell array required', name)
            arrays[name] = a.copy()
        if np.any(arrays['u_r_m_s'] != 0) or initial['velocity_policy'] != 'EXPLICIT_ZERO_GIVEN':
            reject('UNSUPPORTED_INITIAL_STATE_MODE', 'any nonzero initial radial velocity is unsupported; no clipping')
        if initial['source_id'] not in sources:
            reject('SOURCE_UNAVAILABLE', 'profile source missing')
        if sources[initial['source_id']][1] not in ('SYNTHETIC_REFERENCE','EXTERNAL_GIVEN_STATE','DERIVED_EXISTING_DATA'):
            reject('MISSING_PHYSICAL_CLOSURE','a cold native export is not an independent thermal profile')
        supplied = read_json(sources[initial['source_id']][0])
        if any(not np.array_equal(np.asarray(supplied.get(k), float), arrays[k]) for k in arrays):
            reject('PROFILE_SOURCE_MISMATCH', 'profile arrays differ from hash-bound source')
        ec = package['energy_convention']
        air = EquilibriumAirProperties()
        if ec['coefficient_sha256'] != air.coefficient_sha256 or ec['eos_id'] != air.model_status:
            reject('EOS_VERSION_MISMATCH', 'EOS or coefficients differ from the installed native G1')
        for p in (ROOT/'python/streamer_rf/thermal/properties.py', ROOT/'python/streamer_rf/thermal/radial.py'):
            if p.resolve() not in [q for q,_ in sources.values()]:
                reject('EOS_VERSION_MISMATCH', 'native module hash provenance missing', str(p))
        if not ec['energy_zero_id'] or ec['energy_zero_id'] != config['energy_zero_id'] or ec['cold_energy_zero_mapping_id'] is not None:
            reject('ENERGY_REFERENCE_MISMATCH', 'explicit native EOS convention, no cold zero mapping')
        if ec['lte_evidence_id'] != 'LTE_ASSUMED_FOR_REFERENCE' or ec['species_scope_id'] != 'EQUILIBRIUM_AIR' or ec['nonoverlap_partition_id'] != 'EOS_THERMOCHEMICAL_INTERNAL_ENERGY_NO_EXTRA_CHEMICAL_ADDEND':
            reject('ENERGY_REFERENCE_MISMATCH', 'LTE reference and nonoverlapping material convention required')
        props = air.evaluate(arrays['T_K'], arrays['p_Pa'])
        if not np.all(props['valid']):
            reject('INITIAL_PROPERTY_INVALID', 'native EOS reports invalid T/p; no projection or clipping')
        if np.any(props['e_J_kg'] <= 0):
            reject('INITIAL_ENERGY_NONPOSITIVE', 'native initial e_m must be strictly positive')
        b = package['boundary'];options = config['native_options']
        required(options, ('cfl', 'conduction'), 'native_options')
        if type(options['conduction']) is not bool or not 0 < number(options['cfl'], 'cfl') <= .4:
            reject('INVALID_BOUNDARY_OR_OPTIONS', 'invalid native options')
        if b['radial_model'] not in ('ambient', 'closed'):
            reject('INVALID_BOUNDARY_OR_OPTIONS', 'native boundary must be ambient or closed')
        ambient_t=number(b['ambient_T_K'], 'ambient_T_K', positive=True)
        ambient_p=number(b['ambient_p_Pa'], 'ambient_p_Pa', positive=True)
        if not np.all(air.evaluate(ambient_t, ambient_p)['valid']):
            reject('INITIAL_PROPERTY_INVALID', 'native ambient properties invalid')
        if b['radiation_status'] != 'RADIATION_MODEL_NOT_ENABLED' or b['axial_electrode_loss_status'] != 'NOT_IMPLEMENTED_IN_CURRENT_G1':
            reject('UNAVAILABLE_CAPABILITY', 'radiation and axial/electrode losses not modeled')
        tm=package['time'];anchor=number(tm['global_t_h_s'], 'global_t_h_s')
        if anchor < 0 or tm['thermal_tau0_s'] != 0 or tm['cold_time_layer'] != 'NOT_APPLICABLE_STANDALONE' or not tm['new_budget_segment_id']:
            reject('TIME_BASIS_UNDECLARED', 'explicit independent clock anchor and new segment required')
        interval=array(config['run_interval_local_s'],'run_interval_local_s')
        if interval.shape!=(2,) or not np.all(np.isfinite(interval)) or interval[0]!=0 or interval[1]<=0:
            reject('TIME_SUPPORT_INVALID', 'local interval must start at zero and have positive end')
        drive_spec=source_by_uri(current['source_uri'])
        required(drive_spec, ('drive_identity','unit','kind','time_basis','support_s','parameters'), 'drive')
        if drive_spec['drive_identity']!='COLUMN_CONDUCTION_CURRENT' or drive_spec['unit']!='A':
            reject('CURRENT_BRANCH_MISMATCH', 'total port or inductor current is not column conduction')
        basis=drive_spec['time_basis']
        if basis not in ('THERMAL_LOCAL','GLOBAL') or config['current_time_basis']!=basis:
            reject('TIME_BASIS_UNDECLARED', 'explicit consistent THERMAL_LOCAL or GLOBAL required')
        support=array(drive_spec['support_s'],'support_s')
        if support.shape!=(2,) or not np.all(np.isfinite(support)) or support[1]<=support[0] or not np.array_equal(support,np.asarray(current['time_support_s'],float)):
            reject('TIME_SUPPORT_INVALID', 'finite increasing matching signal support required')
        shift=0. if basis=='THERMAL_LOCAL' else anchor
        if support[0]>shift or support[1]<shift+interval[1]:
            reject('SIGNAL_SUPPORT_INSUFFICIENT', 'no signal extrapolation allowed')
        par=drive_spec['parameters']
        if drive_spec['kind'] != 'SIN_SQUARED_OFFSET':
            reject('UNSUPPORTED_CURRENT_SIGNAL', 'only explicit parameterized sine-squared reference signal supported')
        required(par, ('base_A','modulation','period_s','origin_s'), 'signal.parameters')
        base=number(par['base_A'],'base_A');mod=number(par['modulation'],'modulation')
        period=number(par['period_s'],'period_s',positive=True);origin=number(par['origin_s'],'origin_s')
        def drive(tau):
            t=number(tau,'thermal_local_time_s')+shift
            if t<support[0] or t>support[1]:
                reject('SIGNAL_SUPPORT_INSUFFICIENT', 'requested time outside explicit support')
            result=base*(1+mod*np.sin(np.pi*(t-origin)/period)**2)
            if not np.isfinite(result):reject('NONFINITE_CURRENT','current overflow')
            return float(result)
        for tau in (0.,float(interval[1])): drive(tau)
        inv=package['inventories']
        for name, quantity in inv.items():
            required(quantity, schema['$defs']['quantity']['required'], 'inventories.'+name)
            if quantity['availability'] not in schema['$defs']['availability']['enum']:
                reject('INVENTORY_CONTRACT_INVALID','unknown availability',name)
            if quantity['value'] is not None:number(quantity['value'],name)
            if quantity['availability'] in ('MODEL_CLOSURE_REQUIRED','EXTERNAL_INITIAL_STATE_REQUIRED','NOT_IMPLEMENTED','NOT_APPLICABLE') and quantity['value'] is not None:
                reject('INVENTORY_CONTRACT_INVALID','unavailable quantities must be null',name)
        if inv['R_map']['value'] is not None or inv['R_map']['availability']!='NOT_APPLICABLE':
            reject('COLD_HISTORY_MAPPING_UNAVAILABLE','independent reference has no cold energy mapping residual')
        if inv['cold_material_energy']['value'] is not None or inv['net_charge']['value'] is not None:
            reject('COLD_HISTORY_MAPPING_UNAVAILABLE','cold inventory cannot be inherited')
        if current['source_id'] not in sources or sources[current['source_id']][0] != self.path(current['source_uri']):
            reject('SOURCE_UNAVAILABLE','drive source ID and URI disagree')
        if current['i_port']['unit']!='A' or current['i_port']['value'] is not None or current['i_port']['availability']!='NOT_APPLICABLE':
            reject('CURRENT_BRANCH_MISMATCH','this standalone reference provides no physical port current')
        volume=ell*grid.volumes
        mass=float(np.sum(props['rho_kg_m3']*volume))
        material=float(np.sum(props['rho_kg_m3']*props['e_J_kg']*volume))
        inventory=dict(mass0_kg=mass,material_energy0_J=material,kinetic_energy0_J=0.,
            total_material_energy0_J=material,energy_zero_id=ec['energy_zero_id'],
            R_map=None,R_map_status='NOT_APPLICABLE_INDEPENDENT_REFERENCE',
            prior_cold_history=None,prior_cold_history_status='NOT_APPLICABLE',
            radiation_status='MODEL_DISABLED',gas_charge_inventory_status='NOT_MODELED_BY_G1')
        for name,expected,unit in (('mass',mass,'kg'),('thermal_material_energy',material,'J'),('gas_kinetic_energy',0.,'J'),('radial_momentum',0.,'kg*m/s')):
            q=inv[name]
            if q['unit']!=unit:reject('UNIT_MISMATCH','inventory unit mismatch',name)
            if q['value'] is not None and (q['value']!=0 if expected==0 else abs(q['value']-expected)/abs(expected)>1e-12):
                reject('INITIAL_INVENTORY_MISMATCH','supplied inventory inconsistent with explicit state',name)
        audit=dict(status='INPUT_SEMANTICS_VALIDATED',source_count=len(sources),
            schema='FP4_existing_field_names_explicit_supported_subset_checks_NO_DRAFT_META_CERTIFICATION',
            velocity='EXACT_ZERO_INITIAL_ONLY_NO_THRESHOLD',sampling='EXACT_NATIVE_CELL_CENTER_NO_REMAP',
            current_branch='COLUMN_CONDUCTION_CURRENT',time_basis=basis,t_anchor_s=anchor,
            segment_id=tm['new_budget_segment_id'],material='EQUILIBRIUM_AIR',
            lte='LTE_ASSUMED_FOR_REFERENCE',profile_source_type=sources[initial['source_id']][1],
            run_end_local_s=float(interval[1]),initial_inventory=inventory,
            constructor_entry='streamer_rf.thermal.radial.RadialThermalSolver')
        return grid,ell,arrays,air,props,b,options,drive,inventory,audit

    def construct(self, package):
        p=deepcopy(package)
        grid,ell,a,air,props,b,opt,drive,inventory,audit=self.validate(p)
        try:
            solver=RadialThermalSolver(grid,ell,a['T_K'],a['p_Pa'],properties=air,
                current=drive,boundary=b['radial_model'],ambient_T=b['ambient_T_K'],
                ambient_p=b['ambient_p_Pa'],cfl=opt['cfl'],conduction=opt['conduction'])
        except ValueError as e:
            reject('NATIVE_INITIALIZATION_REJECTED',str(e))
        # Independent material inventory above did not read the constructor state.
        for name,actual,expected in [('mass_kg',solver.initial_mass,inventory['mass0_kg']),
                                     ('energy_J',solver.initial_energy,inventory['material_energy0_J'])]:
            if abs(actual-expected)/abs(expected)>1e-12:
                reject('INITIAL_INVENTORY_MISMATCH','native constructor vs independent EOS inventory',name)
        return GivenThermalState(solver,inventory,audit,p['time']['global_t_h_s'],p['time']['new_budget_segment_id'])
