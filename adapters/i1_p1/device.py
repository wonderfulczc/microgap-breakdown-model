"""Provenance-aware device inputs; raw metrology is never silently converted."""
import copy
import math
from .archive import load
from adapters.i1_cold_physics_v3.scenario import CapabilityRefusal

SOURCE_TYPES = ('MEASURED_VERIFIED','USER_TRANSCRIBED_MEASUREMENT','NOMINAL_DESIGN',
                'DOCUMENTED_LITERATURE_PRIOR','DERIVED_CONDITIONAL','SYNTHETIC_ENGINEERING','UNKNOWN')


def need(ok, code, reason=''):
    if not ok: raise CapabilityRefusal(code,reason)


def value(q, unit, positive=False):
    need(isinstance(q,dict) and q.get('unit')==unit and q.get('source_type') in SOURCE_TYPES,
         'QUANTITY_UNIT_OR_SOURCE')
    v=q.get('value')
    need(type(v) in (int,float) and math.isfinite(v),'PARAMETER_VALUE_UNAVAILABLE')
    need(not positive or v>0,'NONPOSITIVE_PARAMETER')
    need(q['source_type']!='UNKNOWN' and q.get('source'),'PARAMETER_SOURCE_MISSING')
    payload=load(q['source'])
    for key in q.get('source_keys',[]): payload=payload[key]
    need(payload==v,'PARAMETER_SOURCE_PAYLOAD_MISMATCH')
    return float(v)


class I1DeviceCasePackage:
    def __init__(self, data): self.data=data

    def validate(self):
        d=self.data
        need(d.get('version')=='I1_DEVICE_CASE_1.0','DEVICE_CASE_VERSION')
        need(d.get('conductor_ids')==['H','Return'],'PERMANENT_CONDUCTOR_IDENTITY')
        need(d.get('geometry_branch') in ('B0_NOMINAL','B0_MEASUREMENT_CANDIDATE'),'GEOMETRY_BRANCH')
        nominal=load(d['b0_nominal_source'])
        for k,unit in [('gap_m','m'),('Rs_ohm','ohm'),('Lext_H','H'),('Cext_F','F')]:
            value(d['parameters'][k],unit,True)
        need(value(d['parameters']['gap_m'],'m')==nominal['parameters']['gap_m'],'NOMINAL_GAP_MISMATCH')
        for key,q in d['parameters'].items():
            need(q.get('source_type') in SOURCE_TYPES and 'unit' in q and 'uncertainty' in q,'PARAMETER_CONTRACT',key)
            if q['value'] is None:
                need(q['source_type']=='UNKNOWN' and bool(q.get('reason')),'UNKNOWN_PARAMETER_REASON',key)
            elif isinstance(q['value'],(int,float)):
                value(q,q['unit'])
        for key in ('electrode_material','insulator_support','second_electrode_symmetry','leads_return',
                    'probes','external_response','surface_charge_mapping','true_source_initial_state'):
            item=d.get(key)
            need(isinstance(item,dict) and item.get('value') is None and item.get('reason'),'UNSUPPORTED_UNVERIFIED_DEVICE_FIELD',key)
        return dict(status='SOURCE_BOUND_CASE_NOT_CALIBRATED',source_initial_state='SOURCE_INITIAL_STATE_UNRESOLVED')

    def geometry(self, assumptions=None):
        self.validate()
        p=copy.deepcopy(load(self.data['b0_nominal_source'])['parameters'])
        if self.data['geometry_branch']=='B0_NOMINAL':
            return dict(branch='B0_NOMINAL',parameters=p,qualification='HISTORICAL_NOMINAL_NOT_MEASURED')
        a=assumptions or {}
        need(a.get('R1_unit')=='mm' and a.get('DAG1_definition')=='CONE_FULL_ANGLE_DEG',
             'MEASUREMENT_DEFINITION_UNRESOLVED')
        p['tip_radius_m']=value(self.data['parameters']['R1_raw'],'UNCONFIRMED')*1e-3
        p['half_angle_rad']=math.radians(value(self.data['parameters']['DAG1_raw'],'UNCONFIRMED')/2)
        return dict(branch='B0_MEASUREMENT_CANDIDATE',parameters=p,
                    qualification='DERIVED_CONDITIONAL_NOT_NATIVE_GEOMETRY_APPROVED',
                    assumptions=a,D2_mapping=None,second_electrode_symmetry=None,
                    reason='Only radius/angle candidates; no CAD, mask, opposite electrode or PDE instance generated')


def compatibility(spec):
    """Optional shared identity on old modes cannot change their actual source."""
    if spec.get('requested_electrical_source') not in (None,
            'PRESCRIBED_DIFFERENTIAL_VOLTAGE' if spec.get('mode')=='COLD_FULL' else 'SYNTHETIC_IDEAL_10V_STEP'):
        need(False,'SOURCE_INITIAL_STATE_OR_PORT_QUALIFICATION_UNAVAILABLE',
             'No implicit substitution of requested TENG/amplifier by legacy constant drive')
    if 'device_case' not in spec: return
    case=I1DeviceCasePackage(load(spec['device_case']));case.validate()
    need(case.data['geometry_branch']=='B0_NOMINAL','CANDIDATE_GEOMETRY_NOT_NATIVE_BACKEND')
    if spec['mode']=='GIVEN_THERMAL_SERIES':
        need(spec.get('thermal_geometry_relation')=='INDEPENDENT_REFERENCE_NOT_DEVICE_GEOMETRY',
             'TRUE_THERMAL_INITIAL_STATE_UNAVAILABLE')
    elif spec['mode']=='COLD_FULL':
        native=load(spec['native_scenario'])
        need(native['geometry_id']==load(case.data['b0_nominal_source'])['geometry_id'],'GEOMETRY_IDENTITY_MISMATCH')
