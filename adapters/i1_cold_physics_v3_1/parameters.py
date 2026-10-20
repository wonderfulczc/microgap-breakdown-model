"""Versioned semantic contract. Source classification is not device calibration."""
import hashlib
import math
from pathlib import Path

VERSION='I1-ION-SURFACE-1.0'
class ParameterRefusal(ValueError):pass

def validate_parameters(c,root=Path('.')):
    if c.get('version')!=VERSION:raise ParameterRefusal('PARAMETER_CONTRACT_VERSION')
    mode=c.get('mode')
    if mode not in ('SYNTHETIC_ENGINEERING','CONDITIONAL_LITERATURE','REAL_DEVICE'):raise ParameterRefusal('PHYSICS_MODE_REQUIRED')
    requested=c.get('requested_domain',{})
    for k in ('E_over_N_Td','temperature_K','pressure_Pa'):
        bounds=requested.get(k)
        if not isinstance(bounds,list) or len(bounds)!=2 or any(type(x) not in (int,float) or not math.isfinite(x) for x in bounds) or bounds[0]<0 or bounds[1]<bounds[0]:
            raise ParameterRefusal('REQUESTED_DOMAIN_REQUIRED:'+k)
    values={};grades=[];contexts=[]
    for key,unit in [('mu_p','m2/(V*s)'),('D_p','m2/s'),('gamma_i','1')]:
        p=c.get('parameters',{}).get(key,{})
        required=('value','unit','classification','source','source_level','species_id','gas_composition',
          'validity_domain','interpolation','uncertainty','out_of_domain','surface_material','surface_state',
          'incident_ion_energy','emission_mechanism')
        if any(k not in p for k in required):raise ParameterRefusal('PARAMETER_FIELDS_MISSING:'+key)
        grade=p['classification'];grades.append(grade)
        if grade not in ('A','B','C','D'):raise ParameterRefusal('SOURCE_CLASSIFICATION_INVALID')
        if grade=='D' or p['value'] is None:raise ParameterRefusal('REAL_DEVICE_PARAMETERS_MISSING:'+key)
        if type(p['value']) not in (float,int) or not math.isfinite(p['value']) or p['value']<0 or (key=='mu_p' and p['value']==0):raise ParameterRefusal('INVALID_VALUE:'+key)
        if p['unit']!=unit or p['out_of_domain']!='REJECT' or p['interpolation']!='CONSTANT_EXPLICIT':raise ParameterRefusal('UNSUPPORTED_UNITS_OR_INTERPOLATION:'+key)
        if not p['species_id'] or not p['gas_composition'] or not p['source'] or not p['source_level']:raise ParameterRefusal('SOURCE_IDENTITY_MISSING')
        if not isinstance(p['uncertainty'],dict) or ('value' not in p['uncertainty']) or (p['uncertainty']['value'] is None and not p['uncertainty'].get('reason')):raise ParameterRefusal('UNCERTAINTY_DECLARATION_MISSING')
        for k,wanted in requested.items():
            bounds=p['validity_domain'].get(k)
            if not isinstance(bounds,list) or len(bounds)!=2 or any(type(x) not in (float,int) or not math.isfinite(x) for x in bounds) or not bounds[0]<=wanted[0]<=wanted[1]<=bounds[1]:raise ParameterRefusal('PARAMETER_DOMAIN_REJECT:'+key+':'+k)
        if grade in ('A','B') or isinstance(p['source'],dict):
            source=p['source']
            if not isinstance(source,dict) or not source.get('path') or not source.get('sha256'):raise ParameterRefusal('TRACEABLE_SOURCE_REQUIRED')
            path=root/source['path']
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=source['sha256']:raise ParameterRefusal('PARAMETER_SOURCE_SHA_MISMATCH')
        if grade=='B' and not p.get('partial_match_reason'):raise ParameterRefusal('PARTIAL_SOURCE_LIMIT_MISSING')
        contexts.append((p['species_id'],p['gas_composition'],p['validity_domain']))
        if key=='gamma_i' and p['emission_mechanism']!='ION_IMPACT_EFFECTIVE_YIELD':raise ParameterRefusal('UNSUPPORTED_EMISSION_MECHANISM')
        values[key]=p['value']
    if contexts[0]!=contexts[1] or contexts[0][:2]!=contexts[2][:2]:raise ParameterRefusal('ION_MOBILITY_DIFFUSION_SURFACE_IDENTITY_MISMATCH')
    if mode=='REAL_DEVICE':
        g=c['parameters']['gamma_i']
        if any(x!='A' for x in grades) or not c.get('device_calibration_evidence') or any(not g[k] or g[k]=='UNKNOWN' for k in ('surface_material','surface_state','incident_ion_energy')):
            raise ParameterRefusal('REAL_DEVICE_PARAMETERS_MISSING')
    elif mode=='CONDITIONAL_LITERATURE' and any(x not in ('A','B') for x in grades):raise ParameterRefusal('LITERATURE_MODE_CANNOT_USE_SYNTHETIC_VALUES')
    return dict(PHYSICS_PARAMETER_STATUS='CONDITIONAL_SYNTHETIC_PHYSICS' if 'C' in grades else 'SOURCE_DOMAIN_PARTIAL' if 'B' in grades else 'TRACEABLE_PARAMETERS',
      PARAMETER_SOURCE_STATUS='SYNTHETIC_TEST_ONLY' if 'C' in grades else 'SOURCE_DOMAIN_PARTIAL' if 'B' in grades else 'SOURCE_TRACEABLE',
      PARAMETER_DOMAIN_STATUS='WITHIN_DECLARED_TEST_DOMAIN' if 'C' in grades else 'WITHIN_DECLARED_DOMAIN',
      REAL_DEVICE_CALIBRATION_STATUS='NOT_CALIBRATED' if mode!='REAL_DEVICE' else 'USER_EVIDENCE_PROVIDED_NOT_INDEPENDENTLY_VALIDATED',
      GLOBAL_PHYSICAL_PORT='GLOBAL_PHYSICAL_PORT_NOT_READY',RF_READINESS='RF_NOT_VALIDATED')

def bind_legacy_scenario(data,contract,root):
    status=validate_parameters(contract,root)
    for key,group,field in [('mu_p','ion_parameters','mobility_m2_Vs'),('D_p','ion_parameters','diffusion_m2_s'),('gamma_i','secondary_parameters','gamma_i')]:
        if contract['parameters'][key]['value']!=data[group][field]:raise ParameterRefusal('NATIVE_PARAMETER_CONTRACT_MISMATCH:'+key)
    p=contract['parameters']['mu_p']['validity_domain'];native=data['ion_parameters']
    for k,a,b in [('E_over_N_Td','min_EN_Td','max_EN_Td'),('temperature_K','min_temperature_K','max_temperature_K'),('pressure_Pa','min_pressure_Pa','max_pressure_Pa')]:
        if [native[a],native[b]]!=p[k]:raise ParameterRefusal('NATIVE_DOMAIN_CONTRACT_MISMATCH')
    return status
