"""Source-bound extension of the accepted V3.1 scenario/recovery protocol."""
import copy,csv,json,math
from pathlib import Path
from adapters.i1_cold_physics_v3_1.scenario import RecoverableI1ScenarioRunner
from adapters.i1_cold_physics_v3.scenario import CapabilityRefusal
from .archive import ROOT,sha

def validate_package(data):
    p=data['ion_field_package'];s=data['surface_yield_package']
    required=('ion_species_id','neutral_gas_composition','field_unit','mobility_unit','diffusion_unit',
      'temperature_K','pressure_Pa','EN_domain_Td','source_doi','database_id','file','file_sha256',
      'evidence_class','aggregate_mapping','interpolation','out_of_domain','diffusion_status','uncertainty')
    if any(k not in p for k in required):raise CapabilityRefusal('ION_TABLE_FIELDS_MISSING','explicit metadata required')
    if p['ion_species_id']!='EFFECTIVE_SINGLE_POSITIVE_ION' or p['neutral_gas_composition']!='SYNTHETIC_AIR_EFFECTIVE_MIXTURE':
        raise CapabilityRefusal('SMOKE_SPECIES_CONTEXT_MISMATCH','this native fixture has a fixed effective-air identity; generic C++ provider supports explicit species separately')
    if p['field_unit']!='Td' or p['diffusion_unit']!='m2/s' or p['mobility_unit'] not in ('m2/(V*s)','K0_m2/(V*s)'):
        raise CapabilityRefusal('ION_TABLE_UNITS','SI only, reduced mobility requires Nref')
    if p['mobility_unit'].startswith('K0') and (not isinstance(p.get('reference_density_m3'),(int,float)) or not math.isfinite(p['reference_density_m3']) or p['reference_density_m3']<=0):
        raise CapabilityRefusal('REDUCED_MOBILITY_REFERENCE_REQUIRED','Nref must be explicit')
    if p['aggregate_mapping']!='EFFECTIVE_SPECIES_APPROXIMATION' or p['out_of_domain']!='REJECT' or p['interpolation']!='LINEAR_EN':
        raise CapabilityRefusal('ION_TABLE_METHOD','no extrapolation or unacknowledged aggregated np')
    if p['evidence_class'] not in ('A','B','C') or p['diffusion_status']!='AVAILABLE':raise CapabilityRefusal('ION_TRANSPORT_DATA_UNAVAILABLE','cannot invent diffusion')
    if not p['source_doi'] or not p['database_id'] or sha(ROOT/p['file'])!=p['file_sha256']:
        raise CapabilityRefusal('ION_TABLE_SHA_OR_SOURCE','file identity mismatch')
    with (ROOT/p['file']).open() as f:rows=list(csv.DictReader(f))
    previous=-1
    for row in rows:
        vals=[float(row[k]) for k in ('EN_Td','mobility','D_parallel_m2_s','D_perpendicular_m2_s')]
        if not all(math.isfinite(v) for v in vals) or vals[0]<=previous or vals[0]<0 or vals[1]<=0 or min(vals[2:])<0:
            raise CapabilityRefusal('ION_TABLE_VALUES','finite monotone field and positive mobility required')
        if vals[2]!=vals[3]:raise CapabilityRefusal('ANISOTROPIC_DIFFUSION_FV_NOT_IMPLEMENTED','provider can store tensor data; native scalar FV cannot use it')
        previous=vals[0]
    if len(rows)<2 or p['EN_domain_Td']!=[float(rows[0]['EN_Td']),float(rows[-1]['EN_Td'])]:
        raise CapabilityRefusal('ION_TABLE_DOMAIN','endpoint mismatch')
    native=data['ion_parameters']
    for key,low,high in [('EN_domain_Td','min_EN_Td','max_EN_Td'),('temperature_K','min_temperature_K','max_temperature_K'),('pressure_Pa','min_pressure_Pa','max_pressure_Pa')]:
        if p[key]!=[native[low],native[high]]:raise CapabilityRefusal('ION_NATIVE_DOMAIN_MISMATCH',key)
    needed=('incident_ion_species','gas','cathode_material','surface_state','incident_energy_assumption',
      'source_file','source_sha256','evidence_class','uncertainty','gamma_i','emission_energy','near_electrode_drop')
    if any(k not in s for k in needed):raise CapabilityRefusal('SURFACE_FIELDS_MISSING','explicit source contract')
    if sha(ROOT/s['source_file'])!=s['source_sha256'] or s['incident_ion_species']!=p['ion_species_id'] or s['gas']!=p['neutral_gas_composition']:
        raise CapabilityRefusal('SURFACE_SOURCE_IDENTITY','source/species/gas mismatch')
    if s['gamma_i']!=data['secondary_parameters']['gamma_i'] or s['emission_energy']!='UNAVAILABLE' or s['near_electrode_drop']!='UNAVAILABLE':
        raise CapabilityRefusal('SURFACE_UNSUPPORTED_PHYSICS','no fabricated energy or drop')
    if p['evidence_class']!='C' or s['evidence_class']!='C' or data.get('physics_mode')!='SYNTHETIC_ENGINEERING':
        raise CapabilityRefusal('REAL_DEVICE_PARAMETERS_MISSING','this bounded fixture is synthetic; source-backed device package not available')
    return p,s

class SourceBoundI1ScenarioRunner(RecoverableI1ScenarioRunner):
    def validate(self,data):
        if data.get('ion_transport_model_id')!='SOURCE_BOUND_EN_TABLE':raise CapabilityRefusal('V32_MODEL_ID','source-bound table required')
        validate_package(data)
        compat=copy.deepcopy(data)
        # Reuse immutable B0, source, resource and legacy scalar-parameter identity checks.
        # Native input below explicitly selects TABLE; compatibility scalars are not used as transport.
        compat['ion_transport_model_id']='CONSTANT_EFFECTIVE_POSITIVE_ION'
        super().validate(compat)
        return data
    def native_extension_fields(self,data):
        p,s=validate_package(data)
        return dict(table_file=str(ROOT/p['file']),table_sha256=p['file_sha256'],table_gas=p['neutral_gas_composition'],
          table_doi=p['source_doi'],table_database=p['database_id'],table_class=p['evidence_class'],
          table_mu_unit=p['mobility_unit'],table_Nref=p.get('reference_density_m3',0),
          surface_material=s['cathode_material'],surface_state=s['surface_state'],surface_energy_assumption=s['incident_energy_assumption'],
          surface_sha256=s['source_sha256'],surface_class=s['evidence_class'],surface_uncertainty=s['uncertainty'])
