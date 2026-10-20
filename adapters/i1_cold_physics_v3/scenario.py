"""Minimal validated scenario entry; no replacement PDE or circuit equations."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]

class CapabilityRefusal(ValueError):
    def __init__(self, code, reason):
        self.code=code
        super().__init__(f'{code}: {reason}')

class I1ScenarioRunner:
    supported={'NATIVE_SPECIES','MECHANISM_BUDGET','INTERNAL_GR','LOCAL_CHARGE_PORT_V2','ELECTRICAL_WORK'}
    def validate(self, data):
        if 'parameter_contract' in data:
            from adapters.i1_cold_physics_v3_1.parameters import bind_legacy_scenario
            if data.get('physics_mode',data['parameter_contract']['mode'])!=data['parameter_contract']['mode']:
                raise CapabilityRefusal('PHYSICS_MODE_CONTRACT_MISMATCH','cannot relabel synthetic parameters as real device input')
            bind_legacy_scenario(data,data['parameter_contract'],ROOT)
        elif data.get('physics_mode','SYNTHETIC_ENGINEERING')!='SYNTHETIC_ENGINEERING':
            raise CapabilityRefusal('REAL_DEVICE_PARAMETERS_MISSING','V3 legacy input is synthetic only; V3.1 requires a parameter contract')
        if data.get('scenario_id')!='COLD_FULL_B0':
            raise CapabilityRefusal('UNAVAILABLE_CAPABILITY','only COLD_FULL_B0 is implemented')
        required=('geometry_id','mesh_id','boundary_id','gas_model_id','reaction_model_id',
                  'electron_transport_model_id','photoionization_mode','secondary_emission_model_id',
                  'ion_transport_model_id','electrical_source_id','requested_output_capabilities',
                  'end_time_s','max_step_s','max_native_steps','wall_budget_s','ion_parameters',
                  'secondary_parameters','source','excitation','source_hashes','build_sha256')
        for field in required:
            if field not in data: raise CapabilityRefusal('MISSING_INPUT',field)
        expected={'geometry_id':'B0_R02_A30_g70_Nominal','mesh_id':'AXISYM_200x400_h2p5um',
                  'boundary_id':'LOCAL_ZERO_NEUMANN','gas_model_id':'AIR_101325Pa_300K',
                  'reaction_model_id':'NATIVE_REACTION_MODEL','electron_transport_model_id':'ELECTRON_TRANSPORT_ISG0',
                  'ion_transport_model_id':'CONSTANT_EFFECTIVE_POSITIVE_ION',
                  'secondary_emission_model_id':'ION_IMPACT_EFFECTIVE_YIELD',
                  'electrical_source_id':'PRESCRIBED_DIFFERENTIAL_VOLTAGE',
                  'photoionization_mode':'OPAQUE_ABSORBING_METAL_PHOTON_BOUNDARY'}
        for key,value in expected.items():
            if data[key]!=value: raise CapabilityRefusal('UNSUPPORTED_MODEL',key)
        if not set(data['requested_output_capabilities'])<=self.supported:
            raise CapabilityRefusal('UNAVAILABLE_CAPABILITY','no gas heat, global port, lifting or validated RF')
        for key,ceiling in [('end_time_s',1e-11),('max_step_s',2e-14),('max_native_steps',1500),('wall_budget_s',2400)]:
            value=data[key]
            if not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<value<=ceiling:
                raise CapabilityRefusal('RESOURCE_OR_TIME_DOMAIN',key)
        if type(data['max_native_steps']) is not int: raise CapabilityRefusal('INVALID_STEP_LIMIT','integer required')
        for group,fields in [('ion_parameters',('mobility_m2_Vs','diffusion_m2_s','min_EN_Td','max_EN_Td',
                  'min_temperature_K','max_temperature_K','min_pressure_Pa','max_pressure_Pa')),
                  ('secondary_parameters',('gamma_i',)),('source',('voltage_V',)),('excitation',('ratio',))]:
            item=data[group]
            for field in (*fields,'source','qualification'):
                if field not in item: raise CapabilityRefusal('MISSING_PARAMETER_PROVENANCE',f'{group}.{field}')
            for field in fields:
                if not isinstance(item[field],(int,float)) or not math.isfinite(item[field]):
                    raise CapabilityRefusal('INVALID_PARAMETER',f'{group}.{field}')
            if not item['source'] or not item['qualification']: raise CapabilityRefusal('MISSING_PARAMETER_PROVENANCE',group)
        p=data['ion_parameters']
        if p['mobility_m2_Vs']<=0 or p['diffusion_m2_s']<0 or data['secondary_parameters']['gamma_i']<0:
            raise CapabilityRefusal('INVALID_PARAMETER','mu/D/gamma')
        for low,high in [('min_EN_Td','max_EN_Td'),('min_temperature_K','max_temperature_K'),('min_pressure_Pa','max_pressure_Pa')]:
            if p[low]<0 or p[high]<p[low]: raise CapabilityRefusal('INVALID_ION_DOMAIN',low)
        if not p['min_temperature_K']<=300<=p['max_temperature_K'] or not p['min_pressure_Pa']<=101325<=p['max_pressure_Pa']:
            raise CapabilityRefusal('ION_OUT_OF_MODEL_DOMAIN','gas conditions')
        if data['source']['voltage_V']!=500 or data['excitation']['ratio']!=.1:
            raise CapabilityRefusal('OUT_OF_AUTHORIZED_FIXTURE','only +500V and original development excitation ratio')
        for path,expected_sha in data['source_hashes'].items():
            actual=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            if actual!=expected_sha: raise CapabilityRefusal('SOURCE_IDENTITY_MISMATCH',path)
        return data

    def run(self, data, binary, output):
        self.validate(data)
        if hashlib.sha256(binary.read_bytes()).hexdigest()!=data['build_sha256']:
            raise CapabilityRefusal('BUILD_IDENTITY_MISMATCH',str(binary))
        output.mkdir(exist_ok=False,parents=True)
        (output/'scenario_input.json').write_text(json.dumps(data,indent=2)+'\n')
        lines={k:data[k] for k in ('end_time_s','max_step_s','max_native_steps','wall_budget_s')}
        lines.update(data['ion_parameters'])
        lines['gamma_i']=data['secondary_parameters']['gamma_i']
        lines['ion_source']=data['ion_parameters']['source']
        lines['surface_source']=data['secondary_parameters']['source']
        lines['excitation_source']=data['excitation']['source']
        lines['qualification']='SYNTHETIC_OR_CONDITIONAL_PHYSICS'
        native_input=output/'native_input.txt'
        native_input.write_text(''.join(f'{k}={v}\n' for k,v in lines.items()))
        started=time.monotonic()
        with (output/'native_stdout.log').open('w') as stdout,(output/'native_stderr.log').open('w') as stderr:
            try:
                proc=subprocess.run([str(binary),str(native_input),str(output)],stdout=stdout,stderr=stderr,
                                    timeout=data['wall_budget_s']+10,cwd=ROOT)
                returncode=proc.returncode
                termination='PROCESS_EXITED'
            except subprocess.TimeoutExpired:
                returncode=None
                termination='WALL_LIMIT_PROCESS_TERMINATED_NO_AUTOMATIC_RESTART'
        record={'command':[str(binary),str(native_input),str(output)],'returncode':returncode,
                'termination':termination,
                'wall_s':time.monotonic()-started,'source_identity':'GENERIC_DEV_PLUS_AUTHORIZED_V3',
                'physical_parameter_calibration':'CONDITIONAL_NOT_CALIBRATED'}
        (output/'execution.json').write_text(json.dumps(record,indent=2)+'\n')
        return record

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--binary',type=Path,required=True)
    args=parser.parse_args()
    result=I1ScenarioRunner().run(json.loads(args.input.read_text()),args.binary.resolve(),args.output.resolve())
    print(json.dumps(result,indent=2))
