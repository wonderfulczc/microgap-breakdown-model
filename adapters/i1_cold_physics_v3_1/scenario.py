"""V3.1 accepted checkpoint branch on the existing scenario contract."""
import json,subprocess,time
from pathlib import Path
from adapters.i1_cold_physics_v3.scenario import I1ScenarioRunner,CapabilityRefusal
from .parameters import bind_legacy_scenario
from .checkpoints import identity,verify,restore_outputs,sha,seal
ROOT=Path(__file__).resolve().parents[2]
def execute_process(cmd,out,data,timeout_s,checkpoint_destination=None):
    """Publish the accepted checkpoint while the child continues its next step."""
    if checkpoint_destination and (checkpoint_destination.exists() or checkpoint_destination.with_name(checkpoint_destination.name+'.incomplete').exists()):
        raise CapabilityRefusal('REFUSE_CHECKPOINT_OVERWRITE',str(checkpoint_destination))
    started=time.monotonic();published=False;reason='NORMAL_PROCESS_EXIT'
    with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
        proc=subprocess.Popen(cmd,cwd=ROOT,stdout=stdout,stderr=stderr)
        try:
            while True:
                meta=out/'checkpoint_metadata.json';payload=out/'accepted_payload.bin'
                if checkpoint_destination and not published and meta.exists() and payload.exists():
                    try:description=json.loads(meta.read_text())
                    except json.JSONDecodeError:description=None
                    if description is not None:
                        seal(payload,description,out,checkpoint_destination,data);published=True
                code=proc.poll()
                if code is not None:
                    if code:reason='NATIVE_REJECTION'
                    break
                if time.monotonic()-started>=timeout_s:
                    proc.kill();proc.wait();code=None;reason='COMPUTATIONAL_BUDGET_REACHED';break
                time.sleep(.05)
        except BaseException:
            proc.kill();proc.wait();raise
    return dict(command=cmd,returncode=code,reason=reason,wall_s=time.monotonic()-started,
                accepted_checkpoint_published=published,no_automatic_retry=True)

class RecoverableI1ScenarioRunner(I1ScenarioRunner):
    def native_extension_fields(self,data):
        return {}
    def preflight(self,data,binary,checkpoint=None):
        if 'parameter_contract' not in data:raise CapabilityRefusal('PARAMETER_CONTRACT_REQUIRED','V3.1 requires explicit source classification')
        self.validate(data)
        status=bind_legacy_scenario(data,data['parameter_contract'],ROOT)
        if sha(binary)!=data['build_sha256']:raise CapabilityRefusal('BUILD_IDENTITY_MISMATCH',str(binary))
        if data.get('capacitance_mode') not in ('TOTAL_CHARGE','SPLIT_GEOMETRIC_CAPACITANCE'):raise CapabilityRefusal('CAPACITANCE_OWNERSHIP_REQUIRED','exclusive mode')
        if data.get('circuit_assembles_Cp')!=(data['capacitance_mode']=='SPLIT_GEOMETRIC_CAPACITANCE'):raise CapabilityRefusal('DOUBLE_CAPACITANCE_OWNERSHIP','incompatible assembly')
        for p,h in data.get('runtime_artifact_hashes',{}).items():
            if sha(ROOT/p)!=h:raise CapabilityRefusal('RUNTIME_IDENTITY_MISMATCH',p)
        status['CHECKPOINT_STATUS']='NEW_ACCEPTED_SEGMENT'
        if checkpoint:
            verify(checkpoint,data);status['CHECKPOINT_STATUS']='CHECKPOINT_RESTORABLE'
        return status
    def run_segment(self,data,binary,out,checkpoint=None,timeout_s=450,checkpoint_destination=None):
        status=self.preflight(data,binary,checkpoint)
        if checkpoint:restore_outputs(checkpoint,data,out)
        else:out.mkdir(parents=True,exist_ok=False)
        (out/'capability_status.json').write_text(json.dumps(status,indent=2)+'\n')
        lines=dict(data['ion_parameters'],gamma_i=data['secondary_parameters']['gamma_i'],
          ion_source=data['ion_parameters']['source'],surface_source=data['secondary_parameters']['source'],
          qualification='CONDITIONAL_SYNTHETIC_PHYSICS',excitation_source=data['excitation']['source'])
        lines['capacitance_mode']=data['capacitance_mode']
        lines.update(self.native_extension_fields(data))
        native=out/'native_input.txt';native.write_text(''.join(f'{k}={v}\n' for k,v in lines.items()))
        payload=Path(checkpoint)/'payload.bin' if checkpoint else out/'accepted_payload.bin'
        cmd=[str(binary),'resume' if checkpoint else 'reference',identity(data),str(native),str(out),str(payload)]
        record=execute_process(cmd,out,data,timeout_s,checkpoint_destination)
        record['subprocess_restart']=bool(checkpoint)
        (out/'process_execution.json').write_text(json.dumps(record,indent=2)+'\n')
        return record
