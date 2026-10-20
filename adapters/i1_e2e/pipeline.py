"""One configurable endpoint; legacy modes delegate rather than copying kernels."""
import json
from pathlib import Path
import time
import numpy as np
from adapters.i1_p1.sources import TENGContactSeparationSource
from adapters.i1_p1.archive import load
from adapters.i1_teng_source.circuit import ReferenceParallelRCPort,TENGSeriesCircuit
from adapters.i1_fp3_port_circuit.circuit import CircuitParameters
from adapters.i1_fp3_port_circuit.ports import CapacitanceOwner
from .artifacts import dump,table,sha,ROOT
from .radiation import propagate


class GeometryProvider:
    def __init__(self,c):
        self.c=c
        for key in ('gap_m','tip_radius_m','cone_full_angle_deg','effective_radius_m'):
            if not np.isfinite(c[key]) or c[key]<=0:raise ValueError('INVALID_GEOMETRY_'+key)
        if c['cone_full_angle_deg']>=180:raise ValueError('CONE_ANGLE_DOMAIN')
        if c['field_mapping']!='PLANAR_AXIAL_1D_EFFECTIVE_AREA':
            raise ValueError('GEOMETRY_FIELD_MAPPING_NOT_IMPLEMENTED')
        self.area=np.pi*c['effective_radius_m']**2
        self.length=c['gap_m']


class I1PredictivePipeline:
    @staticmethod
    def temporal_resolution(s):
        if s['channel']['kind']=='SYNTHETIC_TEST_CHANNEL':
            tau=s['channel']['R_ohm']*s['channel']['Cp_F']
            ratio=s['dt_s']/tau
            return dict(dt_over_RC=ratio,RC_s=tau,
                status='FAST_INITIAL_MODE_UNDERRESOLVED' if ratio>.1 else 'SMALL_RC_STEP_NOT_A_CONVERGENCE_PROOF',
                scope='Diagnostic only; no silent parameter or accepted-step alteration')
        return dict(status='TRANSPORT_POSITIVITY_SUBCYCLING_NOT_TEMPORAL_CONVERGENCE_CERTIFICATION')

    def run(self,s,output,deadline=None):
        if s.get('mode')=='AUTO_COLD_TO_THERMAL':
            raise ValueError('MISSING_PHYSICAL_CLOSURE')
        if s.get('mode') in ('COLD_FULL','GIVEN_THERMAL_SERIES'):
            from adapters.i1_cold_physics_v3_5.scenario import I1MultiphaseScenarioRunner
            from .interfaces import legacy_envelope
            if s['legacy_scenario']['mode']!=s['mode']:raise ValueError('PROVIDER_MODE_MISMATCH')
            result=I1MultiphaseScenarioRunner().run_mode(s['legacy_scenario'],output)
            dump(Path(output)/'channel_state_envelope.json',legacy_envelope(Path(output)/'state_envelope.json'))
            return result
        if s.get('topology')!='TENG_SERIES_L_RS_CEXT_GAP':raise ValueError('TOPOLOGY_NOT_IMPLEMENTED')
        geometry=GeometryProvider(s['geometry'])
        if any(not np.isfinite(x) or x<=0 for x in s['circuit'].values()):raise ValueError('POSITIVE_RLC_REQUIRED')
        source=TENGContactSeparationSource(load(s['source_contract'])).freeze(s['xT'],s['QT0'])
        v0=source.voltage(source.q_T0_C)
        if s['initial']!=dict(i0_A=0.,Vc0_V=0.,port_voltage='SOURCE_OPEN_EQUILIBRIUM',identity='SYNTHETIC_EVENT'):
            raise ValueError('INITIAL_CHARGING_HISTORY_NOT_SUPPORTED')
        if not (0<s['dt_s']<=s['end_s']<=1e-7):raise ValueError('TIME_RESOURCE_DOMAIN')
        out=Path(output);out.mkdir(parents=True,exist_ok=False)
        dump(out/'scenario_inputs.json',s)
        start=time.monotonic()
        deadline=min(deadline or start+2400,start+s['wall_limit_s'])
        params=CircuitParameters(**s['circuit'])
        if s['channel']['kind']=='SYNTHETIC_TEST_CHANNEL':
            port=ReferenceParallelRCPort(s['channel']['Cp_F'],s['channel']['R_ohm'],
                    voltage0_V=v0,owner=CapacitanceOwner.TOTAL)
            circuit=TENGSeriesCircuit(params,port,source,i0_A=0.,Vc0_V=0.,circuit_assembles_Cp=False)
        elif s['channel']['kind']=='PHYSICS_INFORMED_CHANNEL_ROM':
            from .rom import AxialParticlePort,ROMSeriesCircuit
            port=AxialParticlePort(s['channel'],geometry,v0,deadline)
            circuit=ROMSeriesCircuit(params,port,source,i0_A=0.,Vc0_V=0.,circuit_assembles_Cp=False)
        else:raise ValueError('CHANNEL_NOT_IMPLEMENTED')
        rows=[];plasma=[];failure=None
        count=int(np.ceil(s['end_s']/s['dt_s']))
        try:
            for k in range(1,count+1):
                if time.monotonic()>=deadline:raise TimeoutError('BUDGET_REACHED')
                dt=min(k*s['dt_s'],s['end_s'])-circuit.accepted.time_s
                row=circuit.advance(dt);rows.append(row)
                if hasattr(port,'summary'):plasma.append(dict(time_s=row['time_s'],**port.summary()))
                else:plasma.append(dict(time_s=row['time_s'],ne_m3=None,np_m3=None,
                    G_channel_S=1/port.R_test,R_channel_ohm=port.R_test,
                    model='PRESCRIBED_MATH_RESISTOR_NOT_AUTONOMOUS_DISCHARGE'))
        except Exception as e:failure=dict(type=type(e).__name__,message=str(e))
        if rows:
            table(out/'time_series.csv',rows);table(out/'circuit_state_and_energy.csv',rows)
            table(out/'plasma_state_summary.csv',plasma)
        result=dict(id=s['id'],mode=s['channel']['kind'],accepted_steps=len(rows),
            end_time_s=circuit.accepted.time_s,failure=failure,physical_device_validated=False,
            temporal_resolution=self.temporal_resolution(s))
        if hasattr(port,'state'):
            np.savez(out/'final_channel_state.npz',species_m3=port.state,phi_V=port.field(port.state,port.accepted.voltage_V)[0],
                z_m=port.z)
            table(out/'cell_continuity_audit.csv',port.cell_audits)
        if rows and failure is None:
            result['receiver']=propagate(rows,s,out)
            result['peak_current_A']=max(abs(r['i_A']) for r in rows)
            result['max_KVL_V']=max(abs(r['kvl_residual_V']) for r in rows)
            result['max_circuit_energy_residual_J']=max(abs(r['energy_residual_J']) for r in rows)
        result['wall_s']=time.monotonic()-start
        dump(out/'result.json',result)
        if rows:
            from .interfaces import predictive_envelope
            dump(out/'channel_state_envelope.json',predictive_envelope(out))
        dump(out/'provenance.json',dict(source_contract=s['source_contract'],input_sha=sha(out/'scenario_inputs.json'),
             current_code={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'adapters/i1_e2e').glob('*.py')},
             assumptions=s['geometry'],TENG='Frozen position; transferred charge genuinely dynamic',
             actual_geometry='Equivalent mapping, not device CAD',RF='CONDITIONAL_NOT_VALIDATED'))
        return result


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();print(json.dumps(I1PredictivePipeline().run(json.loads(a.input.read_text()),a.output),indent=2))

if __name__=='__main__':main()
