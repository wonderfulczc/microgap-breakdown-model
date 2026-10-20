"""Source-reference mode reusing FP3 transactions and the frozen TENG circuit."""
from pathlib import Path
import time
import numpy as np
from scipy.linalg import expm
from adapters.i1_fp3_port_circuit.circuit import CircuitParameters
from adapters.i1_fp3_port_circuit.ports import CapacitanceOwner
from adapters.i1_teng_source.circuit import ReferenceParallelRCPort, TENGSeriesCircuit
from .amplifier import AmplifierSeriesCircuit
from .archive import load,dump,table
from .device import I1DeviceCasePackage,need,value
from .sources import TENGContactSeparationSource,HVAmplifierSource


def preflight(s):
    need(s.get('version')=='I1_MULTIPHASE_SCENARIO_1.0' and
         s.get('mode')=='SOURCE_NETWORK_REFERENCE','SCENARIO_VERSION_OR_MODE')
    need(s.get('qualification')=='SYNTHETIC_ENGINEERING','REAL_DEVICE_PARAMETERS_MISSING')
    need(s.get('source_type') in ('TENG','HV_AMPLIFIER'),'SOURCE_TYPE_UNAVAILABLE')
    case=I1DeviceCasePackage(load(s['device_case']));case.validate()
    need(s.get('segment_id') and s.get('time_anchor_s')==0.,'INDEPENDENT_SEGMENT_TIME')
    for k,upper in [('end_time_s',1e-7),('dt_s',1e-9),('wall_budget_s',900)]:
        need(type(s.get(k)) in (float,int) and np.isfinite(s[k]) and 0<s[k]<=upper,'RESOURCE_OR_TIME_LIMIT',k)
    need(s['dt_s']<=s['end_time_s'],'TIME_INTERVAL')
    need(s.get('port_identity')=='REFERENCE_FROZEN_R_PARALLEL_SYNTHETIC_CP','PHYSICAL_PORT_NOT_READY')
    for key in ('R_test_ohm','Cp_F'):
        need(type(s.get(key)) in (int,float) and np.isfinite(s[key]) and s[key]>0,'INVALID_MATH_PORT',key)
    need(s.get('capacitance_owner') in [x.value for x in CapacitanceOwner] and
         s.get('circuit_assembles_Cp')==(s['capacitance_owner']==CapacitanceOwner.SPLIT.value),
         'DOUBLE_CAPACITANCE_OWNERSHIP')
    initial=s.get('initial',{})
    need(initial==dict(i0_A=0.,Vc0_V=0.,Vp0_policy='DERIVED_SOURCE_OPEN_CIRCUIT',
         qualification='SYNTHETIC_EVENT_OPEN_EQUILIBRIUM_THEN_CONNECT_R',
         node_charge_policy='EXPLICIT_CAPACITOR_OFFSETS'),'INITIAL_SOURCE_STATE_UNSUPPORTED',
         'Only explicit synthetic open equilibrium; no inference of experimental precharge')
    topology=load(s['topology_source'])[s['source_type']]
    contract=load(s['source_contract'])
    if s['source_type']=='TENG':
        need(topology['ordered_elements']==['TENG+','GAP','Lext','Rs','Cext','TENG-'] and
             topology['capacitor_branch']=='SERIES' and
             all(topology[k] is None for k in ('bleed','rectifier','isolation_switch','extra_ground')),
             'TENG_TOPOLOGY_MISMATCH')
        source=TENGContactSeparationSource(contract).freeze(s['xT'],s['QT0'])
        v0=source.voltage(source.q_T0_C)
    else:
        need(topology['ordered_elements']==['AMP+','GAP','Lext','Rs','Cext||Rbleed','AMP-'] and
             topology['capacitor_branch']=='CEXT_PARALLEL_RBLEED_IN_SERIES_LOOP' and topology['extra_ground'] is None,
             'AMPLIFIER_TOPOLOGY_MISMATCH')
        value(topology['bleed'],'ohm',True)
        source=HVAmplifierSource(contract)
        source.voltage(s['end_time_s']);v0=source.voltage(0.)
    need(v0!=0,'ZERO_REFERENCE_NORMALIZATION_NOT_EVALUABLE')
    return dict(case=case,source=source,voltage0=v0,topology=topology)


def construct(s):
    d=preflight(s);p=d['case'].data['parameters']
    params=CircuitParameters(L_H=value(p['Lext_H'],'H'),Cext_F=value(p['Cext_F'],'F'),Rs_ohm=value(p['Rs_ohm'],'ohm'))
    port=ReferenceParallelRCPort(s['Cp_F'],s['R_test_ohm'],voltage0_V=d['voltage0'],
                                owner=CapacitanceOwner(s['capacitance_owner']))
    kw=dict(i0_A=s['initial']['i0_A'],Vc0_V=s['initial']['Vc0_V'],
            circuit_assembles_Cp=s['circuit_assembles_Cp'])
    if s['source_type']=='TENG': circuit=TENGSeriesCircuit(params,port,d['source'],**kw)
    else:
        circuit=AmplifierSeriesCircuit(params,port,d['source'],
             Rbleed_ohm=value(d['topology']['bleed'],'ohm',True),**kw)
    return circuit,d


def independent_matrix(c,s,v0):
    """Separate continuous state equations, not the midpoint residual/Jacobian."""
    p=c.params;R=s['R_test_ohm'];cp=s['Cp_F']
    if s['source_type']=='TENG':
        a=np.zeros((5,5));a[0]=[-p.Rs_ohm/p.L_H,-1/p.L_H,-1/p.L_H,
                               -c.source.inverse_capacitance/p.L_H,c.source.V_oc_V/p.L_H]
        a[1,0]=1/p.Cext_F;a[2,0]=1/cp;a[2,2]=-1/(R*cp);a[3,0]=1.
        x=np.array([0.,0.,v0,c.source.q_T0_C,1.])
        scale=np.array([abs(v0)/np.sqrt(p.L_H/p.Cext_F),abs(v0),abs(v0),p.Cext_F*abs(v0),1.])
    else:
        need(np.all(c.source.volts==c.source.volts[0]),'REFERENCE_REQUIRES_CONSTANT_WAVEFORM')
        l=p.L_H+c.source.Lout
        a=np.zeros((4,4));a[0]=[-(p.Rs_ohm+c.source.Rout)/l,-1/l,-1/l,c.source.voltage(0)/l]
        a[1,0]=1/p.Cext_F;a[1,1]=-1/(c.Rbleed*p.Cext_F)
        a[2,0]=1/cp;a[2,2]=-1/(R*cp)
        x=np.array([0.,0.,v0,1.]);scale=np.array([abs(v0)/np.sqrt(p.L_H/p.Cext_F),abs(v0),abs(v0),1.])
    return a*scale[None,:]/scale[:,None],x/scale,scale


def run(s,output):
    started=time.monotonic();c,d=construct(s)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    dump(out/'executed_input.json',s)
    v=d['voltage0'];p=c.params
    initial=dict(i0_A=0.,Vc0_V=0.,Vp0_V=v,qT0_C=c.q_T_C if s['source_type']=='TENG' else None,
         Cext_charge_C=0.,Cp_charge_C=c.port.Cp*v,initial_energy_J=c.initial_energy,
         KVL_open_equilibrium_V=0.,node_offsets='Explicit opposite Cp plate charges; not neutral-node assumption',
         event='Close only resistor branch; i, Vc, QT and Cp*Vp continuous; Icol jumps to Vp/R',
         R_test_ohm=s['R_test_ohm'],source_initial_identity=s['initial']['qualification'],
         real_initial_state='SOURCE_INITIAL_STATE_UNRESOLVED')
    dump(out/'initial_state.json',initial)
    failure=None
    try: c.run(s['end_time_s'],s['dt_s'],deadline=started+s['wall_budget_s'])
    except Exception as e: failure=dict(type=type(e).__name__,reason=str(e))
    if c.records: table(out/'accepted.csv',c.records)
    rows=c.records
    metrics=dict(source_type=s['source_type'],mode=s['mode'],accepted_steps=len(rows),
         final_time_s=c.accepted.time_s,trial_evaluations=c.trial_evaluations,failure=failure,
         status='PARTIAL' if failure else 'COMPLETED',native_steps=0,thermal_steps=0)
    if rows:
        a,x,sc=independent_matrix(c,s,v)
        indices=sorted(set(range(min(50,len(rows))))|set(np.linspace(0,len(rows)-1,101,dtype=int)))
        comparisons=[];error=0.
        for k in indices:
            r=rows[k];reference=expm(a*r['time_s'])@x
            measured=np.array([r['i_A'],r['Vc_V'],r['Vp_V']]+([r['q_T_C']] if s['source_type']=='TENG' else [])+[1.])/sc
            err=float(max(abs(measured[:-1]-reference[:-1])));error=max(error,err)
            comparisons.append(dict(time_s=r['time_s'],scaled_state_error=err,
                i_reference_A=float(reference[0]*sc[0]),Vc_reference_V=float(reference[1]*sc[1]),
                Vp_reference_V=float(reference[2]*sc[2]),qT_reference_C=float(reference[3]*sc[3]) if s['source_type']=='TENG' else None))
        table(out/'independent_expm_comparison.csv',comparisons)
        scales=dict(V=abs(v),Q=p.Cext_F*abs(v),W=p.Cext_F*v*v)
        quantities={'kvl_residual_V':'V','port_charge_residual_C':'Q','Cext_charge_residual_C':'Q','energy_residual_J':'W'}
        for key,scale in quantities.items():
            metrics['max_abs_'+key]=max(abs(r[key]) for r in rows)
            metrics['max_scaled_'+key]=metrics['max_abs_'+key]/scales[scale]
        last=rows[-1]
        cumulative=last['E_total_J']-c.initial_energy+last['W_Rs_cumulative_J']+last['W_column_cumulative_J']-last['W_oc_cumulative_J']
        cumulative+=last.get('W_Rout_cumulative_J',0.)+last.get('W_bleed_cumulative_J',0.)
        metrics.update(max_scaled_state_error=error,cumulative_energy_residual_J=cumulative,
            cumulative_energy_residual_scaled=abs(cumulative)/scales['W'],
            initial_voltage_V=v,normalization=scales,peak_abs_i_A=max(abs(r['i_A']) for r in rows),
            W_column_J=last['W_column_cumulative_J'],W_source_coordinate_drive_J=last['W_oc_cumulative_J'])
        metrics['engineering_pass']=not failure and error<=.01 and abs(cumulative)/scales['W']<=1e-8 and all(
             metrics['max_scaled_'+key]<=1e-8 for key in quantities)
        dump(out/'state_envelope.json',dict(global_time_s=c.accepted.time_s,segment_local_time_s=c.accepted.time_s,
            segment_id=s['segment_id'],mode=s['mode'],qualification=s['qualification'],
            device_case=s['device_case'],topology=s['topology_source'],last_accepted=last,
            capacitance_owner=s['capacitance_owner'],material_energy=None,
            material_energy_reason='Mathematical resistor only; electrical dissipation is not gas heat',
            cold_history=None,RF='NOT_VALIDATED',real_source='NOT_CALIBRATED'))
    metrics['wall_s']=time.monotonic()-started
    dump(out/'execution_record.json',metrics)
    return metrics
