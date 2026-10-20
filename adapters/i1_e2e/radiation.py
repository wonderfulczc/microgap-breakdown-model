"""Finite supported equivalent sources -> existing Stage F/H -> receiver voltage."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from streamer_rf.rf.source.global_packet import ConservativeGraph,GlobalSourcePacket,I1GlobalSourceAssembler,save_global_source_bundle
from .source_identity import identity
from streamer_rf.rf.jefimenko.solver import evaluate_waveform
from streamer_rf.rf.jefimenko.observer import Observer
from streamer_rf.rf.jefimenko.constants import C0
from streamer_rf.fullwave.native_receiver import project_field
from .artifacts import dump,table,sha


class ReceiverVoltageModel:
    """Explicit Thevenin effective-height reference, optionally causal sampled FIR."""
    def __init__(self,c):
        self.c=c;self.h=np.asarray(c['effective_height_vector_m'],float)
        if self.h.shape!=(3,) or not np.all(np.isfinite(self.h)) or np.linalg.norm(self.h)==0:
            raise ValueError('RECEIVER_EFFECTIVE_HEIGHT_REQUIRED')
        self.r=c['antenna_resistance_ohm'];self.load=c['load_resistance_ohm']
        if not np.isfinite(self.r+self.load) or min(self.r,self.load)<=0:raise ValueError('RECEIVER_IMPEDANCE_REQUIRED')
        if c['identity']!='CONDITIONAL_EFFECTIVE_HEIGHT_NOT_MEASURED_ANTENNA':raise ValueError('RECEIVER_QUALIFICATION')
        if c.get('Hrx') and c['Hrx'].get('type')!='CAUSAL_FIR':raise ValueError('H_RX_CAUSALITY_UNAVAILABLE')

    def evaluate(self,E,dt):
        direction=self.h/np.linalg.norm(self.h)
        projection=project_field(E,direction)  # Actual Stage H function.
        voc=np.linalg.norm(self.h)*projection
        v=voc*self.load/(self.r+self.load)
        if self.c.get('Hrx'):
            h=self.c['Hrx']
            if not np.isclose(h['dt_s'],dt,rtol=1e-12,atol=0):raise ValueError('RECEIVER_FIR_TIME_SUPPORT')
            taps=np.asarray(h['impulse_response'],float)
            if taps.ndim!=1 or not len(taps) or not np.all(np.isfinite(taps)):raise ValueError('INVALID_FIR')
            v=np.convolve(v,taps,mode='full')[:len(v)]
        return projection,voc,v


def compact_packet(t,I,Q,config,name):
    """I is a declared equivalent conductor branch, never volumetric gas current."""
    spacing=config['cell_length_m'];width=config['cross_section_width_m']
    if min(spacing,width)<=0:raise ValueError('FINITE_SOURCE_SIZE_REQUIRED')
    kind=config['kind']
    if kind=='COMPACT_ELECTRIC_DIPOLE':
        xyz=np.array([[0.,0.,(k-2)*spacing] for k in range(5)])
        widths=np.tile([width,width,spacing],(5,1))
        edges=np.array([[0,1],[1,2],[2,3],[3,4]])
        roles=('LEAD',)*4;electrodes=[4,0]
        q=np.zeros((len(t),5));q[:,0]=-Q;q[:,4]=Q
        returns=dict(type='FINITE_CAPACITIVE_ENDPOINT',node_id=0)
    elif kind=='COMPACT_CLOSED_LOOP':
        if width!=spacing:raise ValueError('FIRST_LOOP_REQUIRES_EQUAL_VOXEL_WIDTH')
        xyz=np.array([[0.,0.,0.],[spacing,0.,0.],[spacing,0.,spacing],[0.,0.,spacing]])
        xyz-=np.mean(xyz,axis=0);widths=np.full((4,3),spacing)
        edges=np.array([[0,1],[1,2],[2,3],[3,0]]);roles=('LEAD','RETURN','RETURN','RETURN')
        electrodes=[1,0];q=np.zeros((len(t),4))
        returns=dict(type='CLOSED_CONDUCTION_RETURN',edge_ids=[1,2,3])
    else:raise ValueError('UNSUPPORTED_COMPACT_SOURCE')
    if config['axis'] not in (0,1,2):raise ValueError('FIRST_SOURCE_CARTESIAN_AXIS_ONLY')
    if config['axis']!=2:
        order=[0,1,2];order[config['axis']],order[2]=order[2],order[config['axis']]
        xyz=xyz[:,order];widths=widths[:,order]
    graph=ConservativeGraph(t,xyz,widths,edges,q,np.tile(I[:,None],(1,len(edges))),
        np.tile(np.diff(Q)[:,None],(1,len(edges))),
        tuple('H' if j==electrodes[0] else 'RETURN' if j==electrodes[1] else 'LEAD_VOLUME' for j in range(len(xyz))),roles)
    m=identity(name);m['config_hash']=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    m['source_hash']=sha(Path(__file__))
    m['build_hash']='PYTHON_EQUIVALENT_SOURCE_ADAPTER:'+m['source_hash']
    m['time_interval']=[float(t[0]),float(t[-1])]
    m['equivalent_source_scope']='CONDITIONAL_CIRCUIT_TO_COMPACT_CONDUCTOR_MAPPING_NOT_ACTUAL_DEVICE'
    m['current_origin']='Declared equivalent external conductor branch from i_loop; not plasma volume or duplicate displacement current'
    rp=m['reference_plane_id']
    m['surface_charge']=dict(reference_plane_id=rp,mapping_id='FINITE_VOLUME_ENDCAP_OR_LOOP_NODE_V1',
        node_ids=list(range(len(xyz))),charge_C=q.tolist())
    m['external_branches']=dict(reference_plane_id=rp,edge_ids=list(range(len(edges))))
    m['external_return']=dict(reference_plane_id=rp,**returns)
    m['electrode_storage']=dict(reference_plane_id=rp,node_ids=electrodes,charge_C=q[:,electrodes].tolist())
    lead=np.zeros((len(t)-1,2))
    for k,(a,b) in enumerate(edges):
        for j,node in enumerate(electrodes):lead[:,j]+=(int(b==node)-int(a==node))*graph.edge_charge_C[:,k]
    m['electrode_exchange']=dict(reference_plane_id=rp,gas_into_metal_C=np.zeros_like(lead).tolist(),
        reason='Equivalent conductor source has no resolved gas interface; not a claim of zero physical gas current')
    m['port_charge']=dict(reference_plane_id=rp,lead_into_metal_C=lead.tolist())
    return GlobalSourcePacket(m,graph=graph)


def waveform_metrics(t,v):
    dt=float(t[1]-t[0]);peak=float(np.max(np.abs(v)))
    if peak==0:return dict(peak_abs=0.,frequency_status='NOT_EVALUABLE_ZERO_SIGNAL')
    prominence=.02*peak
    pos,_=find_peaks(v,prominence=prominence)
    neg,_=find_peaks(-v,prominence=prominence)
    cycles=max(len(pos)-1,len(neg)-1,0)
    signal=np.abs(v)>.02*peak
    active=np.flatnonzero(signal);duration=float(t[active[-1]]-t[active[0]])
    period=float(np.median(np.diff(t[pos]))) if len(pos)>=4 else None
    decay=None
    if len(pos)>=4 and np.all(v[pos]>0):
        slope=np.polyfit(t[pos]-t[pos[0]],np.log(v[pos]),1)[0]
        if slope<0:decay=float(-1/slope)
    window=np.hanning(len(v));freq=np.fft.rfftfreq(len(v),dt)
    amplitudes=2*np.abs(np.fft.rfft((v-np.mean(v))*window))/window.sum()
    amplitudes[0]*=.5
    if len(v)%2==0:amplitudes[-1]*=.5
    band=(freq>=10e6)&(freq<=80e6)
    band_peak=float(freq[np.flatnonzero(band)[np.argmax(amplitudes[band])]]) if np.any(band) else None
    return dict(peak_abs=peak,positive_peaks=len(pos),negative_peaks=len(neg),complete_same_sign_cycles=cycles,
        zero_crossings=int(np.count_nonzero(v[:-1]*v[1:]<0)),prominence=prominence,
        active_duration_s=duration,median_positive_period_s=period,log_decay_time_s=decay,
        frequency_status=('INSUFFICIENT_WINDOW_FOR_10_80_MHZ' if not np.any(band) else
                          'RESOLVABLE_CONDITIONAL_WINDOW' if cycles>=3 else 'INSUFFICIENT_CYCLES'),
        band_peak_Hz=band_peak,df_Hz=1/(len(v)*dt),nyquist_Hz=.5/dt,window='HANN_DEMEANED_NO_ZERO_PADDING',
        spectral_amplitude='one-sided coherent-gain amplitude; PSD separately per Hz',
        samples=len(v),duration_s=float(t[-1]-t[0]))


def propagate(rows,s,out):
    # Sample only accepted states; preserve actual transferred charge from Cext.
    stride=max(1,int(s['radiation']['source_stride']))
    pick=sorted(set(range(0,len(rows),stride))|{len(rows)-1})
    t=np.array([0.]+[rows[k]['time_s'] for k in pick])
    I=np.array([0.]+[rows[k]['i_A'] for k in pick])
    Q=np.array([0.]+[s['circuit']['Cext_F']*rows[k]['Vc_V'] for k in pick])
    dt=float(np.median(np.diff(t)))
    obs=np.asarray(s['receiver']['position_m'],float);delay=np.linalg.norm(obs)/C0
    # Explicit neutral zero-current prehistory, not inferred from absent samples.
    t=np.r_[-max(2*delay,10*dt),-max(delay,5*dt),-dt,t]
    I=np.r_[0.,0.,0.,I];Q=np.r_[0.,0.,0.,Q]
    packet=compact_packet(t,I,Q,s['radiation'],s['id'])
    series=I1GlobalSourceAssembler().assemble(packet,requested_scope='GLOBAL_FREE_SOURCE',purpose='TIME_DOMAIN_FIELD')
    save_global_source_bundle(packet,out/'global_source')
    graph=packet.graph;distance=np.linalg.norm(graph.xyz_m-obs,axis=1)
    first_delay=float(np.min(distance)/C0)
    end=rows[-1]['time_s'];start=max(0.,first_delay-.05*end)
    stop=first_delay+end-2*dt
    times=np.linspace(start,stop,s['radiation']['observer_samples'])
    fields=evaluate_waveform(series,[Observer(s['id']+'_rx',*obs)],times,source_manifest_id=s['id'])
    Eraw=np.array([f.E_total for f in fields]);E=Eraw.copy()
    # Known causal support; centered derivatives otherwise leak across the t=0 knot.
    before=times<first_delay
    leakage=float(np.max(np.abs(Eraw[before]))) if np.any(before) else 0.
    E[before]=0.
    projection,voc,voltage=ReceiverVoltageModel(s['receiver']).evaluate(E,float(times[1]-times[0]))
    table(out/'transmitter_field.csv',[dict(**f.to_row(),causal_support_gate=bool(before[k]),
        Ex_causal=E[k,0],Ey_causal=E[k,1],Ez_causal=E[k,2]) for k,f in enumerate(fields)])
    table(out/'receiver_voltage.csv',[dict(time_s=float(tt),E_projected_V_m=float(ep),
        Voc_rx_V=float(a),Vrx_V=float(b),load_ohm=s['receiver']['load_resistance_ohm'],
        P_load_W=float(b*b/s['receiver']['load_resistance_ohm']),
        source_retarded_center_s=float(tt-delay),valid_fraction=fields[k].retarded_time_valid_fraction)
        for k,(tt,ep,a,b) in enumerate(zip(times,projection,voc,voltage))])
    metrics=waveform_metrics(times,voltage)
    metrics.update(E_peak_V_m=float(np.max(np.linalg.norm(E,axis=1))),first_possible_delay_s=first_delay,
        prearrival_derivative_interpolation_leak_V_m=leakage,
        causal_rule='Exact known zero prehistory; centered derivative boundary leakage removed before earliest source support arrival',
        graph=graph.audit(),source_records=len(series.records),field_samples=len(fields),
        consumer='I1GlobalSourceAssembler -> SourceSeries -> evaluate_waveform -> StageH.project_field -> ReceiverVoltageModel',
        RF_scope='CONDITIONAL_COMPACT_REFERENCE_NOT_EXPERIMENTAL_RECEIVER',
        radiative_backreaction='NEGLECTED_FIRST_ORDER_POSTPROCESS_NO_RF_EFFICIENCY_CLAIM')
    table(out/'spectrum_metrics.csv',[{k:v for k,v in metrics.items() if not isinstance(v,dict)}])
    w=np.hanning(len(voltage));f=np.fft.rfftfreq(len(voltage),times[1]-times[0])
    psd=np.abs(np.fft.rfft((voltage-voltage.mean())*w))**2*(times[1]-times[0])/np.sum(w*w)
    psd[1:]*=2
    if len(voltage)%2==0:psd[-1]*=.5
    table(out/'receiver_spectrum.csv',[dict(frequency_Hz=float(a),PSD_V2_Hz=float(b)) for a,b in zip(f,psd)])
    table(out/'radiation_source_summary.csv',[dict(time_s=float(tt),I_equivalent_branch_A=float(ii),
        transferred_charge_C=float(qq),net_charge_C=float(np.sum(graph.node_charge_C[k])),
        kind=s['radiation']['kind'],source_identity='CONDITIONAL_FINITE_CONDUCTOR_NOT_GAS_J')
        for k,(tt,ii,qq) in enumerate(zip(t,I,Q))])
    dump(out/'receiver_metrics.json',metrics)
    return metrics
