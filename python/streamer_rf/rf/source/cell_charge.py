"""Accepted native per-cell mechanism ledger, independent of source residuals."""
import hashlib
import json
from pathlib import Path
import numpy as np
from .global_packet import require

VERSION='I1_NATIVE_CELL_CHARGE_LEDGER_V1'
E_CHARGE=np.longdouble('1.602176634e-19')


def audit_cell_charge(c, faces, nr):
    """Use recorded raw/pre/post densities, source RHS and real corrections.

    Extended precision is for summing already-rounded native doubles; it does
    not manufacture higher precision states. Net charge defects of reaction
    arithmetic are included once as the actual net source, not physical charge
    production. Pair creation/attachment/recombination remain separate arrays.
    """
    require(len(c)>0 and len(faces)>0,'EMPTY_CELL_LEDGER')
    ids=c['cell_id'].astype(int);require(len(set(ids))==len(ids),'DUPLICATE_CELL')
    require(all(np.all(np.isfinite(c[k])) for k in c.dtype.names),'NONFINITE_CELL_LEDGER')
    require(all(np.all(np.isfinite(faces[k])) for k in faces.dtype.names),'NONFINITE_FACE_LEDGER')
    dt=c['dt_s'].astype(np.longdouble);vol=c['volume_m3'].astype(np.longdouble)
    require(np.all(dt>0) and np.all(vol>0) and np.all(dt==dt[0]) and
            np.allclose(c['post_time_s']-c['pre_time_s'],dt,rtol=1e-12,atol=0),'TIME_LAYER_MISMATCH')
    loc={k:i for i,k in enumerate(ids)};n=len(ids)
    outward=np.zeros((n,3),np.longdouble);unsigned=np.zeros(n,np.longdouble)
    face_groups=np.zeros((n,3),np.longdouble);pairs={};seen=set()
    metal={1:np.zeros(3,np.longdouble),2:np.zeros(3,np.longdouble)}
    for f in faces:
        cell=int(f['j'])*nr+int(f['i']);neighbor=int(f['nj'])*nr+int(f['ni'])
        key=(cell,int(f['direction']));require(key not in seen and cell in loc,'FACE_CELL_MAPPING_MISMATCH');seen.add(key)
        j=loc[cell];area=np.longdouble(f['area_m2']);require(area>=0,'NEGATIVE_FACE_AREA')
        e,p,emit=(np.longdouble(f[k])*area*dt[j] for k in ('gamma_e_net_m2s','gamma_p_m2s','gamma_secondary_m2s'))
        outward[j,0]+=e;outward[j,1]+=p;unsigned[j]+=E_CHARGE*(abs(e)+abs(p))
        role=int(f['role']);group=0 if role==0 else 1 if role==1 else 2
        face_groups[j,group]+=E_CHARGE*(p-e)
        if role==0:
            require(neighbor in loc,'INCOMPLETE_CELL_FACE_COVERAGE');pairs[(cell,neighbor)]=(e,p)
        if role==1:
            conductor=int(f['conductor']);require(conductor in metal,'ELECTRODE_IDENTITY_MISMATCH')
            metal[conductor]+=np.array([e+emit,p,emit])
    pair_error=0.
    for (a,b),flow in pairs.items():
        require((b,a) in pairs,'UNPAIRED_NATIVE_FACE')
        pair_error=max(pair_error,float(E_CHARGE*max(abs(flow[k]+pairs[(b,a)][k]) for k in range(2))))
    delta=[];source=[];correction=[];species_residual=[]
    precheck=0.
    for species in ('e','p','n'):
        raw=c['raw_n'+species+'_m3'].astype(np.longdouble)
        pre=c['pre_n'+species+'_m3'].astype(np.longdouble)
        post=c['post_n'+species+'_m3'].astype(np.longdouble)
        corr=sum(c[prefix+species+'_m3'].astype(np.longdouble) for prefix in ('pre_floor_','clip_','post_floor_'))*vol
        src=c['source_'+species+'_m3s'].astype(np.longdouble)*vol*dt
        d=(post-raw)*vol
        j=('e','p','n').index(species)
        species_residual.append(d+outward[:,j]-src-corr)
        delta.append(d);source.append(src);correction.append(corr)
        precheck=max(precheck,float(np.max(abs(pre-raw-c['pre_floor_'+species+'_m3']))))
    delta=np.array(delta).T;source=np.array(source).T;correction=np.array(correction).T
    signs=np.array([-1,1,-1],dtype=np.longdouble)
    dq=E_CHARGE*(delta@signs);sourceq=E_CHARGE*(source@signs);corrq=E_CHARGE*(correction@signs)
    netout=E_CHARGE*(outward@signs);residual=dq+netout-sourceq-corrq
    scale=unsigned+E_CHARGE*np.sum(abs(source)+abs(correction),axis=1)
    ratio=np.divide(abs(residual),scale,out=np.full(n,np.nan,np.longdouble),where=scale!=0)
    raw_defect=E_CHARGE*c['reaction_defect_m3s'].astype(np.longdouble)*vol*dt
    # Chemistry sums are audited independently from the stored actual RHS.
    impact=c['impact_m3s'].astype(np.longdouble);photo=c['photo_m3s'].astype(np.longdouble)
    attach=c['attach2_m3s'].astype(np.longdouble)+c['attach3_m3s'].astype(np.longdouble)
    ep=c['ep_recomb_m3s'].astype(np.longdouble);pn=c['pn_recomb_m3s'].astype(np.longdouble)
    expected=np.array([impact+photo-attach-ep,impact+photo-ep-pn,attach-pn]).T*vol[:,None]*dt[:,None]
    data=dict(cell_id=ids,delta_gas_C=dq,physical_net_source_C=np.zeros(n),
        actual_reaction_defect_C=sourceq,recorded_native_defect_C=raw_defect,numerical_correction_C=corrq,
        electron_out_C=-E_CHARGE*outward[:,0],ion_out_C=E_CHARGE*outward[:,1],
        internal_out_C=face_groups[:,0],metal_out_C=face_groups[:,1],outer_out_C=face_groups[:,2],
        residual_C=residual,unsigned_exchange_C=scale,residual_relative=ratio)
    for j,s in enumerate(('e','p','n')):data['species_'+s+'_residual_N']=species_residual[j]
    total=float(np.sum(scale));signed=float(np.sum(residual));absolute=float(np.sum(abs(residual)))
    stats=dict(cell_count=n,face_count=len(faces),max_cell_residual_C=float(max(abs(residual))),
        max_cell_exchange_relative=float(np.nanmax(ratio)) if np.any(scale>0) else None,
        zero_exchange_cells=int(np.sum(scale==0)),zero_exchange_nonzero_residual=int(np.sum((scale==0)&(residual!=0))),
        signed_residual_C=signed,L1_residual_C=absolute,unsigned_exchange_C=total,
        signed_relative=signed/total if total else None,L1_relative=absolute/total if total else None,
        internal_pair_max_error_C=pair_error,pre_floor_identity_max_m3=precheck,
        reaction_rhs_roundoff_max_C=float(E_CHARGE*np.max(abs(expected-source))),
        native_defect_arithmetic_difference_C=float(np.max(abs(sourceq-raw_defect))),
        species_residual_max_N=[float(np.max(abs(a))) for a in species_residual],
        electrode_exchange={str(k):dict(electron_absorbed_N=float(v[0]),ion_absorbed_N=float(v[1]),
                            secondary_emitted_N=float(v[2]),gas_into_metal_C=float(E_CHARGE*(v[1]-v[0]+v[2])))
                            for k,v in metal.items()})
    return data,stats


def load_cell_ledger(index_path):
    path=Path(index_path);m=json.loads(path.read_text())
    require(m.get('version')==VERSION and m.get('time_layer')=='RAW_PRE_AND_ACCEPTED_POST_WITH_PRE_RHS',
            'CELL_LEDGER_VERSION_OR_TIME_MISMATCH')
    require(m.get('coverage')=='ALL_GAS_CELLS_ALL_RECORDED_FACES','INCOMPLETE_CELL_LEDGER_FOR_GLOBAL_SOURCE')
    require(m.get('secondary_in_net_face') is True and m.get('negative_ion_transport')=='IMMOBILE_BY_MODEL',
            'SPECIES_CURRENT_IDENTITY_MISMATCH')
    files=m.get('files',{});require(set(files)=={'cells','faces'},'CELL_LEDGER_COMPLETION_MISSING')
    arrays=[]
    for key in ('cells','faces'):
        entry=files[key];p=path.parent/entry['file']
        require(Path(entry['file']).name==entry['file'] and p.is_file() and
                hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'],'CELL_LEDGER_HASH_MISMATCH')
        arrays.append(np.genfromtxt(p,delimiter=',',names=True,dtype=float,encoding='utf-8'))
    c,f=arrays
    require(len(c)==m['gas_cells'] and len(f)==m['faces'],'CELL_LEDGER_COVERAGE_MISMATCH')
    require(np.all(c['pre_time_s']==m['pre_time_s']) and np.all(c['post_time_s']==m['post_time_s']),
            'TIME_LAYER_MISMATCH')
    return m,c,f


def write_selected_cells(path, rows, selected_ids, full_count):
    """Optional sparse *diagnostic* export; never silently a complete source."""
    path=Path(path);require(not path.exists(),'REFUSE_OUTPUT_OVERWRITE')
    selected={int(k) for k in selected_ids};mask=np.array([int(k) in selected for k in rows['cell_id']])
    require(int(mask.sum())==len(selected),'UNKNOWN_SELECTED_CELL')
    np.savetxt(path,rows[mask].tolist(),delimiter=',',header=','.join(rows.dtype.names),comments='')
    return dict(coverage='EXPLICIT_CELL_SELECTION_DIAGNOSTIC_ONLY',included_cell_ids=sorted(selected),
                full_cell_count=full_count,omitted_count=full_count-len(selected),
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),global_continuity_qualified=False)
