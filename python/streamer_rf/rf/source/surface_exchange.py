"""Explicit surface distributions and conductor contact conservation.

Surface-to-volume mapping is supplied by its source, never inferred from Q_H.
The graph remains the authoritative finite-volume free-current representation.
"""
import hashlib
from pathlib import Path
import numpy as np
from .global_packet import require, finite
from .conductor_geometry import I1ConductorGeometryAdapter
from .external_boundary import I1ExternalBoundaryResponse


def validate_conductor_extension(packet):
    m=packet.manifest;g=packet.graph;ext=m.get('conductor_extension')
    require(ext is not None,'SURFACE_CHARGE_DISTRIBUTION_UNAVAILABLE')
    require(ext.get('version')=='I1_SURFACE_CONTACT_EXCHANGE_V1','SURFACE_CONTRACT_VERSION')
    geometry=ext.get('geometry');require(geometry,'EXTERNAL_RETURN_GEOMETRY_MISSING')
    I1ConductorGeometryAdapter().validate(geometry,complete=True)
    require(geometry['reference_plane_id']==m['reference_plane_id'],'REFERENCE_PLANE_MISMATCH')
    require(ext.get('time_layer')==m['time_layer'] and ext.get('time_interval')==m['time_interval'],
            'TIME_LAYER_MISMATCH')
    require(ext.get('units')=={'sigma_s':'C/m2','surface_current':'A/m','area':'m2','charge':'C'},
            'SURFACE_UNIT_MISMATCH')
    source=Path(ext.get('distribution_source_file',''))
    require(source.is_file() and hashlib.sha256(source.read_bytes()).hexdigest()==ext.get('distribution_source_sha256'),
            'SURFACE_DISTRIBUTION_SOURCE_MISMATCH')
    xyz=finite(geometry['terminal_connection_graph']['xyz_m'],g.xyz_m.shape,'graph coordinates')
    T=np.array(geometry['coordinate_transform'])
    require(np.allclose(xyz@T[:3,:3].T+T[:3,3],g.xyz_m,rtol=1e-12,atol=0), 'COORDINATE_MAPPING_MISMATCH')
    segments=geometry['lead_segments']+geometry['return_segments']
    pairs={tuple(s['node_ids']) for s in segments}
    require(pairs=={tuple(e) for e in g.edges},'LEAD_SPATIAL_MAPPING_MISMATCH')
    for s in geometry['return_segments']:
        idx=next(i for i,e in enumerate(g.edges) if tuple(e)==tuple(s['node_ids']))
        require(g.edge_roles[idx]=='RETURN','RETURN_SPATIAL_MAPPING_MISMATCH')
    patches=ext.get('patches',[]);require(patches,'SURFACE_CHARGE_DISTRIBUTION_UNAVAILABLE')
    nt=len(g.times_s);np_=len(patches)
    sigma=finite(ext.get('sigma_s_C_m2'),(nt,np_),'surface density')
    area=finite([p['area_m2'] for p in patches],(np_,),'patch area');require(np.all(area>0),'INVALID_SURFACE_AREA')
    nodes=[p['node_id'] for p in patches]
    require(len(set(p['patch_id'] for p in patches))==np_ and all(isinstance(n,int) and 0<=n<len(g.xyz_m) for n in nodes),
            'SURFACE_MAPPING_INVALID')
    q=np.zeros_like(g.node_charge_C)
    for k,p in enumerate(patches):
        normal=finite(p.get('normal_metal_to_gas'),(3,),'normal')
        require(abs(np.linalg.norm(normal)-1)<=1e-12,'SURFACE_ORIENTATION_MISMATCH')
        center=finite(p.get('center_xyz_m'),(3,),'surface center')
        # Version 1 supports an explicit boundary patch of its finite FV support.
        node=p['node_id'];offset=center-g.xyz_m[node];half=g.widths_m[node]/2
        axis=int(np.argmax(abs(normal)))
        require(abs(normal[axis])==1 and np.count_nonzero(normal)==1
                and abs(abs(offset[axis])-half[axis])<=1e-12*half[axis]
                and np.all(abs(offset)<=half*(1+1e-12)), 'SURFACE_SUPPORT_MAPPING_MISMATCH')
        q[:,p['node_id']]+=sigma[:,k]*area[k]
    scale=float(np.max(np.abs(g.node_charge_C)));err=float(np.max(np.abs(q-g.node_charge_C)))
    require(err==0 if scale==0 else err/scale<=1e-12,'SURFACE_CHARGE_MAPPING_MISMATCH')
    distribution=finite(ext.get('conductor_edge_current_A'),g.edge_current_A.shape,'conductor current')
    require(np.array_equal(distribution,g.edge_current_A),'CONDUCTOR_CURRENT_MAPPING_MISMATCH')
    require(ext.get('surface_current_A_m') is None and ext.get('surface_current_reason')=='VOLUME_CURRENT_DISTRIBUTION_SUPPLIED',
            'UNSUPPORTED_SURFACE_CURRENT_RECONSTRUCTION')
    if geometry['source_class']=='ACTUAL_SOURCE_BOUND_GEOMETRY':
        # Bind the supplied spatial arrays to a numerical payload, not just a citation.
        try:
            with np.load(source,allow_pickle=False) as f:
                matches=np.array_equal(f['sigma_s_C_m2'],sigma) and np.array_equal(f['conductor_edge_current_A'],distribution)
        except (ValueError,KeyError,OSError): require(False,'SURFACE_DISTRIBUTION_PAYLOAD_UNREADABLE')
        require(matches,'SURFACE_DISTRIBUTION_PAYLOAD_MISMATCH')
        op=ext.get('external_operator');require(op,'GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE')
        response=I1ExternalBoundaryResponse(op['contract'],np.array(op['A']),np.array(op['B']),np.array(op['D']),tuple(op['local_indices']))
        response.validate()
        require(response.contract['evidence_class']=='SOURCE_BOUND_OPERATOR' and
                response.contract['reference_plane_id']==m['reference_plane_id'] and
                response.contract['geometry_id']==geometry['geometry_source_id'], 'EXTERNAL_OPERATOR_IDENTITY_MISMATCH')
        # Global source qualification still needs actual endpoint charge data, not a capability flag.
        potentials=finite(ext.get('terminal_voltage_V'),(nt,len(response.D)),'terminal potentials')
        freeq=finite(ext.get('free_node_charge_C'),(nt,len(response.A)),'free charges')
        terminal=finite(m['electrode_storage']['charge_C'],(nt,2),'terminal storage')
        require(len(response.D)==2,'UNSUPPORTED_TERMINAL_MAPPING')
        for k in range(nt):
            predicted=response.evaluate(freeq[k],potentials[k])['terminal_charge_C']
            error=float(np.max(np.abs(predicted-terminal[k])));norm=float(np.max(np.abs(terminal[k])))
            require(error==0 if norm==0 else error/norm<=1e-10,'GLOBAL_TERMINAL_CHARGE_MAPPING_MISMATCH')
    else:
        require(geometry['source_class']=='SYNTHETIC_TEST_GEOMETRY','GLOBAL_SOURCE_INCOMPLETE')
    return dict(surface_mapping_error_C=err,patch_count=np_,node_count=len(g.xyz_m),
                source_scope='SYNTHETIC_COMPLETE_GLOBAL_SOURCE' if geometry['source_class']=='SYNTHETIC_TEST_GEOMETRY'
                else 'SOURCE_BOUND_GLOBAL_CONDUCTOR',physical_validation='NOT_GRANTED')
