"""Native B0 surface mapping and explicit external conductor geometry contract."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .global_packet import require, finite

VERSION='I1_CONDUCTOR_GEOMETRY_V1'
CLASSES=('ACTUAL_SOURCE_BOUND_GEOMETRY','PARAMETRIC_B0_GEOMETRY',
         'SYNTHETIC_TEST_GEOMETRY','UNKNOWN_EXPERIMENTAL_GEOMETRY')


class I1ConductorGeometryAdapter:
    def from_native(self, directory, geometry_source, profile_contract=None):
        directory=Path(directory); source=Path(geometry_source)
        meta=json.loads((directory/'geometry_native.json').read_text())
        profile_evidence=None
        if profile_contract is not None:
            profile_contract=Path(profile_contract);reference=json.loads(profile_contract.read_text());p=reference['parameters']
            require(reference['geometry_id']==meta['geometry_id'],'PROFILE_GEOMETRY_IDENTITY_MISMATCH')
            for k,profile in enumerate(meta['profiles']):
                for native_key,source_key in (('tip_radius_m','tip_radius_m'),('cone_half_angle_rad','half_angle_rad'),
                        ('shank_radius_m','shank_radius_m'),('total_length_m','tip_to_back_length_m')):
                    require(np.isclose(profile[native_key],p[source_key],rtol=1e-14,atol=0),
                            'PROFILE_PARAMETER_MISMATCH',native_key)
                require(profile['extension_direction']==(-1 if k==0 else 1),'PROFILE_ORIENTATION_MISMATCH')
            require(np.isclose(meta['profiles'][1]['tip_z_m']-meta['profiles'][0]['tip_z_m'],p['gap_m'],rtol=1e-14,atol=0),
                    'PROFILE_GAP_MISMATCH')
            profile_evidence=dict(file=str(profile_contract),sha256=hashlib.sha256(profile_contract.read_bytes()).hexdigest(),
                parameters=p,qualification=reference['source_type'],original_CAD=reference['original_editable_CAD_status'],
                global_terminal_selection=reference['physical_final_terminal_selection'])
        mask=np.fromfile(directory/'geometry_mask.bin',dtype=np.uint8).reshape(meta['nz'],meta['nr'])
        faces=list(csv.DictReader((directory/'geometry_metal_faces.csv').open()))
        surfaces={'H':[],'Return':[]}; seen=set()
        for row in faces:
            i,j,ni,nj,d,conductor=(int(row[k]) for k in ('i','j','ni','nj','direction','conductor'))
            require(0<=i<meta['nr'] and 0<=j<meta['nz'] and 0<=ni<meta['nr'] and 0<=nj<meta['nz'],
                    'METAL_CONTACT_OUTSIDE_MASK')
            require(mask[j,i]==0 and mask[nj,ni]==conductor and conductor in (1,2), 'METAL_CONTACT_MISMATCH')
            nr,nz=(( -1,0),(1,0),(0,-1),(0,1))[d]
            require((ni-i,nj-j)==(nr,nz) and float(row['normal_r'])==nr
                    and float(row['normal_z'])==nz, 'SURFACE_ORIENTATION_MISMATCH')
            key=f'{j*meta["nr"]+i}:{d}';require(key not in seen,'DUPLICATE_SURFACE');seen.add(key)
            r=(i+.5)*meta['dr_m']+.5*nr*meta['dr_m']
            z=meta['z0_center_m']+j*meta['dz_m']+.5*nz*meta['dz_m']
            area=2*np.pi*r*meta['dz_m'] if nr else np.pi*((i+1)**2-i**2)*meta['dr_m']**2
            require(np.isclose(float(row['area_m2']),area,rtol=1e-12,atol=0)
                    and abs(float(row['r_m'])-r)<=1e-15 and abs(float(row['z_m'])-z)<=1e-15,
                    'SURFACE_GEOMETRY_MISMATCH')
            surfaces['H' if conductor==1 else 'Return'].append(dict(face_id=key,gas_cell=j*meta['nr']+i,
                metal_cell=nj*meta['nr']+ni,r_m=r,z_m=z,area_m2=float(row['area_m2']),
                normal_gas_to_metal_rz=[nr,nz]))
        require(all(surfaces.values()),'ELECTRODE_IDENTITY_MISSING')
        # Local stair-step faces, not the missing whole electrode CAD or contacts.
        return dict(version=VERSION,geometry_source_id=meta['geometry_id'],
            geometry_source_file=str(source),geometry_source_sha=hashlib.sha256(source.read_bytes()).hexdigest(),
            source_class='PARAMETRIC_B0_GEOMETRY',coordinate_system='AXISYMMETRIC_R_Z_SI',
            coordinate_transform=np.eye(4).tolist(),profiles=meta['profiles'],profile_contract=profile_evidence,
            electrode_H_surfaces=surfaces['H'],electrode_Return_surfaces=surfaces['Return'],
            reference_plane_H=meta['reference_plane_id']+':H',reference_plane_Return=meta['reference_plane_id']+':Return',
            reference_plane_id=meta['reference_plane_id'],lead_segments=None,return_segments=None,material_id=None,
            surface_orientation='GAS_TO_METAL',mesh_or_surface_mapping=dict(method='NATIVE_GEOMETRY2DVIEW_FACES',
                nr=meta['nr'],nz=meta['nz'],files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (directory/'geometry_native.json',directory/'geometry_mask.bin',directory/'geometry_metal_faces.csv')}),
            external_boundary='LOCAL_ZERO_NEUMANN_TRUNCATION',terminal_connection_graph=None,
            status='B0_LOCAL_CONDUCTOR_GEOMETRY_READY',unavailable=dict(
                lead_segments='UNKNOWN_EXPERIMENTAL_GEOMETRY',return_segments='UNKNOWN_EXPERIMENTAL_GEOMETRY',
                material_id='ACTUAL_MATERIAL_UNKNOWN',external_operator='GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE'))

    def validate(self, c, *, complete=False):
        require(c.get('version')==VERSION and c.get('source_class') in CLASSES,'GEOMETRY_SOURCE_CLASS')
        for k in ('geometry_source_id','geometry_source_sha','reference_plane_H','reference_plane_Return',
                  'reference_plane_id','coordinate_system','surface_orientation','external_boundary','mesh_or_surface_mapping'):
            require(c.get(k),'GEOMETRY_IDENTITY_MISSING',k)
        source=Path(c.get('geometry_source_file',''))
        require(source.is_file() and hashlib.sha256(source.read_bytes()).hexdigest()==c['geometry_source_sha'],
                'GEOMETRY_SOURCE_HASH_MISMATCH')
        if c['source_class']=='ACTUAL_SOURCE_BOUND_GEOMETRY':
            try: payload=json.loads(source.read_text())
            except (ValueError,OSError): require(False,'GEOMETRY_PAYLOAD_UNREADABLE')
            for key in ('coordinate_system','coordinate_transform','electrode_H_surfaces','electrode_Return_surfaces',
                        'reference_plane_id','lead_segments','return_segments','terminal_connection_graph',
                        'material_id','mesh_or_surface_mapping','external_boundary','surface_orientation'):
                require(payload.get(key)==c.get(key),'GEOMETRY_PAYLOAD_MISMATCH',key)
        T=finite(c.get('coordinate_transform'),(4,4),'transform')
        require(np.allclose(T[3],[0,0,0,1],rtol=0,atol=1e-14)
                and np.allclose(T[:3,:3].T@T[:3,:3],np.eye(3),rtol=0,atol=1e-12)
                and np.linalg.det(T[:3,:3])>0, 'INVALID_COORDINATE_TRANSFORM')
        require(c['surface_orientation']=='GAS_TO_METAL','SURFACE_ORIENTATION_MISMATCH')
        for key in ('electrode_H_surfaces','electrode_Return_surfaces'):
            require(c.get(key),'SURFACE_GEOMETRY_MISSING',key)
            for patch in c[key]:
                require(patch.get('face_id') and np.isfinite(patch.get('area_m2',np.nan)) and patch['area_m2']>0,
                        'SURFACE_GEOMETRY_MISMATCH')
                normal=patch.get('normal_gas_to_metal',patch.get('normal_gas_to_metal_rz'))
                require(normal is not None and np.all(np.isfinite(normal)) and abs(np.linalg.norm(normal)-1)<1e-12,
                        'SURFACE_ORIENTATION_MISMATCH')
        if not complete: return 'B0_LOCAL_CONDUCTOR_GEOMETRY_READY'
        require(c.get('lead_segments') is not None,'EXTERNAL_LEAD_GEOMETRY_MISSING')
        require(c.get('return_segments') is not None,'EXTERNAL_RETURN_GEOMETRY_MISSING')
        graph=c.get('terminal_connection_graph');require(graph,'TERMINAL_CONNECTION_GRAPH_MISSING')
        nodes=finite(graph.get('xyz_m'),(len(graph.get('xyz_m',[])),3),'connection nodes')
        segments=c['lead_segments']+c['return_segments'];require(segments,'EXTERNAL_RETURN_GEOMETRY_MISSING')
        edges=[]
        for segment in segments:
            a,b=segment['node_ids'];require(isinstance(a,int) and isinstance(b,int) and a!=b
                and 0<=a<len(nodes) and 0<=b<len(nodes),'INVALID_CONNECTIVITY')
            require(segment.get('reference_plane_id')==c['reference_plane_id'],'REFERENCE_PLANE_MISMATCH')
            require(np.isfinite(segment.get('area_m2',np.nan)) and segment['area_m2']>0,'FINITE_SOURCE_SIZE_REQUIRED')
            require(np.linalg.norm(nodes[a]-nodes[b])>0,'ZERO_LENGTH_EDGE');edges.append((a,b))
        require(len(set(tuple(sorted(e)) for e in edges))==len(edges),'DUPLICATE_FACE')
        contacts=graph.get('electrode_contacts',{})
        require(set(contacts)=={'H','Return'},'ELECTRODE_CONTACT_MISSING')
        for electrode,contact in contacts.items():
            available={s['face_id'] for s in c['electrode_'+electrode+'_surfaces']}
            require(contact.get('face_id') in available and contact.get('node_id') in range(len(nodes)),
                    'ELECTRODE_CONTACT_MISMATCH')
        # A closed conduction graph must have no dangling or disconnected pieces.
        degree=np.zeros(len(nodes),int); adjacency=[set() for _ in nodes]
        for a,b in edges: degree[a]+=1;degree[b]+=1;adjacency[a].add(b);adjacency[b].add(a)
        visited=set();pending=[contacts['H']['node_id']]
        while pending:
            a=pending.pop()
            if a not in visited: visited.add(a);pending.extend(adjacency[a]-visited)
        require(np.all(degree>=2) and len(visited)==len(nodes),'RETURN_CONNECTION_NOT_CLOSED')
        require(c.get('material_id'),'CONDUCTOR_MATERIAL_MISSING')
        return 'COMPLETE_GEOMETRY_CONTRACT_READY_NOT_PHYSICAL_VALIDATION'
