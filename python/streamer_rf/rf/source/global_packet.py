"""Conservative free-source handoff. No field, circuit, or discharge solver.

The graph is a finite-volume charge/current representation, not a netlist
solver. Endpoint charge and every return edge must be supplied by the caller.
The existing Jefimenko consumer receives a versioned volume reconstruction.
"""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np

from .schema import SourceMetadata, SourceRecord, SourceSeries

VERSION = 'I1-GLOBAL-SOURCE-1.0'
UNITS = {'length': 'm', 'time': 's', 'charge': 'C', 'current': 'A',
         'J': 'A/m2', 'rho': 'C/m3', 'volume': 'm3'}
IDENTITIES = ('geometry_id', 'mesh_id', 'operator_id', 'boundary_id',
              'reference_plane_id', 'accepted_state_id', 'species_model_id',
              'parameter_source_class', 'source_hash', 'build_hash', 'config_hash')


class SourceRefusal(ValueError):
    def __init__(self, code, details=None):
        self.code, self.details = code, details
        super().__init__(code + (': ' + str(details) if details is not None else ''))


def require(condition, code, details=None):
    if not condition:
        raise SourceRefusal(code, details)


def finite(value, shape, name):
    a = np.asarray(value, dtype=float)
    require(a.shape == shape and np.all(np.isfinite(a)), 'INVALID_ARRAY', name)
    return a


@dataclass(frozen=True)
class ConservativeGraph:
    """Outward incidence is + at edge tail, - at edge head.

    node_charge_C is independently supplied at samples; edge_charge_C is
    the actual interval integral, not endpoint current times dt. Zero-radius
    point/line sources are forbidden. Each node has a positive finite volume.
    """
    times_s: np.ndarray
    xyz_m: np.ndarray
    widths_m: np.ndarray
    edges: np.ndarray
    node_charge_C: np.ndarray
    edge_current_A: np.ndarray
    edge_charge_C: np.ndarray
    node_roles: tuple
    edge_roles: tuple

    def audit(self):
        t = np.asarray(self.times_s, dtype=float)
        require(t.ndim == 1 and len(t) >= 3 and np.all(np.isfinite(t))
                and np.all(np.diff(t) > 0), 'TIME_LAYER_MISMATCH')
        n, nt = len(self.xyz_m), len(t)
        x = finite(self.xyz_m, (n,3), 'xyz')
        widths = finite(self.widths_m, (n,3), 'widths')
        require(np.all(widths > 0), 'FINITE_SOURCE_SIZE_REQUIRED')
        edges = np.asarray(self.edges)
        require(edges.ndim == 2 and edges.shape[1] == 2
                and np.issubdtype(edges.dtype, np.integer), 'INVALID_CONNECTIVITY')
        m = len(edges)
        require(np.all((edges >= 0) & (edges < n))
                and np.all(edges[:,0] != edges[:,1]), 'INVALID_CONNECTIVITY')
        require(len(self.node_roles) == n and len(self.edge_roles) == m,
                'SPATIAL_SUPPORT_IDENTITY_MISSING')
        require(len({tuple(sorted(e)) for e in edges}) == m, 'DUPLICATE_FACE')
        length = np.linalg.norm(x[edges[:,1]]-x[edges[:,0]], axis=1)
        require(np.all(length > 0), 'ZERO_LENGTH_EDGE')
        require(len(np.unique(x,axis=0)) == n, 'DUPLICATE_SPATIAL_NODE')
        # This version reconstructs adjoining Cartesian finite-volume cells,
        # not arbitrary zero-radius wire segments or disconnected volume lumps.
        for a,b in edges:
            delta = x[b]-x[a]
            axes = np.flatnonzero(delta != 0)
            require(len(axes) == 1, 'FV_FACE_GEOMETRY_MISMATCH')
            axis = axes[0]
            transverse = [k for k in range(3) if k != axis]
            require(np.isclose(abs(delta[axis]),(widths[a,axis]+widths[b,axis])/2,
                               rtol=1e-12,atol=0)
                    and np.array_equal(widths[a,transverse],widths[b,transverse]),
                    'FV_FACE_GEOMETRY_MISMATCH')
        q = finite(self.node_charge_C, (nt,n), 'node_charge')
        current = finite(self.edge_current_A, (nt,m), 'edge_current')
        transfer = finite(self.edge_charge_C, (nt-1,m), 'edge_charge')
        divergence = np.zeros((nt-1,n))
        for k,(a,b) in enumerate(edges):
            divergence[:,a] += transfer[:,k]
            divergence[:,b] -= transfer[:,k]
        residual = np.diff(q, axis=0) + divergence
        scale = max(float(np.max(np.abs(q))), float(np.max(np.abs(transfer))))
        absolute = float(np.max(np.abs(residual)))
        relative = absolute / scale if scale else None
        require(absolute == 0 if scale == 0 else relative <= 1e-10,
                'FREE_SOURCE_CONTINUITY_VIOLATION', relative)
        # Reject sample currents unrelated to the supplied charge integral.
        # This is a quadrature consistency check, not a continuity correction.
        trapezoid = .5*(current[:-1]+current[1:])*np.diff(t)[:,None]
        integral_scale = float(np.max(np.abs(transfer)))
        integral_error = float(np.max(np.abs(trapezoid-transfer)))
        require(integral_error == 0 if integral_scale == 0
                else integral_error/integral_scale <= 1e-3,
                'CURRENT_INTEGRAL_TIME_MISMATCH')
        return dict(max_cell_residual_C=absolute, charge_scale_C=scale,
                    relative_residual=relative,
                    max_global_charge_change_C=float(np.max(np.abs(q.sum(axis=1)-q[0].sum()))),
                    current_quadrature_error_C=integral_error,
                    reconstruction_id='GRAPH_HALF_EDGE_VOLUME_V1',
                    continuum_limit='FINITE_VOLUME_QUADRATURE_NOT_POINT_SOURCE')

    def to_series(self, identity):
        self.audit()
        x, edges = np.asarray(self.xyz_m), np.asarray(self.edges)
        widths = np.asarray(self.widths_m)
        volume = widths.prod(axis=1)
        current = np.asarray(self.edge_current_A)
        nt, n = len(self.times_s), len(x)
        moment_density = np.zeros((nt,n,3))
        for k,(a,b) in enumerate(edges):
            contribution = current[:,k,None]*(x[b]-x[a])[None,:]/2
            moment_density[:,a] += contribution
            moment_density[:,b] += contribution
        records = []
        for k,t in enumerate(self.times_s):
            meta = SourceMetadata(
                case_id=identity['geometry_id'], solver='I1GlobalSourceAssembler',
                solver_version=VERSION, time_s=float(t), coordinate_system='cartesian_3d',
                pressure_Pa=0, temperature_K=0, geometry_id=identity['geometry_id'],
                voltage_state='NOT_APPLICABLE_MANUFACTURED_FREE_SOURCE',
                photoionization='NOT_APPLICABLE', source_definition='FREE_CHARGE_AND_CONDUCTION',
                units=UNITS, extra={**identity, 'thermodynamic_fields':'NOT_APPLICABLE',
                                   'reconstruction':'GRAPH_HALF_EDGE_VOLUME_V1'})
            J = moment_density[k]/volume[:,None]
            records.append(SourceRecord(meta, dict(level=np.zeros(n,dtype=int),
                x_center=x[:,0], y_center=x[:,1], z_center=x[:,2],
                dx=widths[:,0], dy=widths[:,1], dz=widths[:,2], cell_volume=volume,
                rho=self.node_charge_C[k]/volume, Jx=J[:,0], Jy=J[:,1], Jz=J[:,2])))
        return SourceSeries(tuple(records))


@dataclass(frozen=True)
class GlobalSourcePacket:
    manifest: dict
    graph: ConservativeGraph | None = None
    local_record: SourceRecord | None = None


class I1GlobalSourceAssembler:
    """Explicit local diagnostic or complete, charge-conserving free source."""
    def assemble(self, packet, *, requested_scope='LOCAL_GAS_ONLY', purpose='CURRENT_MOMENT'):
        m = packet.manifest
        require(m.get('version') in (VERSION,'I1-GLOBAL-SOURCE-1.1'), 'SOURCE_VERSION_MISMATCH')
        require(all(m.get(k) for k in IDENTITIES), 'SOURCE_IDENTITY_MISSING')
        require(m.get('parameter_source_class') in ('A','B','C','D'), 'PARAMETER_SOURCE_CLASS_INVALID')
        require(m.get('units') == UNITS, 'SOURCE_UNIT_MISMATCH')
        require(m.get('includes_displacement_current') is False, 'DISPLACEMENT_DOUBLE_COUNT')
        require(m.get('secondary_emission_in_net_electron_flux') is True
                and m.get('secondary_emission_added_as_source') is False, 'EMISSION_DOUBLE_COUNT')
        require(m.get('negative_ion_transport') in ('IMMOBILE_BY_MODEL','NOT_APPLICABLE'),
                'UNSUPPORTED_NEGATIVE_ION_SOURCE')
        require(m.get('capacitance_owner') in ('TOTAL_CHARGE','SPLIT_GEOMETRIC_CAPACITANCE'),
                'CAPACITANCE_OWNERSHIP_MISMATCH')
        require(m.get('capacitance_copies') == 1, 'CAPACITANCE_DOUBLE_COUNT')
        interval = m.get('time_interval', [])
        require(len(interval) == 2 and np.all(np.isfinite(interval))
                and interval[0] < interval[1], 'TIME_LAYER_MISMATCH')
        require(purpose in ('CURRENT_MOMENT','TIME_DOMAIN_FIELD'), 'UNSUPPORTED_SPECTRUM_REQUEST')
        if requested_scope == 'LOCAL_GAS_ONLY':
            require(purpose == 'CURRENT_MOMENT', 'GLOBAL_SOURCE_INCOMPLETE')
            require(packet.local_record is not None, 'LOCAL_RECORD_MISSING')
            record = packet.local_record
            require(record.metadata.extra.get('reference_plane_id') == m['reference_plane_id'],
                    'REFERENCE_PLANE_MISMATCH')
            require(m.get('time_layer') == 'ACCEPTED_INTERVAL_PRE_FLUX'
                    and record.metadata.time_s == interval[0]
                    and record.metadata.extra.get('post_time_s') == interval[1], 'TIME_LAYER_MISMATCH')
            for key in ('geometry_id','mesh_id','operator_id','source_hash'):
                require(record.metadata.extra.get(key) == m[key], 'SOURCE_IDENTITY_MISMATCH', key)
            if m.get('cell_charge_ledger') is not None:
                from .cell_charge import load_cell_ledger,audit_cell_charge
                entry=m['cell_charge_ledger'];path=Path(entry['index_path'])
                require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==entry.get('index_sha256'),
                        'CELL_LEDGER_HASH_MISMATCH')
                ledger,cells,faces=load_cell_ledger(path)
                require(ledger['geometry_id']==m['geometry_id'] and ledger['reference_plane_id']==m['reference_plane_id']
                        and ledger['mesh_id']==record.metadata.extra.get('native_mesh_id')
                        and ledger['operator_id']==record.metadata.extra.get('native_operator_id')
                        and [ledger['pre_time_s'],ledger['post_time_s']]==interval
                        and ledger['accepted_step_id']==entry.get('accepted_step_id'), 'CELL_LEDGER_IDENTITY_MISMATCH')
                require(np.array_equal(cells['cell_id'].astype(int),record.columns['cell_id']),
                        'CELL_LEDGER_COVERAGE_MISMATCH')
                _,audit=audit_cell_charge(cells,faces,entry['nr'])
                require(audit['zero_exchange_nonzero_residual']==0 and audit['internal_pair_max_error_C']==0
                        and audit['max_cell_exchange_relative'] is not None
                        and audit['max_cell_exchange_relative']<=1e-6 and audit['L1_relative']<=1e-8,
                        'CELL_CHARGE_CONTINUITY_NOT_QUALIFIED',audit)
            return record
        require(requested_scope == 'GLOBAL_FREE_SOURCE', 'UNSUPPORTED_SOURCE_SCOPE')
        required = {'external_return':'EXTERNAL_RETURN_GEOMETRY_MISSING',
                    'surface_charge':'SURFACE_CHARGE_DISTRIBUTION_UNAVAILABLE',
                    'external_branches':'EXTERNAL_LEAD_GEOMETRY_MISSING',
                    'electrode_exchange':'ELECTRODE_EXCHANGE_UNAVAILABLE',
                    'electrode_storage':'ELECTRODE_STORAGE_UNAVAILABLE',
                    'port_charge':'PORT_CHARGE_UNAVAILABLE'}
        missing = [code for key,code in required.items() if m.get(key) is None]
        require(not missing, 'GLOBAL_SOURCE_INCOMPLETE', missing)
        require(packet.graph is not None, 'GLOBAL_CONSERVATIVE_GRAPH_MISSING')
        graph = packet.graph
        require(m.get('time_layer') in ('ACCEPTED_SYNTHETIC_TIME_SERIES','ACCEPTED_SOURCE_BOUND_TIME_SERIES')
                and interval == [float(graph.times_s[0]),float(graph.times_s[-1])],
                'TIME_LAYER_MISMATCH')
        if m.get('version') == 'I1-GLOBAL-SOURCE-1.1':
            from .surface_exchange import validate_conductor_extension
            audit = validate_conductor_extension(packet)
            expected = ('SYNTHETIC_ENGINEERING_REFERENCE' if audit['source_scope']=='SYNTHETIC_COMPLETE_GLOBAL_SOURCE'
                        else 'SOURCE_BOUND_GLOBAL_CONDUCTOR')
            require(m.get('evidence_class') == expected, 'REAL_GLOBAL_SOURCE_NOT_QUALIFIED')
            require(m['time_layer'] == ('ACCEPTED_SYNTHETIC_TIME_SERIES' if expected=='SYNTHETIC_ENGINEERING_REFERENCE'
                                      else 'ACCEPTED_SOURCE_BOUND_TIME_SERIES'), 'TIME_LAYER_MISMATCH')
        else:
            require(m.get('evidence_class') == 'SYNTHETIC_ENGINEERING_REFERENCE'
                    and m.get('time_layer')=='ACCEPTED_SYNTHETIC_TIME_SERIES',
                    'REAL_GLOBAL_SOURCE_NOT_QUALIFIED')
        for key in required:
            require(m[key].get('reference_plane_id') == m['reference_plane_id'],
                    'REFERENCE_PLANE_MISMATCH', key)
        # Surface storage is a finite-volume regularization with explicit mapping.
        surface = m['surface_charge']
        require(surface.get('mapping_id') == 'FINITE_VOLUME_ENDCAP_OR_LOOP_NODE_V1',
                'SURFACE_CHARGE_MAPPING_UNVERIFIED')
        nodes = surface.get('node_ids', [])
        require(nodes and all(isinstance(i,int) and 0 <= i < len(graph.xyz_m) for i in nodes),
                'SURFACE_CHARGE_DISTRIBUTION_UNAVAILABLE')
        sq = finite(surface.get('charge_C'), (len(graph.times_s),len(nodes)), 'surface_charge')
        require(np.array_equal(sq, graph.node_charge_C[:,nodes]), 'SURFACE_CHARGE_MAPPING_MISMATCH')
        leads = m['external_branches'].get('edge_ids', [])
        require(leads and all(isinstance(i,int) and 0 <= i < len(graph.edges) for i in leads),
                'EXTERNAL_LEAD_GEOMETRY_MISSING')
        returns = m['external_return']
        if returns.get('type') == 'CLOSED_CONDUCTION_RETURN':
            ids = returns.get('edge_ids', [])
            require(ids and all(i in leads and graph.edge_roles[i] == 'RETURN' for i in ids),
                    'EXTERNAL_RETURN_GEOMETRY_MISSING')
        else:
            require(returns.get('type') == 'FINITE_CAPACITIVE_ENDPOINT'
                    and returns.get('node_id') in nodes,
                    'EXTERNAL_RETURN_GEOMETRY_MISSING')
        electrodes = m['electrode_storage'].get('node_ids', [])
        require(len(electrodes) == 2 and len(set(electrodes)) == 2
                and all(i in nodes for i in electrodes), 'ELECTRODE_IDENTITY_MISMATCH')
        storage = finite(m['electrode_storage'].get('charge_C'),
                         (len(graph.times_s),2), 'electrode_storage')
        require(np.array_equal(storage, graph.node_charge_C[:,electrodes]),
                'ELECTRODE_STORAGE_MISMATCH')
        gas_in = np.zeros((len(graph.times_s)-1,2))
        lead_in = np.zeros_like(gas_in)
        for edge,(a,b) in enumerate(graph.edges):
            for j,node in enumerate(electrodes):
                sign = int(b == node)-int(a == node)
                target = gas_in if graph.edge_roles[edge] == 'GAS' else lead_in
                target[:,j] += sign*graph.edge_charge_C[:,edge]
        supplied_gas = finite(m['electrode_exchange'].get('gas_into_metal_C'),
                              gas_in.shape, 'electrode_exchange')
        supplied_port = finite(m['port_charge'].get('lead_into_metal_C'),
                               lead_in.shape, 'port_charge')
        require(np.array_equal(gas_in,supplied_gas) and np.array_equal(lead_in,supplied_port),
                'ELECTRODE_FACE_PORT_MISMATCH')
        balance = np.diff(storage,axis=0)-gas_in-lead_in
        scale = max(float(np.max(np.abs(storage))),float(np.max(np.abs(graph.edge_charge_C))))
        require(np.max(np.abs(balance)) == 0 if scale == 0
                else np.max(np.abs(balance))/scale <= 1e-10, 'ELECTRODE_CHARGE_BALANCE')
        require(m.get('outer_exchange_C') == 0 and m.get('numerical_correction_C') == 0
                and m.get('physical_net_source_C') == 0, 'UNBALANCED_EXTERNAL_FREE_SOURCE')
        graph.audit()
        return graph.to_series({k:m[k] for k in IDENTITIES})


def audit_native_faces(faces, cell_ids, pre_charge_C, post_charge_C, dt_s, nr,
                       physical_source_C=None, numerical_correction_C=None):
    """Keep unknown corrections unknown; never infer them from residuals."""
    require(np.isfinite(dt_s) and dt_s > 0, 'TIME_LAYER_MISMATCH')
    ids = np.asarray(cell_ids, dtype=int)
    index = {int(k): j for j,k in enumerate(ids)}
    n = len(ids)
    pre = finite(pre_charge_C, (n,), 'pre_charge')
    post = finite(post_charge_C, (n,), 'post_charge')
    out_e, out_p = np.zeros(n), np.zeros(n)
    face_map = {}
    boundary = {}
    qe = 1.602176634e-19
    for f in faces:
        k = int(f['j'])*nr+int(f['i'])
        neighbor = int(f['nj'])*nr+int(f['ni'])
        area = float(f['area_m2'])
        ge,gp = float(f['gamma_e_net_m2s']),float(f['gamma_p_m2s'])
        require(area >= 0 and np.isfinite(area+ge+gp), 'INVALID_FACE_DATA')
        de,dp = -qe*ge*area*dt_s, qe*gp*area*dt_s
        out_e[index[k]] += de
        out_p[index[k]] += dp
        if int(f['role']) == 0:
            face_map[(k,neighbor)] = de+dp
        else:
            label = f"role_{int(f['role'])}_conductor_{int(f['conductor'])}"
            boundary[label] = boundary.get(label,0.)+de+dp
    pairs = []
    for (a,b),v in face_map.items():
        require((b,a) in face_map, 'UNPAIRED_NATIVE_FACE')
        pairs.append(abs(v+face_map[(b,a)]))
    raw = post-pre+out_e+out_p
    qualified = physical_source_C is not None and numerical_correction_C is not None
    corrected = (raw-finite(physical_source_C,(n,),'physical_source')
                 -finite(numerical_correction_C,(n,),'numerical_correction')) if qualified else None
    return dict(cell_ids=ids, delta_gas_C=post-pre, electron_out_C=out_e, ion_out_C=out_p,
                raw_balance_C=raw, qualified_residual_C=corrected,
                source_correction_status='AVAILABLE' if qualified else 'NOT_AVAILABLE_PER_CELL_IN_ARCHIVE',
                max_internal_pair_error_C=max(pairs,default=0.), boundary_C=boundary)


def save_global_source_bundle(packet, directory):
    """Publish a new finite-source bundle; never overwrite an existing directory."""
    I1GlobalSourceAssembler().assemble(packet, requested_scope='GLOBAL_FREE_SOURCE')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    g = packet.graph
    data_path = directory/'graph.npz'
    np.savez_compressed(data_path, times_s=g.times_s, xyz_m=g.xyz_m, widths_m=g.widths_m,
                        edges=g.edges, node_charge_C=g.node_charge_C,
                        edge_current_A=g.edge_current_A, edge_charge_C=g.edge_charge_C)
    manifest = dict(packet.manifest)
    manifest['graph_payload'] = dict(file='graph.npz',sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),
                                    node_roles=list(g.node_roles),edge_roles=list(g.edge_roles))
    path = directory/'manifest.json'
    path.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    (directory/'COMPLETE').write_text(hashlib.sha256(path.read_bytes()).hexdigest()+'\n')
    return path


def load_global_source_bundle(manifest_path):
    """Verify completion and byte identity before materializing any source arrays."""
    path = Path(manifest_path)
    marker = path.parent/'COMPLETE'
    require(marker.is_file(), 'INCOMPLETE_SOURCE_BUNDLE')
    require(marker.read_text().strip() == hashlib.sha256(path.read_bytes()).hexdigest(),
            'SOURCE_MANIFEST_HASH_MISMATCH')
    m = json.loads(path.read_text()); payload = m.get('graph_payload',{})
    require(payload.get('file') == 'graph.npz', 'SOURCE_PAYLOAD_PATH_INVALID')
    data_path = path.parent/payload['file']
    require(data_path.is_file() and hashlib.sha256(data_path.read_bytes()).hexdigest() == payload.get('sha256'),
            'SOURCE_PAYLOAD_HASH_MISMATCH')
    with np.load(data_path,allow_pickle=False) as f:
        arrays = {key:f[key].copy() for key in ('times_s','xyz_m','widths_m','edges',
                  'node_charge_C','edge_current_A','edge_charge_C')}
    graph = ConservativeGraph(**arrays,node_roles=tuple(payload['node_roles']),edge_roles=tuple(payload['edge_roles']))
    packet = GlobalSourcePacket(m,graph=graph)
    I1GlobalSourceAssembler().assemble(packet,requested_scope='GLOBAL_FREE_SOURCE')
    return packet
