"""Physical assembly inputs extend the V3.4 geometry, not its physical claims."""
import json
from pathlib import Path
import numpy as np
from streamer_rf.rf.source.conductor_geometry import I1ConductorGeometryAdapter
from adapters.i1_cold_physics_v3.scenario import CapabilityRefusal
from .artifacts import ROOT, sha

VERSION = 'I1_PHYSICAL_ASSEMBLY_1.0'
CLASSES = ('PARAMETRIC_B0_LOCAL', 'SOURCE_BOUND_MEASURED_ASSEMBLY',
           'SYNTHETIC_ENGINEERING_ASSEMBLY')
SECTIONS = ('external_geometry', 'circuit_topology', 'component_locations',
            'third_conductors_and_ground', 'measurement_planes_and_loads',
            'electrostatic_response', 'dynamic_surface_mapping', 'device_calibration')


def require(condition, code, reason=''):
    if not condition:
        raise CapabilityRefusal(code, reason)


def bound_json(ref):
    require(isinstance(ref, dict) and set(ref) >= {'file', 'sha256'}, 'SOURCE_REFERENCE_REQUIRED')
    p = (ROOT/ref['file']).resolve()
    require(p.is_relative_to(ROOT) and p.is_file() and sha(p) == ref['sha256'],
            'SOURCE_IDENTITY_MISMATCH', str(p))
    return json.loads(p.read_text())


def unavailable(reason):
    return dict(value=None, reason=reason, source=None, uncertainty=None)


class I1PhysicalAssemblyPackage:
    """Validate geometry/connection evidence separately from field and source data."""
    def __init__(self, data):
        self.data = data

    def validate(self, *, require_global_geometry=False):
        d = self.data
        require(d.get('version') == VERSION and d.get('identity') in CLASSES, 'ASSEMBLY_IDENTITY')
        require(d.get('conductor_identities') == ['H', 'Return'], 'PERMANENT_CONDUCTOR_IDENTITY')
        require(d.get('units') == {'length': 'm', 'coordinates': 'm', 'angle': 'rad'}, 'ASSEMBLY_UNITS')
        geometry = bound_json(d.get('local_geometry'))
        I1ConductorGeometryAdapter().validate(geometry)
        require(d.get('reference_plane_id') == geometry['reference_plane_id'], 'REFERENCE_PLANE_MISMATCH')
        mapping = {'PARAMETRIC_B0_LOCAL': 'PARAMETRIC_B0_GEOMETRY',
                   'SOURCE_BOUND_MEASURED_ASSEMBLY': 'ACTUAL_SOURCE_BOUND_GEOMETRY',
                   'SYNTHETIC_ENGINEERING_ASSEMBLY': 'SYNTHETIC_TEST_GEOMETRY'}
        require(geometry['source_class'] == mapping[d['identity']], 'ASSEMBLY_SOURCE_CLASS_MISMATCH')
        gap = d.get('nominal_gap')
        require(isinstance(gap, dict) and gap.get('unit') == 'm' and
                isinstance(gap.get('value'), (int, float)) and np.isfinite(gap['value'])
                and gap['value'] > 0 and gap.get('provenance'), 'GAP_PROVENANCE_REQUIRED')
        if d['identity'] == 'PARAMETRIC_B0_LOCAL':
            reference = geometry.get('profile_contract')
            require(reference and gap['value'] == reference['parameters']['gap_m'], 'B0_GAP_IDENTITY_MISMATCH')
        for section in SECTIONS:
            entry = d.get(section)
            require(isinstance(entry, dict) and set(entry) >= {'value', 'reason', 'source', 'uncertainty'},
                    'ASSEMBLY_FIELD_MISSING', section)
            if entry['value'] is None:
                require(bool(entry['reason']), 'MISSING_VALUE_REASON', section)
            else:
                payload = bound_json(entry['source'])
                require(payload == entry['value'], 'ASSEMBLY_PAYLOAD_MISMATCH', section)
        topology = d['circuit_topology']['value']
        if topology is not None:
            self._topology(topology, d['reference_plane_id'])
        if require_global_geometry:
            require(d['external_geometry']['value'] is not None, 'EXTERNAL_RETURN_GEOMETRY_MISSING')
            external = d['external_geometry']['value']
            I1ConductorGeometryAdapter().validate(external, complete=True)
            require(external['reference_plane_id'] == d['reference_plane_id'], 'REFERENCE_PLANE_MISMATCH')
            require(external['geometry_source_id'] == geometry['geometry_source_id'], 'GEOMETRY_IDENTITY_MISMATCH')
            require(topology is not None, 'TERMINAL_CONNECTION_GRAPH_MISSING')
        return dict(status='ASSEMBLY_INPUT_VALIDATED_NOT_PHYSICAL_CERTIFICATION',
                    identity=d['identity'], local_geometry=geometry['status'] if 'status' in geometry else 'GEOMETRY_CONTRACT_READY',
                    real_global_geometry='WAITING_DATA' if d['identity'] != 'SOURCE_BOUND_MEASURED_ASSEMBLY'
                    else 'SOURCE_BOUND_INPUT_NOT_CALIBRATION',
                    missing=[s for s in SECTIONS if d[s]['value'] is None])

    @staticmethod
    def _topology(t, ref):
        require(t.get('reference_plane_id') == ref, 'REFERENCE_PLANE_MISMATCH')
        require(t.get('topology_id') == 'SINGLE_SERIES_SOURCE_RS_L_CEXT_PORT', 'UNSUPPORTED_TOPOLOGY')
        nodes = t.get('nodes', [])
        edges = t.get('branches', [])
        require(len(nodes) == len(set(nodes)) == 5 and len(edges) == 5, 'SERIES_CONNECTION_REQUIRED')
        require({e.get('role') for e in edges} == {'SOURCE', 'RS', 'L', 'CEXT', 'PORT'}, 'SERIES_COMPONENT_IDENTITY')
        degree = {n: 0 for n in nodes}
        adjacency = {n: set() for n in nodes}
        seen = set()
        for e in edges:
            a, b = e.get('nodes', (None, None))
            require(a in degree and b in degree and a != b, 'INVALID_CONNECTION')
            pair = tuple(sorted((a, b)))
            require(pair not in seen, 'DUPLICATE_BRANCH')
            seen.add(pair)
            degree[a] += 1; degree[b] += 1
            adjacency[a].add(b); adjacency[b].add(a)
        visited, queue = set(), [nodes[0]]
        while queue:
            node = queue.pop()
            if node not in visited:
                visited.add(node); queue.extend(adjacency[node]-visited)
        require(set(degree.values()) == {2} and len(visited) == 5, 'RETURN_CONNECTION_NOT_CLOSED')
        require(t.get('return_is_voltage_gauge_only') is True, 'UNDECLARED_GROUND_BRANCH')
        require(t.get('terminal_nodes') == {'H': 'H', 'Return': 'Return'}, 'PERMANENT_CONDUCTOR_IDENTITY')
        port = next(e for e in edges if e['role'] == 'PORT')
        require(port['nodes'] == ['H', 'Return'], 'PORT_POLARITY_MISMATCH')

    def require_physical_source(self):
        self.validate(require_global_geometry=True)
        require(self.data['identity'] == 'SOURCE_BOUND_MEASURED_ASSEMBLY', 'REAL_DEVICE_PARAMETERS_MISSING')
        require(self.data['electrostatic_response']['value'] is not None, 'GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE')
        require(self.data['dynamic_surface_mapping']['value'] is not None, 'SURFACE_CHARGE_DISTRIBUTION_UNAVAILABLE')
        # Even complete input metadata cannot confer a physical validation Gate.
        raise CapabilityRefusal('GLOBAL_PHYSICAL_PORT_NOT_READY', 'requires existing assembler conservation and device validation')
