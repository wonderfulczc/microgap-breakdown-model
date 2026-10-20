"""Existing native scenario plus evidence-checked spatial source capabilities."""
from adapters.i1_cold_physics_v3_2.scenario import SourceBoundI1ScenarioRunner
from streamer_rf.rf.source.global_packet import I1GlobalSourceAssembler, require
from streamer_rf.rf.source.conductor_geometry import I1ConductorGeometryAdapter
from streamer_rf.rf.source.surface_exchange import validate_conductor_extension

class I1ScenarioRunner(SourceBoundI1ScenarioRunner):
    def validate(self,data):
        super().validate(data)
        require(data.get('source_output_scope','LOCAL_GAS_SOURCE')=='LOCAL_GAS_SOURCE',
                'NATIVE_RUN_CANNOT_INVENT_GLOBAL_SOURCE')
        return data

    def source(self,scope,packet,geometry=None):
        if scope=='FULL_GLOBAL_PHYSICAL_SOURCE':
            require(False,'GLOBAL_PHYSICAL_PORT_NOT_READY')
        if scope in ('LOCAL_GAS_SOURCE','B0_LOCAL_CONDUCTOR_SOURCE'):
            record=I1GlobalSourceAssembler().assemble(packet,requested_scope='LOCAL_GAS_ONLY')
            if scope=='B0_LOCAL_CONDUCTOR_SOURCE':
                require(geometry is not None,'LOCAL_CONDUCTOR_GEOMETRY_MISSING')
                I1ConductorGeometryAdapter().validate(geometry)
                require(geometry['geometry_source_id']==packet.manifest['geometry_id'] and
                        geometry['reference_plane_id']==packet.manifest['reference_plane_id'],
                        'REFERENCE_PLANE_MISMATCH')
            return record
        require(scope in ('SYNTHETIC_COMPLETE_GLOBAL_SOURCE','SOURCE_BOUND_GLOBAL_CONDUCTOR'),
                'UNSUPPORTED_SOURCE_SCOPE')
        require(validate_conductor_extension(packet)['source_scope']==scope,'SOURCE_CLASS_MISMATCH')
        return I1GlobalSourceAssembler().assemble(packet,requested_scope='GLOBAL_FREE_SOURCE')
