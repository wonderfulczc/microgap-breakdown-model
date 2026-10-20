"""Metadata for conditional compact sources, independent of local test fixtures."""
import hashlib
from pathlib import Path

from streamer_rf.rf.source.global_packet import VERSION, UNITS


def identity(name):
    source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(
        version=VERSION, geometry_id=name, mesh_id=name + ':finite_voxels',
        operator_id='ORIENTED_FV_GRAPH_V1', boundary_id='ISOLATED_FREE_SPACE',
        reference_plane_id='SYNTHETIC_FIXED_H_RETURN', accepted_state_id=name + ':series',
        species_model_id='PRESCRIBED_CONDUCTION_NO_PLASMA', parameter_source_class='C',
        source_hash=source_sha, build_hash='PYTHON_NO_NATIVE_BUILD:' + source_sha,
        config_hash=None, units=UNITS, time_interval=None,
        time_layer='ACCEPTED_SYNTHETIC_TIME_SERIES',
        evidence_class='SYNTHETIC_ENGINEERING_REFERENCE',
        includes_displacement_current=False,
        secondary_emission_in_net_electron_flux=True,
        secondary_emission_added_as_source=False,
        negative_ion_transport='NOT_APPLICABLE',
        capacitance_owner='TOTAL_CHARGE', capacitance_copies=1,
        capacitance_value_F=None,
        capacitance_reason='No lumped capacitor inferred from this equivalent free source',
        outer_exchange_C=0, numerical_correction_C=0, physical_net_source_C=0,
        physical_capabilities=['SYNTHETIC_CONSERVATIVE_FREE_SOURCE'],
        unavailable={'receiver': 'NO_ANTENNA_OR_INSTRUMENT_TRANSFER',
                     'real_geometry': 'NOT_A_DEVICE_MODEL'},
    )
