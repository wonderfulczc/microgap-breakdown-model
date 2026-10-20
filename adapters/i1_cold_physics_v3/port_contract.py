"""Reuse FP3 accounting ownership; this is not a calibrated physical port."""
from adapters.i1_fp3_port_circuit.ports import (
    PortDescriptor,CapacitanceOwner,ChargeScope,EnergyScope,PortError,
    require_capacitance_assembly,
)

def descriptor(scenario, operator_id, owner=CapacitanceOwner.TOTAL):
    total=owner is CapacitanceOwner.TOTAL
    return PortDescriptor(
        provider_id='NATIVE_FULL_COLD_V3',geometry_id=scenario['geometry_id'],
        operator_id=operator_id,reference_plane_id='LOCAL_HV_RETURN',
        capacitance_owner=owner,
        charge_scope=ChargeScope.TOTAL if total else ChargeScope.REMAINING,
        energy_scope=EnergyScope.TOTAL if total else EnergyScope.REMAINING,
        model_identity=('ELECTRON_TRANSPORT_ISG0','POSITIVE_ION_TRANSPORT',
                        'CATHODE_ION_SECONDARY_EMISSION','GEOMETRY_AWARE_SP3','LOCAL_CHARGE_PORT_V2'),
        schema_version='I1-COLD-PHYSICS-V3',
        capabilities=('provide_local_numerical_port_charge','provide_enabled_species_face_current',
                      'native_serial_transactions','internal_GR','electrical_field_work'),
        validity_domain='LOCAL_NUMERICAL_TEST_SERIAL_FIXED_GEOMETRY_PRESCRIBED_VOLTAGE',
        source_identity='SYNTHETIC_OR_CONDITIONAL_PHYSICS_NOT_DEVICE_CALIBRATION')

def require_operation(desc, capability, *, operation_scope):
    if operation_scope!='LOCAL_NUMERICAL_TEST' or capability not in desc.capabilities:
        raise PortError('UNAVAILABLE_CAPABILITY: no global port, gas heat, automatic thermal or RF')

__all__=['descriptor','require_operation','require_capacitance_assembly']
