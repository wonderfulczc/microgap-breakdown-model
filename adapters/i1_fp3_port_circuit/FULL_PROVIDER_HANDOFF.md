# FULL provider handoff: protocol only, native connection not implemented

Contract: [I1-PORT-1.0](../../docs/i1_part_vi_port_circuit_energy_contract.md).
Implementation: [ports.py](ports.py), [circuit.py](circuit.py), [guards.py](guards.py).
All accepted engineering results in this directory are **SYNTHETIC** R/C tests.
`FullProviderUnavailable` rejects evaluation with `FULL_PROVIDER_NOT_IMPLEMENTED`;
lifting from G_R rejects with `STATE_LIFTING_UNAVAILABLE`. It supplies no fake zero
charge/current/heat. `PortProvider` is the future structural protocol, not a C++ binding.

## Minimum native requirements (not authorized here)

1. Bind the existing Geometry2DView, actual boundary/operator identity, mesh and
   physical reference plane to a provider instance. Keep q_A, candidate local q_B
   and verified physical Q_H separate; never silently relabel one as another.
2. Snapshot **all** state and histories: ne/np/nn, time and field constraints,
   electrode-charge/head histories, SP3 warm state where supported, external
   observers and cumulative budgets, as well as input/cache identities. Existing
   native file checkpoint reconstructs fields and resets electrode history; it
   does not yet satisfy this transactional protocol.
3. Add trial voltage evaluation from the same accepted state without permanent
   observer output. Validate token, geometry/operator/reference plane and time
   layer; commit PDE and circuit once together. A failed trial must not retain
   species updates, charge histories, source/field state or accumulated work.
4. Export actual native electron faces, implemented-species sources and all
   floor/clip changes. Record DeltaQ_H and integrated transfer charge separately;
   no numerical correction may be invented from a remaining residual.
5. Select TOTAL_CHARGE or SPLIT_GEOMETRIC_CAPACITANCE once. Global C_p, local
   diagnostic C_local, series C_ext, source capacitance and environment parasitics
   must retain distinct ownership. The current memory x is a synthetic C_d voltage,
   **not** a proposed additional I1 plasma state.
6. Match local-to-global external response and charge/energy extraction. Weighted
   particle flux F_o, electrostatic Q_o and outer power P_o,in are not the unweighted
   native outer electron count. Local zero-Neumann algebra is not global validation.
7. Define time integration/energy quadrature for the actual coupled FULL/DAE
   problem separately. The midpoint ODE tests here do not select a general DAE
   backend. Reuse the prescribed-voltage native kernel behind a current-residual
   adapter rather than changing transport or reaction physics.

## Scope and independent pending gates

- Native SP3-OFF, B0 local-A, 500 V / 10 ps references retain their original
  conditional and case-scoped evidence. This engineering work extends none of it.
- C1-L0 predicts only internal FULL_GAS G_R at the existing supported centered
  points. It supplies no terminal current, gas heat, validated feedback or lifting.
  k_star remains 9134636280.4 s^-1. The wrapper is a capability guard, not a physical
  port. Absolute elapsed time must retain the same origin; recreating a wrapper
  with a caller-falsified new origin cannot be prevented by a standalone library
  and is expressly outside its contract.
- Synthetic transactions PASS does not test native StreamerSolver rollback.
- Physical global port, B0 independent validation, low-field material validity,
  TENG source/mechanical energy, cold-hot transfer and RF remain unchanged.
- Historical G2/G3 parallel Cext topology stays frozen; the new series-Cext
  reference is implemented independently and is not an arbitrary netlist solver.

Next action: scientific review may authorize a native FULL adapter and its own
rollback/interval budgets on an agreed reference plane. No such connection,
PDE/Poisson run, source calibration or subsequent stage is automatically started.
