# Stage 2 Method

## Scope, symbols, and units

Stage 2 contains only ISG-0 flux verification, an axisymmetric PETSc Poisson building block, and three-group SP3 photoionization. Internal quantities are SI. `Gamma_e=W n_e-D_e grad(n_e)` is a particle flux and `W=-mu_e E`; it is not conventional current.

`AxisymmetricGrid` is cell centred and uses exact annular volumes `pi(r_+^2-r_-^2) dz`. `ScalarField2D`, `SolverTolerances`, boundary metadata, PETSc context, and unit metadata are shared.

## ISG-0 and n_ref

The ordinary SG implementation uses a stable Bernoulli function. ISG-0 follows Kulikovsky's virtual-width construction and records ordinary-SG fallback and zero-width branches. Interpolation is applied to `u=n/n_ref`; callers must explicitly supply positive `n_ref`. Stage 2 uses `n_ref=1` only for nondimensional benchmarks. A physical choice is deferred.

## Poisson and open boundaries

The C++ production framework uses a 2-D PETSc DMDA, conservative radial face coefficients, KSPBCGS, and PCGAMG. The electrostatic boundary reference uses the azimuthally integrated ring kernel with the complete elliptic integral; the axis limit is analytic. Manufactured quadratic and quartic fields exercise the axis stencil and second-order trend.

## SP3 equations and boundaries

Bourdon Table 3 coefficients retain their original units in configuration and are converted before matrix assembly. `kappa` and `gamma` are evaluated from their radical formulas. Two shifted elliptic equations per group give six PETSc solves. Liu's no-reflection/no-emission coupled Robin formula is the governing boundary reference; zero Dirichlet is a sensitivity control only.

## PETSc, MPI, and tolerances

CMake discovers PETSc through pkg-config and links MPI explicitly because the distribution PETSc `.pc` file does not expose `mpi.h`. Candidate tolerances are scanned. The frozen Poisson and SP3 settings are `rtol=1e-10`, `atol=1e-14`, with changes below 0.1% relative to the strict run. MPI checks cover one and two ranks, with four ranks when resources permit.

## Validation and exclusions

The tests cover local flux identities, mass conservation, manufactured convergence, ring and Gaussian references, Table-3 conversion, the SP3 kernel fit, C++/Python consistency, and MPI execution. This stage contains no chemistry, coupled electron/ion stepping, streamer propagation, collision, Case I, current moment, or radiation output.

## Stage 3 interfaces and transport decision

The generic `TransportTable` rejects out-of-range input. Before Stage 3, Liu and Pasko curves are to be audibly digitized to CSV with original and SI units and interpolated by PCHIP. High-order spline extrapolation is forbidden; Morrow-Lowke is sensitivity-only. This decision is recorded here rather than modifying the frozen Stage 1 `decision_log.md`.

