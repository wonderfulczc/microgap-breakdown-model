# Stage 3 Method

## Scope and model identity

Stage 3 develops and validates a coupled axisymmetric three-species fluid solver. It uses Morrow–Lowke (1997) Appendix A transport functions and therefore is not a quantitative Shi et al. (2019) reconstruction. Liu–Pasko Figure 1 calibration is deferred to Stage 4.

## Equations and coupling

The state contains electron, positive-ion and negative-ion number density. Ions are immobile. Electron finite-volume fluxes use the Stage 2 ISG-0 interface with the left and right cell drift velocities, `v_e=-mu_e E`. Charge density is passed to the existing Stage 2 OpenCharge Poisson solver every accepted step. The impact-ionization emission field is passed to the existing three-group SP3 solver with its coupled Robin boundary.

The first-order step order is charge/Poisson, electric field, Morrow–Lowke coefficients, emission/SP3, ISG-0 fluxes, reactions, density update, positivity check and diagnostics. A rejected step leaves all three densities and time unchanged and is retried with half the time step.

## Units and material functions

Internal units are SI. Morrow–Lowke reduced fields printed in V cm² are converted from SI `E/N` by multiplying by `1e4`; densities are converted to cm⁻³ only inside the Appendix formula transcription. Townsend coefficients are converted to m⁻¹ and multiplied by drift speed to obtain s⁻¹ reaction frequencies. Details and the documented A10 junction regularization are in `stage3_formula_audit.md`.

## Time step

The selected time step is the minimum of drift, diffusion, ionization, dielectric relaxation and particle-loss limits with factors 0.5, 0.5, 0.05, 0.2 and 0.2. Every history row records the actual accepted step and controller.

## Boundaries and conservation

At the axis, radial electron flux and radial electric field are zero. Outer electron boundaries have zero diffusive gradient, permit outward upwind drift and prohibit inflow. Ions have zero normal flux. The potential uses the updated background-plus-ring-charge boundary. Conservation diagnostics close the electron balance against integrated reaction and boundary flux, and close net charge against electron boundary flux; the latter is normalized by the total charged-particle inventory rather than nearly zero net charge.

## Evidence

Executables produce raw CSV fields and scalar histories. Python converts raw fields to HDF5, derives trajectories and sensitivities, and registers executable/config/result hashes. Development runs invalidated by implementation defects are isolated under `results/stage3/forensic` and are excluded from acceptance.

## Stage 4 handoff

Before target reconstruction, Stage 4 must use either an audited Liu–Pasko solid-curve digitization or an independent Boltzmann transport table checked graphically against Figure 1. Stage 3 Morrow–Lowke values cannot be relabeled as Shi et al. target parameters.
