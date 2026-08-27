#pragma once
#include "streamer_rf/electrode.hpp"
#include "streamer_rf/petsc_solver.hpp"
namespace streamer_rf {
struct PoissonBoundaryConfig { BoundaryCondition r_outer,z_lower,z_upper; };
enum class RingQuadratureMode { CellCenterRing, GaussianCellQuadrature };
struct OpenBoundaryOptions { RingQuadratureMode mode{RingQuadratureMode::GaussianCellQuadrature}; double threshold{0.0}; };
void solve_potential(const ScalarField2D& charge_density,double background_field,
 const PoissonBoundaryConfig&,ScalarField2D& potential,const SolverTolerances& = {},const OpenBoundaryOptions& = {});
void solve_potential_with_electrodes(const ScalarField2D& charge_density,
 const AxisymmetricNeedlePlaneGeometry& geometry,double applied_voltage,
 const PoissonBoundaryConfig&,ScalarField2D& potential,const SolverTolerances& = {},const OpenBoundaryOptions& = {});
double ring_potential_on_axis(double charge,double ring_radius,double axial_separation);
double ring_potential(double charge,double source_radius,double observation_radius,double axial_separation);
double complete_elliptic_k(double parameter);
ScalarField2D open_charge_boundary(const ScalarField2D&,const OpenBoundaryOptions& = {});
void solve_axisymmetric_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig&,
 ScalarField2D& solution,double shift,const SolverTolerances& = {});
void solve_axisymmetric_dirichlet(const ScalarField2D& rhs,const ScalarField2D& boundary,
 ScalarField2D& solution,double shift,const SolverTolerances& = {});
void solve_axisymmetric_robin(const ScalarField2D& rhs,ScalarField2D& solution,double shift,
 double self_coefficient,double cross_coefficient,const ScalarField2D& partner,const SolverTolerances& = {});
int last_elliptic_iterations();
int last_poisson_iterations();
}
