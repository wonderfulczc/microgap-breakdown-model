#pragma once
#include "streamer_rf/electrode.hpp"
#include "streamer_rf/geometry2d.hpp"
#include "streamer_rf/petsc_solver.hpp"
#include <memory>
namespace streamer_rf {
struct PoissonBoundaryConfig { BoundaryCondition r_outer,z_lower,z_upper; };
enum class RingQuadratureMode { CellCenterRing, GaussianCellQuadrature };
struct OpenBoundaryOptions { RingQuadratureMode mode{RingQuadratureMode::GaussianCellQuadrature}; double threshold{0.0}; };
struct PoissonTiming {double assembly_s{},setup_s{},solve_s{},transfer_s{},residual_s{},true_relative_residual{};int iterations{},reason{};bool reused{};};
// Serial transaction diagnostics only; no matrix, PETSc object, or field state.
struct PoissonDiagnosticSnapshot {int elliptic_iterations{},poisson_iterations{};bool profiling{};std::vector<PoissonTiming> records;};
PoissonDiagnosticSnapshot snapshot_poisson_diagnostics();
void restore_poisson_diagnostics(const PoissonDiagnosticSnapshot&);
void set_poisson_profiling(bool);
const std::vector<PoissonTiming>& poisson_timing_records();
void clear_poisson_timing_records();
std::string electrode_operator_identity(const Geometry2DView&,const PoissonBoundaryConfig&,const SolverTolerances&);
class FixedElectrodePoisson {
 public:
  FixedElectrodePoisson(const Geometry2DView&,const PoissonBoundaryConfig&,const SolverTolerances&);
  ~FixedElectrodePoisson();
  FixedElectrodePoisson(const FixedElectrodePoisson&)=delete;
  void solve(const ScalarField2D&,double,ScalarField2D&);
  bool matches_fresh_matrix(const ScalarField2D&,double)const;
 private: struct Impl;std::unique_ptr<Impl> impl_;
};
void solve_potential(const ScalarField2D& charge_density,double background_field,
 const PoissonBoundaryConfig&,ScalarField2D& potential,const SolverTolerances& = {},const OpenBoundaryOptions& = {});
void solve_potential_with_electrodes(const ScalarField2D& charge_density,
 const AxisymmetricNeedlePlaneGeometry& geometry,double applied_voltage,
 const PoissonBoundaryConfig&,ScalarField2D& potential,const SolverTolerances& = {},const OpenBoundaryOptions& = {});
void solve_potential_with_electrodes(const ScalarField2D& charge_density,
 const Geometry2DView& geometry,double applied_voltage,
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
