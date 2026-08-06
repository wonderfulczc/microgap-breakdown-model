#pragma once
#include "streamer_rf/poisson.hpp"
#include <array>
#include <vector>
namespace streamer_rf {
struct Sp3Group { double A_original,lambda_original,a_si,k_si; };
struct Sp3Constants { double kappa1_sq,kappa2_sq,gamma1,gamma2; std::array<Sp3Group,3> groups; };
struct Sp3BoundaryOptions { double rtol{1e-8},atol{1e-14},relaxation{0.5};int max_iterations{500};bool zero_dirichlet{false}; };
struct Sp3GroupDiagnostics { int iterations{};double initial_residual{},final_residual{},sph_change{};bool converged{};std::vector<double> residual_history; };
struct Sp3WarmStart { std::array<std::vector<double>,6> components; };
Sp3Constants sp3_constants(const PressureConfig&);
void solve_photoionization(const ScalarField2D& emission_source,const PressureConfig&,
 const PoissonBoundaryConfig&,ScalarField2D& photoionization_rate,const SolverTolerances& = {},
 const Sp3BoundaryOptions& = {},std::array<Sp3GroupDiagnostics,3>* = nullptr,
 Sp3WarmStart* = nullptr);
int last_sp3_ksp_iterations();
int last_sp3_boundary_iterations();
}
