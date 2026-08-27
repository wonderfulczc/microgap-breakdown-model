#include "streamer_rf/poisson.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace streamer_rf;

namespace {
void electric_field(const ScalarField2D& phi, ScalarField2D& er, ScalarField2D& ez) {
  const auto& g = phi.grid();
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      auto dr = [&](int a, int b) { return (phi(b, j) - phi(a, j)) / ((b - a) * g.dr()); };
      auto dz = [&](int a, int b) { return (phi(i, b) - phi(i, a)) / ((b - a) * g.dz()); };
      er(i, j) = i == 0 ? 0.0 : -(i == g.nr() - 1 ? dr(i - 1, i) : dr(i - 1, i + 1));
      ez(i, j) = -(j == 0 ? dz(0, 1) : j == g.nz() - 1 ? dz(j - 1, j) : dz(j - 1, j + 1));
    }
  }
}
}  // namespace

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank = 0, size = 1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  AxisymmetricGrid g(18, 36, 0.9e-3, 0.0, 1.8e-3);
  AxisymmetricNeedlePlaneGeometry geom("stage-c1-smoke-needle-plane", 0.0, 1.55e-3, 0.10e-3, 0.05e-3, 0.10e-3);
  ConstantVoltage voltage(1500.0);
  ScalarField2D rho(g), phi(g), er(g), ez(g);
  PoissonBoundaryConfig bc;
  bc.r_outer.kind = bc.z_lower.kind = bc.z_upper.kind = BoundaryKind::Dirichlet;
  solve_potential_with_electrodes(rho, geom, voltage.value(0.0), bc, phi, {1e-12, 1e-14, 20000});
  electric_field(phi, er, ez);
  const auto diag = evaluate_electrode_diagnostics(phi, er, ez, geom, 0.0, voltage.value(0.0), last_poisson_iterations());
  double checksum = 0.0;
  for (std::size_t k = 0; k < phi.values().size(); ++k) checksum += (k + 1) * phi.values()[k];
  double lo = 0.0, hi = 0.0;
  MPI_Allreduce(&checksum, &lo, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
  MPI_Allreduce(&checksum, &hi, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
  const double rel = std::abs(hi - lo) / std::max(1.0, std::abs(hi));
  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "stage_c1_smoke ranks=" << size
              << " voltage=" << diag.applied_voltage_V
              << " Emax=" << diag.emax_V_m
              << " phi_HV_residual=" << diag.phi_hv_residual_V
              << " phi_ground_residual=" << diag.phi_ground_residual_V
              << " poisson_iterations=" << diag.ksp_iterations
              << " checksum=" << checksum
              << " relative_difference=" << rel << '\n';
  }
  const bool ok = std::isfinite(diag.emax_V_m) && diag.emax_V_m > 0.0 &&
                  diag.phi_hv_residual_V < 1e-8 && diag.phi_ground_residual_V < 1e-8 &&
                  rel < 1e-11;
  PetscFinalize();
  return ok ? 0 : 1;
}
