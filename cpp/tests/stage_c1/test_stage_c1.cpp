#include "streamer_rf/poisson.hpp"
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
int checks = 0;

void check(bool ok, const char* name) {
  ++checks;
  if (!ok) throw std::runtime_error(name);
  std::cout << "ok " << checks << " - " << name << '\n';
}

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
  try {
    AxisymmetricGrid g(16, 32, 0.8e-3, 0.0, 1.6e-3);
    AxisymmetricNeedlePlaneGeometry geom("stage-c1-test-needle-plane", 0.0, 1.35e-3, 0.08e-3, 0.04e-3, 0.08e-3);

    check(geom.classify(0.02e-3, 1.35e-3) == ElectrodeCellType::HighVoltageElectrode, "HV tip classification");
    check(geom.classify(0.4e-3, 0.025e-3) == ElectrodeCellType::GroundElectrode, "ground plane classification");
    check(geom.classify(0.4e-3, 0.8e-3) == ElectrodeCellType::Gas, "gas classification");
    check(std::abs(geom.metadata().gap_m - 1.35e-3) < 1e-15, "geometry metadata gap");

    ConstantVoltage cv(1200.0);
    check(cv.value(1.0) == 1200.0, "constant voltage");
    SampledVoltage sv({0.0, 1.0e-9, 2.0e-9}, {0.0, 1000.0, 500.0});
    check(sv.value(0.0) == 0.0 && sv.value(1.0e-9) == 1000.0, "sampled voltage knot values");
    check(std::abs(sv.value(0.5e-9) - 500.0) < 1e-12, "sampled voltage linear midpoint");
    bool rejected = false;
    try {
      SampledVoltage bad({0.0, 0.0}, {0.0, 1.0});
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    check(rejected, "sampled voltage rejects nonmonotone time");
    rejected = false;
    try {
      (void)sv.value(3.0e-9);
    } catch (const std::out_of_range&) {
      rejected = true;
    }
    check(rejected, "sampled voltage rejects default extrapolation");

    ScalarField2D rho(g), phi(g), legacy(g);
    PoissonBoundaryConfig bc;
    solve_potential(rho, 2.0e5, bc, legacy);
    double legacy_error = 0.0;
    for (int j = 0; j < g.nz(); ++j) {
      for (int i = 0; i < g.nr(); ++i) legacy_error = std::max(legacy_error, std::abs(legacy(i, j) + 2.0e5 * g.z(j)));
    }
    check(legacy_error < 1e-9, "legacy solve_potential background field path");

    bc.r_outer.kind = bc.z_lower.kind = bc.z_upper.kind = BoundaryKind::Dirichlet;
    solve_potential_with_electrodes(rho, geom, cv.value(0.0), bc, phi, {1e-12, 1e-14, 20000});
    ScalarField2D er(g), ez(g);
    electric_field(phi, er, ez);
    auto diag = evaluate_electrode_diagnostics(phi, er, ez, geom, 0.0, 1200.0, last_poisson_iterations());
    check(diag.phi_hv_residual_V < 1e-8, "HV electrode Dirichlet residual");
    check(diag.phi_ground_residual_V < 1e-8, "ground electrode Dirichlet residual");
    check(std::isfinite(diag.emax_V_m) && diag.emax_V_m > 0.0, "electrostatic smoke finite Emax");
    check(last_poisson_iterations() >= 0, "electrostatic smoke KSP iterations recorded");

    StreamerConfig cfg;
    cfg.electrode_geometry = &geom;
    cfg.voltage_waveform = &cv;
    cfg.photoionization = false;
    StreamerSolver solver(g, cfg);
    solver.initialize_gaussian_at_tip_offset(1e14, 0.05e-3, -0.12e-3);
    const double expected_z = geom.tip_z_m() - 0.12e-3;
    double weighted_z = 0.0, weight = 0.0;
    for (int j = 0; j < g.nz(); ++j) {
      for (int i = 0; i < g.nr(); ++i) {
        const double w = solver.state().ne(i, j);
        weighted_z += w * g.z(j);
        weight += w;
      }
    }
    check(std::abs(weighted_z / weight - expected_z) < 0.35 * g.dz(), "GaussianSeed tip-relative positioning");
    ScalarField2D ser(g), sez(g);
    electric_field(solver.state().phi, ser, sez);
    auto sdiag = evaluate_electrode_diagnostics(solver.state().phi, ser, sez, geom, 0.0, 1200.0, last_poisson_iterations());
    check(sdiag.geometry_id == "stage-c1-test-needle-plane" && sdiag.applied_voltage_V == 1200.0, "Stage C1 diagnostics metadata");
    check(sdiag.phi_hv_residual_V < 1e-8 && sdiag.phi_ground_residual_V < 1e-8, "Stage C1 diagnostics conductor residuals");

    std::cout << "passed " << checks << " Stage C1 C++ checks\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
  PetscFinalize();
  return 0;
}
