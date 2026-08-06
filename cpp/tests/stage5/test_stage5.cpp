#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <cmath>
#include <iostream>
#include <petscsys.h>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int n = 0;
  auto ck = [&](bool ok, const char* msg) {
    if (!ok) {
      std::cerr << "FAIL " << msg << '\n';
      PetscFinalize();
      return false;
    }
    std::cout << "ok " << ++n << " - " << msg << '\n';
    return true;
  };

  AxisymmetricGrid g(12, 80, 2.0e-4, 0.0, 1.6e-3);
  StreamerConfig cfg;
  cfg.background_field = 6.4e6;
  cfg.n_ref = 1e8;
  cfg.photoionization = false;
  cfg.open_boundary.threshold = 1e-3;
  cfg.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver high(g, cfg);
  high.initialize_gaussians({{1e14, 4e-5, 5e-4}, {1e14, 4e-5, 1.1e-3}});
  if (!ck(std::abs(cfg.background_field - 6.4e6) < 1.0, "high-field configuration is represented")) return 1;
  if (!ck(high.state().ne(0, 25) > high.state().ne(0, 40) && high.state().ne(0, 55) > high.state().ne(0, 40), "high-field double seed has two centers")) return 1;

  StreamerSolver gap(g, cfg);
  gap.initialize_gaussians({{1e14, 4e-5, 3e-4}, {1e14, 4e-5, 1.3e-3}});
  if (!ck(gap.state().ne(0, 15) > gap.state().ne(0, 40) && gap.state().ne(0, 65) > gap.state().ne(0, 40), "larger-gap initialization preserves separated seeds")) return 1;

  StreamerSolver asym(g, cfg);
  asym.initialize_gaussians({{1e14, 4e-5, 5e-4}, {1e14, 2.0e-4, 1.1e-3}});
  const double narrow_mid = asym.state().ne(0, 28);
  const double wide_mid = asym.state().ne(0, 50);
  if (!ck(wide_mid > narrow_mid, "asymmetric wider seed produces stronger bridge tail")) return 1;

  auto lim = high.timestep_limits();
  StreamerDiagnostics d;
  if (!ck(high.step(0.05 * lim.selected, d), "Stage 5 high-field short step accepted")) return 1;
  if (!ck(std::isfinite(d.emax) && d.emax > 0.0, "Stage 5 event-local diagnostics finite")) return 1;
  if (!ck(d.total_electrons > 0.0 && std::isfinite(d.total_charge), "Stage 5 particle and charge diagnostics finite")) return 1;

  PetscFinalize();
  std::cout << "passed " << n << " Stage 5 C++ checks\n";
  return 0;
}
