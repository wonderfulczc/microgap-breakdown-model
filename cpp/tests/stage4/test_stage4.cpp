#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
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

  AxisymmetricGrid g(8, 24, 2.0e-4, 0.0, 6.0e-4);
  StreamerConfig cfg;
  cfg.background_field = 4.8e6;
  cfg.n_ref = 1e8;
  cfg.photoionization = false;
  cfg.open_boundary.threshold = 1e-3;
  cfg.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver collision(g, cfg);
  collision.initialize_gaussians({{1e14, 4e-5, 2e-4}, {1e14, 4e-5, 4e-4}});
  const double left = collision.state().ne(0, 8);
  const double mid = collision.state().ne(0, 12);
  const double right = collision.state().ne(0, 16);
  if (!ck(left > mid && right > mid, "double-seed initialization produces two separated peaks")) return 1;

  StreamerSolver left_only(g, cfg), right_only(g, cfg);
  left_only.initialize_gaussians({{1e14, 4e-5, 2e-4}});
  right_only.initialize_gaussians({{1e14, 4e-5, 4e-4}});
  if (!ck(left_only.state().ne(0, 8) > left_only.state().ne(0, 16), "left isolated keeps only left seed")) return 1;
  if (!ck(right_only.state().ne(0, 16) > right_only.state().ne(0, 8), "right isolated keeps only right seed")) return 1;

  auto lim = collision.timestep_limits();
  StreamerDiagnostics d;
  if (!ck(collision.step(0.05 * lim.selected, d), "double-seed short coupled step accepted")) return 1;
  if (!ck(std::isfinite(d.emax) && d.total_electrons > 0.0, "double-seed diagnostics finite")) return 1;

  const auto cp = std::filesystem::temp_directory_path() / "streamer_rf_stage4_checkpoint.bin";
  collision.save_checkpoint(cp);
  StreamerSolver fork(g, cfg);
  fork.load_checkpoint(cp);
  std::filesystem::remove(cp);
  if (!ck(fork.state().time == collision.state().time && fork.state().ne.values() == collision.state().ne.values(), "checkpoint fork preserves Stage 4 state")) return 1;

  PetscFinalize();
  std::cout << "passed " << n << " Stage 4 C++ checks\n";
  return 0;
}
