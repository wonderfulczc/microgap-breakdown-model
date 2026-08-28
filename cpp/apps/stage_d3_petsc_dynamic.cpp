#include "streamer_rf/electrode.hpp"
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/types.hpp"
#include "streamer_rf/voltage_waveform.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Options {
  double voltage{500.0};
  std::filesystem::path out{"solver3d/afivo_reference/common_benchmark/d3/petsc/output/d3_petsc_dynamic.csv"};
};

Options parse(int argc, char** argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--voltage" && i + 1 < argc) {
      o.voltage = std::stod(argv[++i]);
    } else if (a == "--out" && i + 1 < argc) {
      o.out = argv[++i];
    } else {
      throw std::invalid_argument("usage: stage_d3_petsc_dynamic [--voltage V] [--out file.csv]");
    }
  }
  return o;
}

struct CompactDiag {
  double time{}, emax{}, ne_max{}, total_electrons{}, head_z{}, head_v{};
  bool bridge{};
};

CompactDiag sample(const streamer_rf::AxisymmetricGrid& g,
                   const streamer_rf::AxisymmetricNeedlePlaneGeometry& geom,
                   const streamer_rf::streamer::StreamerSolver& solver,
                   double head_threshold, double previous_head, double previous_time) {
  const auto& s = solver.state();
  CompactDiag d;
  d.time = s.time;
  double best_ne = -1.0;
  double best_z = 0.0;
  bool threshold_hit = false;
  bool in_gap = false;
  d.bridge = true;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (geom.classify(g, i, j) != streamer_rf::ElectrodeCellType::Gas) continue;
      const double vol = g.cell_volume(i);
      d.total_electrons += s.ne(i, j) * vol;
      d.ne_max = std::max(d.ne_max, s.ne(i, j));
      d.emax = std::max(d.emax, s.emag(i, j));
    }
    if (geom.classify(g, 0, j) == streamer_rf::ElectrodeCellType::Gas) {
      const double n = s.ne(0, j);
      if (n >= head_threshold) {
        if (!threshold_hit) {
          best_z = g.z(j);
          threshold_hit = true;
        } else {
          best_z = std::min(best_z, g.z(j));
        }
      }
      if (!threshold_hit && n > best_ne) {
        best_ne = n;
        best_z = g.z(j);
      }
      if (g.z(j) > 0.0 && g.z(j) < 70e-6) {
        in_gap = true;
        if (n < 1e18) d.bridge = false;
      }
    }
  }
  if (!in_gap) d.bridge = false;
  d.head_z = best_z;
  d.head_v = std::isfinite(previous_head) && d.time > previous_time ? (d.head_z - previous_head) / (d.time - previous_time) : 0.0;
  return d;
}

void write_axis_profile(const std::filesystem::path& path, const streamer_rf::AxisymmetricGrid& g,
                        const streamer_rf::AxisymmetricNeedlePlaneGeometry& geom,
                        const streamer_rf::streamer::StreamerSolver& solver) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream f(path);
  f << std::setprecision(17) << "z_m,ne_m3,phi_V,Eabs_Vpm,cell_class\n";
  const auto& s = solver.state();
  for (int j = 0; j < g.nz(); ++j) {
    const int i = 0;
    const auto c = geom.classify(g, i, j);
    const char* cls = c == streamer_rf::ElectrodeCellType::Gas
                          ? "gas"
                          : (c == streamer_rf::ElectrodeCellType::HighVoltageElectrode ? "hv" : "ground");
    f << g.z(j) << ',' << s.ne(i, j) << ',' << s.phi(i, j) << ',' << s.emag(i, j) << ',' << cls << '\n';
  }
}
}  // namespace

int main(int argc, char** argv) {
  bool petsc_initialized = false;
  try {
    const auto opt = parse(argc, argv);
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    petsc_initialized = true;

    streamer_rf::AxisymmetricGrid grid(64, 148, 80e-6, -2.5e-6, 90e-6);
    streamer_rf::AxisymmetricNeedlePlaneGeometry geom(
        "D3-axisymmetric-dynamic", -2.5e-6, 75e-6, 5e-6, 5e-6, 2.5e-6);
    streamer_rf::ConstantVoltage voltage(opt.voltage);
    streamer_rf::streamer::StreamerConfig cfg;
    cfg.electrode_geometry = &geom;
    cfg.voltage_waveform = &voltage;
    cfg.photoionization = false;
    cfg.n_ref = 1e12;
    cfg.head_ne_threshold = 1e14;
    cfg.bridge_ne_threshold = 1e18;
    cfg.elliptic = {1e-10, 1e-14, 20000};
    cfg.poisson_boundary.r_outer = {streamer_rf::BoundaryKind::Neumann, 0.0};
    cfg.poisson_boundary.z_lower = {streamer_rf::BoundaryKind::Dirichlet, 0.0};
    cfg.poisson_boundary.z_upper = {streamer_rf::BoundaryKind::Dirichlet, opt.voltage};

    streamer_rf::streamer::StreamerSolver solver(grid, cfg);
    solver.initialize_gaussian(1e16, std::sqrt(2.0) * 3e-6, 65e-6);

    const std::vector<double> targets{0.0, 0.5e-12, 1e-12, 2e-12, 3e-12, 5e-12};
    std::filesystem::create_directories(opt.out.parent_path());
    std::ofstream f(opt.out);
    f << std::setprecision(17)
      << "time_s,voltage_V,Emax_Vpm,ne_max_m3,total_electrons,head_z_m,head_velocity_mps,bridge_flag\n";

    double prev_head = std::numeric_limits<double>::quiet_NaN();
    double prev_time = std::numeric_limits<double>::quiet_NaN();
    int output_index = 0;
    int steps = 0;
    bool ok = true;
    for (double target : targets) {
      while (solver.state().time < target - 1e-18) {
        auto lim = solver.timestep_limits();
        double dt = std::min({0.5 * lim.selected, 2.0e-14, target - solver.state().time});
        streamer_rf::streamer::StreamerDiagnostics step_diag;
        ok = false;
        for (int retry = 0; retry < 12; ++retry) {
          ok = solver.step(dt, step_diag);
          if (ok) break;
          dt *= 0.5;
        }
        if (!ok) throw std::runtime_error("PETSc D3 dynamic step failed");
        ++steps;
      }
      auto d = sample(grid, geom, solver, cfg.head_ne_threshold, prev_head, prev_time);
      f << d.time << ',' << opt.voltage << ',' << d.emax << ',' << d.ne_max << ','
        << d.total_electrons << ',' << d.head_z << ',' << d.head_v << ',' << d.bridge << '\n';
      write_axis_profile(opt.out.parent_path() / ("axis_profile_" + std::to_string(output_index) + ".csv"), grid, geom, solver);
      prev_head = d.head_z;
      prev_time = d.time;
      ++output_index;
    }
    std::ofstream meta(opt.out.parent_path() / "d3_petsc_metadata.txt");
    meta << std::setprecision(17)
         << "voltage_V=" << opt.voltage << "\nfinal_time_s=" << solver.state().time
         << "\nsteps=" << steps << "\nhead_threshold_m3=" << cfg.head_ne_threshold
         << "\nseed_n0_m3=1e16\nseed_sigma_common_m=3e-6\nseed_internal_width_m=" << std::sqrt(2.0) * 3e-6
         << "\nphotoionization=OFF\n";
    std::cout << "stage_d3_petsc_dynamic status=PASS voltage=" << opt.voltage
              << " final_time_s=" << solver.state().time << " steps=" << steps
              << " out=" << opt.out << '\n';
  } catch (const std::exception& e) {
    std::cerr << "stage_d3_petsc_dynamic status=FAIL error=" << e.what() << '\n';
    if (petsc_initialized) PetscFinalize();
    return 2;
  }
  PetscFinalize();
  return 0;
}
