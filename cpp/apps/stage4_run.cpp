#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <petscsys.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
constexpr double qe = 1.602176634e-19;
constexpr double pi = 3.14159265358979323846;

void write_fields(const std::filesystem::path& p, const AxisymmetricGrid& g, const StreamerSolver& solver) {
  const auto& s = solver.state();
  const auto source = solver.electron_transport_current_source();
  std::ofstream f(p);
  f << "i,j,time_s,r_m,z_m,ne_m_3,np_m_3,nn_m_3,rho_C_m_3,phi_V,Er_V_m,Ez_V_m,E_V_m,Sph_m_3_s_1,Jr_RF_A_m2,Jz_RF_A_m2\n"
    << std::setprecision(17);
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      f << i << ',' << j << ',' << s.time << ',' << g.r(i) << ',' << g.z(j) << ',' << s.ne(i, j) << ','
        << s.np(i, j) << ',' << s.nn(i, j) << ',' << s.rho(i, j) << ',' << s.phi(i, j) << ','
        << s.er(i, j) << ',' << s.ez(i, j) << ',' << s.emag(i, j) << ',' << s.sph(i, j) << ','
        << source.jr(i, j) << ',' << source.jz(i, j) << '\n';
    }
  }
}

struct CurrentMoment {
  double drift{}, electron{}, integral_abs_jz{}, max_abs_jz{};
};

CurrentMoment current_moment(const AxisymmetricGrid& g, const StreamerSolver& solver, const StreamerConfig& c) {
  const auto& s = solver.state();
  const auto source = solver.electron_transport_current_source();
  CurrentMoment out;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      const auto q = evaluate_morrow_lowke(s.emag(i, j), c.neutral_density, c.pressure, c.temperature);
      double grad = 0.0;
      if (j == 0) grad = (s.ne(i, 1) - s.ne(i, 0)) / g.dz();
      else if (j == g.nz() - 1) grad = (s.ne(i, j) - s.ne(i, j - 1)) / g.dz();
      else grad = (s.ne(i, j + 1) - s.ne(i, j - 1)) / (2.0 * g.dz());
      const double j_drift = qe * s.ne(i, j) * q.mobility * s.ez(i, j);
      const double vol = g.cell_volume(i);
      out.drift += j_drift * vol;
      (void)grad;
    }
  }
  out.electron = source.current_moment_z;
  out.integral_abs_jz = source.integral_abs_jz;
  out.max_abs_jz = source.max_abs_jz;
  return out;
}

struct CollisionMetrics {
  double left_inner_z{}, right_inner_z{}, d_head{};
  double bridge_min_ne{}, bridge_mean_ne{}, gap_emax{}, gap_mean_e{}, inner_rho_left{}, inner_rho_right{};
  int bridge_cells{};
};

CollisionMetrics collision_metrics(const AxisymmetricGrid& g, const StreamerState& s, double z1, double z2) {
  const double mid = 0.5 * (z1 + z2);
  int jl = 0, jr = g.nz() - 1;
  double best_l = -1.0, best_r = -1.0;
  for (int j = 1; j < g.nz() - 1; ++j) {
    const double z = g.z(j);
    if (z > z1 && z < mid) {
      const double score = std::abs(s.rho(0, j)) * (1.0 + s.emag(0, j) / 1e7);
      if (score > best_l) { best_l = score; jl = j; }
    } else if (z >= mid && z < z2) {
      const double score = std::abs(s.rho(0, j)) * (1.0 + s.emag(0, j) / 1e7);
      if (score > best_r) { best_r = score; jr = j; }
    }
  }
  if (jl > jr) std::swap(jl, jr);
  CollisionMetrics m;
  m.left_inner_z = g.z(jl);
  m.right_inner_z = g.z(jr);
  m.d_head = std::max(0.0, m.right_inner_z - m.left_inner_z);
  m.bridge_min_ne = std::numeric_limits<double>::infinity();
  double ne_sum = 0.0, e_sum = 0.0;
  m.gap_emax = 0.0;
  for (int j = jl; j <= jr; ++j) {
    m.bridge_min_ne = std::min(m.bridge_min_ne, s.ne(0, j));
    ne_sum += s.ne(0, j);
    e_sum += s.emag(0, j);
    m.gap_emax = std::max(m.gap_emax, s.emag(0, j));
    ++m.bridge_cells;
  }
  m.bridge_mean_ne = ne_sum / std::max(m.bridge_cells, 1);
  m.gap_mean_e = e_sum / std::max(m.bridge_cells, 1);
  m.inner_rho_left = s.rho(0, jl);
  m.inner_rho_right = s.rho(0, jr);
  if (!std::isfinite(m.bridge_min_ne)) m.bridge_min_ne = 0.0;
  return m;
}

double output_interval(double d_head) {
  if (d_head > 1.0e-3) return 5.0e-11;
  if (d_head > 2.0e-4) return 1.0e-11;
  return 5.0e-12;
}
}

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank = 0, size = 1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  const auto wall0 = std::chrono::steady_clock::now();

  if (argc < 19) {
    if (!rank) {
      std::cerr << "usage: output nr nz tend E0 nref sp3 eta rmax zmax dt_scale n0 sigma z1 z2 mode max_steps minimum_dt_s post_event_s [--sigma2 value] [--field-output-interval seconds] [resume_checkpoint] [initial_step]\n";
    }
    PetscFinalize();
    return 2;
  }

  const std::filesystem::path out = argv[1];
  const int nr = std::stoi(argv[2]), nz = std::stoi(argv[3]);
  const double tend = std::stod(argv[4]), E0 = std::stod(argv[5]), nref = std::stod(argv[6]);
  const bool sp3 = std::stoi(argv[7]);
  const double eta = std::stod(argv[8]), rmax = std::stod(argv[9]), zmax = std::stod(argv[10]);
  const double dt_scale = std::stod(argv[11]), n0 = std::stod(argv[12]), sigma = std::stod(argv[13]);
  const double z1 = std::stod(argv[14]), z2 = std::stod(argv[15]);
  const std::string mode = argv[16];
  const int max_steps = std::stoi(argv[17]);
  const double minimum_dt = std::stod(argv[18]);
  const double post_event = argc > 19 ? std::stod(argv[19]) : 0.2e-9;
  double sigma2 = sigma;
  double fixed_field_output_interval = -1.0;
  std::filesystem::path resume_checkpoint;
  int initial_step = 0;
  for (int a = 20; a < argc; ++a) {
    const std::string opt = argv[a];
    if (opt == "--sigma2" && a + 1 < argc) {
      sigma2 = std::stod(argv[++a]);
    } else if (opt == "--field-output-interval" && a + 1 < argc) {
      fixed_field_output_interval = std::stod(argv[++a]);
    } else if (resume_checkpoint.empty()) {
      resume_checkpoint = std::filesystem::path(argv[a]);
    } else {
      initial_step = std::stoi(argv[a]);
    }
  }

  AxisymmetricGrid g(nr, nz, rmax, 0.0, zmax);
  StreamerConfig c;
  c.background_field = E0;
  c.n_ref = nref;
  c.photoionization = sp3;
  c.open_boundary.mode = RingQuadratureMode::CellCenterRing;
  c.open_boundary.threshold = eta;
  c.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver solver(g, c);
  if (mode == "collision") solver.initialize_gaussians({{n0, sigma, z1}, {n0, sigma2, z2}});
  else if (mode == "left") solver.initialize_gaussians({{n0, sigma, z1}});
  else if (mode == "right") solver.initialize_gaussians({{n0, sigma2, z2}});
  else {
    if (!rank) std::cerr << "invalid mode " << mode << '\n';
    PetscFinalize();
    return 2;
  }
  const bool resuming = !resume_checkpoint.empty();
  if (resuming) solver.load_checkpoint(resume_checkpoint);

  if (!rank) {
    std::filesystem::create_directories(out);
    if (!resuming) {
      write_fields(out / ("fields_" + std::to_string(initial_step) + ".csv"), g, solver);
      solver.save_checkpoint(out / "checkpoint.bin");
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  std::ofstream hist, cm, coll;
  if (!rank) {
    hist.open(out / "scalar_history.csv", resuming ? std::ios::app : std::ios::out);
    if (!resuming) hist << "step,time_s,dt_s,controller,E_max_V_m,ne_max_m_3,total_electrons,total_charge_C,rejected_steps,poisson_iterations,sp3_ksp_iterations,sp3_boundary_iterations,conservation_residual,current_continuity_residual,outer_boundary_current_A,plasma_charge_derivative_A\n";
    hist << std::setprecision(17);
    cm.open(out / "current_moment.csv", resuming ? std::ios::app : std::ios::out);
    if (!resuming) cm << "time,I_CM_drift,I_CM_electron,integral_abs_Jz,max_abs_Jz\n";
    cm << std::setprecision(17);
    coll.open(out / "collision_metrics.csv", resuming ? std::ios::app : std::ios::out);
    if (!resuming) coll << "time,left_inner_z,right_inner_z,d_head,bridge_min_ne,bridge_mean_ne,E_gap_max,E_gap_mean,inner_rho_left,inner_rho_right,bridge_cells\n";
    coll << std::setprecision(17);
  }

  int step = initial_step, rejected_steps = 0;
  double last_dt = 0.0, event_time = -1.0, next_output = 0.0;
  std::string termination = "unknown";
  StreamerDiagnostics last_diag{};
  CollisionMetrics last_metrics = collision_metrics(g, solver.state(), z1, z2);
  const double initial_bridge = std::max(last_metrics.bridge_mean_ne, 1.0);
  double gap_peak = std::max(last_metrics.gap_emax, 1.0);
  next_output = solver.state().time + (fixed_field_output_interval > 0.0 ? fixed_field_output_interval : output_interval(last_metrics.d_head));

  while (solver.state().time < tend && step < max_steps) {
    auto lim = solver.timestep_limits();
    double dt = std::min(dt_scale * lim.selected, tend - solver.state().time);
    if (dt < minimum_dt) { termination = "minimum_dt"; last_dt = dt; break; }
    StreamerDiagnostics d;
    bool accepted = false;
    int local_rejects = 0;
    for (int retry = 0; retry <= 20 && !accepted; ++retry) {
      accepted = solver.step(dt, d);
      if (!accepted) { dt *= 0.5; ++local_rejects; }
    }
    rejected_steps += local_rejects;
    if (!accepted) { termination = "negativity_retry_exhausted"; break; }
    ++step;
    last_dt = dt;
    last_diag = d;
    const auto& s = solver.state();
    last_metrics = collision_metrics(g, s, z1, z2);
    gap_peak = std::max(gap_peak, last_metrics.gap_emax);
    const auto I = current_moment(g, solver, c);
    const bool distance_event = last_metrics.d_head <= std::max(2.0 * g.dz(), 2.0e-4);
    const bool bridge_field_event = s.time > 0.5e-9 &&
      last_metrics.bridge_mean_ne >= 50.0 * initial_bridge &&
      last_metrics.gap_emax <= 0.8 * gap_peak;
    const bool event_candidate = (mode == "collision" && (distance_event || bridge_field_event));
    if (event_candidate && event_time < 0.0) {
      event_time = s.time;
      if (!rank) solver.save_checkpoint(out / "checkpoint_pre_collision.bin");
    }

    if (!rank) {
      hist << step << ',' << d.time << ',' << d.dt << ',' << lim.controller << ',' << d.emax << ','
           << d.ne_max << ',' << d.total_electrons << ',' << d.total_charge << ',' << rejected_steps << ','
           << d.poisson_iterations << ',' << d.sp3_ksp_iterations << ',' << d.sp3_boundary_iterations << ','
           << d.conservation_residual << ',' << d.current_continuity_residual << ','
           << d.outer_boundary_current << ',' << d.plasma_charge_derivative << '\n';
      cm << d.time << ',' << I.drift << ',' << I.electron << ',' << I.integral_abs_jz << ',' << I.max_abs_jz << '\n';
      coll << d.time << ',' << last_metrics.left_inner_z << ',' << last_metrics.right_inner_z << ','
           << last_metrics.d_head << ',' << last_metrics.bridge_min_ne << ',' << last_metrics.bridge_mean_ne << ','
           << last_metrics.gap_emax << ',' << last_metrics.gap_mean_e << ',' << last_metrics.inner_rho_left << ','
           << last_metrics.inner_rho_right << ',' << last_metrics.bridge_cells << '\n';
    }

    const bool due = s.time + 1e-18 >= next_output || s.time >= tend || (event_time > 0.0 && s.time - event_time < 2.0e-11);
    if (!rank && due) {
      write_fields(out / ("fields_" + std::to_string(step) + ".csv"), g, solver);
      solver.save_checkpoint(out / "checkpoint.bin");
      hist.flush(); cm.flush(); coll.flush();
      next_output = s.time + (fixed_field_output_interval > 0.0 ? fixed_field_output_interval : output_interval(last_metrics.d_head));
    }
    if (event_time > 0.0 && post_event > 0.0 && s.time - event_time >= post_event) {
      termination = "reached_collision_post_window";
      if (!rank) solver.save_checkpoint(out / "checkpoint_post_collision.bin");
      break;
    }
    if (!std::isfinite(s.time) || !std::isfinite(d.emax) || !std::isfinite(d.ne_max)) { termination = "invalid_state"; break; }
  }

  if (termination == "unknown") {
    if (solver.state().time >= tend) termination = "reached_end_time";
    else if (step >= max_steps) termination = "max_steps";
  }
  const auto wall1 = std::chrono::steady_clock::now();
  const double wall_time = std::chrono::duration_cast<std::chrono::duration<double>>(wall1 - wall0).count();
  if (!rank) {
    write_fields(out / "fields_final.csv", g, solver);
    solver.save_checkpoint(out / "checkpoint.bin");
    std::ofstream meta(out / "termination.csv");
    meta << "termination_reason,requested_end_time,actual_end_time,completed_steps,accepted_steps,rejected_steps,last_dt,wall_time,checkpoint_path,event_time,mpi_ranks\n";
    meta << termination << ',' << tend << ',' << solver.state().time << ',' << (step - initial_step) << ','
         << (step - initial_step) << ',' << rejected_steps << ',' << last_dt << ',' << wall_time << ','
         << (out / "checkpoint.bin").string() << ',' << event_time << ',' << size << '\n';
    std::ofstream reason(out / "termination_reason.txt");
    reason << termination << '\n';
    std::cout << "steps=" << step << " final_time_s=" << solver.state().time << " mode=" << mode
              << " event_time_s=" << event_time << " termination_reason=" << termination << " ranks=" << size << '\n';
  }

  PetscFinalize();
  if (termination == "negativity_retry_exhausted") return 3;
  if (termination == "invalid_state") return 5;
  if (termination == "minimum_dt") return 6;
  return 0;
}
