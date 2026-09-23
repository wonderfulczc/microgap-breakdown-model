#include "streamer_rf/streamer/HeadTracker.hpp"
#include "streamer_rf/streamer/UltrafastEventDiagnostics.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <sys/resource.h>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
struct Options {
  std::filesystem::path out{"rf/c_r5/development"};
  int steps{24};
  double voltage_V{500.0};
  double dt_cap_s{1e-15};
};

Options parse(int argc, char** argv) {
  Options options;
  if (argc > 1) options.out = argv[1];
  for (int index = 2; index < argc; ++index) {
    const std::string arg = argv[index];
    auto need = [&](const char* name) {
      if (index + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
      return argv[++index];
    };
    if (arg == "--steps") options.steps = std::stoi(need("--steps"));
    else if (arg == "--voltage") options.voltage_V = std::stod(need("--voltage"));
    else if (arg == "--dt-cap") options.dt_cap_s = std::stod(need("--dt-cap"));
    else throw std::runtime_error("unknown option " + arg);
  }
  if (options.steps < 1 || !(options.dt_cap_s > 0.0)) throw std::runtime_error("invalid dry-run options");
  return options;
}

void write_trace_header(std::ostream& stream) {
  stream << "step,time_s,dt_s,topology,topology_confidence,topology_evidence,topology_status,head_count,"
            "head_position_m,head_velocity_m_s,bridge_status,roi_definition,roi_status,roi_cell_count,"
            "roi_volume_m3,roi_centroid_r_m,roi_centroid_z_m,roi_r_min_m,roi_r_max_m,roi_z_min_m,roi_z_max_m,"
            "electrode_distance_m,Ek_V_m,gas,pressure_Pa,temperature_K,Ek_semantics,E0_source,"
            "E0_max_V_m,E0_mean_V_m,E0_p95_V_m,eta_0_max,eta_0_mean,eta_0_p95,"
            "E_peak_V_m,E_peak_time_s,E_peak_r_m,E_peak_z_m,eta_peak,ne_at_E_peak_m3,ne_peak_m3,"
            "nu_i_at_E_peak_s_1,mu_e_at_E_peak_m2_V_s,tau_i_at_E_peak_s,tau_i_at_ne_peak_s,tau_i_min_s,"
            "tau_i_median_s,tau_i_p95_s,sigma_e_at_E_peak_S_m,sigma_e_at_ne_peak_S_m,sigma_e_peak_S_m,"
            "sigma_e_median_S_m,sigma_e_p95_S_m,tau_M_at_E_peak_s,tau_M_at_ne_peak_s,tau_M_min_s,"
            "tau_M_median_s,tau_M_p95_s,Pi_RF_at_E_peak,Pi_RF_at_ne_peak,Pi_RF_min,Pi_RF_median,Pi_RF_p95,"
            "N_tau_i_at_E_peak,N_tau_M_at_E_peak,tau_i_resolution_status,tau_M_resolution_status,"
            "current_moment_r_A_m,current_moment_z_A_m\n";
}

void write_trace_row(std::ostream& stream, const UltrafastEventSample& value) {
  stream << value.step << ',' << value.time_s << ',' << value.dt_s << ',' << event_topology_name(value.topology) << ','
         << value.topology_confidence << ',' << value.topology_evidence << ',' << value.topology_status << ','
         << value.head_count << ',' << value.head_position_m << ',' << value.head_velocity_m_s << ','
         << value.bridge_status << ',' << value.roi_definition << ',' << value.roi_status << ','
         << value.roi_cell_count << ',' << value.roi_volume_m3 << ',' << value.roi_centroid_r_m << ','
         << value.roi_centroid_z_m << ',' << value.roi_r_min_m << ',' << value.roi_r_max_m << ','
         << value.roi_z_min_m << ',' << value.roi_z_max_m << ',' << value.electrode_distance_m << ','
         << value.Ek_V_m << ',' << value.gas << ',' << value.pressure_Pa << ',' << value.temperature_K << ','
         << value.Ek_semantics << ',' << value.E0_source << ',' << value.E0_max_V_m << ','
         << value.E0_mean_V_m << ',' << value.E0_p95_V_m << ','
         << value.eta_0_max << ',' << value.eta_0_mean << ',' << value.eta_0_p95 << ',' << value.E_peak_V_m << ','
         << value.E_peak_time_s << ',' << value.E_peak_r_m << ',' << value.E_peak_z_m << ',' << value.eta_peak << ','
         << value.ne_at_E_peak_m3 << ',' << value.ne_peak_m3 << ',' << value.nu_i_at_E_peak_s_1 << ','
         << value.mu_e_at_E_peak_m2_V_s << ',' << value.tau_i_at_E_peak_s << ',' << value.tau_i_at_ne_peak_s << ','
         << value.tau_i_min_s << ',' << value.tau_i_median_s << ',' << value.tau_i_p95_s << ','
         << value.sigma_e_at_E_peak_S_m << ',' << value.sigma_e_at_ne_peak_S_m << ',' << value.sigma_e_peak_S_m << ','
         << value.sigma_e_median_S_m << ',' << value.sigma_e_p95_S_m << ',' << value.tau_M_at_E_peak_s << ','
         << value.tau_M_at_ne_peak_s << ',' << value.tau_M_min_s << ',' << value.tau_M_median_s << ','
         << value.tau_M_p95_s << ',' << value.Pi_RF_at_E_peak << ',' << value.Pi_RF_at_ne_peak << ','
         << value.Pi_RF_min << ',' << value.Pi_RF_median << ',' << value.Pi_RF_p95 << ','
         << value.N_tau_i_at_E_peak << ',' << value.N_tau_M_at_E_peak << ','
         << value.tau_i_resolution_status << ',' << value.tau_M_resolution_status << ','
         << value.current_moment_r_A_m << ',' << value.current_moment_z_A_m << '\n';
}
}  // namespace

int main(int argc, char** argv) {
  Options options;
  try { options = parse(argc, argv); }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
  int petsc_argc = 1;
  char** petsc_argv = argv;
  PetscInitialize(&petsc_argc, &petsc_argv, nullptr, nullptr);
  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  const auto total_start = std::chrono::steady_clock::now();

  AxisymmetricGrid grid(12, 24, 80e-6, 0.0, 90e-6);
  AxisymmetricNeedlePlaneGeometry geometry("stage-c-r5-development-needle-plane", 0.0, 75e-6,
                                           5e-6, 2.5e-6, 5e-6);
  ConstantVoltage voltage(options.voltage_V);
  StreamerConfig config;
  config.electrode_geometry = &geometry;
  config.voltage_waveform = &voltage;
  config.photoionization = false;
  config.n_ref = 1e12;
  config.head_ne_threshold = 1e15;
  config.bridge_ne_threshold = 1e15;
  config.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver electrostatic(grid, config);
  electrostatic.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
  StreamerSolver solver(grid, config);
  solver.initialize_gaussian_at_tip_offset(1e16, 3e-6, -10e-6);
  StreamerHeadTracker head_tracker(StreamerHeadSegmentationConfig{0.2, 1});
  UltrafastEventDiagnostics diagnostics(grid, config, electrostatic.state().emag,
                                        UltrafastEventConfig{0.2, 1, 2, 3.0, 8.0, 1.5, 2});

  std::ofstream trace;
  if (!rank) {
    std::filesystem::create_directories(options.out);
    trace.open(options.out / "ultrafast_event_trace.csv");
    trace << std::setprecision(17);
    write_trace_header(trace);
  }
  int accepted = 0;
  double solver_seconds = 0.0, diagnostic_seconds = 0.0, io_seconds = 0.0;
  std::vector<double> accepted_dt;
  bool okay = true;
  for (int step = 1; step <= options.steps; ++step) {
    const auto limits = solver.timestep_limits();
    double dt = std::min(0.05 * limits.selected, options.dt_cap_s);
    StreamerDiagnostics solver_diagnostics;
    bool accepted_step = false;
    const auto solver_start = std::chrono::steady_clock::now();
    for (int retry = 0; retry <= 12; ++retry) {
      if (solver.step(dt, solver_diagnostics)) {
        accepted_step = true;
        break;
      }
      dt *= 0.5;
    }
    solver_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - solver_start).count();
    if (!accepted_step) {
      okay = false;
      break;
    }
    const auto diagnostic_start = std::chrono::steady_clock::now();
    const auto head = head_tracker.sample(grid, solver.state(), config);
    const auto current = solver.electron_transport_current_source();
    const auto sample = diagnostics.sample(step, solver.state(), head, current,
                                           solver_diagnostics.bridge_flag, dt);
    diagnostic_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - diagnostic_start).count();
    if (!rank) {
      const auto io_start = std::chrono::steady_clock::now();
      write_trace_row(trace, sample);
      trace.flush();
      io_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - io_start).count();
    }
    ++accepted;
    accepted_dt.push_back(dt);
  }
  if (!rank) trace.close();
  const double total_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
  struct rusage usage {};
  getrusage(RUSAGE_SELF, &usage);
  const auto trace_path = options.out / "ultrafast_event_trace.csv";
  const auto trace_size = !rank && std::filesystem::exists(trace_path) ? std::filesystem::file_size(trace_path) : 0;
  if (!rank) {
    const double dt_min = accepted_dt.empty() ? 0.0 : *std::min_element(accepted_dt.begin(), accepted_dt.end());
    const double dt_max = accepted_dt.empty() ? 0.0 : *std::max_element(accepted_dt.begin(), accepted_dt.end());
    std::vector<double> sorted_dt = accepted_dt;
    std::sort(sorted_dt.begin(), sorted_dt.end());
    const double dt_median = sorted_dt.empty() ? 0.0 : sorted_dt[sorted_dt.size() / 2];
    const std::string overhead_status = accepted >= 100 && solver_seconds >= 0.1
                                            ? ((diagnostic_seconds + io_seconds) / solver_seconds < 0.1
                                                   ? "PASS"
                                                   : "ABOVE_TARGET")
                                            : "OVERHEAD_NOT_RESOLVED_SHORT_RUN";
    std::ofstream resource(options.out / "resource_report.json");
    resource << std::setprecision(17)
             << "{\n  \"schema_version\": \"1.0\",\n  \"mode\": \"LIGHTWEIGHT_STAGE_C_DEVELOPMENT_REFERENCE\",\n"
             << "  \"steps_requested\": " << options.steps << ",\n  \"steps_accepted\": " << accepted
             << ",\n  \"physical_duration_s\": " << solver.state().time
             << ",\n  \"runtime_total_s\": " << total_seconds
             << ",\n  \"runtime_solver_s\": " << solver_seconds
             << ",\n  \"runtime_diagnostic_s\": " << diagnostic_seconds
             << ",\n  \"runtime_trace_io_s\": " << io_seconds
             << ",\n  \"diagnostic_plus_io_over_solver\": " << (diagnostic_seconds + io_seconds) / std::max(solver_seconds, 1e-300)
             << ",\n  \"overhead_status\": \"" << overhead_status
             << "\",\n  \"peak_RSS_kB\": " << usage.ru_maxrss
             << ",\n  \"trace_size_bytes\": " << trace_size
             << ",\n  \"history_2d_fields_retained\": false,\n  \"status\": \"" << (okay ? "PASS" : "FAILED_STEP") << "\"\n}\n";
    std::ofstream timing(options.out / "temporal_metadata.json");
    timing << std::setprecision(17)
           << "{\n  \"schema_version\": \"1.0\",\n  \"SOLVER_INTERNAL_DT\": \"adaptive timestep limit\",\n"
           << "  \"DIAGNOSTIC_TRACE_DT\": \"accepted solver-step spacing\",\n"
           << "  \"FIELD_SNAPSHOT_DT\": \"not changed by C-R5\",\n"
           << "  \"dt_min_s\": " << dt_min << ",\n  \"dt_median_s\": " << dt_median
           << ",\n  \"dt_max_s\": " << dt_max
           << ",\n  \"interpolation_used\": false,\n  \"pulse_width_resolution_claimed\": false\n}\n";
    std::cout << "stage_c_r5_diagnostics status=" << (okay ? "PASS" : "FAILED_STEP")
              << " accepted=" << accepted << " physical_duration=" << solver.state().time
              << " runtime=" << total_seconds << " trace_bytes=" << trace_size << '\n';
  }
  PetscFinalize();
  return okay && accepted == options.steps ? 0 : 1;
}
