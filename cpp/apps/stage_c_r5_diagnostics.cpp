#include "streamer_rf/streamer/HeadTracker.hpp"
#include "streamer_rf/streamer/UltrafastEventDiagnostics.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <sys/resource.h>
#include <algorithm>
#include <chrono>
#include <deque>
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
  double dt_scale{0.05};
  std::string waveform{"constant"};
  double pre_hold_s{2e-12};
  double ramp_time_s{5e-12};
  double trigger_settle_duration_s{};
  std::string trace_mode{"full"};
  bool stop_on_event_completion{};
  double pre_buffer_duration_s{20e-12};
  double post_tail_duration_s{2e-12};
  double tail_fraction{0.1};
  int post_tail_samples{64};
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
    else if (arg == "--dt-scale") options.dt_scale = std::stod(need("--dt-scale"));
    else if (arg == "--waveform") options.waveform = need("--waveform");
    else if (arg == "--pre-hold") options.pre_hold_s = std::stod(need("--pre-hold"));
    else if (arg == "--ramp-time") options.ramp_time_s = std::stod(need("--ramp-time"));
    else if (arg == "--trigger-settle-duration") options.trigger_settle_duration_s = std::stod(need("--trigger-settle-duration"));
    else if (arg == "--trace-mode") options.trace_mode = need("--trace-mode");
    else if (arg == "--stop-on-event-completion") options.stop_on_event_completion = std::stoi(need("--stop-on-event-completion")) != 0;
    else if (arg == "--pre-buffer-duration") options.pre_buffer_duration_s = std::stod(need("--pre-buffer-duration"));
    else if (arg == "--post-tail-duration") options.post_tail_duration_s = std::stod(need("--post-tail-duration"));
    else if (arg == "--tail-fraction") options.tail_fraction = std::stod(need("--tail-fraction"));
    else if (arg == "--post-tail-samples") options.post_tail_samples = std::stoi(need("--post-tail-samples"));
    else throw std::runtime_error("unknown option " + arg);
  }
  if (options.steps < 1 || !(options.dt_cap_s > 0.0) || !(options.dt_scale > 0.0) ||
      (options.waveform != "constant" && options.waveform != "ramp") ||
      (options.trace_mode != "full" && options.trace_mode != "targeted") ||
      !(options.pre_hold_s > 0.0) || !(options.ramp_time_s > 0.0) || options.trigger_settle_duration_s < 0.0 ||
      !(options.pre_buffer_duration_s > 0.0) || !(options.post_tail_duration_s > 0.0) ||
      !(options.tail_fraction > 0.0 && options.tail_fraction < 0.5) || options.post_tail_samples < 1)
    throw std::runtime_error("invalid dry-run options");
  return options;
}

void write_compact_header(std::ostream& stream) {
  stream << "step,time_s,dt_s,topology,topology_status,head_position_m,E_peak_V_m,eta_peak,"
            "ne_peak_m3,sigma_e_peak_S_m,tau_i_at_E_peak_s,tau_M_at_E_peak_s,Pi_RF_at_E_peak,"
            "K_ion_z_A_m_s,current_moment_z_A_m\n";
}

void write_compact_row(std::ostream& stream, const UltrafastEventSample& value) {
  stream << value.step << ',' << value.time_s << ',' << value.dt_s << ','
         << event_topology_name(value.topology) << ',' << value.topology_status << ','
         << value.head_position_m << ',' << value.E_peak_V_m << ',' << value.eta_peak << ','
         << value.ne_peak_m3 << ',' << value.sigma_e_peak_S_m << ','
         << value.tau_i_at_E_peak_s << ',' << value.tau_M_at_E_peak_s << ','
         << value.Pi_RF_at_E_peak << ',' << value.K_ion_z_A_m_s << ','
         << value.current_moment_z_A_m << '\n';
}

void write_trace_header(std::ostream& stream) {
  stream << "step,time_s,dt_s,topology,topology_confidence,topology_evidence,topology_status,head_count,"
            "head_position_m,head_velocity_m_s,bridge_status,roi_definition,roi_status,roi_cell_count,"
            "roi_volume_m3,roi_centroid_r_m,roi_centroid_z_m,roi_r_min_m,roi_r_max_m,roi_z_min_m,roi_z_max_m,"
            "electrode_distance_m,Ek_V_m,gas,pressure_Pa,temperature_K,Ek_semantics,E0_source,"
            "E0_max_V_m,E0_mean_V_m,E0_p95_V_m,eta_0_max,eta_0_mean,eta_0_p95,"
            "E_peak_V_m,E_peak_time_s,E_peak_r_m,E_peak_z_m,eta_peak,ne_at_E_peak_m3,ne_peak_m3,"
            "E_at_ne_peak_V_m,nu_i_at_E_peak_s_1,mu_e_at_E_peak_m2_V_s,nu_i_at_ne_peak_s_1,"
            "mu_e_at_ne_peak_m2_V_s,tau_i_at_E_peak_s,tau_i_at_ne_peak_s,tau_i_min_s,"
            "tau_i_median_s,tau_i_p95_s,sigma_e_at_E_peak_S_m,sigma_e_at_ne_peak_S_m,sigma_e_peak_S_m,"
            "sigma_e_median_S_m,sigma_e_p95_S_m,tau_M_at_E_peak_s,tau_M_at_ne_peak_s,tau_M_min_s,"
            "tau_M_median_s,tau_M_p95_s,Pi_RF_at_E_peak,Pi_RF_at_ne_peak,Pi_RF_min,Pi_RF_median,Pi_RF_p95,"
            "N_tau_i_at_E_peak,N_tau_M_at_E_peak,tau_i_resolution_status,tau_M_resolution_status,"
            "sigma_e_peak_r_m,sigma_e_peak_z_m,E_at_sigma_peak_V_m,ne_at_sigma_peak_m3,"
            "mu_e_at_sigma_peak_m2_V_s,nu_i_at_sigma_peak_s_1,tau_i_at_sigma_peak_s,"
            "tau_M_at_sigma_peak_s,Pi_RF_at_sigma_peak,E_roi_median_V_m,ne_roi_median_m3,"
            "mu_e_roi_median_m2_V_s,nu_i_roi_median_s_1,K_ion_z_A_m_s,K_ion_abs_A_m_s,"
            "signed_proxy_status,mechanism_proxy_role,"
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
         << value.ne_at_E_peak_m3 << ',' << value.ne_peak_m3 << ',' << value.E_at_ne_peak_V_m << ','
         << value.nu_i_at_E_peak_s_1 << ',' << value.mu_e_at_E_peak_m2_V_s << ','
         << value.nu_i_at_ne_peak_s_1 << ',' << value.mu_e_at_ne_peak_m2_V_s << ','
         << value.tau_i_at_E_peak_s << ',' << value.tau_i_at_ne_peak_s << ','
         << value.tau_i_min_s << ',' << value.tau_i_median_s << ',' << value.tau_i_p95_s << ','
         << value.sigma_e_at_E_peak_S_m << ',' << value.sigma_e_at_ne_peak_S_m << ',' << value.sigma_e_peak_S_m << ','
         << value.sigma_e_median_S_m << ',' << value.sigma_e_p95_S_m << ',' << value.tau_M_at_E_peak_s << ','
         << value.tau_M_at_ne_peak_s << ',' << value.tau_M_min_s << ',' << value.tau_M_median_s << ','
         << value.tau_M_p95_s << ',' << value.Pi_RF_at_E_peak << ',' << value.Pi_RF_at_ne_peak << ','
         << value.Pi_RF_min << ',' << value.Pi_RF_median << ',' << value.Pi_RF_p95 << ','
         << value.N_tau_i_at_E_peak << ',' << value.N_tau_M_at_E_peak << ','
         << value.tau_i_resolution_status << ',' << value.tau_M_resolution_status << ','
         << value.sigma_e_peak_r_m << ',' << value.sigma_e_peak_z_m << ','
         << value.E_at_sigma_peak_V_m << ',' << value.ne_at_sigma_peak_m3 << ','
         << value.mu_e_at_sigma_peak_m2_V_s << ',' << value.nu_i_at_sigma_peak_s_1 << ','
         << value.tau_i_at_sigma_peak_s << ',' << value.tau_M_at_sigma_peak_s << ','
         << value.Pi_RF_at_sigma_peak << ',' << value.E_roi_median_V_m << ','
         << value.ne_roi_median_m3 << ',' << value.mu_e_roi_median_m2_V_s << ','
         << value.nu_i_roi_median_s_1 << ',' << value.K_ion_z_A_m_s << ','
         << value.K_ion_abs_A_m_s << ',' << value.signed_proxy_status << ','
         << value.mechanism_proxy_role << ','
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
  SampledVoltage ramp_voltage(
      {0.0, options.pre_hold_s, options.pre_hold_s + options.ramp_time_s},
      {0.0, 0.0, options.voltage_V}, WaveformOutOfRangePolicy::HoldEndpoint);
  StreamerConfig config;
  config.electrode_geometry = &geometry;
  config.voltage_waveform = options.waveform == "ramp"
                                ? static_cast<const VoltageWaveform*>(&ramp_voltage)
                                : static_cast<const VoltageWaveform*>(&voltage);
  config.photoionization = false;
  config.n_ref = 1e12;
  config.head_ne_threshold = 1e15;
  config.bridge_ne_threshold = 1e15;
  config.elliptic = {1e-10, 1e-14, 20000};

  StreamerConfig electrostatic_config = config;
  electrostatic_config.voltage_waveform = &voltage;
  StreamerSolver electrostatic(grid, electrostatic_config);
  electrostatic.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
  StreamerSolver solver(grid, config);
  solver.initialize_gaussian_at_tip_offset(1e16, 3e-6, -10e-6);
  StreamerHeadTracker head_tracker(StreamerHeadSegmentationConfig{0.2, 1});
  UltrafastEventDiagnostics diagnostics(grid, config, electrostatic.state().emag,
                                        UltrafastEventConfig{0.2, 1, 2, 3.0, 8.0, 1.5, 2});

  std::ofstream trace, compact, event_trace;
  if (!rank) {
    std::filesystem::create_directories(options.out);
    if (options.trace_mode == "full") {
      trace.open(options.out / "ultrafast_event_trace.csv");
      trace << std::setprecision(17);
      write_trace_header(trace);
    } else {
      event_trace.open(options.out / "event_detail_trace.csv");
      event_trace << std::setprecision(17);
      write_trace_header(event_trace);
    }
    compact.open(options.out / "compact_trace.csv");
    compact << std::setprecision(17);
    write_compact_header(compact);
  }
  int accepted = 0;
  double solver_seconds = 0.0, diagnostic_seconds = 0.0, io_seconds = 0.0;
  std::vector<double> accepted_dt;
  std::deque<UltrafastEventSample> pre_event_buffer;
  std::vector<double> derivative_history;
  std::vector<double> derivative_time_history;
  double previous_time = std::numeric_limits<double>::quiet_NaN();
  double previous_moment = std::numeric_limits<double>::quiet_NaN();
  double baseline_sum = 0.0;
  int baseline_count = 0;
  bool event_triggered = false, event_complete = false, left_half = false, left_one_over_e = false;
  int primary_peak_derivative_index = -1, post_tail_count = 0;
  double primary_peak_abs = 0.0, primary_peak_time = std::numeric_limits<double>::quiet_NaN();
  double tail_start_time = std::numeric_limits<double>::quiet_NaN();
  std::string trigger_reason{"NOT_TRIGGERED"};
  bool okay = true;
  for (int step = 1; step <= options.steps; ++step) {
    const auto limits = solver.timestep_limits();
    double dt = std::min(options.dt_scale * limits.selected, options.dt_cap_s);
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
      write_compact_row(compact, sample);
      if (options.trace_mode == "full") write_trace_row(trace, sample);
      pre_event_buffer.push_back(sample);
      while (!pre_event_buffer.empty() &&
             sample.time_s - pre_event_buffer.front().time_s > options.pre_buffer_duration_s)
        pre_event_buffer.pop_front();

      if (std::isfinite(previous_time) && sample.time_s > previous_time) {
        const double derivative = (sample.current_moment_z_A_m - previous_moment) /
                                  (sample.time_s - previous_time);
        const double magnitude = std::abs(derivative);
        derivative_history.push_back(magnitude);
        derivative_time_history.push_back(sample.time_s);
        if (baseline_count < 32) {
          baseline_sum += magnitude;
          ++baseline_count;
        }
        if (!event_triggered && derivative_history.size() >= 3 && baseline_count >= 8) {
          const std::size_t n = derivative_history.size();
          const double before = derivative_history[n - 3];
          const double candidate = derivative_history[n - 2];
          const double after = derivative_history[n - 1];
          const double baseline = baseline_sum / baseline_count;
          const bool drive_transition_complete = options.waveform != "ramp" ||
                                                 sample.time_s > options.pre_hold_s + options.ramp_time_s +
                                                                     options.trigger_settle_duration_s;
          const bool kinetic_activity = sample.topology_status == "CONFIRMED" &&
                                        sample.roi_status == "VALID" &&
                                        std::isfinite(sample.ne_peak_m3) && sample.ne_peak_m3 > 0.0 &&
                                        std::isfinite(sample.sigma_e_peak_S_m) && sample.sigma_e_peak_S_m > 0.0;
          if (drive_transition_complete && kinetic_activity && candidate > before && candidate >= after &&
              candidate > 3.0 * std::max(baseline, 1e-300)) {
            event_triggered = true;
            primary_peak_derivative_index = static_cast<int>(n - 2);
            primary_peak_abs = candidate;
            primary_peak_time = derivative_time_history[n - 2];
            for (std::size_t index = 0; index + 1 < n - 2; ++index) {
              left_half = left_half || derivative_history[index] <= 0.5 * primary_peak_abs;
              left_one_over_e = left_one_over_e || derivative_history[index] <= primary_peak_abs / std::exp(1.0);
            }
            trigger_reason = "INTERIOR_DMDT_PEAK_WITH_KINETIC_ACTIVITY";
            if (options.trace_mode == "targeted")
              for (const auto& buffered : pre_event_buffer) write_trace_row(event_trace, buffered);
          }
        } else if (event_triggered) {
          if (options.trace_mode == "targeted") write_trace_row(event_trace, sample);
          if (magnitude > primary_peak_abs) {
            primary_peak_abs = magnitude;
            primary_peak_time = sample.time_s;
            primary_peak_derivative_index = static_cast<int>(derivative_history.size() - 1);
            left_half = left_one_over_e = false;
            for (std::size_t index = 0; index + 1 < derivative_history.size(); ++index) {
              left_half = left_half || derivative_history[index] <= 0.5 * primary_peak_abs;
              left_one_over_e = left_one_over_e || derivative_history[index] <= primary_peak_abs / std::exp(1.0);
            }
            post_tail_count = 0;
            tail_start_time = std::numeric_limits<double>::quiet_NaN();
          }
          if (magnitude < options.tail_fraction * primary_peak_abs) {
            if (!std::isfinite(tail_start_time)) tail_start_time = sample.time_s;
            ++post_tail_count;
          } else {
            post_tail_count = 0;
            tail_start_time = std::numeric_limits<double>::quiet_NaN();
          }
          const bool right_half = magnitude <= 0.5 * primary_peak_abs;
          const bool right_one_over_e = magnitude <= primary_peak_abs / std::exp(1.0);
          const bool sustained_tail = post_tail_count >= options.post_tail_samples &&
                                      std::isfinite(tail_start_time) &&
                                      sample.time_s - tail_start_time >= options.post_tail_duration_s;
          event_complete = primary_peak_derivative_index > 0 && left_half && left_one_over_e &&
                           right_half && right_one_over_e && sustained_tail;
        }
      }
      previous_time = sample.time_s;
      previous_moment = sample.current_moment_z_A_m;
      io_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - io_start).count();
    }
    ++accepted;
    accepted_dt.push_back(dt);
    if (options.stop_on_event_completion && event_complete) break;
  }
  if (!rank) {
    if (trace.is_open()) trace.close();
    compact.close();
    if (event_trace.is_open()) event_trace.close();
  }
  const double total_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
  struct rusage usage {};
  getrusage(RUSAGE_SELF, &usage);
  const auto trace_path = options.trace_mode == "full" ? options.out / "ultrafast_event_trace.csv"
                                                        : options.out / "event_detail_trace.csv";
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
             << "  \"waveform\": \"" << options.waveform << "\",\n"
             << "  \"trigger_settle_duration_s\": " << options.trigger_settle_duration_s << ",\n"
             << "  \"dt_scale\": " << options.dt_scale << ",\n  \"dt_cap_s\": " << options.dt_cap_s << ",\n"
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
    std::ofstream completion(options.out / "event_completion.json");
    completion << std::setprecision(17)
               << "{\n  \"schema_version\": \"1.0\",\n"
               << "  \"stop_policy\": \"PHYSICAL_EVENT_COMPLETION_WITH_MAX_STEPS_SAFETY_CAP\",\n"
               << "  \"triggered\": " << (event_triggered ? "true" : "false") << ",\n"
               << "  \"trigger_reason\": \"" << trigger_reason << "\",\n"
               << "  \"primary_peak_time_s\": ";
    if (std::isfinite(primary_peak_time)) completion << primary_peak_time;
    else completion << "null";
    completion << ",\n  \"primary_peak_abs_dMdt_A_m_s\": " << primary_peak_abs << ",\n"
               << "  \"left_half_crossing_available\": " << (left_half ? "true" : "false") << ",\n"
               << "  \"left_one_over_e_crossing_available\": " << (left_one_over_e ? "true" : "false") << ",\n"
               << "  \"post_tail_fraction\": " << options.tail_fraction << ",\n"
               << "  \"post_tail_samples_required\": " << options.post_tail_samples << ",\n"
               << "  \"post_tail_duration_s_required\": " << options.post_tail_duration_s << ",\n"
               << "  \"post_tail_samples_observed\": " << post_tail_count << ",\n"
               << "  \"event_complete\": " << (event_complete ? "true" : "false") << "\n}\n";
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
  return okay && (accepted == options.steps || (options.stop_on_event_completion && event_complete)) ? 0 : 1;
}
