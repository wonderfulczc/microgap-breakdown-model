#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/streamer/LfaAudit.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
constexpr double kb = 1.380649e-23;

struct Options {
  std::filesystem::path out{"results/stage_c2/development"};
  std::string case_id{"stage-c2-avalanche-control"};
  double voltage{3000.0};
  bool sp3{false};
  int max_steps{120};
  double dt_scale{0.001};
  double dt_cap{1e-16};
  double n0{1e16};
  double sigma{3.0e-6};
  double seed_offset{-10.0e-6};
  double head_threshold{1e15};
  double bridge_threshold{1e15};
  bool screen_only{false};
  bool debug_reaction_limit{false};
  bool lfa_audit{false};
  bool source_decomposition{false};
};

Options parse(int argc, char** argv) {
  Options o;
  if (argc > 1) o.out = argv[1];
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    auto need = [&](const char* name) -> const char* {
      if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
      return argv[++i];
    };
    if (a == "--case") o.case_id = need("--case");
    else if (a == "--voltage") o.voltage = std::stod(need("--voltage"));
    else if (a == "--sp3") o.sp3 = std::stoi(need("--sp3")) != 0;
    else if (a == "--steps") o.max_steps = std::stoi(need("--steps"));
    else if (a == "--dt-scale") o.dt_scale = std::stod(need("--dt-scale"));
    else if (a == "--dt-cap") o.dt_cap = std::stod(need("--dt-cap"));
    else if (a == "--n0") o.n0 = std::stod(need("--n0"));
    else if (a == "--sigma") o.sigma = std::stod(need("--sigma"));
    else if (a == "--seed-offset") o.seed_offset = std::stod(need("--seed-offset"));
    else if (a == "--head-threshold") o.head_threshold = std::stod(need("--head-threshold"));
    else if (a == "--bridge-threshold") o.bridge_threshold = std::stod(need("--bridge-threshold"));
    else if (a == "--screen-only") o.screen_only = true;
    else if (a == "--debug-reaction-limit") o.debug_reaction_limit = true;
    else if (a == "--lfa-audit") o.lfa_audit = true;
    else if (a == "--source-decomposition") o.source_decomposition = true;
    else throw std::runtime_error("unknown option " + a);
  }
  return o;
}

void write_fields(const std::filesystem::path& p, const AxisymmetricGrid& g, const StreamerState& s,
                  const AxisymmetricNeedlePlaneGeometry& geom) {
  std::ofstream f(p);
  f << "i,j,r_m,z_m,cell_type,ne_m_3,np_m_3,nn_m_3,rho_C_m_3,phi_V,Er_V_m,Ez_V_m,E_V_m,Sph_m_3_s_1\n"
    << std::setprecision(17);
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      const auto t = geom.classify(g, i, j);
      const char* name = t == ElectrodeCellType::Gas ? "gas" : (t == ElectrodeCellType::HighVoltageElectrode ? "hv" : "ground");
      f << i << ',' << j << ',' << g.r(i) << ',' << g.z(j) << ',' << name << ',' << s.ne(i, j) << ','
        << s.np(i, j) << ',' << s.nn(i, j) << ',' << s.rho(i, j) << ',' << s.phi(i, j) << ','
        << s.er(i, j) << ',' << s.ez(i, j) << ',' << s.emag(i, j) << ',' << s.sph(i, j) << '\n';
    }
  }
}

void write_axis_profile(const std::filesystem::path& p, const AxisymmetricGrid& g, const StreamerState& s) {
  std::ofstream f(p);
  f << "j,z_m,ne_m_3,np_m_3,nn_m_3,rho_C_m_3,phi_V,Ez_V_m,E_V_m,Sph_m_3_s_1\n" << std::setprecision(17);
  for (int j = 0; j < g.nz(); ++j) {
    f << j << ',' << g.z(j) << ',' << s.ne(0, j) << ',' << s.np(0, j) << ',' << s.nn(0, j)
      << ',' << s.rho(0, j) << ',' << s.phi(0, j) << ',' << s.ez(0, j) << ',' << s.emag(0, j)
      << ',' << s.sph(0, j) << '\n';
  }
}

double median(std::vector<double> v) {
  if (v.empty()) return 0.0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}
}

int main(int argc, char** argv) {
  Options opt;
  try {
    opt = parse(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 2;
  }
  int petsc_argc = 1;
  char** petsc_argv = argv;
  PetscInitialize(&petsc_argc, &petsc_argv, nullptr, nullptr);
  int rank = 0, size = 1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  const auto wall0 = std::chrono::steady_clock::now();

  AxisymmetricGrid g(20, 48, 80e-6, 0.0, 90e-6);
  AxisymmetricNeedlePlaneGeometry geom("stage-c2-dev-70um-needle-plane", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
  ConstantVoltage voltage(opt.voltage);
  StreamerConfig cfg;
  cfg.electrode_geometry = &geom;
  cfg.voltage_waveform = &voltage;
  cfg.photoionization = opt.sp3;
  cfg.n_ref = 1e12;
  cfg.head_ne_threshold = opt.head_threshold;
  cfg.bridge_ne_threshold = opt.bridge_threshold;
  cfg.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver electrostatic(g, cfg);
  electrostatic.initialize_gaussian_at_tip_offset(0.0, opt.sigma, opt.seed_offset);
  const double Ek = morrow_lowke_breakdown_field(cfg.neutral_density, cfg.pressure, cfg.temperature);
  const double gas_gap = geom.tip_z_m() - (geom.ground_z_m() + geom.ground_thickness_m());
  const double eavg = opt.voltage / gas_gap;
  const double initial_electrostatic_emax = *std::max_element(electrostatic.state().emag.values().begin(), electrostatic.state().emag.values().end());
  const double initial_eover_n_max = initial_electrostatic_emax / cfg.neutral_density * 1e21;

  StreamerSolver solver(g, cfg);
  solver.initialize_gaussian_at_tip_offset(opt.n0, opt.sigma, opt.seed_offset);

  if (!rank) {
    std::filesystem::create_directories(opt.out);
    std::ofstream meta(opt.out / "case_metadata.txt");
    meta << "case_id=" << opt.case_id << "\n"
         << "geometry_id=" << geom.metadata().geometry_id << "\n"
         << "gap_m=" << geom.metadata().gap_m << "\n"
         << "gas_gap_m=" << gas_gap << "\n"
         << "tip_radius_m=" << geom.metadata().tip_radius_m << "\n"
         << "polarity=positive\nhv_electrode_z=75e-6\nground_electrode_z=0\npropagation_direction=-z\n"
         << "gas=air_morrow_lowke_baseline\npressure_Pa=" << cfg.pressure << "\ntemperature_K=" << cfg.temperature << "\n"
         << "voltage_V=" << opt.voltage << "\nseed_n0_m_3=" << opt.n0 << "\nseed_sigma_m=" << opt.sigma
         << "\nseed_tip_offset_m=" << opt.seed_offset << "\nseed_center_z_m=" << geom.seed_z_from_tip_offset(opt.seed_offset) << "\n"
         << "seed_in_gas=" << (geom.classify(0.0, geom.seed_z_from_tip_offset(opt.seed_offset)) == ElectrodeCellType::Gas) << "\n"
         << "grid_nr=20\ngrid_nz=48\nrmax_m=8e-5\nzmax_m=9e-5\nphotoionization=" << opt.sp3 << "\n"
         << "head_ne_threshold_m_3=" << cfg.head_ne_threshold << "\nbridge_ne_threshold_m_3=" << cfg.bridge_ne_threshold << "\n"
         << "max_steps=" << opt.max_steps << "\ndt_scale=" << opt.dt_scale << "\ndt_cap_s=" << opt.dt_cap << "\n";
    if (opt.lfa_audit) meta << "lfa_audit=1\n";
    if (opt.source_decomposition) meta << "source_decomposition=1\n";
    meta << "Ek_V_m=" << Ek << "\nEavg_V_m=" << eavg << "\nzero_charge_Emax_V_m=" << initial_electrostatic_emax
         << "\nEavg_over_Ek=" << eavg / Ek << "\nzero_charge_Emax_over_Ek=" << initial_electrostatic_emax / Ek
         << "\nzero_charge_EoverN_max_Td=" << initial_eover_n_max << "\n";
    write_fields(opt.out / "fields_initial.csv", g, solver.state(), geom);
    write_axis_profile(opt.out / "axis_profile_initial.csv", g, solver.state());
  }
  if (opt.screen_only) {
    if (!rank) {
      std::ofstream summary(opt.out / "summary.txt");
      summary << std::setprecision(17) << "ranks=" << size << "\nstatus=SCREEN_ONLY\nEk=" << Ek
              << "\nEavg=" << eavg << "\nzero_charge_Emax=" << initial_electrostatic_emax
              << "\nEavg_over_Ek=" << eavg / Ek << "\nzero_charge_Emax_over_Ek=" << initial_electrostatic_emax / Ek << "\n";
      std::cout << "stage_c2_dynamic ranks=" << size << " status=SCREEN_ONLY voltage=" << opt.voltage
                << " Ek=" << Ek << " Eavg=" << eavg << " Emax=" << initial_electrostatic_emax
                << " Eavg_over_Ek=" << eavg / Ek << " Emax_over_Ek=" << initial_electrostatic_emax / Ek << '\n';
    }
    PetscFinalize();
    return 0;
  }

  std::ofstream diag;
  std::ofstream reaction_debug;
  std::ofstream lfa_csv;
  std::ofstream source_csv;
  if (!rank) {
    diag.open(opt.out / "diagnostics.csv");
    diag << "step,time,dt,dt_controller,voltage,Emax,EoverN_max_Td,ne_max,np_max,nn_max,total_electrons,total_charge,conservation_residual,sigma_max,head_position,head_velocity,bridge_flag,absorbed_electron_hv,absorbed_electron_ground,poisson_iterations,rejected_retries\n"
         << std::setprecision(17);
    if (opt.debug_reaction_limit) {
      reaction_debug.open(opt.out / "reaction_timestep_debug.csv");
      reaction_debug << "step,global_controller,i,j,r_m,z_m,cell_type,E_V_m,EoverN_Td,ne_m_3,np_m_3,nn_m_3,nu_ion_s_1,nu_att2_s_1,nu_att3_s_1,recombination_frequency_s_1,reaction_dt_s,ne_max_m_3,numerical_density_tolerance_m_3\n"
                     << std::setprecision(17);
    }
    if (opt.lfa_audit) {
      lfa_csv.open(opt.out / "lfa_audit.csv");
      lfa_csv << "step,time_s,dt_s,region,EoverN_max_Td,LE_min_valid_m,LE_p05_m,LE_median_m,LE_valid_fraction,"
                 "tauE_min_valid_s,tauE_p05_s,tauE_median_s,tauE_valid_fraction,gas_cells,selected_gas_cells,"
                 "LE_valid_cells,tauE_valid_cells,insufficient_stencil_cells,near_zero_gradE_cells,"
                 "near_zero_dEdt_cells,nonfinite_input_cells,chiL_p95,chiL_median,chiT_p95,chiT_median,"
                 "relaxation_coverage_fraction,relaxation_table_min_Td,relaxation_table_max_Td,"
                 "fraction_outside_relaxation_table,temporal_status,lfa_relaxation_data_status,"
                 "lfa_applicability\n"
              << std::setprecision(17);
    }
    if (opt.source_decomposition) {
      source_csv.open(opt.out / "reaction_source_diagnostics.csv");
      source_csv << "step,time_s,dt_s,impact_rate_s_1,photo_rate_s_1,attach2_rate_s_1,attach3_rate_s_1,"
                    "recomb_e_rate_s_1,net_electron_reaction_rate_s_1,source_closure_abs_s_1,"
                    "source_closure_rel,local_source_closure_max_abs_m_3_s_1,gas_cells,source_stage\n"
                 << std::setprecision(17);
    }
  }

  StreamerDiagnostics d;
  const auto lfa_regions = default_lfa_audit_regions();
  std::vector<LfaAuditHistory> lfa_histories(lfa_regions.size());
  std::vector<LfaAuditSummary> last_lfa_summaries(lfa_regions.size());
  ReactionSourceDecompositionDiagnostics last_source_decomp;
  double max_head_velocity = 0.0;
  double bridge_time = -1.0;
  double max_emax = initial_electrostatic_emax;
  double dt_min = std::numeric_limits<double>::infinity(), dt_max = 0.0;
  std::vector<double> accepted_dt;
  std::map<std::string, int> controller_counts;
  int snapshots = 0;
  bool ok = true;
  int accepted_steps = 0, total_retries = 0;
  for (int step = 1; step <= opt.max_steps; ++step) {
    auto lim = solver.timestep_limits();
    auto rlim = solver.reaction_timestep_diagnostic();
    if (!rank && opt.debug_reaction_limit) {
      reaction_debug << step << ',' << lim.controller << ',' << rlim.i << ',' << rlim.j << ',' << rlim.r << ',' << rlim.z << ','
                     << rlim.cell_classification << ',' << rlim.E << ',' << rlim.E_over_N_Td << ',' << rlim.ne << ','
                     << rlim.np << ',' << rlim.nn << ',' << rlim.ionization_frequency << ',' << rlim.attachment_two_body_frequency
                     << ',' << rlim.attachment_three_body_frequency << ',' << rlim.recombination_frequency << ','
                     << rlim.reaction_dt << ',' << rlim.ne_max << ',' << rlim.numerical_density_tolerance << '\n';
    }
    double dt = std::min(opt.dt_scale * lim.selected, opt.dt_cap);
    int retries = 0;
    ok = false;
    for (; retries <= 16; ++retries) {
      ok = solver.step(dt, d);
      if (ok) break;
      dt *= 0.5;
    }
    if (!ok || !std::isfinite(d.emax) || !std::isfinite(d.ne_max)) break;
    if (opt.lfa_audit) {
      for (std::size_t n = 0; n < lfa_regions.size(); ++n) {
        last_lfa_summaries[n] = lfa_histories[n].sample(g, solver.state(), cfg, d.dt, nullptr, nullptr, lfa_regions[n]);
      }
    }
    if (opt.source_decomposition) last_source_decomp = solver.reaction_source_decomposition_diagnostics();
    ++accepted_steps;
    total_retries += retries;
    controller_counts[lim.controller]++;
    dt_min = std::min(dt_min, dt);
    dt_max = std::max(dt_max, dt);
    accepted_dt.push_back(dt);
    max_emax = std::max(max_emax, d.emax);
    max_head_velocity = std::max(max_head_velocity, std::abs(d.head_velocity));
    if (d.bridge_flag && bridge_time < 0.0) bridge_time = d.time;
    if (!rank) {
      diag << step << ',' << d.time << ',' << d.dt << ',' << lim.controller << ',' << d.applied_voltage << ',' << d.emax << ','
           << d.emax / cfg.neutral_density * 1e21 << ','
           << d.ne_max << ','
           << d.np_max << ',' << d.nn_max << ',' << d.total_electrons << ',' << d.total_charge << ','
           << d.conservation_residual << ',' << d.sigma_max << ',' << d.head_position << ',' << d.head_velocity << ','
           << d.bridge_flag << ',' << d.absorbed_electron_hv << ',' << d.absorbed_electron_ground << ','
           << d.poisson_iterations << ',' << retries << '\n';
      if (opt.lfa_audit) {
        for (const auto& last_lfa : last_lfa_summaries) {
          lfa_csv << step << ',' << d.time << ',' << d.dt << ',' << last_lfa.region_name << ','
                  << last_lfa.eovern_max_Td << ',' << last_lfa.le_min_valid_m << ',' << last_lfa.le_p05_m << ','
                  << last_lfa.le_median_m << ',' << last_lfa.le_valid_fraction << ',' << last_lfa.tauE_min_valid_s
                  << ',' << last_lfa.tauE_p05_s << ',' << last_lfa.tauE_median_s << ','
                  << last_lfa.tauE_valid_fraction << ',' << last_lfa.gas_cells << ',' << last_lfa.selected_gas_cells
                  << ',' << last_lfa.le_valid_cells << ',' << last_lfa.tauE_valid_cells << ','
                  << last_lfa.insufficient_stencil_cells << ',' << last_lfa.near_zero_gradE_cells << ','
                  << last_lfa.near_zero_dEdt_cells << ',' << last_lfa.nonfinite_input_cells << ','
                  << last_lfa.chiL_p95 << ',' << last_lfa.chiL_median << ',' << last_lfa.chiT_p95 << ','
                  << last_lfa.chiT_median << ',' << last_lfa.relaxation_coverage_fraction << ','
                  << last_lfa.relaxation_table_min_Td << ',' << last_lfa.relaxation_table_max_Td << ','
                  << last_lfa.fraction_outside_relaxation_table << ',' << last_lfa.temporal_status << ','
                  << last_lfa.lfa_relaxation_data_status << ',' << last_lfa.lfa_applicability << '\n';
        }
      }
      if (opt.source_decomposition) {
        source_csv << step << ',' << d.time << ',' << d.dt << ',' << last_source_decomp.impact_rate_s1 << ','
                   << last_source_decomp.photo_rate_s1 << ',' << last_source_decomp.attach2_rate_s1 << ','
                   << last_source_decomp.attach3_rate_s1 << ',' << last_source_decomp.recomb_e_rate_s1 << ','
                   << last_source_decomp.net_electron_reaction_rate_s1 << ',' << last_source_decomp.source_closure_abs_s1
                   << ',' << last_source_decomp.source_closure_rel << ','
                   << last_source_decomp.local_source_closure_max_abs_m3s << ',' << last_source_decomp.gas_cells
                   << ',' << last_source_decomp.source_stage << '\n';
      }
      if ((step == 5 || step == 30 || step == 80 || d.bridge_flag) && snapshots < 4) {
        write_fields(opt.out / ("fields_snapshot_" + std::to_string(step) + ".csv"), g, solver.state(), geom);
        write_axis_profile(opt.out / ("axis_profile_" + std::to_string(step) + ".csv"), g, solver.state());
        ++snapshots;
      }
    }
    if (d.bridge_flag && step > 5) break;
  }
  if (!rank) {
    write_fields(opt.out / "fields_final.csv", g, solver.state(), geom);
    write_axis_profile(opt.out / "axis_profile_final.csv", g, solver.state());
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    std::ofstream summary(opt.out / "summary.txt");
    summary << std::setprecision(17)
            << "ranks=" << size << "\nstatus=" << (ok ? "PASS" : "FAILED_STEP") << "\n"
            << "case_id=" << opt.case_id << "\nfinal_physical_time=" << solver.state().time << "\ntotal_steps=" << accepted_steps
            << "\ndt_min=" << (accepted_dt.empty() ? 0.0 : dt_min) << "\ndt_max=" << dt_max
            << "\ndt_median=" << median(accepted_dt) << "\ntotal_retries=" << total_retries
            << "\ncontroller_drift=" << controller_counts["drift"] << "\ncontroller_diffusion=" << controller_counts["diffusion"]
            << "\ncontroller_ionization=" << controller_counts["ionization"] << "\ncontroller_dielectric=" << controller_counts["dielectric"]
            << "\ncontroller_reaction=" << controller_counts["reaction"]
            << "\ncontroller_drift_fraction=" << (accepted_steps ? static_cast<double>(controller_counts["drift"]) / accepted_steps : 0.0)
            << "\ncontroller_diffusion_fraction=" << (accepted_steps ? static_cast<double>(controller_counts["diffusion"]) / accepted_steps : 0.0)
            << "\ncontroller_ionization_fraction=" << (accepted_steps ? static_cast<double>(controller_counts["ionization"]) / accepted_steps : 0.0)
            << "\ncontroller_dielectric_fraction=" << (accepted_steps ? static_cast<double>(controller_counts["dielectric"]) / accepted_steps : 0.0)
            << "\ncontroller_reaction_fraction=" << (accepted_steps ? static_cast<double>(controller_counts["reaction"]) / accepted_steps : 0.0)
            << "\ninitial_Emax=" << initial_electrostatic_emax
            << "\nmaximum_Emax=" << max_emax << "\ninitial_EoverN_max_Td=" << initial_eover_n_max
            << "\nfinal_EoverN_max_Td=" << d.emax / cfg.neutral_density * 1e21
            << "\nEk=" << Ek << "\nEavg=" << eavg << "\nEavg_over_Ek=" << eavg / Ek
            << "\nzero_charge_Emax_over_Ek=" << initial_electrostatic_emax / Ek
            << "\nEmax=" << d.emax << "\nne_max=" << d.ne_max << "\ntotal_electrons=" << d.total_electrons
            << "\nsigma_max=" << d.sigma_max << "\n"
            << "head_position=" << d.head_position << "\nmaximum_head_velocity=" << max_head_velocity << "\n"
            << "bridge_flag=" << d.bridge_flag << "\nbridging_time=" << bridge_time << "\n"
            << "conservation_residual=" << d.conservation_residual << "\nelapsed_s=" << elapsed << "\n";
    if (opt.lfa_audit && accepted_steps > 0 && !last_lfa_summaries.empty()) {
      const auto& all = last_lfa_summaries[0];
      summary << "lfa_EoverN_max_Td=" << all.eovern_max_Td
              << "\nlfa_LE_p05_m=" << all.le_p05_m
              << "\nlfa_LE_median_m=" << all.le_median_m
              << "\nlfa_LE_valid_fraction=" << all.le_valid_fraction
              << "\nlfa_tauE_p05_s=" << all.tauE_p05_s
              << "\nlfa_tauE_median_s=" << all.tauE_median_s
              << "\nlfa_tauE_valid_fraction=" << all.tauE_valid_fraction
              << "\nlfa_relaxation_data_status=" << all.lfa_relaxation_data_status
              << "\nLFA_APPLICABILITY=" << all.lfa_applicability << "\n";
      for (const auto& region_summary : last_lfa_summaries) {
        summary << "lfa_region_" << region_summary.region_name << "_LE_p05_m=" << region_summary.le_p05_m
                << "\nlfa_region_" << region_summary.region_name << "_LE_median_m=" << region_summary.le_median_m
                << "\nlfa_region_" << region_summary.region_name << "_tauE_p05_s=" << region_summary.tauE_p05_s
                << "\nlfa_region_" << region_summary.region_name << "_tauE_median_s=" << region_summary.tauE_median_s
                << "\nlfa_region_" << region_summary.region_name << "_coverage_fraction="
                << region_summary.relaxation_coverage_fraction << "\n";
      }
    }
    if (opt.source_decomposition && accepted_steps > 0) {
      summary << "source_decomposition_stage=" << last_source_decomp.source_stage
              << "\nimpact_rate_s_1=" << last_source_decomp.impact_rate_s1
              << "\nphoto_rate_s_1=" << last_source_decomp.photo_rate_s1
              << "\nattach2_rate_s_1=" << last_source_decomp.attach2_rate_s1
              << "\nattach3_rate_s_1=" << last_source_decomp.attach3_rate_s1
              << "\nrecomb_e_rate_s_1=" << last_source_decomp.recomb_e_rate_s1
              << "\nnet_electron_reaction_rate_s_1=" << last_source_decomp.net_electron_reaction_rate_s1
              << "\nsource_closure_abs_s_1=" << last_source_decomp.source_closure_abs_s1
              << "\nsource_closure_rel=" << last_source_decomp.source_closure_rel
              << "\nlocal_source_closure_max_abs_m_3_s_1=" << last_source_decomp.local_source_closure_max_abs_m3s
              << "\nsource_decomposition_gas_cells=" << last_source_decomp.gas_cells << "\n";
    }
    std::cout << "stage_c2_dynamic ranks=" << size << " status=" << (ok ? "PASS" : "FAILED_STEP")
              << " case=" << opt.case_id << " final_time=" << solver.state().time
              << " steps=" << accepted_steps << " dt_min=" << (accepted_dt.empty() ? 0.0 : dt_min)
              << " Emax=" << d.emax << " ne_max=" << d.ne_max << " sigma_max=" << d.sigma_max
              << " head_position=" << d.head_position << " max_head_velocity=" << max_head_velocity
              << " bridge=" << d.bridge_flag << " bridge_time=" << bridge_time
              << " conservation_residual=" << d.conservation_residual << '\n';
  }
  PetscFinalize();
  return ok ? 0 : 1;
}
