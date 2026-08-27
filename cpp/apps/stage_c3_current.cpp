#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
struct Options {
  std::filesystem::path out{"results/stage_c3/current"};
  std::string case_id{"stage-c3-case-b"};
  std::string mode{"case-b"};
  double voltage{500.0};
  double ramp_time{5e-12};
  bool sp3{true};
  int max_steps{500};
  double dt_scale{0.1};
  double dt_cap{1e-14};
  double n0{1e16};
  double sigma{3e-6};
  double seed_offset{-10e-6};
};

Options parse(int argc, char** argv) {
  Options o;
  if (argc > 1) o.out = argv[1];
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* name) {
      if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
      return argv[++i];
    };
    if (a == "--case") o.case_id = need("--case");
    else if (a == "--mode") o.mode = need("--mode");
    else if (a == "--voltage") o.voltage = std::stod(need("--voltage"));
    else if (a == "--ramp-time") o.ramp_time = std::stod(need("--ramp-time"));
    else if (a == "--sp3") o.sp3 = std::stoi(need("--sp3")) != 0;
    else if (a == "--steps") o.max_steps = std::stoi(need("--steps"));
    else if (a == "--dt-scale") o.dt_scale = std::stod(need("--dt-scale"));
    else if (a == "--dt-cap") o.dt_cap = std::stod(need("--dt-cap"));
    else if (a == "--n0") o.n0 = std::stod(need("--n0"));
    else if (a == "--sigma") o.sigma = std::stod(need("--sigma"));
    else if (a == "--seed-offset") o.seed_offset = std::stod(need("--seed-offset"));
    else throw std::runtime_error("unknown option " + a);
  }
  if (o.mode == "vacuum-ramp" || o.mode == "vacuum-constant" || o.mode == "fixed-rho") {
    o.sp3 = false;
    o.n0 = 0.0;
  }
  return o;
}

double range_min(const std::vector<double>& v) { return v.empty() ? 0.0 : *std::min_element(v.begin(), v.end()); }
double range_max(const std::vector<double>& v) { return v.empty() ? 0.0 : *std::max_element(v.begin(), v.end()); }
}

int main(int argc, char** argv) {
  Options opt;
  try { opt = parse(argc, argv); } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
  int petsc_argc = 1;
  char** petsc_argv = argv;
  PetscInitialize(&petsc_argc, &petsc_argv, nullptr, nullptr);
  int rank = 0, size = 1;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  const auto wall0 = std::chrono::steady_clock::now();

  if (opt.mode == "uniform-conductance") {
    AxisymmetricGrid ug(24, 20, 120e-6, 0.0, 100e-6);
    StreamerConfig ucfg;
    ucfg.photoionization = false;
    StreamerSolver solver(ug, ucfg);
    solver.initialize_gaussian(0.0, 1e-6, 50e-6);
    const double ne = 1e15, efield = 2e6, length = ug.nz() * ug.dz();
    for (int j = 0; j < ug.nz(); ++j) {
      for (int i = 0; i < ug.nr(); ++i) {
        solver.state().ne(i, j) = ne;
        solver.state().np(i, j) = ne;
        solver.state().nn(i, j) = 0.0;
        solver.state().er(i, j) = 0.0;
        solver.state().ez(i, j) = efield;
        solver.state().emag(i, j) = efield;
      }
    }
    const double voltage = efield * length;
    auto q = evaluate_morrow_lowke(efield, ucfg.neutral_density, ucfg.pressure, ucfg.temperature);
    const double sigma = 1.602176634e-19 * q.mobility * ne;
    const double theoretical = sigma * 3.14159265358979323846 * std::pow(ug.nr() * ug.dr(), 2) / length;
    auto gd = solver.conductance_diagnostics(voltage);
    if (!rank) {
      std::filesystem::create_directories(opt.out);
      std::ofstream summary(opt.out / "summary.txt");
      summary << std::setprecision(17) << "ranks=" << size << "\nstatus=PASS\ncase_id=" << opt.case_id
              << "\nmode=uniform-conductance\nuniform_ne_m3=" << ne << "\nuniform_E_Vpm=" << efield
              << "\nvoltage_V=" << voltage << "\nsigma_Spm=" << sigma
              << "\nGb_theoretical_S=" << theoretical << "\nGb_simulated_S=" << gd.gb
              << "\nGb_relative_error=" << std::abs(gd.gb - theoretical) / theoretical
              << "\nRb_ohm=" << gd.rb << "\n";
      std::cout << "stage_c3_current ranks=" << size << " status=PASS case=" << opt.case_id
                << " mode=uniform-conductance Gb=" << gd.gb
                << " error=" << std::abs(gd.gb - theoretical) / theoretical << '\n';
    }
    PetscFinalize();
    return 0;
  }

  AxisymmetricGrid g(20, 48, 80e-6, 0.0, 90e-6);
  AxisymmetricNeedlePlaneGeometry geom("stage-c3-dev-70um-needle-plane", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
  ConstantVoltage constant_voltage(opt.voltage);
  SampledVoltage ramp_voltage({0.0, opt.ramp_time}, {0.0, opt.voltage}, WaveformOutOfRangePolicy::HoldEndpoint);
  const VoltageWaveform* voltage = opt.mode == "vacuum-ramp" ? static_cast<const VoltageWaveform*>(&ramp_voltage) : static_cast<const VoltageWaveform*>(&constant_voltage);

  StreamerConfig cfg;
  cfg.electrode_geometry = &geom;
  cfg.voltage_waveform = voltage;
  cfg.photoionization = opt.sp3;
  cfg.n_ref = 1e12;
  cfg.head_ne_threshold = 1e15;
  cfg.bridge_ne_threshold = 1e15;
  cfg.elliptic = {1e-10, 1e-14, 20000};

  StreamerSolver solver(g, cfg);
  solver.initialize_gaussian_at_tip_offset(opt.n0, opt.sigma, opt.seed_offset);
  const double cgap = solver.vacuum_gap_capacitance();
  const double dvdt = opt.mode == "vacuum-ramp" ? opt.voltage / opt.ramp_time : 0.0;

  if (opt.mode == "fixed-rho") {
    auto set_positive_gaussian = [&](double z0) {
      for (int j = 0; j < g.nz(); ++j) {
        for (int i = 0; i < g.nr(); ++i) {
          const bool gas = geom.classify(g, i, j) == ElectrodeCellType::Gas;
          const double rr = g.r(i);
          const double zz = g.z(j) - z0;
          const double n = gas ? 5e20 * std::exp(-(rr * rr + zz * zz) / std::pow(8e-6, 2)) : 0.0;
          solver.state().ne(i, j) = 0.0;
          solver.state().np(i, j) = n;
          solver.state().nn(i, j) = 0.0;
        }
      }
      solver.refresh_electrostatic_fields();
    };
    const double dt = opt.dt_cap > 0.0 ? opt.dt_cap : 1e-12;
    set_positive_gaussian(38e-6);
    const auto qa = solver.electrode_surface_diagnostics();
    solver.reset_electrode_history();
    set_positive_gaussian(54e-6);
    const auto qb = solver.electrode_surface_diagnostics();
    const StreamerDiagnostics sd = solver.sample_terminal_diagnostics_from_history(dt);
    const double expected = (qb.q_hv - qa.q_hv) / dt;
    const double err = std::abs(sd.i_disp_hv - expected);
    const double rel = err / std::max(std::abs(expected), 1e-300);
    const bool ok = std::isfinite(sd.i_disp_hv) && std::abs(qb.q_hv - qa.q_hv) > 0.0 && rel < 1e-12;
    if (!rank) {
      std::filesystem::create_directories(opt.out);
      std::ofstream summary(opt.out / "summary.txt");
      summary << std::setprecision(17)
              << "ranks=" << size << "\nstatus=" << (ok ? "PASS" : "FAIL") << "\ncase_id=" << opt.case_id
              << "\nmode=fixed-rho\nvoltage_V=" << opt.voltage << "\ndt_s=" << dt
              << "\nQ_HV_A_C=" << qa.q_hv << "\nQ_HV_B_C=" << qb.q_hv
              << "\nQ_ground_A_C=" << qa.q_ground << "\nQ_ground_B_C=" << qb.q_ground
              << "\nexpected_I_disp_HV_A=" << expected << "\ncomputed_I_disp_HV_A=" << sd.i_disp_hv
              << "\nabsolute_error_A=" << err << "\nrelative_error=" << rel
              << "\ntotal_plasma_charge_C=" << sd.total_charge << "\n";
      std::ofstream csv(opt.out / "fixed_rho_diagnostics.csv");
      csv << std::setprecision(17)
          << "state,time_s,Q_HV_C,Q_ground_C,I_disp_HV_A,I_disp_ground_A,total_plasma_charge_C\n"
          << "A," << solver.state().time << ',' << qa.q_hv << ',' << qa.q_ground << ",nan,nan,nan\n"
          << "B," << solver.state().time << ',' << qb.q_hv << ',' << qb.q_ground << ',' << sd.i_disp_hv << ','
          << sd.i_disp_ground << ',' << sd.total_charge << '\n';
      std::cout << "stage_c3_current ranks=" << size << " status=" << (ok ? "PASS" : "FAIL")
                << " case=" << opt.case_id << " mode=fixed-rho qdiff=" << (qb.q_hv - qa.q_hv)
                << " I_disp=" << sd.i_disp_hv << " error=" << rel << '\n';
    }
    PetscFinalize();
    return ok ? 0 : 1;
  }

  std::ofstream csv;
  if (!rank) {
    std::filesystem::create_directories(opt.out);
    csv.open(opt.out / "current_diagnostics.csv");
    csv << "time_s,voltage_V,I_cond_HV_A,I_disp_HV_A,I_total_HV_A,I_cond_ground_A,I_disp_ground_A,I_total_ground_A,Q_HV_C,Q_ground_C,total_plasma_charge_C,Emax_Vpm,ne_max_m3,sigma_max_Spm,P_cond_W,Gb_S,Rb_ohm,Rb_valid,head_position_m,bridge_flag,current_continuity_residual\n"
        << std::setprecision(17);
  }

  std::vector<double> i_cond_hv, i_disp_hv, i_total_hv, p_cond, gb, rb, continuity;
  StreamerDiagnostics d;
  bool ok = true;
  int accepted = 0, retries_total = 0;
  for (int step = 1; step <= opt.max_steps; ++step) {
    auto lim = solver.timestep_limits();
    double dt = std::min(opt.dt_scale * lim.selected, opt.dt_cap);
    if (opt.mode == "vacuum-ramp") dt = std::min(dt, opt.ramp_time / opt.max_steps);
    int retries = 0;
    ok = false;
    for (; retries <= 16; ++retries) {
      ok = solver.step(dt, d);
      if (ok) break;
      dt *= 0.5;
    }
    if (!ok || !std::isfinite(d.emax) || !std::isfinite(d.ne_max)) break;
    ++accepted;
    retries_total += retries;
    i_cond_hv.push_back(d.i_cond_hv);
    i_disp_hv.push_back(d.i_disp_hv);
    i_total_hv.push_back(d.i_total_hv);
    p_cond.push_back(d.p_cond);
    gb.push_back(d.gb);
    if (d.rb_valid && std::isfinite(d.rb)) rb.push_back(d.rb);
    continuity.push_back(d.current_continuity_residual);
    if (!rank) {
      csv << d.time << ',' << d.applied_voltage << ',' << d.i_cond_hv << ',' << d.i_disp_hv << ',' << d.i_total_hv << ','
          << d.i_cond_ground << ',' << d.i_disp_ground << ',' << d.i_total_ground << ',' << d.q_hv << ',' << d.q_ground << ','
          << d.total_charge << ',' << d.emax << ',' << d.ne_max << ',' << d.sigma_max << ',' << d.p_cond << ',' << d.gb << ',' << d.rb << ','
          << d.rb_valid << ',' << d.head_position << ',' << d.bridge_flag << ',' << d.current_continuity_residual << '\n';
    }
  }

  if (!rank) {
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    const double pred = cgap * dvdt;
    const double sim_disp = i_disp_hv.empty() ? 0.0 : i_disp_hv.back();
    std::ofstream summary(opt.out / "summary.txt");
    summary << std::setprecision(17)
            << "ranks=" << size << "\nstatus=" << (ok ? "PASS" : "FAILED_STEP")
            << "\ncase_id=" << opt.case_id << "\nmode=" << opt.mode << "\nsteps=" << accepted
            << "\nfinal_time_s=" << solver.state().time << "\nvoltage_V=" << d.applied_voltage
            << "\nC_gap_vacuum_F=" << cgap << "\ndVdt_V_s=" << dvdt
            << "\npredicted_CdVdt_A=" << pred << "\nlast_I_disp_HV_A=" << sim_disp
            << "\nramp_relative_error=" << (pred != 0.0 ? std::abs(sim_disp - pred) / std::abs(pred) : 0.0)
            << "\nI_cond_HV_min_A=" << range_min(i_cond_hv) << "\nI_cond_HV_max_A=" << range_max(i_cond_hv)
            << "\nI_disp_HV_min_A=" << range_min(i_disp_hv) << "\nI_disp_HV_max_A=" << range_max(i_disp_hv)
            << "\nI_total_HV_min_A=" << range_min(i_total_hv) << "\nI_total_HV_max_A=" << range_max(i_total_hv)
            << "\nP_cond_min_W=" << range_min(p_cond) << "\nP_cond_max_W=" << range_max(p_cond)
            << "\nGb_min_S=" << range_min(gb) << "\nGb_max_S=" << range_max(gb)
            << "\nRb_min_ohm=" << range_min(rb) << "\nRb_max_ohm=" << range_max(rb)
            << "\nEmax_Vpm=" << d.emax << "\nne_max_m3=" << d.ne_max << "\nsigma_max_Spm=" << d.sigma_max
            << "\nbridge_flag=" << d.bridge_flag << "\nconservation_residual=" << d.conservation_residual
            << "\ncurrent_continuity_residual_max=" << range_max(continuity)
            << "\ntotal_retries=" << retries_total << "\nelapsed_s=" << elapsed << "\n";
    std::cout << "stage_c3_current ranks=" << size << " status=" << (ok ? "PASS" : "FAILED_STEP")
              << " case=" << opt.case_id << " mode=" << opt.mode << " final_time=" << solver.state().time
              << " steps=" << accepted << " I_cond_HV_max=" << range_max(i_cond_hv)
              << " I_disp_HV_max=" << range_max(i_disp_hv) << " Gb_max=" << range_max(gb)
              << " bridge=" << d.bridge_flag << " continuity=" << range_max(continuity) << '\n';
  }
  PetscFinalize();
  return ok ? 0 : 1;
}
