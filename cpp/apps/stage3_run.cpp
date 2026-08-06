#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <petscsys.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
constexpr double qe = 1.602176634e-19;
constexpr double kb = 1.380649e-23;
constexpr double pi = 3.14159265358979323846;

void write_fields(const std::filesystem::path& p, const AxisymmetricGrid& g, const StreamerState& s) {
  std::ofstream f(p);
  f << "i,j,r_m,z_m,ne_m_3,np_m_3,nn_m_3,rho_C_m_3,phi_V,Er_V_m,Ez_V_m,E_V_m,Sph_m_3_s_1\n"
    << std::setprecision(17);
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      f << i << ',' << j << ',' << g.r(i) << ',' << g.z(j) << ',' << s.ne(i, j) << ','
        << s.np(i, j) << ',' << s.nn(i, j) << ',' << s.rho(i, j) << ',' << s.phi(i, j) << ','
        << s.er(i, j) << ',' << s.ez(i, j) << ',' << s.emag(i, j) << ',' << s.sph(i, j) << '\n';
    }
  }
}

double median(std::vector<double> v) {
  if (v.empty()) return 0.0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}

struct AxisHeads {
  int lower_charge{1}, upper_charge{1}, lower_field{1}, upper_field{1};
};

AxisHeads axis_heads(const AxisymmetricGrid& g, const StreamerState& s) {
  const int mid = g.nz() / 2;
  AxisHeads h;
  double alo = 0, ahi = 0, elo = 0, ehi = 0;
  for (int j = 1; j < mid; ++j) {
    if (std::abs(s.rho(0, j)) > alo) {
      alo = std::abs(s.rho(0, j));
      h.lower_charge = j;
    }
    if (s.emag(0, j) > elo) {
      elo = s.emag(0, j);
      h.lower_field = j;
    }
  }
  for (int j = mid; j < g.nz() - 1; ++j) {
    if (std::abs(s.rho(0, j)) > ahi) {
      ahi = std::abs(s.rho(0, j));
      h.upper_charge = j;
    }
    if (s.emag(0, j) > ehi) {
      ehi = s.emag(0, j);
      h.upper_field = j;
    }
  }
  return h;
}

void write_ledger_row(std::ofstream& ledger, const AxisymmetricGrid& g, const StreamerState& s,
                      const StreamerConfig& c, const StreamerDiagnostics& d) {
  std::vector<double> td;
  td.reserve(s.emag.values().size());
  double td_min = std::numeric_limits<double>::infinity();
  double td_max = 0.0;
  double mu_min = std::numeric_limits<double>::infinity(), mu_max = 0.0;
  double D_min = std::numeric_limits<double>::infinity(), D_max = 0.0;
  double nui_max = 0.0, nua2_max = 0.0, nua3_max = 0.0;
  double Si = 0.0, Sa2 = 0.0, Sa3 = 0.0, Sph = 0.0, Sep = 0.0, Snp = 0.0;
  double Re = 0.0, Rp = 0.0, Rn = 0.0, charge = 0.0;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      const double vol = g.cell_volume(i);
      const auto q = evaluate_morrow_lowke(s.emag(i, j), c.neutral_density, c.pressure, c.temperature);
      const auto r = evaluate_reactions(s.ne(i, j), s.np(i, j), s.nn(i, j), s.sph(i, j), q, c.temperature);
      const double x = s.emag(i, j) / c.neutral_density * 1e21;
      td.push_back(x);
      td_min = std::min(td_min, x);
      td_max = std::max(td_max, x);
      mu_min = std::min(mu_min, q.mobility);
      mu_max = std::max(mu_max, q.mobility);
      D_min = std::min(D_min, q.diffusion);
      D_max = std::max(D_max, q.diffusion);
      nui_max = std::max(nui_max, q.ionization_frequency);
      nua2_max = std::max(nua2_max, q.attachment_two_body_frequency);
      nua3_max = std::max(nua3_max, q.attachment_three_body_frequency);
      Si += r.ionization * vol;
      Sa2 += r.attachment_two_body * vol;
      Sa3 += r.attachment_three_body * vol;
      Sph += s.sph(i, j) * vol;
      Sep += r.electron_positive_recombination * vol;
      Snp += r.ion_recombination * vol;
      Re += r.electron * vol;
      Rp += r.positive_ion * vol;
      Rn += r.negative_ion * vol;
      charge += qe * (r.positive_ion - r.electron - r.negative_ion) * vol;
    }
  }
  ledger << d.time << ',' << td_min << ',' << median(td) << ',' << td_max << ',' << mu_min << ','
         << mu_max << ',' << D_min << ',' << D_max << ',' << nui_max << ',' << nua2_max << ','
         << nua3_max << ',' << Si << ',' << Sa2 << ',' << Sa3 << ',' << Sph << ',' << Sep << ','
         << Snp << ',' << 0.0 << ',' << Re << ',' << Sph << ',' << Rp << ',' << Rn << ','
         << charge << ',' << d.poisson_iterations << ',' << d.sp3_ksp_iterations << ','
         << d.sp3_boundary_iterations << '\n';
}
}

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rank, size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  const auto wall0 = std::chrono::steady_clock::now();

  if (argc < 9) {
    if (!rank) {
      std::cerr << "usage: output nr nz tend E0 nref sp3 eta [rmax] [zmax] [dt_scale] [n0] [sigma] [max_steps] [physical_stop_after_s] [resume_checkpoint] [initial_step] [minimum_dt_s] [output_interval_s]\n";
    }
    PetscFinalize();
    return 2;
  }

  const std::filesystem::path out = argv[1];
  const int nr = std::stoi(argv[2]), nz = std::stoi(argv[3]);
  const double tend = std::stod(argv[4]), E0 = std::stod(argv[5]), nref = std::stod(argv[6]);
  const bool sp3 = std::stoi(argv[7]);
  const double eta = std::stod(argv[8]);
  const double rmax = argc > 9 ? std::stod(argv[9]) : 1e-3;
  const double zmax = argc > 10 ? std::stod(argv[10]) : 1e-2;
  const double dt_scale = argc > 11 ? std::stod(argv[11]) : 1.0;
  const double n0 = argc > 12 ? std::stod(argv[12]) : 1e20;
  const double sigma = argc > 13 ? std::stod(argv[13]) : 1e-4;
  const int max_steps = argc > 14 ? std::stoi(argv[14]) : 50000;
  const double physical_stop_after = argc > 15 ? std::stod(argv[15]) : 0.0;
  const std::filesystem::path resume_checkpoint = argc > 16 ? std::filesystem::path(argv[16]) : std::filesystem::path();
  const int initial_step = argc > 17 ? std::stoi(argv[17]) : 0;
  const double minimum_dt = argc > 18 ? std::stod(argv[18]) : 1e-18;
  const double output_interval = argc > 19 ? std::stod(argv[19]) : 1e-10;

  AxisymmetricGrid g(nr, nz, rmax, 0, zmax);
  StreamerConfig c;
  c.background_field = E0;
  c.n_ref = nref;
  c.photoionization = sp3;
  c.open_boundary.mode = RingQuadratureMode::CellCenterRing;
  c.open_boundary.threshold = eta;
  c.elliptic = {1e-10, 1e-14, 20000};
  StreamerSolver solver(g, c);
  solver.initialize_gaussian(n0, sigma, .5 * zmax);
  if (!resume_checkpoint.empty()) {
    solver.load_checkpoint(resume_checkpoint);
  }

  if (!rank) {
    std::filesystem::create_directories(out);
    write_fields(out / ("fields_" + std::to_string(initial_step) + ".csv"), g, solver.state());
    solver.save_checkpoint(out / "checkpoint.bin");
  }
  MPI_Barrier(PETSC_COMM_WORLD);

  std::ofstream hist, ledger;
  if (!rank) {
    hist.open(out / "scalar_history.csv");
    hist << "step,time_s,dt_s,controller,E_max_V_m,ne_max_m_3,np_max_m_3,nn_max_m_3,total_electrons,total_charge_C,charge_source_residual,poisson_iterations,sp3_ksp_iterations,sp3_boundary_iterations,lower_head_z_m,upper_head_z_m,lower_head_rho_C_m_3,upper_head_rho_C_m_3,lower_field_head_z_m,upper_field_head_z_m,channel_ne_m_3,channel_charge_fraction,sph_ahead_m_3_s_1,formation_candidate,rejected_steps,negative_density_corrections,particle_balance_residual,charge_balance_residual\n"
         << std::setprecision(17);
    ledger.open(out / "coupling_ledger.csv");
    ledger << "time,E_over_N_min,E_over_N_median,E_over_N_max,mu_min,mu_max,D_min,D_max,nu_i_max,nu_a2_max,nu_a3_max,integral_S_i,integral_S_a2,integral_S_a3,integral_S_ph,integral_S_ep,integral_S_np,electron_transport_change,electron_reaction_change,electron_photoionization_change,positive_ion_change,negative_ion_change,total_charge_change,Poisson_residual,SP3_residual,SP3_boundary_iterations\n"
           << std::setprecision(17);
  }

  int step = initial_step, sustained = 0, rejected_steps = 0, negative_corrections = 0;
  bool formed = false;
  double formed_time = -1.0, last_dt = 0.0;
  std::string termination = "unknown";
  const auto initial_heads = axis_heads(g, solver.state());
  double initial_lower = g.z(initial_heads.lower_charge);
  double initial_upper = g.z(initial_heads.upper_charge);
  double initial_field_lower = g.z(initial_heads.lower_field);
  double initial_field_upper = g.z(initial_heads.upper_field);
  double next_output_time = (std::floor(solver.state().time / output_interval) + 1.0) * output_interval;
  StreamerDiagnostics last_diag{};

  while (solver.state().time < tend && step < max_steps) {
    auto lim = solver.timestep_limits();
    double dt = std::min(dt_scale * lim.selected, tend - solver.state().time);
    if (dt < minimum_dt) {
      termination = "minimum_dt";
      last_dt = dt;
      break;
    }
    StreamerDiagnostics d;
    bool accepted = false;
    int local_rejects = 0;
    for (int retry = 0; retry <= 20 && !accepted; ++retry) {
      accepted = solver.step(dt, d);
      if (!accepted) {
        dt *= .5;
        ++local_rejects;
      }
    }
    rejected_steps += local_rejects;
    negative_corrections += local_rejects;
    if (!accepted) {
      termination = "negativity_retry_exhausted";
      if (!rank) std::cerr << "negative-density retry exhausted\n";
      break;
    }
    ++step;
    last_dt = dt;
    last_diag = d;
    if (dt < minimum_dt) {
      termination = "minimum_dt";
    }

    const auto& s = solver.state();
    const int mid = nz / 2;
    const auto h = axis_heads(g, s);
    const double channel_ne = s.ne(0, mid);
    const double channel_fraction =
        std::abs(s.rho(0, mid)) / (qe * std::max(s.ne(0, mid) + s.np(0, mid) + s.nn(0, mid), 1.0));
    const double sph_ahead = s.sph(0, std::min(nz - 1, h.upper_charge + 2));
    const bool candidate =
        s.rho(0, h.lower_charge) * s.rho(0, h.upper_charge) < 0 &&
        initial_lower - g.z(h.lower_charge) > 3 * g.dz() &&
        g.z(h.upper_charge) - initial_upper > 3 * g.dz() &&
        initial_field_lower - g.z(h.lower_field) > 3 * g.dz() &&
        g.z(h.upper_field) - initial_field_upper > 3 * g.dz() &&
        channel_ne > 1e18 && channel_fraction < .1 && (!sp3 || sph_ahead > 0);
    sustained = candidate ? sustained + 1 : 0;
    if (!formed && sustained >= 25) {
      formed = true;
      formed_time = d.time;
    }

    if (!rank) {
      hist << step << ',' << d.time << ',' << d.dt << ',' << lim.controller << ',' << d.emax << ','
           << d.ne_max << ',' << d.np_max << ',' << d.nn_max << ',' << d.total_electrons << ','
           << d.total_charge << ',' << d.conservation_residual << ',' << d.poisson_iterations << ','
           << d.sp3_ksp_iterations << ',' << d.sp3_boundary_iterations << ',' << g.z(h.lower_charge)
           << ',' << g.z(h.upper_charge) << ',' << s.rho(0, h.lower_charge) << ','
           << s.rho(0, h.upper_charge) << ',' << g.z(h.lower_field) << ',' << g.z(h.upper_field)
           << ',' << channel_ne << ',' << channel_fraction << ',' << sph_ahead << ',' << candidate
           << ',' << rejected_steps << ',' << negative_corrections << ',' << d.conservation_residual
           << ',' << d.conservation_residual << '\n';
    }

    const bool output_due = solver.state().time + 1e-18 >= next_output_time || solver.state().time >= tend;
    if (!rank && output_due) {
      const std::string tag = "fields_" + std::to_string(step);
      write_fields(out / (tag + ".csv"), g, s);
      solver.save_checkpoint(out / "checkpoint.bin");
      write_ledger_row(ledger, g, s, c, d);
      ledger.flush();
      hist.flush();
      while (next_output_time <= solver.state().time + 1e-18) next_output_time += output_interval;
    }
    if (formed && physical_stop_after > 0.0 && solver.state().time - formed_time >= physical_stop_after) {
      termination = "reached_physical_stop";
      break;
    }
    if (!std::isfinite(solver.state().time) || !std::isfinite(d.emax) || !std::isfinite(d.ne_max)) {
      termination = "invalid_state";
      break;
    }
    if (termination == "minimum_dt") break;
  }

  if (termination == "unknown") {
    if (solver.state().time >= tend) termination = "reached_end_time";
    else if (step >= max_steps) termination = "max_steps";
  }
  const auto wall1 = std::chrono::steady_clock::now();
  const double wall_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(wall1 - wall0).count();

  if (!rank) {
    write_fields(out / "fields_final.csv", g, solver.state());
    solver.save_checkpoint(out / "checkpoint.bin");
    std::ofstream reason(out / "termination_reason.txt");
    reason << termination << '\n';
    std::ofstream meta(out / "termination.csv");
    meta << "termination_reason,requested_end_time,actual_end_time,completed_steps,accepted_steps,rejected_steps,last_dt,wall_time,checkpoint_path,formed,formation_time\n";
    meta << termination << ',' << tend << ',' << solver.state().time << ',' << (step - initial_step) << ',' << (step - initial_step) << ','
         << rejected_steps << ',' << last_dt << ',' << wall_time << ',' << (out / "checkpoint.bin").string()
         << ',' << formed << ',' << formed_time << '\n';
    if (ledger.is_open()) write_ledger_row(ledger, g, solver.state(), c, last_diag);
    std::cout << "steps=" << step << " final_time_s=" << solver.state().time << " ranks=" << size
              << " formed=" << formed << " termination_reason=" << termination << '\n';
  }

  PetscFinalize();
  if (termination == "negativity_retry_exhausted") return 3;
  if (termination == "invalid_state") return 5;
  if (termination == "minimum_dt") return 6;
  return (termination == "reached_end_time" || termination == "reached_physical_stop") ? 0 : 4;
}
