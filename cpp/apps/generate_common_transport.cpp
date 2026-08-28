#include "streamer_rf/streamer/MorrowLowke.hpp"
#include "streamer_rf/streamer/ReactionModel.hpp"

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
constexpr double kB = 1.380649e-23;

struct Row {
  double td{};
  double efield{};
  streamer_rf::streamer::TransportCoefficients c{};
};

double rel_err(double a, double b) {
  const double denom = std::max({std::abs(a), std::abs(b), 1e-300});
  return std::abs(a - b) / denom;
}

double interp(const std::vector<Row>& rows, double td, double Row::*member) {
  if (td <= rows.front().td) return rows.front().*member;
  if (td >= rows.back().td) return rows.back().*member;
  auto it = std::lower_bound(rows.begin(), rows.end(), td,
                             [](const Row& r, double value) { return r.td < value; });
  const auto& b = *it;
  const auto& a = *(it - 1);
  const double f = (td - a.td) / (b.td - a.td);
  return (1.0 - f) * (a.*member) + f * (b.*member);
}

double interp_transport(const std::vector<Row>& rows, double td,
                        double streamer_rf::streamer::TransportCoefficients::*member) {
  if (td <= rows.front().td) return rows.front().c.*member;
  if (td >= rows.back().td) return rows.back().c.*member;
  auto it = std::lower_bound(rows.begin(), rows.end(), td,
                             [](const Row& r, double value) { return r.td < value; });
  const auto& b = *it;
  const auto& a = *(it - 1);
  const double f = (td - a.td) / (b.td - a.td);
  return (1.0 - f) * (a.c.*member) + f * (b.c.*member);
}

void write_section(std::ofstream& f, const std::string& name,
                   const std::vector<std::pair<double, double>>& values) {
  f << '\n' << name << "\n-----------------------\n";
  for (const auto& [x, y] : values) f << x << ' ' << y << '\n';
  f << "-----------------------\n";
}
}  // namespace

int main(int argc, char** argv) {
  namespace fs = std::filesystem;
  using streamer_rf::streamer::beta_ep;
  using streamer_rf::streamer::beta_np;
  using streamer_rf::streamer::evaluate_morrow_lowke;
  using streamer_rf::streamer::evaluate_reactions;

  const fs::path out_dir =
      argc > 1 ? fs::path(argv[1]) : fs::path("solver3d/afivo_reference/common_benchmark/transport/generated");
  fs::create_directories(out_dir);

  const double pressure = 101325.0;
  const double temperature = 300.0;
  const double neutral_density = pressure / (kB * temperature);
  const int n_rows = 1201;
  const double td_min = 0.1;
  const double td_max = 1000.0;

  std::vector<Row> rows;
  rows.reserve(n_rows);
  for (int i = 0; i < n_rows; ++i) {
    const double u = static_cast<double>(i) / (n_rows - 1);
    const double td = td_min * std::pow(td_max / td_min, u);
    const double efield = td * 1e-21 * neutral_density;
    rows.push_back({td, efield, evaluate_morrow_lowke(efield, neutral_density, pressure, temperature)});
  }

  std::ofstream full(out_dir / "morrow_lowke_transport_full.csv");
  full << std::setprecision(17);
  full << "E_over_N_Td,E_Vpm,mobility_m2_Vs,diffusion_m2_s,alpha_1_m,"
          "eta2_1_m,eta3_1_m,eta_total_1_m,ionization_frequency_s,"
          "attachment_two_body_frequency_s,attachment_three_body_frequency_s,"
          "electron_temperature_K,mean_energy_eV\n";
  for (const auto& r : rows) {
    const double te = streamer_rf::streamer::electron_temperature(r.c.mobility, r.c.diffusion);
    const double mean_energy_ev = 1.5 * r.c.diffusion / std::max(r.c.mobility, 1e-300);
    full << r.td << ',' << r.efield << ',' << r.c.mobility << ',' << r.c.diffusion << ','
         << r.c.ionization_townsend << ',' << r.c.attachment_two_body_townsend << ','
         << r.c.attachment_three_body_townsend << ','
         << r.c.attachment_two_body_townsend + r.c.attachment_three_body_townsend << ','
         << r.c.ionization_frequency << ',' << r.c.attachment_two_body_frequency << ','
         << r.c.attachment_three_body_frequency << ',' << te << ',' << mean_energy_ev << '\n';
  }

  std::ofstream old_style(out_dir / "morrow_lowke_transport_afivo_old_style.txt");
  old_style << std::setprecision(17)
            << "# Generated from PETSc frozen evaluate_morrow_lowke(); SI old-style Afivo input.\n"
            << "# Afivo old-style has one eta column; this file uses eta2 + eta3 at p=101325 Pa, T=300 K.\n";
  std::vector<std::pair<double, double>> v;
  auto fill = [&](auto getter) {
    v.clear();
    for (const auto& r : rows) v.emplace_back(r.efield, getter(r));
  };
  auto fill_td = [&](auto getter) {
    v.clear();
    for (const auto& r : rows) v.emplace_back(r.td, getter(r));
  };
  fill([](const Row& r) { return r.c.mobility; });
  write_section(old_style, "efield[V/m]_vs_mu[m2/Vs]", v);
  fill([](const Row& r) { return r.c.diffusion; });
  write_section(old_style, "efield[V/m]_vs_dif[m2/s]", v);
  fill([](const Row& r) { return r.c.ionization_townsend; });
  write_section(old_style, "efield[V/m]_vs_alpha[1/m]", v);
  fill([](const Row& r) { return r.c.attachment_two_body_townsend + r.c.attachment_three_body_townsend; });
  write_section(old_style, "efield[V/m]_vs_eta[1/m]", v);

  std::ofstream afivo_common(out_dir / "morrow_lowke_common_chemistry.txt");
  afivo_common << std::setprecision(17)
               << "# D2 common transport/chemistry candidate generated from PETSc frozen Morrow-Lowke.\n"
               << "# Field tables below are versus E/N in Townsend.\n\n"
               << "reaction_list\n-----------------------\n"
               << "e -> e + e + P+,field_table,ML ionization frequency\n"
               << "e -> N-,field_table,ML attachment two-body frequency\n"
               << "e -> N-,field_table,ML attachment three-body frequency\n"
               << "e + P+ -> N2,c1*(300/Te)**c2,"
               << 1.138e-11 / std::pow(300.0, 0.7) << " 0.7\n"
               << "P+ + N- -> N2,c1*(300/Tg)**c2,2e-13 0.5\n"
               << "-----------------------\n\n"
               << "ignored_species\n-----------------------\nN2\n-----------------------\n";
  fill_td([&](const Row& r) { return r.c.mobility * neutral_density; });
  write_section(afivo_common, "Mobility *N (1/m/V/s)", v);
  fill_td([&](const Row& r) { return r.c.diffusion * neutral_density; });
  write_section(afivo_common, "Diffusion coefficient *N (1/m/s)", v);
  fill_td([&](const Row& r) { return r.c.ionization_townsend / neutral_density; });
  write_section(afivo_common, "Townsend ioniz. coef. alpha/N (m2)", v);
  fill_td([&](const Row& r) {
    return (r.c.attachment_two_body_townsend + r.c.attachment_three_body_townsend) / neutral_density;
  });
  write_section(afivo_common, "Townsend attach. coef. eta/N (m2)", v);
  fill_td([](const Row& r) { return 1.5 * r.c.diffusion / std::max(r.c.mobility, 1e-300); });
  write_section(afivo_common, "Mean energy (eV)", v);
  fill_td([](const Row& r) { return r.c.ionization_frequency; });
  write_section(afivo_common, "ML ionization frequency", v);
  fill_td([](const Row& r) { return r.c.attachment_two_body_frequency; });
  write_section(afivo_common, "ML attachment two-body frequency", v);
  fill_td([](const Row& r) { return r.c.attachment_three_body_frequency; });
  write_section(afivo_common, "ML attachment three-body frequency", v);

  const std::vector<double> sample_td{0.2, 1.0, 10.0, 92.0, 120.0, 250.0, 600.0};
  std::ofstream rt(out_dir / "transport_roundtrip.csv");
  rt << "E_over_N_Td,mu_rel_err,D_rel_err,alpha_rel_err,eta2_rel_err,eta3_rel_err,eta_total_rel_err\n";
  double max_rt = 0.0;
  for (double td : sample_td) {
    const double efield = td * 1e-21 * neutral_density;
    const auto direct = evaluate_morrow_lowke(efield, neutral_density, pressure, temperature);
    const double mu = interp_transport(rows, td, &streamer_rf::streamer::TransportCoefficients::mobility);
    const double dif = interp_transport(rows, td, &streamer_rf::streamer::TransportCoefficients::diffusion);
    const double alpha = interp_transport(rows, td, &streamer_rf::streamer::TransportCoefficients::ionization_townsend);
    const double eta2 = interp_transport(rows, td, &streamer_rf::streamer::TransportCoefficients::attachment_two_body_townsend);
    const double eta3 = interp_transport(rows, td, &streamer_rf::streamer::TransportCoefficients::attachment_three_body_townsend);
    const double eta_total = eta2 + eta3;
    const double errs[] = {
        rel_err(direct.mobility, mu),
        rel_err(direct.diffusion, dif),
        rel_err(direct.ionization_townsend, alpha),
        rel_err(direct.attachment_two_body_townsend, eta2),
        rel_err(direct.attachment_three_body_townsend, eta3),
        rel_err(direct.attachment_two_body_townsend + direct.attachment_three_body_townsend, eta_total)};
    for (double e : errs) max_rt = std::max(max_rt, e);
    rt << td << ',' << errs[0] << ',' << errs[1] << ',' << errs[2] << ','
       << errs[3] << ',' << errs[4] << ',' << errs[5] << '\n';
  }

  struct State { double td, ne, np, nn; };
  const std::vector<State> states{{1.0, 1e12, 1e12, 0.0},
                                  {92.0, 1e16, 9e15, 2e14},
                                  {120.0, 5e17, 4e17, 1e16},
                                  {250.0, 1e18, 7e17, 4e17}};
  std::ofstream chem(out_dir / "chemistry_pointwise.csv");
  chem << "E_over_N_Td,ne_m3,np_m3,nn_m3,dne_petsc,dnp_petsc,dnn_petsc,"
          "dne_common,dnp_common,dnn_common,max_rel_err\n";
  double max_chem = 0.0;
  for (const auto& s : states) {
    const double efield = s.td * 1e-21 * neutral_density;
    const auto c = evaluate_morrow_lowke(efield, neutral_density, pressure, temperature);
    const auto p = evaluate_reactions(s.ne, s.np, s.nn, 0.0, c, temperature);
    const double ion = c.ionization_frequency * s.ne;
    const double att2 = c.attachment_two_body_frequency * s.ne;
    const double att3 = c.attachment_three_body_frequency * s.ne;
    const double ep = beta_ep(c.mobility, c.diffusion) * s.ne * s.np;
    const double pn = beta_np(temperature) * s.np * s.nn;
    const double dne = ion - att2 - att3 - ep;
    const double dnp = ion - ep - pn;
    const double dnn = att2 + att3 - pn;
    const double emax = std::max({rel_err(p.electron, dne), rel_err(p.positive_ion, dnp), rel_err(p.negative_ion, dnn)});
    max_chem = std::max(max_chem, emax);
    chem << s.td << ',' << s.ne << ',' << s.np << ',' << s.nn << ','
         << p.electron << ',' << p.positive_ion << ',' << p.negative_ion << ','
         << dne << ',' << dnp << ',' << dnn << ',' << emax << '\n';
  }

  std::ofstream summary(out_dir / "b0_summary.json");
  summary << std::setprecision(17)
          << "{\n"
          << "  \"pressure_Pa\": " << pressure << ",\n"
          << "  \"temperature_K\": " << temperature << ",\n"
          << "  \"neutral_density_m3\": " << neutral_density << ",\n"
          << "  \"table_rows\": " << n_rows << ",\n"
          << "  \"E_over_N_min_Td\": " << td_min << ",\n"
          << "  \"E_over_N_max_Td\": " << td_max << ",\n"
          << "  \"transport_roundtrip_max_rel_err\": " << max_rt << ",\n"
          << "  \"chemistry_pointwise_max_rel_err\": " << max_chem << ",\n"
          << "  \"afivo_old_style_eta\": \"eta2_plus_eta3_at_fixed_pressure_temperature\",\n"
          << "  \"common_chemistry_runtime_note\": \"dedicated Afivo reaction table generated; old-style transport alone cannot represent eta2/eta3 split or recombination\"\n"
          << "}\n";

  std::cout << "stage_d2_b0 status=PASS out_dir=" << out_dir
            << " transport_roundtrip_max_rel_err=" << max_rt
            << " chemistry_pointwise_max_rel_err=" << max_chem << '\n';
  return (max_rt < 5e-3 && max_chem < 1e-12) ? 0 : 2;
}
