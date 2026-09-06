#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/streamer/LfaAudit.hpp"
#include "streamer_rf/transport.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
int checks = 0;
void check(bool ok, const char* name) {
  ++checks;
  if (!ok) throw std::runtime_error(name);
  std::cout << "ok " << checks << " - " << name << '\n';
}

AxisymmetricNeedlePlaneGeometry tiny_geometry() {
  return AxisymmetricNeedlePlaneGeometry("stage-c2-test-needle-plane", 0.0, 1.35e-3, 0.08e-3, 0.04e-3, 0.08e-3);
}

StreamerConfig electrode_config(const AxisymmetricNeedlePlaneGeometry& geom, const VoltageWaveform& voltage, bool sp3=false) {
  StreamerConfig cfg;
  cfg.electrode_geometry = &geom;
  cfg.voltage_waveform = &voltage;
  cfg.photoionization = sp3;
  cfg.n_ref = 1e10;
  cfg.head_ne_threshold = 1e9;
  cfg.bridge_ne_threshold = 1e9;
  cfg.elliptic = {1e-10, 1e-14, 20000};
  return cfg;
}

void set_uniform_field(StreamerSolver& solver, const AxisymmetricGrid& g, double e) {
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      solver.state().er(i, j) = 0.0;
      solver.state().ez(i, j) = e;
      solver.state().emag(i, j) = std::abs(e);
      solver.state().sph(i, j) = 0.0;
    }
  }
}

void fill_exponential_z(StreamerState& state, const AxisymmetricGrid& g, double e0, double scale) {
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      state.er(i, j) = 0.0;
      state.ez(i, j) = e0 * std::exp(g.z(j) / scale);
      state.emag(i, j) = std::abs(state.ez(i, j));
    }
  }
}

double relerr(double a, double b) {
  return std::abs(a - b) / std::max(std::abs(b), 1e-300);
}

ElectronRelaxationMetadata synthetic_relaxation_metadata() {
  ElectronRelaxationMetadata m;
  m.source_id = "synthetic-unit-test";
  m.source_reference = "C-R1c deterministic synthetic table";
  m.gas_composition = "synthetic air";
  m.pressure_Pa = 101325.0;
  m.temperature_K = 300.0;
  m.neutral_density_m3 = 101325.0 / (1.380649e-23 * 300.0);
  return m;
}

bool throws_invalid_relaxation_table(std::vector<double> eovern, std::optional<std::vector<double>> tau,
                                     std::optional<std::vector<double>> lambda) {
  try {
    ElectronRelaxationTable bad(synthetic_relaxation_metadata(), std::move(eovern), std::move(tau), std::move(lambda));
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}
}  // namespace

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  try {
    AxisymmetricGrid g(16, 32, 0.8e-3, 0.0, 1.6e-3);
    auto geom = tiny_geometry();
    ConstantVoltage voltage(100.0);
    auto cfg = electrode_config(geom, voltage, false);

    check(std::abs(isg0_flux(2.0, 1.0, 3.0, 3.0, 1.0, 0.1, 0.01, 1.0) - sg_flux(2.0, 1.0, 3.0, 1.0, 0.1)) < 1e-14,
          "gas-gas ISG-0 regression");

    const double f_plus = absorbing_electrode_flux(10.0, -2.0, 0.5, 0.1, true);
    const double f_minus = absorbing_electrode_flux(10.0, 2.0, 0.5, 0.1, false);
    check(f_plus > 0.0 && f_minus < 0.0, "absorbing electrode diffusion loss without injection");
    check(absorbing_electrode_flux(0.0, 1.0, 0.5, 0.1, true) == 0.0, "absorbing electrode zero density gives zero flux");
    check(absorbing_electrode_flux(10.0, 3.0, 0.0, 0.1, true) > 0.0, "absorbing electrode drift outflow direction");

    StreamerSolver solver(g, cfg);
    solver.initialize_gaussian_at_tip_offset(1e11, 0.04e-3, -0.16e-3);
    for (int j = 0; j < g.nz(); ++j) {
      for (int i = 0; i < g.nr(); ++i) {
        if (geom.classify(g, i, j) != ElectrodeCellType::Gas) {
          check(solver.state().ne(i, j) == 0.0 && solver.state().np(i, j) == 0.0 && solver.state().nn(i, j) == 0.0,
                "electrode cells plasma-free after seed");
          j = g.nz();
          break;
        }
      }
    }

    StreamerDiagnostics d;
    check(solver.step(1e-14, d), "electrode-aware tiny dynamic step accepted");
    check(std::isfinite(d.emax) && std::isfinite(d.ne_max) && d.conservation_residual < 1e-9, "electrode-aware conservation with absorption");
    check(d.absorbed_electron_hv >= 0.0 && d.absorbed_electron_ground >= 0.0, "electrode absorption diagnostics nonnegative");
    for (int j = 0; j < g.nz(); ++j) {
      for (int i = 0; i < g.nr(); ++i) {
        if (geom.classify(g, i, j) != ElectrodeCellType::Gas) {
          check(solver.state().ne(i, j) == 0.0 && solver.state().np(i, j) == 0.0 && solver.state().nn(i, j) == 0.0,
                "electrode cells plasma-free after step");
          j = g.nz();
          break;
        }
      }
    }

    SampledVoltage waveform({0.0, 1e-12}, {100.0, 200.0});
    auto wcfg = electrode_config(geom, waveform, false);
    StreamerSolver wsolver(g, wcfg);
    wsolver.initialize_gaussian_at_tip_offset(1e9, 0.04e-3, -0.16e-3);
    StreamerDiagnostics wd;
    check(wsolver.step(1e-13, wd), "waveform dynamic Poisson step accepted");
    check(std::abs(wd.applied_voltage - 110.0) < 1e-9 && wd.phi_hv_residual < 1e-8, "waveform-dynamic Poisson follows V(t)");

    auto pcfg = electrode_config(geom, voltage, true);
    StreamerSolver psolver(g, pcfg);
    psolver.initialize_gaussian_at_tip_offset(1e10, 0.04e-3, -0.16e-3);
    StreamerDiagnostics pd;
    check(psolver.step(1e-15, pd), "SP3 gas-mask step accepted");
    for (int j = 0; j < g.nz(); ++j) {
      for (int i = 0; i < g.nr(); ++i) {
        if (geom.classify(g, i, j) != ElectrodeCellType::Gas) {
          check(psolver.state().sph(i, j) == 0.0 && psolver.state().ne(i, j) == 0.0, "SP3 source masked in electrodes");
          j = g.nz();
          break;
        }
      }
    }

    ConstantVoltage high_voltage(800.0);
    auto scfg = electrode_config(geom, high_voltage, false);
    scfg.head_ne_threshold = 1e8;
    StreamerSolver smoke(g, scfg);
    smoke.initialize_gaussian_at_tip_offset(1e10, 0.04e-3, -0.16e-3);
    bool ok = true;
    StreamerDiagnostics sd;
    for (int n = 0; n < 5; ++n) ok = ok && smoke.step(1e-15, sd);
    check(ok && std::isfinite(sd.ne_max) && std::isfinite(sd.sigma_max) && std::isfinite(sd.head_position), "short dynamic streamer smoke finite");

    AxisymmetricGrid rg(8, 16, 80e-6, 0.0, 90e-6);
    AxisymmetricNeedlePlaneGeometry rgeom("stage-c2-reaction-test", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
    ConstantVoltage rvoltage(500.0);
    auto rcfg = electrode_config(rgeom, rvoltage, false);
    rcfg.n_ref = 1e12;

    StreamerSolver inactive(rg, rcfg);
    inactive.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_uniform_field(inactive, rg, 1e-6);
    inactive.state().ne(0, 8) = inactive.state().np(0, 8) = 1e14;
    inactive.state().emag(0, 8) = 7e7;
    inactive.state().ez(0, 8) = 7e7;
    inactive.state().ne(7, 15) = inactive.state().np(7, 15) = 500.0 * inactive.numerical_density_tolerance();
    inactive.state().emag(7, 15) = 1e-6;
    auto ilim = inactive.timestep_limits();
    check(ilim.reaction > 1e-13 && ilim.controller != "reaction", "inactive electron cell does not dominate reaction dt");

    StreamerSolver active_decay(rg, rcfg);
    active_decay.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_uniform_field(active_decay, rg, 1e-6);
    active_decay.state().ne(4, 8) = active_decay.state().np(4, 8) = 1e14;
    auto alim = active_decay.timestep_limits();
    check(alim.controller == "reaction" && alim.reaction < 1e-18, "active attachment decay still limits reaction dt");
    StreamerDiagnostics ad;
    check(active_decay.step(0.5 * alim.reaction, ad), "reaction positivity timestep accepts protected substep");
    StreamerSolver rejected_decay(rg, rcfg);
    rejected_decay.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_uniform_field(rejected_decay, rg, 1e-6);
    rejected_decay.state().ne(4, 8) = rejected_decay.state().np(4, 8) = 1e14;
    StreamerDiagnostics rd;
    check(!rejected_decay.step(6.0 * alim.reaction, rd), "reaction positivity timestep rejects overlarge decay step");

    StreamerSolver conductor_base(rg, rcfg);
    conductor_base.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_uniform_field(conductor_base, rg, 7e7);
    conductor_base.state().ne(0, 8) = conductor_base.state().np(0, 8) = 1e14;
    StreamerSolver conductor_polluted(rg, rcfg);
    conductor_polluted.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_uniform_field(conductor_polluted, rg, 7e7);
    conductor_polluted.state().ne(0, 8) = conductor_polluted.state().np(0, 8) = 1e14;
    bool marked_conductor = false;
    for (int j = 0; j < rg.nz(); ++j) {
      for (int i = 0; i < rg.nr(); ++i) {
        if (rgeom.classify(rg, i, j) != ElectrodeCellType::Gas) {
          conductor_polluted.state().ne(i, j) = conductor_polluted.state().np(i, j) = 1e30;
          conductor_polluted.state().emag(i, j) = 1e-6;
          marked_conductor = true;
        }
      }
    }
    check(marked_conductor, "reaction test geometry contains conductor cells");
    auto base_lim = conductor_base.timestep_limits();
    auto polluted_lim = conductor_polluted.timestep_limits();
    check(base_lim.controller == polluted_lim.controller && base_lim.selected == polluted_lim.selected &&
          base_lim.dielectric == polluted_lim.dielectric && base_lim.reaction == polluted_lim.reaction,
          "conductor excluded from timestep scan");

    StreamerConfig legacy_cfg;
    legacy_cfg.photoionization = false;
    legacy_cfg.n_ref = 1e12;
    StreamerSolver legacy(rg, legacy_cfg);
    legacy.initialize_gaussian(0.0, 3e-6, 45e-6);
    set_uniform_field(legacy, rg, 1e-6);
    legacy.state().ne(7, 15) = legacy.state().np(7, 15) = 50.0 * legacy.numerical_density_tolerance();
    auto legacy_lim = legacy.timestep_limits();
    check(legacy_lim.controller == "reaction" && legacy_lim.reaction < 1e-18, "legacy timestep regression unchanged");

    {
      AxisymmetricGrid lg(12, 80, 0.6e-3, 0.0, 2.0e-3);
      StreamerConfig lcfg;
      lcfg.photoionization = false;
      StreamerState ls(lg);
      const double L0 = 0.45e-3;
      fill_exponential_z(ls, lg, 2.0e5, L0);
      auto audit = compute_lfa_audit(lg, ls, lcfg, 0.0);
      check(audit.le_valid_fraction > 0.95 && relerr(audit.le_median_m, L0) < 2e-3,
            "LFA spatial field scale recovers exponential scale");
      check(audit.tauE_valid_fraction == 0.0 && audit.temporal_status == "INVALID_INITIAL_SAMPLE",
            "LFA initial temporal sample invalid");
    }

    {
      AxisymmetricGrid tg(10, 16, 0.5e-3, 0.0, 0.8e-3);
      StreamerConfig tcfg;
      tcfg.photoionization = false;
      StreamerState old_s(tg), new_s(tg);
      const double tau0 = 2.5e-12;
      const double dt1 = 3.0e-14;
      fill_exponential_z(old_s, tg, 1.0e6, 0.3e-3);
      fill_exponential_z(new_s, tg, 1.0e6 * std::exp(dt1 / tau0), 0.3e-3);
      auto ta = compute_lfa_audit(tg, new_s, tcfg, dt1, &old_s.emag);
      check(ta.temporal_status == "VALID" && ta.tauE_valid_fraction > 0.99 && relerr(ta.tauE_median_s, tau0) < 0.01,
            "LFA temporal field scale recovers exponential scale");
      const double dt2 = 7.0e-14;
      fill_exponential_z(new_s, tg, 1.0e6 * std::exp(dt2 / tau0), 0.3e-3);
      auto tb = compute_lfa_audit(tg, new_s, tcfg, dt2, &old_s.emag);
      check(tb.temporal_status == "VALID" && tb.tauE_valid_fraction > 0.99 && relerr(tb.tauE_median_s, tau0) < 0.02,
            "LFA temporal field scale handles variable accepted dt");
    }

    {
      AxisymmetricGrid zg(8, 10, 0.4e-3, 0.0, 0.5e-3);
      StreamerConfig zcfg;
      StreamerState zs(zg);
      for (int j = 0; j < zg.nz(); ++j) {
        for (int i = 0; i < zg.nr(); ++i) zs.emag(i, j) = 3.0e6;
      }
      auto za = compute_lfa_audit(zg, zs, zcfg, 0.0);
      check(za.le_valid_fraction == 0.0 && za.near_zero_gradE_cells == za.selected_gas_cells,
            "LFA near-zero spatial gradient is explicitly invalid");
      StreamerState zn(zg);
      fill_exponential_z(zn, zg, 2.0e5, 0.2e-3);
      zn.emag(3, 4) = std::numeric_limits<double>::quiet_NaN();
      auto nf = compute_lfa_audit(zg, zn, zcfg, 0.0);
      check(nf.nonfinite_input_cells > 0 && std::isfinite(nf.eovern_max_Td), "LFA nonfinite field input is counted");
    }

    {
      AxisymmetricGrid mg(8, 16, 80e-6, 0.0, 90e-6);
      AxisymmetricNeedlePlaneGeometry mgeom("stage-c2-lfa-mask-test", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
      ConstantVoltage mvoltage(500.0);
      auto mcfg = electrode_config(mgeom, mvoltage, false);
      StreamerState ms(mg);
      fill_exponential_z(ms, mg, 1.0e6, 25e-6);
      int conductor_cells = 0;
      for (int j = 0; j < mg.nz(); ++j) {
        for (int i = 0; i < mg.nr(); ++i) {
          if (mgeom.classify(mg, i, j) != ElectrodeCellType::Gas) {
            ++conductor_cells;
            ms.emag(i, j) = 1.0e99;
          }
        }
      }
      auto ma = compute_lfa_audit(mg, ms, mcfg, 0.0);
      check(conductor_cells > 0 && ma.gas_cells + conductor_cells == mg.nr() * mg.nz() && ma.eovern_max_Td < 1e10,
            "LFA audit excludes conductor cells");
      std::vector<unsigned char> mask(mg.size(), 0);
      int selected = 0;
      for (int j = 0; j < mg.nz(); ++j) {
        for (int i = 0; i < mg.nr(); ++i) {
          if (mgeom.classify(mg, i, j) == ElectrodeCellType::Gas && mg.z(j) < 50e-6) {
            mask[static_cast<std::size_t>(j) * mg.nr() + i] = 1;
            ++selected;
          }
        }
      }
      auto masked = compute_lfa_audit(mg, ms, mcfg, 0.0, nullptr, &mask);
      check(masked.selected_gas_cells == selected && masked.selected_gas_cells < masked.gas_cells,
            "LFA optional gas-cell mask restricts summary region");
    }

    {
      auto meta = synthetic_relaxation_metadata();
      ElectronRelaxationTable table(meta, {10.0, 20.0, 40.0}, std::vector<double>{1e-12, 2e-12, 4e-12},
                                    std::vector<double>{1e-6, 2e-6, 4e-6});
      auto mid = table.lookup(15.0);
      check(mid.status == "IN_RANGE" && mid.has_tau_epsilon && mid.has_lambda_epsilon &&
                relerr(mid.tau_epsilon_s, 1.5e-12) < 1e-14 && relerr(mid.lambda_epsilon_m, 1.5e-6) < 1e-14,
            "electron relaxation table interpolates linearly");
      auto lo = table.lookup(10.0);
      auto hi = table.lookup(40.0);
      check(lo.status == "IN_RANGE" && hi.status == "IN_RANGE" && lo.tau_epsilon_s == 1e-12 &&
                hi.lambda_epsilon_m == 4e-6,
            "electron relaxation table handles exact bounds");
      auto outside = table.lookup(5.0);
      check(outside.status == "OUTSIDE_RELAXATION_TABLE_RANGE" && !outside.has_tau_epsilon,
            "electron relaxation table rejects silent extrapolation");
      check(throws_invalid_relaxation_table({10.0, 10.0, 40.0}, std::vector<double>{1e-12, 2e-12, 4e-12},
                                            std::vector<double>{1e-6, 2e-6, 4e-6}),
            "electron relaxation table rejects duplicate E/N");
      check(throws_invalid_relaxation_table({20.0, 10.0, 40.0}, std::vector<double>{1e-12, 2e-12, 4e-12},
                                            std::vector<double>{1e-6, 2e-6, 4e-6}),
            "electron relaxation table rejects nonmonotonic E/N");
      check(throws_invalid_relaxation_table({10.0, 20.0}, std::vector<double>{1e-12, -2e-12},
                                            std::vector<double>{1e-6, 2e-6}),
            "electron relaxation table rejects nonpositive relaxation quantities");
      ElectronRelaxationTable tau_only(meta, {10.0, 20.0}, std::vector<double>{1e-12, 2e-12}, std::nullopt);
      ElectronRelaxationTable lambda_only(meta, {10.0, 20.0}, std::nullopt, std::vector<double>{1e-6, 2e-6});
      check(tau_only.has_tau_epsilon() && !tau_only.has_lambda_epsilon() &&
                lambda_only.has_lambda_epsilon() && !lambda_only.has_tau_epsilon(),
            "electron relaxation table supports missing tau or lambda");
    }

    {
      AxisymmetricGrid cg(10, 64, 0.5e-3, 0.0, 1.0e-3);
      StreamerConfig ccfg;
      ccfg.photoionization = false;
      StreamerState cold(cg), cnew(cg);
      const double L0 = 0.8e-3;
      const double tau0 = 3.0e-12;
      const double dt = 4.0e-14;
      fill_exponential_z(cold, cg, 1.0e6, L0);
      fill_exponential_z(cnew, cg, 1.0e6 * std::exp(dt / tau0), L0);
      ElectronRelaxationTable table(synthetic_relaxation_metadata(), {1.0, 100.0, 1000.0},
                                    std::vector<double>{6e-12, 6e-12, 6e-12},
                                    std::vector<double>{2e-6, 2e-6, 2e-6});
      auto chi = compute_lfa_audit(cg, cnew, ccfg, dt, &cold.emag, nullptr, &table);
      check(chi.lfa_relaxation_data_status == "AVAILABLE_TAU_AND_LAMBDA" &&
                chi.lfa_applicability == "UNRESOLVED_NO_APPROVED_CRITERION",
            "LFA chi audit keeps applicability unresolved without approved criterion");
      check(chi.chiL_valid_cells > 0 && relerr(chi.chiL_median, 2e-6 / L0) < 0.01,
            "LFA chi_L reconstructs lambda over L_E");
      check(chi.chiT_valid_cells > 0 && relerr(chi.chiT_median, 6e-12 / tau0) < 0.02,
            "LFA chi_t reconstructs tau over tau_E");
      auto no_data = compute_lfa_audit(cg, cnew, ccfg, dt, &cold.emag);
      check(no_data.lfa_relaxation_data_status == "NOT_AVAILABLE" && no_data.chiL_valid_cells == 0 &&
                no_data.chiT_valid_cells == 0,
            "LFA no-relaxation-data path produces no fake chi values");
      ElectronRelaxationTable narrow(synthetic_relaxation_metadata(), {1.0, 2.0}, std::vector<double>{1e-12, 2e-12},
                                     std::vector<double>{1e-6, 2e-6});
      auto outside = compute_lfa_audit(cg, cnew, ccfg, dt, &cold.emag, nullptr, &narrow);
      check(outside.fraction_outside_relaxation_table > 0.9 && outside.relaxation_coverage_fraction < 0.1,
            "LFA audit reports outside relaxation table coverage");
      StreamerState flat_old(cg), flat_new(cg);
      for (int j = 0; j < cg.nz(); ++j) {
        for (int i = 0; i < cg.nr(); ++i) flat_old.emag(i, j) = flat_new.emag(i, j) = 1.0e6;
      }
      auto invalid = compute_lfa_audit(cg, flat_new, ccfg, dt, &flat_old.emag, nullptr, &table);
      check(invalid.chiL_valid_cells == 0 && invalid.chiT_valid_cells == 0,
            "LFA chi audit excludes invalid L_E and tau_E cells");
    }

    {
      AxisymmetricGrid rg2(10, 30, 0.5e-3, 0.0, 1.5e-3);
      StreamerConfig rcfg2;
      rcfg2.photoionization = false;
      rcfg2.n_ref = 1e8;
      StreamerState rs(rg2);
      fill_exponential_z(rs, rg2, 1.0e5, 0.3e-3);
      for (int j = 0; j < rg2.nz(); ++j) {
        for (int i = 0; i < rg2.nr(); ++i) rs.ne(i, j) = j > rg2.nz() / 2 ? 1e12 : 0.0;
      }
      auto regions = default_lfa_audit_regions();
      auto all = compute_lfa_audit(rg2, rs, rcfg2, 0.0, nullptr, nullptr, nullptr, regions[0]);
      auto active = compute_lfa_audit(rg2, rs, rcfg2, 0.0, nullptr, nullptr, nullptr, regions[1]);
      auto high = compute_lfa_audit(rg2, rs, rcfg2, 0.0, nullptr, nullptr, nullptr, regions[2]);
      check(all.region_name == "ALL_GAS" && active.region_name == "ACTIVE_ELECTRON" &&
                high.region_name == "HIGH_FIELD",
            "LFA default region names are stable");
      check(active.selected_gas_cells < all.selected_gas_cells && active.selected_gas_cells > 0,
            "LFA active-electron region uses numerical electron mask");
      check(high.selected_gas_cells < all.selected_gas_cells && high.selected_gas_cells > 0,
            "LFA high-field region uses diagnostic relative field mask");
    }

    std::cout << "passed " << checks << " Stage C2 C++ checks\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
  PetscFinalize();
  return 0;
}
