#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/streamer/HeadTracker.hpp"
#include "streamer_rf/streamer/JouleHandoff.hpp"
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

void set_single_charge_head(StreamerState& state, const AxisymmetricGrid& g, int i0, int j0, double rho) {
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      state.rho(i, j) = 0.0;
      state.emag(i, j) = 1.0e6;
    }
  }
  state.rho(i0, j0) = rho;
  state.emag(i0, j0) = 2.0e6;
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
      TransportCoefficients tc;
      tc.mobility = 1.0;
      tc.diffusion = 0.5;
      tc.ionization_frequency = 7.0;
      auto impact = electron_reaction_source_components(2.0, 0.0, 0.0, 0.0, tc, 300.0);
      check(impact.impact_source_m3s == 14.0 && impact.net_electron_reaction_source_m3s == 14.0 &&
                impact.algebraic_closure_residual_m3s == 0.0,
            "reaction source decomposition impact-only closure");

      tc.ionization_frequency = 0.0;
      tc.attachment_two_body_frequency = 11.0;
      tc.attachment_three_body_frequency = 13.0;
      auto attach = electron_reaction_source_components(3.0, 0.0, 0.0, 0.0, tc, 300.0);
      check(attach.attachment2_loss_m3s == 33.0 && attach.attachment3_loss_m3s == 39.0 &&
                attach.net_electron_reaction_source_m3s == -72.0,
            "reaction source decomposition separates 2-body and 3-body attachment");

      tc.attachment_two_body_frequency = 0.0;
      tc.attachment_three_body_frequency = 0.0;
      auto recomb = electron_reaction_source_components(4.0, 5.0, 6.0, 0.0, tc, 300.0);
      check(recomb.electron_recombination_loss_m3s > 0.0 &&
                recomb.net_electron_reaction_source_m3s == -recomb.electron_recombination_loss_m3s,
            "reaction source decomposition includes electron-positive recombination");

      auto photo = electron_reaction_source_components(0.0, 0.0, 0.0, 9.0, tc, 300.0);
      auto no_photo = electron_reaction_source_components(0.0, 0.0, 0.0, 0.0, tc, 300.0);
      check(photo.photo_source_m3s == 9.0 && photo.net_electron_reaction_source_m3s == 9.0 &&
                no_photo.photo_source_m3s == 0.0 && no_photo.net_electron_reaction_source_m3s == 0.0,
            "reaction source decomposition maps SP3 source and SP3-off zero source");

      check(impact.impact_source_m3s >= 0.0 && attach.attachment2_loss_m3s >= 0.0 &&
                attach.attachment3_loss_m3s >= 0.0 && recomb.electron_recombination_loss_m3s >= 0.0,
            "reaction source decomposition component magnitudes are nonnegative");
    }

    {
      AxisymmetricGrid sg(8, 16, 80e-6, 0.0, 90e-6);
      AxisymmetricNeedlePlaneGeometry sgeom("stage-c2-source-decomp-test", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
      ConstantVoltage svoltage(500.0);
      auto scfg2 = electrode_config(sgeom, svoltage, false);
      scfg2.n_ref = 1e12;
      StreamerSolver src_solver(sg, scfg2);
      src_solver.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
      set_uniform_field(src_solver, sg, 3.0e6);
      for (int j = 0; j < sg.nz(); ++j) {
        for (int i = 0; i < sg.nr(); ++i) {
          if (sgeom.classify(sg, i, j) == ElectrodeCellType::Gas) {
            src_solver.state().ne(i, j) = 1.0e12 * (1 + i + j);
            src_solver.state().np(i, j) = 0.5e12 * (1 + i);
            src_solver.state().nn(i, j) = 0.25e12 * (1 + j);
            src_solver.state().sph(i, j) = 1.0e18;
          } else {
            src_solver.state().ne(i, j) = src_solver.state().np(i, j) = 1e30;
            src_solver.state().nn(i, j) = 1e30;
            src_solver.state().sph(i, j) = 1e40;
          }
        }
      }
      auto sdg = src_solver.reaction_source_decomposition_diagnostics();
      double manual_impact = 0.0, manual_photo = 0.0, manual_net = 0.0;
      int manual_gas = 0;
      for (int j = 0; j < sg.nz(); ++j) {
        for (int i = 0; i < sg.nr(); ++i) {
          if (sgeom.classify(sg, i, j) != ElectrodeCellType::Gas) continue;
          ++manual_gas;
          auto q = evaluate_morrow_lowke(src_solver.state().emag(i, j), scfg2.neutral_density, scfg2.pressure, scfg2.temperature);
          auto s = evaluate_reactions(src_solver.state().ne(i, j), src_solver.state().np(i, j), src_solver.state().nn(i, j),
                                      src_solver.state().sph(i, j), q, scfg2.temperature);
          manual_impact += s.ionization * sg.cell_volume(i);
          manual_photo += src_solver.state().sph(i, j) * sg.cell_volume(i);
          manual_net += s.electron * sg.cell_volume(i);
        }
      }
      check(sdg.gas_cells == manual_gas && relerr(sdg.impact_rate_s1, manual_impact) < 1e-14 &&
                relerr(sdg.photo_rate_s1, manual_photo) < 1e-14 && relerr(sdg.net_electron_reaction_rate_s1, manual_net) < 1e-14,
            "reaction source decomposition uses gas mask and axisymmetric volume integration");
      check(sdg.source_closure_abs_s1 <= 1e-6 * std::max(std::abs(sdg.net_electron_reaction_rate_s1), 1.0) &&
                sdg.source_closure_rel < 1e-14 && sdg.source_stage == "PRE_LIMITER",
            "reaction source decomposition algebraic closure is exact before positivity limiter");
    }

    {
      AxisymmetricGrid hg(8, 20, 80e-6, 0.0, 100e-6);
      StreamerConfig hcfg;
      hcfg.photoionization = false;
      StreamerState hs(hg);
      StreamerHeadTracker tracker(StreamerHeadSegmentationConfig{0.2});
      set_single_charge_head(hs, hg, 2, 8, 1.0);
      hs.time = 0.0;
      auto h0 = tracker.sample(hg, hs, hcfg);
      hs.time = 1.0e-12;
      auto h1 = tracker.sample(hg, hs, hcfg);
      hs.time = 2.0e-12;
      auto h2 = tracker.sample(hg, hs, hcfg);
      check(h0.head_valid && h1.velocity_valid && h2.acceleration_valid &&
                h1.head_velocity_z_m_s == 0.0 && h2.head_acceleration_z_m_s2 == 0.0,
            "head tracker stationary compact charge has zero velocity and acceleration");
      check(h0.velocity_status == "INSUFFICIENT_HISTORY_FOR_VELOCITY" &&
                h1.acceleration_status == "INSUFFICIENT_HISTORY_FOR_ACCELERATION",
            "head tracker reports insufficient kinematic history");

      StreamerHeadTracker vtracker(StreamerHeadSegmentationConfig{0.2});
      set_single_charge_head(hs, hg, 2, 5, 1.0);
      hs.time = 0.0;
      vtracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 7, 1.0);
      hs.time = 1.0e-12;
      auto v1 = vtracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 9, 1.0);
      hs.time = 2.0e-12;
      auto v2 = vtracker.sample(hg, hs, hcfg);
      const double v_expected = 2.0 * hg.dz() / 1.0e-12;
      check(v1.velocity_valid && v2.acceleration_valid && relerr(v2.head_velocity_z_m_s, v_expected) < 1e-14 &&
                std::abs(v2.head_acceleration_z_m_s2) < 1e-9 * std::abs(v_expected) / 1.0e-12,
            "head tracker recovers constant-velocity moving head");

      StreamerHeadTracker atracker(StreamerHeadSegmentationConfig{0.2});
      set_single_charge_head(hs, hg, 2, 5, 1.0);
      hs.time = 0.0;
      atracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 6, 1.0);
      hs.time = 1.0e-12;
      atracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 9, 1.0);
      hs.time = 2.0e-12;
      auto a2 = atracker.sample(hg, hs, hcfg);
      const double a_expected = 2.0 * hg.dz() / (1.0e-12 * 1.0e-12);
      check(a2.acceleration_valid && relerr(a2.head_acceleration_z_m_s2, a_expected) < 1e-14,
            "head tracker recovers constant acceleration with uniform dt");

      StreamerHeadTracker ndt_tracker(StreamerHeadSegmentationConfig{0.2});
      set_single_charge_head(hs, hg, 2, 5, 1.0);
      hs.time = 0.0;
      ndt_tracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 6, 1.0);
      hs.time = 1.0e-12;
      ndt_tracker.sample(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 2, 8, 1.0);
      hs.time = 3.0e-12;
      auto ndt = ndt_tracker.sample(hg, hs, hcfg);
      check(ndt.acceleration_valid && std::abs(ndt.head_acceleration_z_m_s2) < 1e-9 * std::abs(v_expected) / 1.0e-12,
            "head tracker handles variable accepted dt");

      set_single_charge_head(hs, hg, 3, 8, 2.0);
      auto hp = compute_streamer_head_diagnostics(hg, hs, hcfg);
      set_single_charge_head(hs, hg, 3, 8, -2.0);
      auto hn = compute_streamer_head_diagnostics(hg, hs, hcfg);
      check(hp.head_valid && hp.head_polarity == 1 && hp.head_charge_C > 0.0 &&
                hn.head_valid && hn.head_polarity == -1 && hn.head_charge_C < 0.0,
            "head tracker supports positive and negative polarity");
      set_single_charge_head(hs, hg, 2, 7, 1.0);
      hs.rho(6, 15) = -5.0;
      auto auto_polarity = compute_streamer_head_diagnostics(hg, hs, hcfg);
      auto requested_positive = compute_streamer_head_diagnostics(hg, hs, hcfg, StreamerHeadSegmentationConfig{0.2, 1});
      check(auto_polarity.head_polarity == -1 && requested_positive.head_polarity == 1 &&
                relerr(requested_positive.head_z_m, hg.z(7)) < 1e-14,
            "head tracker honors explicit case polarity");

      set_single_charge_head(hs, hg, 4, 10, 3.0);
      auto qh = compute_streamer_head_diagnostics(hg, hs, hcfg);
      check(qh.head_cell_count == 1 && relerr(qh.head_charge_C, 3.0 * hg.cell_volume(4)) < 1e-14 &&
                qh.head_peak_E_V_m == 2.0e6,
            "head tracker integrates axisymmetric head charge");

      for (int j = 0; j < hg.nz(); ++j) for (int i = 0; i < hg.nr(); ++i) hs.rho(i, j) = 0.0;
      hs.emag(0, 0) = 1.0;
      auto none = compute_streamer_head_diagnostics(hg, hs, hcfg);
      check(!none.head_valid && none.head_status == "INSUFFICIENT_CHARGE",
            "head tracker reports no-head invalid state");

      set_single_charge_head(hs, hg, 1, 3, 10.0);
      hs.rho(6, 16) = 9.0;
      auto island = compute_streamer_head_diagnostics(hg, hs, hcfg);
      check(island.head_valid && island.head_cell_count == 1,
            "head tracker excludes disconnected same-polarity background island");
    }

    {
      AxisymmetricGrid cg(8, 16, 80e-6, 0.0, 90e-6);
      AxisymmetricNeedlePlaneGeometry cgeom("stage-c2-head-conductor-test", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
      ConstantVoltage cvoltage(500.0);
      auto ccfg = electrode_config(cgeom, cvoltage, false);
      StreamerState cs(cg);
      int gas_i = -1, gas_j = -1;
      for (int j = 0; j < cg.nz(); ++j) {
        for (int i = 0; i < cg.nr(); ++i) {
          cs.rho(i, j) = 0.0;
          cs.emag(i, j) = 1.0e6;
          if (gas_i < 0 && cgeom.classify(cg, i, j) == ElectrodeCellType::Gas) {
            gas_i = i;
            gas_j = j;
          }
          if (cgeom.classify(cg, i, j) != ElectrodeCellType::Gas) cs.rho(i, j) = 1.0e9;
        }
      }
      cs.rho(gas_i, gas_j) = 1.0;
      auto masked_head = compute_streamer_head_diagnostics(cg, cs, ccfg);
      check(masked_head.head_valid && masked_head.head_charge_C < 1e-12 &&
                relerr(masked_head.head_z_m, cg.z(gas_j)) < 1e-14,
            "head tracker masks conductor charge");
    }

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
      AxisymmetricGrid jg(4, 4, 4.0e-6, 0.0, 4.0e-6);
      StreamerConfig jcfg;
      jcfg.photoionization = false;
      StreamerState js(jg);
      ElectronTransportCurrentSource jc(jg);
      const double E = 2.0e6;
      const double sigma = 3.0e-4;
      double volume = 0.0;
      for (int j = 0; j < jg.nz(); ++j) {
        for (int i = 0; i < jg.nr(); ++i) {
          js.er(i, j) = 0.0;
          js.ez(i, j) = E;
          js.emag(i, j) = E;
          jc.jr(i, j) = 0.0;
          jc.jz(i, j) = sigma * E;
          volume += jg.cell_volume(i);
        }
      }
      auto jd = compute_joule_handoff_diagnostics(jg, js, jcfg, jc, 1e-12, 1e12, true, false,
                                                  JouleHandoffConfig{0.0});
      check(relerr(jd.PJ_gas_W, sigma * E * E * volume) < 1e-14 &&
                jd.PJ_negative_W == 0.0 && jd.PJ_positive_W > 0.0,
            "Joule handoff integrates uniform J dot E over axisymmetric gas volume");

      for (auto& v : jc.jz.values()) v = 0.0;
      auto zero_current = compute_joule_handoff_diagnostics(jg, js, jcfg, jc, 1e-12, 1e12, true, false);
      check(zero_current.PJ_gas_W == 0.0 && zero_current.PJ_positive_W == 0.0 &&
                zero_current.PJ_negative_W == 0.0,
            "Joule handoff excludes displacement current from changing electric field");

      jc.jz(0, 0) = E;
      jc.jz(1, 0) = -0.25 * E;
      auto signed_power = compute_joule_handoff_diagnostics(jg, js, jcfg, jc, 1e-12, 1e12, true, false);
      check(signed_power.PJ_positive_W > 0.0 && signed_power.PJ_negative_W < 0.0 &&
                relerr(signed_power.PJ_gas_W, signed_power.PJ_positive_W + signed_power.PJ_negative_W) < 1e-14,
            "Joule handoff reports positive and negative J dot E without absolute-value clamping");
    }

    {
      AxisymmetricGrid ag(2, 3, 2.0e-6, 0.0, 3.0e-6);
      StreamerConfig acfg;
      acfg.photoionization = false;
      StreamerState as(ag);
      ElectronTransportCurrentSource ac(ag);
      const double E = 1.0e6;
      const double J0 = 4.0;
      for (int j = 0; j < ag.nz(); ++j) {
        for (int i = 0; i < ag.nr(); ++i) {
          as.er(i, j) = 0.0;
          as.ez(i, j) = E;
          as.emag(i, j) = E;
          ac.jr(i, j) = 0.0;
          ac.jz(i, j) = J0;
        }
      }
      JouleHandoffAccumulator acc(JouleHandoffConfig{0.0});
      as.time = 0.0;
      auto a0 = acc.sample(ag, as, acfg, ac, 1.0e-9, 1.0e9, true, false);
      as.time = 1.0e-12;
      auto a1 = acc.sample(ag, as, acfg, ac, 2.0e-9, 5.0e8, true, false);
      as.time = 3.0e-12;
      auto a2 = acc.sample(ag, as, acfg, ac, 5.0e-9, 2.0e8, true, false);
      check(a0.energy_accumulator_status == "INITIAL_SAMPLE" &&
                relerr(a1.QJ_gas_J, a1.PJ_gas_W * 1.0e-12) < 1e-14 &&
                relerr(a2.QJ_gas_J, a2.PJ_gas_W * 3.0e-12) < 1e-14,
            "Joule accumulator handles constant power and variable accepted dt");
      const double expected_dgbdt = (5.0e-9 - 2.0e-9) / (2.0e-12);
      check(a2.tau_evolution_status == "VALID" &&
                relerr(a2.dGb_dt_S_s, expected_dgbdt) < 1e-14 &&
                relerr(a2.tau_evolution_s, std::abs(5.0e-9 / expected_dgbdt)) < 1e-14,
            "Joule handoff computes tau_evolution from accepted-state Gb history");
      check(a2.thermal_energy_reference_status == "NOT_AVAILABLE" && std::isnan(a2.Pi_H) &&
                a2.handoff_status == "UNRESOLVED_CALIBRATION",
            "Joule handoff keeps Pi_H unresolved without Q_required");
    }

    {
      AxisymmetricGrid cg2(4, 5, 4.0e-6, 0.0, 5.0e-6);
      StreamerConfig ccfg2;
      ccfg2.photoionization = false;
      StreamerState cs2(cg2);
      ElectronTransportCurrentSource current(cg2);
      for (int j = 0; j < cg2.nz(); ++j) {
        for (int i = 0; i < cg2.nr(); ++i) {
          cs2.er(i, j) = 0.0;
          cs2.ez(i, j) = 1.5e6;
          cs2.emag(i, j) = 1.5e6;
          cs2.ne(i, j) = (j == 2 && i < 2) ? 1.0e16 : 1.0e10;
          current.jz(i, j) = 2.0;
        }
      }
      auto cd = compute_joule_handoff_diagnostics(cg2, cs2, ccfg2, current, 1e-12, 1e12, true, true,
                                                  JouleHandoffConfig{0.5});
      check(cd.channel_valid && cd.bridge_flag && cd.channel_volume_m3 > 0.0 &&
                relerr(cd.channel_length_m, cg2.dz()) < 1e-14 && cd.channel_effective_radius_m > 0.0 &&
                cd.tau_sigma_status == "VALID" && cd.tau_sigma_s > 0.0,
            "Joule handoff reports channel mask, geometry, effective radius and tau_sigma");

      JouleHandoffAccumulator invalid_acc;
      cs2.time = 0.0;
      invalid_acc.sample(cg2, cs2, ccfg2, current, 0.0, std::numeric_limits<double>::infinity(), true, false);
      cs2.time = 1.0e-12;
      auto invalid = invalid_acc.sample(cg2, cs2, ccfg2, current, 0.0, std::numeric_limits<double>::infinity(), true, false);
      check(invalid.tau_evolution_status == "ZERO_GB_OR_DGBDT" && invalid.Xi_sigma_status == "UNAVAILABLE",
            "Joule handoff marks invalid zero-Gb evolution cases");
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
