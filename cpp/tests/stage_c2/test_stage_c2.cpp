#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/transport.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

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

    std::cout << "passed " << checks << " Stage C2 C++ checks\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
  PetscFinalize();
  return 0;
}
