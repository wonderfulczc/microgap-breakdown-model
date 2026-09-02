#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
constexpr double qe = 1.602176634e-19;
int checks = 0;

void check(bool ok, const char* name) {
  ++checks;
  if (!ok) throw std::runtime_error(name);
  std::cout << "ok " << checks << " - " << name << '\n';
}

AxisymmetricNeedlePlaneGeometry geom() {
  return AxisymmetricNeedlePlaneGeometry("stage-c3-test-needle-plane", 0.0, 75e-6, 5e-6, 2.5e-6, 5e-6);
}

StreamerConfig cfg_for(const AxisymmetricNeedlePlaneGeometry& g, const VoltageWaveform& v) {
  StreamerConfig cfg;
  cfg.electrode_geometry = &g;
  cfg.voltage_waveform = &v;
  cfg.photoionization = false;
  cfg.n_ref = 1e12;
  cfg.elliptic = {1e-10, 1e-14, 20000};
  return cfg;
}

void set_uniform_state(StreamerSolver& solver, const AxisymmetricGrid& g, double ne, double e) {
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      solver.state().ne(i, j) = ne;
      solver.state().np(i, j) = ne;
      solver.state().nn(i, j) = 0.0;
      solver.state().er(i, j) = 0.0;
      solver.state().ez(i, j) = e;
      solver.state().emag(i, j) = std::abs(e);
      solver.state().sph(i, j) = 0.0;
    }
  }
}

void set_positive_gaussian_space_charge(StreamerSolver& solver, const AxisymmetricGrid& g,
                                        const AxisymmetricNeedlePlaneGeometry& eg, double z0) {
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      const bool gas = eg.classify(g, i, j) == ElectrodeCellType::Gas;
      const double rr = g.r(i);
      const double zz = g.z(j) - z0;
      const double n = gas ? 5e20 * std::exp(-(rr * rr + zz * zz) / std::pow(8e-6, 2)) : 0.0;
      solver.state().ne(i, j) = 0.0;
      solver.state().np(i, j) = n;
      solver.state().nn(i, j) = 0.0;
      solver.state().sph(i, j) = 0.0;
    }
  }
  solver.refresh_electrostatic_fields();
}
}

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  try {
    AxisymmetricGrid area_grid(4, 5, 2.0, 0.0, 5.0);
    check(std::abs(axisymmetric_radial_face_area(area_grid, 2) - 2.0 * M_PI) < 1e-14, "electrode radial face area");
    check(std::abs(axisymmetric_axial_face_area(area_grid, 1) - 0.75 * M_PI) < 1e-14, "electrode axial face area");

    AxisymmetricGrid g(20, 48, 80e-6, 0.0, 90e-6);
    auto eg = geom();
    ConstantVoltage cv(100.0);
    auto cfg = cfg_for(eg, cv);

    StreamerSolver vacuum(g, cfg);
    vacuum.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    const double cgap = vacuum.vacuum_gap_capacitance();
    auto q0 = vacuum.electrode_surface_diagnostics();
    check(cgap > 0.0 && std::abs(std::abs(q0.q_hv / cv.value(0.0)) - cgap) / cgap < 1e-10,
          "electrode charge and vacuum capacitance");
    StreamerDiagnostics vd;
    check(vacuum.step(1e-12, vd), "constant-voltage vacuum step accepted");
    check(std::abs(vd.i_cond_hv) < 1e-30 && std::abs(vd.i_cond_ground) < 1e-30, "vacuum constant-V conduction current zero");
    check(std::abs(vd.i_disp_hv) < 1e-12 && std::abs(vd.i_disp_ground) < 1e-12, "vacuum constant-V displacement current zero");

    const double ramp_time = 1e-9, ramp_voltage = 100.0, dt = 1e-10;
    SampledVoltage ramp({0.0, ramp_time}, {0.0, ramp_voltage}, WaveformOutOfRangePolicy::HoldEndpoint);
    auto rcfg = cfg_for(eg, ramp);
    StreamerSolver ramp_solver(g, rcfg);
    ramp_solver.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    const double ramp_cgap = ramp_solver.vacuum_gap_capacitance();
    StreamerDiagnostics rd;
    check(ramp_solver.step(dt, rd), "linear-ramp vacuum step accepted");
    const double predicted = ramp_cgap * ramp_voltage / ramp_time;
    check(predicted > 0.0 && std::abs(rd.i_disp_hv - predicted) / predicted < 5e-5,
          "vacuum linear ramp displacement current matches C dV/dt");
    check(rd.i_disp_hv > 0.0 && rd.i_disp_ground < 0.0, "displacement current sign convention");

    StreamerSolver fixed_rho(g, cfg);
    fixed_rho.initialize_gaussian_at_tip_offset(0.0, 3e-6, -10e-6);
    set_positive_gaussian_space_charge(fixed_rho, g, eg, 38e-6);
    const auto qa = fixed_rho.electrode_surface_diagnostics();
    fixed_rho.reset_electrode_history();
    set_positive_gaussian_space_charge(fixed_rho, g, eg, 54e-6);
    const auto qb = fixed_rho.electrode_surface_diagnostics();
    const double synthetic_dt = 1e-12;
    const auto sd = fixed_rho.sample_terminal_diagnostics_from_history(synthetic_dt);
    const double expected_disp = (qb.q_hv - qa.q_hv) / synthetic_dt;
    check(std::abs(qb.q_hv - qa.q_hv) > 0.0, "fixed-V changing rho changes induced HV charge");
    check(std::abs(sd.i_disp_hv - expected_disp) / std::max(std::abs(expected_disp), 1e-300) < 1e-12,
          "fixed-V changing rho displacement current equals dQ/dt");

    StreamerSolver absorption(g, cfg);
    absorption.initialize_gaussian_at_tip_offset(1e16, 3e-6, -1e-6);
    StreamerDiagnostics ad;
    check(absorption.step(1e-16, ad), "absorption current step accepted");
    check(std::abs(ad.i_cond_hv - qe * ad.absorbed_electron_hv / ad.dt) /
              std::max(std::abs(ad.i_cond_hv), 1e-30) < 1e-14,
          "electron absorption flux matches HV conduction current");

    AxisymmetricGrid ug(24, 20, 120e-6, 0.0, 100e-6);
    StreamerConfig ucfg;
    ucfg.photoionization = false;
    StreamerSolver uniform(ug, ucfg);
    uniform.initialize_gaussian(0.0, 1e-6, 50e-6);
    const double ne = 1e15, efield = 2e6, length = ug.nz() * ug.dz();
    set_uniform_state(uniform, ug, ne, efield);
    const double voltage = efield * length;
    auto q = evaluate_morrow_lowke(efield, ucfg.neutral_density, ucfg.pressure, ucfg.temperature);
    const double sigma = qe * q.mobility * ne;
    auto rf_source = uniform.electron_transport_current_source();
    check(std::isfinite(rf_source.current_moment_z) && rf_source.integral_abs_jz > 0.0,
          "PETSc flux-derived RF current source finite");
    check(std::abs(rf_source.jz(5, 5) - sigma * efield) / std::max(std::abs(sigma * efield), 1e-300) < 1e-12,
          "PETSc RF current source matches uniform drift face flux");
    double mz_sum = 0.0;
    for (int jj = 0; jj < ug.nz(); ++jj) {
      for (int ii = 0; ii < ug.nr(); ++ii) mz_sum += rf_source.jz(ii, jj) * ug.cell_volume(ii);
    }
    check(std::abs(mz_sum - rf_source.current_moment_z) / std::max(std::abs(mz_sum), 1e-300) < 1e-14,
          "PETSc RF current moment equals cell-volume integral");
    const double analytic_g = sigma * M_PI * std::pow(ug.nr() * ug.dr(), 2) / length;
    auto gd = uniform.conductance_diagnostics(voltage);
    check(gd.valid && std::abs(gd.gb - analytic_g) / analytic_g < 1e-12, "uniform sigma Gb equals sigma A over d");
    check(std::abs(gd.rb - 1.0 / gd.gb) / gd.rb < 1e-14, "Rb is inverse Gb");
    auto invalid = uniform.conductance_diagnostics(0.0);
    check(!invalid.valid && std::isnan(invalid.gb) && std::isnan(invalid.rb), "near-zero voltage invalidates Gb/Rb");

    StreamerSolver continuity(g, cfg);
    continuity.initialize_gaussian_at_tip_offset(1e14, 3e-6, -10e-6);
    StreamerDiagnostics cd;
    check(continuity.step(1e-16, cd), "current continuity step accepted");
    check(std::isfinite(cd.current_continuity_residual) && cd.current_continuity_residual < 1e-2,
          "current continuity residual finite");

    StreamerConfig legacy_cfg;
    legacy_cfg.photoionization = false;
    StreamerSolver legacy(ug, legacy_cfg);
    legacy.initialize_gaussian(1e12, 10e-6, 50e-6);
    auto lim = legacy.timestep_limits();
    StreamerDiagnostics ld;
    check(legacy.step(0.01 * lim.selected, ld) && std::isfinite(ld.ne_max), "legacy Stage 1-5 path still steps");

    std::cout << "passed " << checks << " Stage C3 C++ checks\n";
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
  PetscFinalize();
  return 0;
}
