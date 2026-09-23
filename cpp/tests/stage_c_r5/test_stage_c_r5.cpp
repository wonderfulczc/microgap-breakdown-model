#include "streamer_rf/streamer/ReactionModel.hpp"
#include "streamer_rf/streamer/UltrafastEventDiagnostics.hpp"
#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace streamer_rf;
using namespace streamer_rf::streamer;

namespace {
int checks = 0;
void check(bool condition, const char* name) {
  ++checks;
  if (!condition) throw std::runtime_error(name);
  std::cout << "ok " << checks << " - " << name << '\n';
}
double relerr(double value, double reference) {
  return std::abs(value - reference) / std::max(std::abs(reference), 1e-300);
}
void clear(StreamerState& state, const AxisymmetricGrid& grid, double field = 3e6) {
  for (int j = 0; j < grid.nz(); ++j) for (int i = 0; i < grid.nr(); ++i) {
    state.ne(i, j) = state.np(i, j) = state.nn(i, j) = state.rho(i, j) = 0.0;
    state.er(i, j) = 0.0;
    state.ez(i, j) = state.emag(i, j) = field;
  }
}
void blob(StreamerState& state, int i, int j, double rho, double ne) {
  state.rho(i, j) = rho;
  state.rho(i + 1, j) = rho;
  state.ne(i, j) = state.np(i, j) = ne;
  state.ne(i + 1, j) = state.np(i + 1, j) = ne;
}
}  // namespace

int main(int argc, char** argv) {
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  try {
    AxisymmetricGrid grid(12, 24, 120e-6, 0.0, 240e-6);
    StreamerConfig config;
    config.photoionization = false;
    StreamerState state(grid);
    ScalarField2D E0(grid);
    for (double& value : E0.values()) value = 2e6;
    UltrafastEventConfig event_config;
    event_config.confirmation_samples = 1;
    event_config.requested_polarity = 1;
    event_config.collision_distance_cells = 4.0;
    event_config.interaction_distance_cells = 10.0;
    UltrafastEventDiagnostics diagnostic(grid, config, E0, event_config);
    check(diagnostic.Ek_V_m() == morrow_lowke_breakdown_field(config.neutral_density, config.pressure, config.temperature),
          "Ek reuses frozen Morrow-Lowke implementation");

    const double E = 5e6, ne = 2e17;
    const auto transport = evaluate_morrow_lowke(E, config.neutral_density, config.pressure, config.temperature);
    const auto sources = electron_reaction_source_components(ne, 0.0, 0.0, 0.0, transport, config.temperature);
    check(sources.impact_source_m3s / ne == transport.ionization_frequency,
          "diagnostic ionization frequency is source-term consistent");
    const double sigma = 1.602176634e-19 * ne * transport.mobility;
    const double tau_M = 8.8541878128e-12 / sigma;
    const double tau_i = 1.0 / transport.ionization_frequency;
    check(sigma > 0.0 && tau_M > 0.0 && tau_i > 0.0 &&
              relerr((tau_M / tau_i), tau_M * transport.ionization_frequency) < 1e-15,
          "conductivity tauM and same-cell Pi_RF arithmetic");

    ElectronTransportCurrentSource current(grid);
    for (int j = 0; j < grid.nz(); ++j) for (int i = 0; i < grid.nr(); ++i) {
      current.jz(i, j) = 2.0 + i + 0.5 * j;
      current.current_moment_z += current.jz(i, j) * grid.cell_volume(i);
    }
    StreamerHeadDiagnostics head;
    head.head_valid = true;
    head.head_z_m = grid.z(10);
    head.head_velocity_z_m_s = 1e6;
    state.time = 1e-12;
    clear(state, grid, E);
    blob(state, 2, 10, 1e-3, ne);
    auto one = diagnostic.sample(1, state, head, current, false, 0.5e-12);
    check(one.topology == EventTopology::Propagation && one.head_count == 1,
          "single head conservatively classified as propagation");
    check(one.roi_status == "VALID" && one.roi_definition == "HEAD_COMPONENT_ROI" &&
              one.roi_cell_count == 2,
          "single-head component ROI");
    check(relerr(one.nu_i_at_E_peak_s_1, transport.ionization_frequency) < 1e-15 &&
              relerr(one.mu_e_at_E_peak_m2_V_s, transport.mobility) < 1e-15,
          "nu_i and mobility reuse frozen transport coefficients");
    check(relerr(one.sigma_e_at_E_peak_S_m, sigma) < 1e-15 &&
              relerr(one.tau_M_at_E_peak_s, tau_M) < 1e-15 &&
              relerr(one.tau_i_at_E_peak_s, tau_i) < 1e-15 &&
              relerr(one.Pi_RF_at_E_peak, tau_M / tau_i) < 1e-15,
          "local kinetics diagnostics are same-cell values");
    const double expected_proxy = 1.602176634e-19 * (-transport.mobility * E) *
                                  transport.ionization_frequency * ne *
                                  (grid.cell_volume(2) + grid.cell_volume(3));
    check(relerr(one.K_ion_z_A_m_s, expected_proxy) < 1e-15 &&
              relerr(one.K_ion_abs_A_m_s, std::abs(expected_proxy)) < 1e-15 &&
              one.signed_proxy_status == "DEFINED_FROM_FROZEN_ELECTRON_DRIFT_CONVENTION",
          "Koile ionization proxy preserves electron drift sign and ROI integral");
    check(one.current_moment_z_A_m == current.current_moment_z &&
              one.current_definition == "FROZEN_FINITE_VOLUME_ELECTRON_TRANSPORT_CURRENT",
          "current moment handoff is identical to frozen Stage-F source value");
    check(one.eta_0_mean == 2e6 / diagnostic.Ek_V_m() &&
              one.eta_peak == E / diagnostic.Ek_V_m(),
          "eta0 and eta_peak use reference-only Ek semantics");

    clear(state, grid, E);
    blob(state, 2, 10, 1e-3, ne);
    blob(state, 2, 13, 0.9e-3, ne);
    auto two = diagnostic.sample(2, state, head, current, false, 0.5e-12);
    check(two.topology == EventTopology::StreamerCollision && two.head_count == 2 &&
              two.roi_definition == "EVENT_INTERACTION_ROI",
          "collision requires two distinct close components");

    AxisymmetricNeedlePlaneGeometry geometry("c-r5-attachment-test", 0.0, 220e-6,
                                              10e-6, 5e-6, 10e-6);
    StreamerConfig attachment_config = config;
    attachment_config.electrode_geometry = &geometry;
    StreamerState attachment_state(grid);
    clear(attachment_state, grid, E);
    attachment_state.time = state.time;
    blob(attachment_state, 2, 1, 1e-3, ne);
    UltrafastEventDiagnostics attachment_diagnostic(grid, attachment_config, E0, event_config);
    auto attachment = attachment_diagnostic.sample(1, attachment_state, head, current, false, 0.5e-12);
    check(attachment.topology == EventTopology::ElectrodeAttachment &&
              attachment.roi_definition == "ATTACHMENT_ROI" &&
              attachment.electrode_distance_m <= 1.5 * std::max(grid.dr(), grid.dz()),
          "head near electrode selects attachment ROI");

    clear(state, grid, 0.0);
    blob(state, 2, 10, 1e-3, 0.0);
    auto invalid = diagnostic.sample(3, state, head, current, false, 0.5e-12);
    check(!std::isfinite(invalid.tau_i_at_E_peak_s) &&
              !std::isfinite(invalid.tau_M_at_E_peak_s) &&
              !std::isfinite(invalid.Pi_RF_at_E_peak),
          "zero nu_i and conductivity remain invalid instead of clipped");

    auto bridged = diagnostic.sample(4, state, head, current, true, 0.5e-12);
    check(bridged.topology == EventTopology::Bridging,
          "frozen bridge condition has explicit priority");
    check(UltrafastEventDiagnostics::kinetic_resolution_status(10e-12, 0.5e-12) == "RESOLVED_PREFERRED" &&
              UltrafastEventDiagnostics::kinetic_resolution_status(10e-12, 1.5e-12) == "RESOLVED_MINIMUM" &&
              UltrafastEventDiagnostics::kinetic_resolution_status(10e-12, 3e-12) == "NOT_RESOLVED_KINETIC_TIMESCALE",
          "kinetic temporal resolution thresholds");
    check(std::string(event_topology_name(EventTopology::Unresolved)) == "UNRESOLVED" &&
              std::string(event_topology_name(EventTopology::ElectrodeAttachment)) == "ELECTRODE_ATTACHMENT",
          "frozen event topology enum is complete");
    std::cout << "stage_c_r5 checks=" << checks << " status=PASS\n";
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    PetscFinalize();
    return 1;
  }
  PetscFinalize();
  return 0;
}
