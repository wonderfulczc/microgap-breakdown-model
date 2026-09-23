#include "streamer_rf/streamer/UltrafastEventDiagnostics.hpp"
#include "streamer_rf/streamer/MorrowLowke.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace streamer_rf::streamer {
namespace {
constexpr double elementary_charge = 1.602176634e-19;
constexpr double epsilon_0 = 8.8541878128e-12;

struct Component {
  std::vector<std::pair<int, int>> cells;
  double weight{};
  double centroid_r{};
  double centroid_z{};
  double peak_abs_rho{};
};

std::size_t flat_index(const AxisymmetricGrid& grid, int i, int j) {
  return static_cast<std::size_t>(j) * grid.nr() + i;
}

bool is_gas(const AxisymmetricGrid& grid, const StreamerConfig& config, int i, int j) {
  return !config.electrode_geometry ||
         config.electrode_geometry->classify(grid, i, j) == ElectrodeCellType::Gas;
}

double quantile(std::vector<double> values, double q) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double position = q * static_cast<double>(values.size() - 1);
  const auto lower = static_cast<std::size_t>(std::floor(position));
  const auto upper = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - static_cast<double>(lower);
  return values[lower] + fraction * (values[upper] - values[lower]);
}

std::vector<Component> detect_components(const AxisymmetricGrid& grid, const StreamerState& state,
                                         const StreamerConfig& config,
                                         const UltrafastEventConfig& event_config) {
  double peak = 0.0;
  for (int j = 0; j < grid.nz(); ++j) {
    for (int i = 0; i < grid.nr(); ++i) {
      if (!is_gas(grid, config, i, j)) continue;
      const double rho = state.rho(i, j);
      if (event_config.requested_polarity == 0 || event_config.requested_polarity * rho > 0.0)
        peak = std::max(peak, std::abs(rho));
    }
  }
  if (!(peak > 0.0)) return {};
  const double threshold = event_config.relative_rho_threshold * peak;
  auto selected = [&](int i, int j) {
    if (!is_gas(grid, config, i, j)) return false;
    const double rho = state.rho(i, j);
    return std::abs(rho) >= threshold &&
           (event_config.requested_polarity == 0 || event_config.requested_polarity * rho > 0.0);
  };
  std::vector<unsigned char> visited(grid.size(), 0);
  std::vector<Component> components;
  for (int seed_j = 0; seed_j < grid.nz(); ++seed_j) {
    for (int seed_i = 0; seed_i < grid.nr(); ++seed_i) {
      const auto seed_k = flat_index(grid, seed_i, seed_j);
      if (visited[seed_k] || !selected(seed_i, seed_j)) continue;
      Component component;
      std::vector<std::pair<int, int>> stack{{seed_i, seed_j}};
      visited[seed_k] = 1;
      double rsum = 0.0, zsum = 0.0;
      while (!stack.empty()) {
        const auto [i, j] = stack.back();
        stack.pop_back();
        component.cells.emplace_back(i, j);
        const double weight = std::abs(state.rho(i, j)) * grid.cell_volume(i);
        component.weight += weight;
        rsum += weight * grid.r(i);
        zsum += weight * grid.z(j);
        component.peak_abs_rho = std::max(component.peak_abs_rho, std::abs(state.rho(i, j)));
        constexpr int di[4] = {-1, 1, 0, 0};
        constexpr int dj[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
          const int ni = i + di[n], nj = j + dj[n];
          if (ni < 0 || ni >= grid.nr() || nj < 0 || nj >= grid.nz()) continue;
          const auto k = flat_index(grid, ni, nj);
          if (!visited[k] && selected(ni, nj)) {
            visited[k] = 1;
            stack.emplace_back(ni, nj);
          }
        }
      }
      if (static_cast<int>(component.cells.size()) >= event_config.minimum_component_cells &&
          component.weight > 0.0) {
        component.centroid_r = rsum / component.weight;
        component.centroid_z = zsum / component.weight;
        components.push_back(std::move(component));
      }
    }
  }
  std::sort(components.begin(), components.end(), [](const Component& a, const Component& b) {
    return a.peak_abs_rho > b.peak_abs_rho;
  });
  return components;
}

double component_distance(const Component& a, const Component& b) {
  return std::hypot(a.centroid_r - b.centroid_r, a.centroid_z - b.centroid_z);
}

double electrode_distance(const AxisymmetricGrid& grid, const StreamerConfig& config,
                          const std::vector<std::pair<int, int>>& roi) {
  if (!config.electrode_geometry || roi.empty()) return std::numeric_limits<double>::quiet_NaN();
  double minimum = std::numeric_limits<double>::infinity();
  for (const auto& [ri, rj] : roi) {
    for (int j = 0; j < grid.nz(); ++j) {
      for (int i = 0; i < grid.nr(); ++i) {
        if (config.electrode_geometry->classify(grid, i, j) == ElectrodeCellType::Gas) continue;
        minimum = std::min(minimum, std::hypot(grid.r(ri) - grid.r(i), grid.z(rj) - grid.z(j)));
      }
    }
  }
  return minimum;
}

double valid_inverse(double value) {
  return std::isfinite(value) && value > 0.0 ? 1.0 / value
                                             : std::numeric_limits<double>::quiet_NaN();
}
}  // namespace

const char* event_topology_name(EventTopology topology) {
  switch (topology) {
    case EventTopology::Propagation: return "PROPAGATION";
    case EventTopology::HeadInteraction: return "HEAD_INTERACTION";
    case EventTopology::StreamerCollision: return "STREAMER_COLLISION";
    case EventTopology::ElectrodeAttachment: return "ELECTRODE_ATTACHMENT";
    case EventTopology::Bridging: return "BRIDGING";
    case EventTopology::Unresolved: return "UNRESOLVED";
  }
  return "UNRESOLVED";
}

UltrafastEventDiagnostics::UltrafastEventDiagnostics(
    const AxisymmetricGrid& grid, const StreamerConfig& config,
    const ScalarField2D& electrostatic_reference_E, UltrafastEventConfig event_config)
    : grid_(grid), config_(config), E0_(electrostatic_reference_E), event_config_(event_config),
      Ek_V_m_(morrow_lowke_breakdown_field(config.neutral_density, config.pressure,
                                           config.temperature)) {
  if (&E0_.grid() != &grid_ || E0_.values().size() != grid_.size())
    throw std::invalid_argument("electrostatic reference grid mismatch");
  if (!(event_config_.relative_rho_threshold > 0.0 && event_config_.relative_rho_threshold <= 1.0) ||
      event_config_.minimum_component_cells < 1 || event_config_.confirmation_samples < 1)
    throw std::invalid_argument("invalid ultrafast event configuration");
}

std::string UltrafastEventDiagnostics::kinetic_resolution_status(double timescale_s,
                                                                 double local_dt_s) {
  if (!(std::isfinite(timescale_s) && timescale_s > 0.0 &&
        std::isfinite(local_dt_s) && local_dt_s > 0.0))
    return "NOT_RESOLVED";
  const double samples = timescale_s / local_dt_s;
  if (samples >= 10.0) return "RESOLVED_PREFERRED";
  if (samples >= 5.0) return "RESOLVED_MINIMUM";
  return "NOT_RESOLVED_KINETIC_TIMESCALE";
}

UltrafastEventSample UltrafastEventDiagnostics::sample(
    int step, const StreamerState& state, const StreamerHeadDiagnostics& primary_head,
    const ElectronTransportCurrentSource& current_source, bool bridge_status,
    double accepted_dt_s) {
  if (!(std::isfinite(accepted_dt_s) && accepted_dt_s > 0.0))
    throw std::invalid_argument("accepted diagnostic dt must be finite and positive");
  UltrafastEventSample out;
  out.step = step;
  out.time_s = state.time;
  out.dt_s = accepted_dt_s;
  out.Ek_V_m = Ek_V_m_;
  out.pressure_Pa = config_.pressure;
  out.temperature_K = config_.temperature;
  out.head_position_m = primary_head.head_z_m;
  out.head_velocity_m_s = primary_head.head_velocity_z_m_s;
  out.bridge_status = bridge_status;
  out.current_moment_r_A_m = current_source.current_moment_r;
  out.current_moment_z_A_m = current_source.current_moment_z;

  const auto components = detect_components(grid_, state, config_, event_config_);
  out.head_count = static_cast<int>(components.size());
  const double cell_scale = std::max(grid_.dr(), grid_.dz());
  const double collision_threshold = event_config_.collision_distance_cells * cell_scale;
  const double interaction_threshold = event_config_.interaction_distance_cells * cell_scale;
  double closest_component_distance = std::numeric_limits<double>::infinity();
  int second_component = -1;
  if (components.size() >= 2) {
    for (std::size_t index = 1; index < components.size(); ++index) {
      const double distance = component_distance(components[0], components[index]);
      if (distance < closest_component_distance) {
        closest_component_distance = distance;
        second_component = static_cast<int>(index);
      }
    }
  }

  std::vector<std::pair<int, int>> roi;
  if (!components.empty()) roi = components[0].cells;
  const double primary_electrode_distance = electrode_distance(grid_, config_, roi);
  const double attachment_threshold = event_config_.attachment_distance_cells * cell_scale;

  EventTopology candidate = EventTopology::Unresolved;
  if (bridge_status && !components.empty()) {
    candidate = EventTopology::Bridging;
    out.topology_confidence = "HIGH";
    out.topology_evidence = "FROZEN_STAGE_C_BRIDGE_CONDITION";
  } else if (components.size() >= 2 && closest_component_distance <= collision_threshold) {
    candidate = EventTopology::StreamerCollision;
    out.topology_confidence = "HIGH";
    out.topology_evidence = "TWO_DISTINCT_COMPONENTS_WITHIN_COLLISION_DISTANCE";
    roi.insert(roi.end(), components[second_component].cells.begin(), components[second_component].cells.end());
    out.roi_definition = "EVENT_INTERACTION_ROI";
  } else if (!components.empty() && std::isfinite(primary_electrode_distance) &&
             primary_electrode_distance <= attachment_threshold) {
    candidate = EventTopology::ElectrodeAttachment;
    out.topology_confidence = "HIGH";
    out.topology_evidence = "HEAD_COMPONENT_WITHIN_ELECTRODE_DISTANCE";
    out.roi_definition = "ATTACHMENT_ROI";
  } else if (components.size() >= 2 && closest_component_distance <= interaction_threshold) {
    candidate = EventTopology::HeadInteraction;
    out.topology_confidence = "MEDIUM";
    out.topology_evidence = "TWO_DISTINCT_COMPONENTS_WITHIN_INTERACTION_DISTANCE";
    roi.insert(roi.end(), components[second_component].cells.begin(), components[second_component].cells.end());
    out.roi_definition = "EVENT_INTERACTION_ROI";
  } else if (components.size() == 1) {
    candidate = EventTopology::Propagation;
    out.topology_confidence = "MEDIUM";
    out.topology_evidence = "ONE_RESOLVED_HEAD_COMPONENT";
  } else if (components.size() >= 2) {
    out.topology_confidence = "LOW";
    out.topology_evidence = "MULTIPLE_COMPONENTS_WITHOUT_SPATIAL_EVENT_EVIDENCE";
  }
  if (out.roi_definition == "NOT_AVAILABLE" && !roi.empty()) out.roi_definition = "HEAD_COMPONENT_ROI";
  if (candidate == pending_topology_) ++pending_count_;
  else {
    pending_topology_ = candidate;
    pending_count_ = 1;
  }
  out.topology = candidate;
  if (candidate == EventTopology::Unresolved) out.topology_status = "UNRESOLVED";
  else if (pending_count_ >= event_config_.confirmation_samples) out.topology_status = "CONFIRMED";
  else out.topology_status = "PENDING_CONTIGUOUS_CONFIRMATION";

  if (roi.empty()) return out;
  out.roi_status = "VALID";
  out.roi_cell_count = static_cast<int>(roi.size());
  out.electrode_distance_m = electrode_distance(grid_, config_, roi);
  out.roi_r_min_m = out.roi_z_min_m = std::numeric_limits<double>::infinity();
  out.roi_r_max_m = out.roi_z_max_m = -std::numeric_limits<double>::infinity();
  double volume_r = 0.0, volume_z = 0.0;
  std::vector<double> E0_values, tau_i_values, sigma_values, tau_M_values, pi_values;
  int e_peak_i = -1, e_peak_j = -1, ne_peak_i = -1, ne_peak_j = -1;
  for (const auto& [i, j] : roi) {
    const double volume = grid_.cell_volume(i);
    out.roi_volume_m3 += volume;
    volume_r += volume * grid_.r(i);
    volume_z += volume * grid_.z(j);
    out.roi_r_min_m = std::min(out.roi_r_min_m, grid_.r(i));
    out.roi_r_max_m = std::max(out.roi_r_max_m, grid_.r(i));
    out.roi_z_min_m = std::min(out.roi_z_min_m, grid_.z(j));
    out.roi_z_max_m = std::max(out.roi_z_max_m, grid_.z(j));
    E0_values.push_back(E0_(i, j));
    if (e_peak_i < 0 || state.emag(i, j) > state.emag(e_peak_i, e_peak_j)) {
      e_peak_i = i;
      e_peak_j = j;
    }
    if (ne_peak_i < 0 || state.ne(i, j) > state.ne(ne_peak_i, ne_peak_j)) {
      ne_peak_i = i;
      ne_peak_j = j;
    }
    const auto transport = evaluate_morrow_lowke(state.emag(i, j), config_.neutral_density,
                                                  config_.pressure, config_.temperature);
    const double nu_i = transport.ionization_frequency;
    const double sigma_e = elementary_charge * state.ne(i, j) * transport.mobility;
    const double tau_i = valid_inverse(nu_i);
    const double tau_M = std::isfinite(sigma_e) && sigma_e > 0.0
                             ? epsilon_0 / sigma_e
                             : std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(tau_i)) tau_i_values.push_back(tau_i);
    if (std::isfinite(sigma_e) && sigma_e > 0.0) sigma_values.push_back(sigma_e);
    if (std::isfinite(tau_M)) tau_M_values.push_back(tau_M);
    if (std::isfinite(tau_i) && std::isfinite(tau_M)) pi_values.push_back(tau_M / tau_i);
  }
  out.roi_centroid_r_m = volume_r / out.roi_volume_m3;
  out.roi_centroid_z_m = volume_z / out.roi_volume_m3;
  out.E0_max_V_m = *std::max_element(E0_values.begin(), E0_values.end());
  out.E0_mean_V_m = std::accumulate(E0_values.begin(), E0_values.end(), 0.0) / E0_values.size();
  out.E0_p95_V_m = quantile(E0_values, 0.95);
  out.eta_0_max = out.E0_max_V_m / Ek_V_m_;
  out.eta_0_mean = out.E0_mean_V_m / Ek_V_m_;
  out.eta_0_p95 = out.E0_p95_V_m / Ek_V_m_;

  out.E_peak_V_m = state.emag(e_peak_i, e_peak_j);
  out.E_peak_time_s = state.time;
  out.E_peak_r_m = grid_.r(e_peak_i);
  out.E_peak_z_m = grid_.z(e_peak_j);
  out.ne_at_E_peak_m3 = state.ne(e_peak_i, e_peak_j);
  out.eta_peak = out.E_peak_V_m / Ek_V_m_;
  out.ne_peak_m3 = state.ne(ne_peak_i, ne_peak_j);
  out.ne_peak_r_m = grid_.r(ne_peak_i);
  out.ne_peak_z_m = grid_.z(ne_peak_j);

  auto at_cell = [&](int i, int j, double& mu, double& nu, double& tau_i, double& sigma,
                     double& tau_M, double& pi) {
    const auto transport = evaluate_morrow_lowke(state.emag(i, j), config_.neutral_density,
                                                  config_.pressure, config_.temperature);
    mu = transport.mobility;
    nu = transport.ionization_frequency;
    tau_i = valid_inverse(nu);
    sigma = elementary_charge * state.ne(i, j) * mu;
    tau_M = std::isfinite(sigma) && sigma > 0.0
                ? epsilon_0 / sigma
                : std::numeric_limits<double>::quiet_NaN();
    pi = std::isfinite(tau_i) && std::isfinite(tau_M)
             ? tau_M / tau_i
             : std::numeric_limits<double>::quiet_NaN();
  };
  double unused_mu{}, unused_nu{};
  at_cell(e_peak_i, e_peak_j, out.mu_e_at_E_peak_m2_V_s, out.nu_i_at_E_peak_s_1,
          out.tau_i_at_E_peak_s, out.sigma_e_at_E_peak_S_m, out.tau_M_at_E_peak_s,
          out.Pi_RF_at_E_peak);
  at_cell(ne_peak_i, ne_peak_j, unused_mu, unused_nu, out.tau_i_at_ne_peak_s,
          out.sigma_e_at_ne_peak_S_m, out.tau_M_at_ne_peak_s, out.Pi_RF_at_ne_peak);
  out.tau_i_min_s = quantile(tau_i_values, 0.0);
  out.tau_i_median_s = quantile(tau_i_values, 0.5);
  out.tau_i_p95_s = quantile(tau_i_values, 0.95);
  out.sigma_e_peak_S_m = quantile(sigma_values, 1.0);
  out.sigma_e_median_S_m = quantile(sigma_values, 0.5);
  out.sigma_e_p95_S_m = quantile(sigma_values, 0.95);
  out.tau_M_min_s = quantile(tau_M_values, 0.0);
  out.tau_M_median_s = quantile(tau_M_values, 0.5);
  out.tau_M_p95_s = quantile(tau_M_values, 0.95);
  out.Pi_RF_min = quantile(pi_values, 0.0);
  out.Pi_RF_median = quantile(pi_values, 0.5);
  out.Pi_RF_p95 = quantile(pi_values, 0.95);
  out.N_tau_i_at_E_peak = out.tau_i_at_E_peak_s / accepted_dt_s;
  out.N_tau_M_at_E_peak = out.tau_M_at_E_peak_s / accepted_dt_s;
  out.tau_i_resolution_status = kinetic_resolution_status(out.tau_i_at_E_peak_s, accepted_dt_s);
  out.tau_M_resolution_status = kinetic_resolution_status(out.tau_M_at_E_peak_s, accepted_dt_s);
  return out;
}

}  // namespace streamer_rf::streamer
