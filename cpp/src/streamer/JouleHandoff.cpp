#include "streamer_rf/streamer/JouleHandoff.hpp"
#include "streamer_rf/transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace streamer_rf::streamer {
namespace {
constexpr double qe = 1.602176634e-19;
constexpr double eps0 = 8.8541878128e-12;
constexpr double pi = 3.14159265358979323846;

bool is_gas_cell(const AxisymmetricGrid& grid, const StreamerConfig& config, int i, int j) {
  return !config.electrode_geometry ||
         config.electrode_geometry->classify(grid, i, j) == ElectrodeCellType::Gas;
}

double conductivity(const StreamerState& state, const StreamerConfig& config, int i, int j) {
  const auto q = evaluate_morrow_lowke(state.emag(i, j), config.neutral_density,
                                      config.pressure, config.temperature);
  return qe * q.mobility * state.ne(i, j);
}
}  // namespace

JouleHandoffDiagnostics compute_joule_handoff_diagnostics(
    const AxisymmetricGrid& grid, const StreamerState& state, const StreamerConfig& config,
    const ElectronTransportCurrentSource& current_source, double Gb_S, double Rb_ohm,
    bool Gb_valid, bool bridge_flag, const JouleHandoffConfig& handoff_config) {
  JouleHandoffDiagnostics d;
  d.time_s = state.time;
  d.bridge_flag = bridge_flag;
  d.Gb_S = Gb_S;
  d.Rb_ohm = Rb_ohm;
  d.Gb_valid = Gb_valid && std::isfinite(Gb_S) && Gb_S >= 0.0;
  d.channel_region_definition = handoff_config.channel_region_definition;

  double sigma_max = 0.0;
  for (int j = 0; j < grid.nz(); ++j) {
    for (int i = 0; i < grid.nr(); ++i) {
      if (!is_gas_cell(grid, config, i, j)) continue;
      const double sigma = conductivity(state, config, i, j);
      if (std::isfinite(sigma)) sigma_max = std::max(sigma_max, sigma);
    }
  }

  const double channel_threshold = handoff_config.channel_sigma_relative_threshold * sigma_max;
  bool has_channel = sigma_max > 0.0 && std::isfinite(channel_threshold) &&
                     handoff_config.channel_sigma_relative_threshold >= 0.0;
  double channel_sigma_vol = 0.0, channel_E_vol = 0.0, channel_ne_vol = 0.0;
  double channel_z_min = std::numeric_limits<double>::infinity();
  double channel_z_max = -std::numeric_limits<double>::infinity();

  for (int j = 0; j < grid.nz(); ++j) {
    for (int i = 0; i < grid.nr(); ++i) {
      if (!is_gas_cell(grid, config, i, j)) continue;
      const double jr = current_source.jr(i, j);
      const double jz = current_source.jz(i, j);
      const double p = jr * state.er(i, j) + jz * state.ez(i, j);
      const double vol = grid.cell_volume(i);
      if (!std::isfinite(p) || !std::isfinite(vol)) continue;
      const double p_int = p * vol;
      d.PJ_gas_W += p_int;
      if (p_int >= 0.0) d.PJ_positive_W += p_int;
      else d.PJ_negative_W += p_int;

      const double sigma = conductivity(state, config, i, j);
      if (has_channel && std::isfinite(sigma) && sigma >= channel_threshold) {
        d.PJ_channel_W = std::isfinite(d.PJ_channel_W) ? d.PJ_channel_W + p_int : p_int;
        d.channel_volume_m3 = std::isfinite(d.channel_volume_m3) ? d.channel_volume_m3 + vol : vol;
        channel_sigma_vol += sigma * vol;
        channel_E_vol += state.emag(i, j) * vol;
        channel_ne_vol += state.ne(i, j) * vol;
        d.ne_channel_max_m3 = std::isfinite(d.ne_channel_max_m3)
                                  ? std::max(d.ne_channel_max_m3, state.ne(i, j))
                                  : state.ne(i, j);
        channel_z_min = std::min(channel_z_min, grid.z(j) - 0.5 * grid.dz());
        channel_z_max = std::max(channel_z_max, grid.z(j) + 0.5 * grid.dz());
      }
    }
  }

  if (std::isfinite(d.channel_volume_m3) && d.channel_volume_m3 > 0.0 &&
      channel_z_max > channel_z_min) {
    d.channel_valid = true;
    d.channel_status = "VALID";
    d.channel_length_m = channel_z_max - channel_z_min;
    d.channel_effective_radius_m = std::sqrt(d.channel_volume_m3 / (pi * d.channel_length_m));
    d.sigma_eff_S_m = channel_sigma_vol / d.channel_volume_m3;
    d.E_channel_mean_V_m = channel_E_vol / d.channel_volume_m3;
    d.ne_channel_mean_m3 = channel_ne_vol / d.channel_volume_m3;
    if (d.sigma_eff_S_m > 0.0 && std::isfinite(d.sigma_eff_S_m)) {
      d.tau_sigma_s = eps0 / d.sigma_eff_S_m;
      d.tau_sigma_status = "VALID";
    } else {
      d.tau_sigma_status = "INVALID_SIGMA_EFF";
    }
  } else if (has_channel) {
    d.channel_status = "EMPTY_CHANNEL_REGION";
  } else {
    d.channel_status = "NO_CONDUCTIVITY";
  }

  return d;
}

JouleHandoffAccumulator::JouleHandoffAccumulator(JouleHandoffConfig config)
    : handoff_config_(std::move(config)) {}

JouleHandoffDiagnostics JouleHandoffAccumulator::sample(
    const AxisymmetricGrid& grid, const StreamerState& state, const StreamerConfig& config,
    const ElectronTransportCurrentSource& current_source, double Gb_S, double Rb_ohm,
    bool Gb_valid, bool bridge_flag) {
  auto d = compute_joule_handoff_diagnostics(grid, state, config, current_source, Gb_S, Rb_ohm,
                                             Gb_valid, bridge_flag, handoff_config_);
  if (has_previous_power_ && state.time > previous_time_s_ && std::isfinite(d.PJ_gas_W) &&
      std::isfinite(previous_PJ_gas_W_)) {
    const double dt = state.time - previous_time_s_;
    QJ_gas_J_ += 0.5 * (d.PJ_gas_W + previous_PJ_gas_W_) * dt;
    if (d.channel_valid && previous_channel_valid_ && std::isfinite(d.PJ_channel_W) &&
        std::isfinite(previous_PJ_channel_W_)) {
      QJ_channel_J_ += 0.5 * (d.PJ_channel_W + previous_PJ_channel_W_) * dt;
      d.QJ_channel_J = QJ_channel_J_;
    } else {
      d.QJ_channel_J = std::numeric_limits<double>::quiet_NaN();
    }
    d.energy_accumulator_status = "VALID";
  } else {
    d.energy_accumulator_status = "INITIAL_SAMPLE";
    if (d.channel_valid) d.QJ_channel_J = QJ_channel_J_;
  }
  d.QJ_gas_J = QJ_gas_J_;

  if (has_previous_gb_ && state.time > previous_time_s_ && d.Gb_valid &&
      std::isfinite(previous_Gb_S_)) {
    const double dt = state.time - previous_time_s_;
    d.dGb_dt_S_s = (Gb_S - previous_Gb_S_) / dt;
    if (Gb_S > 0.0 && std::isfinite(d.dGb_dt_S_s) && d.dGb_dt_S_s != 0.0) {
      d.tau_evolution_s = std::abs(Gb_S / d.dGb_dt_S_s);
      d.tau_evolution_status = "VALID";
    } else {
      d.tau_evolution_status = "ZERO_GB_OR_DGBDT";
    }
  }

  if (d.tau_sigma_status == "VALID" && d.tau_evolution_status == "VALID") {
    d.Xi_sigma = d.tau_sigma_s / d.tau_evolution_s;
    d.Xi_sigma_status = "VALID";
  }

  previous_time_s_ = state.time;
  previous_PJ_gas_W_ = d.PJ_gas_W;
  previous_PJ_channel_W_ = d.PJ_channel_W;
  previous_channel_valid_ = d.channel_valid;
  if (d.Gb_valid) {
    previous_Gb_S_ = Gb_S;
    has_previous_gb_ = true;
  }
  has_previous_power_ = true;
  return d;
}

void JouleHandoffAccumulator::reset() {
  has_previous_power_ = false;
  has_previous_gb_ = false;
  previous_time_s_ = 0.0;
  previous_PJ_gas_W_ = 0.0;
  previous_PJ_channel_W_ = 0.0;
  previous_channel_valid_ = false;
  previous_Gb_S_ = 0.0;
  QJ_gas_J_ = 0.0;
  QJ_channel_J_ = 0.0;
}

}  // namespace streamer_rf::streamer
