#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <limits>
#include <string>

namespace streamer_rf::streamer {

struct JouleHandoffConfig {
  double channel_sigma_relative_threshold{0.1};
  std::string channel_region_definition{
      "NUMERICAL_CHANNEL_DIAGNOSTIC_REGION:sigma>=channel_sigma_relative_threshold*sigma_max"};
};

struct JouleHandoffDiagnostics {
  double time_s{};
  bool bridge_flag{};
  double PJ_gas_W{};
  double PJ_channel_W{std::numeric_limits<double>::quiet_NaN()};
  double PJ_positive_W{};
  double PJ_negative_W{};
  double QJ_gas_J{};
  double QJ_channel_J{std::numeric_limits<double>::quiet_NaN()};
  bool channel_valid{};
  std::string channel_status{"NO_CHANNEL_REGION"};
  double channel_volume_m3{std::numeric_limits<double>::quiet_NaN()};
  double channel_length_m{std::numeric_limits<double>::quiet_NaN()};
  double channel_effective_radius_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_eff_S_m{std::numeric_limits<double>::quiet_NaN()};
  double E_channel_mean_V_m{std::numeric_limits<double>::quiet_NaN()};
  double ne_channel_mean_m3{std::numeric_limits<double>::quiet_NaN()};
  double ne_channel_max_m3{std::numeric_limits<double>::quiet_NaN()};
  double Gb_S{std::numeric_limits<double>::quiet_NaN()};
  double Rb_ohm{std::numeric_limits<double>::quiet_NaN()};
  bool Gb_valid{};
  double dGb_dt_S_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_sigma_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_evolution_s{std::numeric_limits<double>::quiet_NaN()};
  double Xi_sigma{std::numeric_limits<double>::quiet_NaN()};
  std::string tau_sigma_status{"INVALID_CHANNEL"};
  std::string tau_evolution_status{"INSUFFICIENT_GB_HISTORY"};
  std::string Xi_sigma_status{"UNAVAILABLE"};
  std::string thermal_energy_reference_status{"NOT_AVAILABLE"};
  double Pi_H{std::numeric_limits<double>::quiet_NaN()};
  std::string handoff_status{"UNRESOLVED_CALIBRATION"};
  std::string energy_accumulator_status{"INITIAL_SAMPLE"};
  std::string joule_current_definition{"J_cond=-e*Gamma_e finite-volume electron transport flux"};
  std::string channel_region_definition;
};

JouleHandoffDiagnostics compute_joule_handoff_diagnostics(
    const AxisymmetricGrid& grid, const StreamerState& state, const StreamerConfig& config,
    const ElectronTransportCurrentSource& current_source, double Gb_S, double Rb_ohm,
    bool Gb_valid, bool bridge_flag, const JouleHandoffConfig& handoff_config = {});

class JouleHandoffAccumulator {
 public:
  explicit JouleHandoffAccumulator(JouleHandoffConfig config = {});
  JouleHandoffDiagnostics sample(const AxisymmetricGrid& grid, const StreamerState& state,
                                 const StreamerConfig& config,
                                 const ElectronTransportCurrentSource& current_source,
                                 double Gb_S, double Rb_ohm, bool Gb_valid,
                                 bool bridge_flag);
  void reset();

 private:
  JouleHandoffConfig handoff_config_;
  bool has_previous_power_{};
  bool has_previous_gb_{};
  double previous_time_s_{};
  double previous_PJ_gas_W_{};
  double previous_PJ_channel_W_{};
  bool previous_channel_valid_{};
  double previous_Gb_S_{};
  double QJ_gas_J_{};
  double QJ_channel_J_{};
};

}  // namespace streamer_rf::streamer
