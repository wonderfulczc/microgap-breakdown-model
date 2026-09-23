#pragma once

#include "streamer_rf/streamer/HeadTracker.hpp"
#include <limits>
#include <string>
#include <vector>

namespace streamer_rf::streamer {

enum class EventTopology {
  Propagation,
  HeadInteraction,
  StreamerCollision,
  ElectrodeAttachment,
  Bridging,
  Unresolved,
};

const char* event_topology_name(EventTopology topology);

struct UltrafastEventConfig {
  double relative_rho_threshold{0.2};
  int requested_polarity{};
  int minimum_component_cells{2};
  double collision_distance_cells{3.0};
  double interaction_distance_cells{8.0};
  double attachment_distance_cells{1.5};
  int confirmation_samples{2};
};

struct UltrafastEventSample {
  int step{};
  double time_s{};
  double dt_s{};
  EventTopology topology{EventTopology::Unresolved};
  std::string topology_confidence{"NONE"};
  std::string topology_evidence{"NO_EVENT_COMPONENT"};
  std::string topology_status{"UNRESOLVED"};
  int head_count{};
  std::string secondary_head_status{"SECONDARY_HEAD_DIAGNOSTIC_EXPERIMENTAL"};

  std::string roi_definition{"NOT_AVAILABLE"};
  std::string roi_status{"NOT_RESOLVED"};
  int roi_cell_count{};
  double roi_volume_m3{};
  double roi_centroid_r_m{std::numeric_limits<double>::quiet_NaN()};
  double roi_centroid_z_m{std::numeric_limits<double>::quiet_NaN()};
  double roi_r_min_m{std::numeric_limits<double>::quiet_NaN()};
  double roi_r_max_m{std::numeric_limits<double>::quiet_NaN()};
  double roi_z_min_m{std::numeric_limits<double>::quiet_NaN()};
  double roi_z_max_m{std::numeric_limits<double>::quiet_NaN()};
  double electrode_distance_m{std::numeric_limits<double>::quiet_NaN()};

  double Ek_V_m{};
  std::string gas{"AIR_MORROW_LOWKE_BASELINE"};
  double pressure_Pa{};
  double temperature_K{};
  std::string Ek_semantics{"REFERENCE_NORMALIZATION_ONLY"};
  std::string E0_source{"STAGE_C_2D_ELECTROSTATIC_REFERENCE"};
  double E0_max_V_m{std::numeric_limits<double>::quiet_NaN()};
  double E0_mean_V_m{std::numeric_limits<double>::quiet_NaN()};
  double E0_p95_V_m{std::numeric_limits<double>::quiet_NaN()};
  double eta_0_max{std::numeric_limits<double>::quiet_NaN()};
  double eta_0_mean{std::numeric_limits<double>::quiet_NaN()};
  double eta_0_p95{std::numeric_limits<double>::quiet_NaN()};

  double E_peak_V_m{std::numeric_limits<double>::quiet_NaN()};
  double E_peak_time_s{std::numeric_limits<double>::quiet_NaN()};
  double E_peak_r_m{std::numeric_limits<double>::quiet_NaN()};
  double E_peak_z_m{std::numeric_limits<double>::quiet_NaN()};
  double ne_at_E_peak_m3{std::numeric_limits<double>::quiet_NaN()};
  double eta_peak{std::numeric_limits<double>::quiet_NaN()};
  double ne_peak_m3{std::numeric_limits<double>::quiet_NaN()};
  double ne_peak_r_m{std::numeric_limits<double>::quiet_NaN()};
  double ne_peak_z_m{std::numeric_limits<double>::quiet_NaN()};

  double E_at_ne_peak_V_m{std::numeric_limits<double>::quiet_NaN()};
  double mu_e_at_ne_peak_m2_V_s{std::numeric_limits<double>::quiet_NaN()};
  double nu_i_at_ne_peak_s_1{std::numeric_limits<double>::quiet_NaN()};

  double mu_e_at_E_peak_m2_V_s{std::numeric_limits<double>::quiet_NaN()};
  double nu_i_at_E_peak_s_1{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_at_E_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_at_ne_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_min_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_median_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_p95_s{std::numeric_limits<double>::quiet_NaN()};

  double sigma_e_at_E_peak_S_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_at_ne_peak_S_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_peak_S_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_peak_r_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_peak_z_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_median_S_m{std::numeric_limits<double>::quiet_NaN()};
  double sigma_e_p95_S_m{std::numeric_limits<double>::quiet_NaN()};
  std::string sigma_e_role{"ELECTRON_CONDUCTIVITY_DIAGNOSTIC"};

  double tau_M_at_E_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_M_at_ne_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_M_min_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_M_median_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_M_p95_s{std::numeric_limits<double>::quiet_NaN()};

  double Pi_RF_at_E_peak{std::numeric_limits<double>::quiet_NaN()};
  double Pi_RF_at_ne_peak{std::numeric_limits<double>::quiet_NaN()};
  double Pi_RF_min{std::numeric_limits<double>::quiet_NaN()};
  double Pi_RF_median{std::numeric_limits<double>::quiet_NaN()};
  double Pi_RF_p95{std::numeric_limits<double>::quiet_NaN()};
  std::string Pi_RF_role{"DIAGNOSTIC_CANDIDATE"};

  double E_at_sigma_peak_V_m{std::numeric_limits<double>::quiet_NaN()};
  double ne_at_sigma_peak_m3{std::numeric_limits<double>::quiet_NaN()};
  double mu_e_at_sigma_peak_m2_V_s{std::numeric_limits<double>::quiet_NaN()};
  double nu_i_at_sigma_peak_s_1{std::numeric_limits<double>::quiet_NaN()};
  double tau_i_at_sigma_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double tau_M_at_sigma_peak_s{std::numeric_limits<double>::quiet_NaN()};
  double Pi_RF_at_sigma_peak{std::numeric_limits<double>::quiet_NaN()};

  double E_roi_median_V_m{std::numeric_limits<double>::quiet_NaN()};
  double ne_roi_median_m3{std::numeric_limits<double>::quiet_NaN()};
  double mu_e_roi_median_m2_V_s{std::numeric_limits<double>::quiet_NaN()};
  double nu_i_roi_median_s_1{std::numeric_limits<double>::quiet_NaN()};

  double K_ion_z_A_m_s{std::numeric_limits<double>::quiet_NaN()};
  double K_ion_abs_A_m_s{std::numeric_limits<double>::quiet_NaN()};
  std::string signed_proxy_status{"DEFINED_FROM_FROZEN_ELECTRON_DRIFT_CONVENTION"};
  std::string mechanism_proxy_role{"MECHANISM_DIAGNOSTIC_PROXY_NOT_CLOSURE_RELATION"};

  double N_tau_i_at_E_peak{std::numeric_limits<double>::quiet_NaN()};
  double N_tau_M_at_E_peak{std::numeric_limits<double>::quiet_NaN()};
  std::string tau_i_resolution_status{"NOT_RESOLVED"};
  std::string tau_M_resolution_status{"NOT_RESOLVED"};
  std::string temporal_semantics{"SOLVER_ACCEPTED_STEP_NO_INTERPOLATION"};

  double head_position_m{std::numeric_limits<double>::quiet_NaN()};
  double head_velocity_m_s{std::numeric_limits<double>::quiet_NaN()};
  bool bridge_status{};
  double current_moment_r_A_m{};
  double current_moment_z_A_m{};
  std::string current_definition{"FROZEN_FINITE_VOLUME_ELECTRON_TRANSPORT_CURRENT"};
};

class UltrafastEventDiagnostics {
 public:
  UltrafastEventDiagnostics(const AxisymmetricGrid& grid, const StreamerConfig& config,
                            const ScalarField2D& electrostatic_reference_E,
                            UltrafastEventConfig event_config = {});

  UltrafastEventSample sample(int step, const StreamerState& state,
                              const StreamerHeadDiagnostics& primary_head,
                              const ElectronTransportCurrentSource& current_source,
                              bool bridge_status, double accepted_dt_s);

  double Ek_V_m() const { return Ek_V_m_; }
  static std::string kinetic_resolution_status(double timescale_s, double local_dt_s);

 private:
  const AxisymmetricGrid& grid_;
  const StreamerConfig& config_;
  const ScalarField2D& E0_;
  UltrafastEventConfig event_config_;
  double Ek_V_m_{};
  EventTopology pending_topology_{EventTopology::Unresolved};
  int pending_count_{};
};

}  // namespace streamer_rf::streamer
