#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <array>
#include <string>

namespace streamer_rf::streamer {

struct StreamerHeadSegmentationConfig {
  double relative_rho_threshold{0.2};
  int requested_polarity{};
  std::string position_method{"same_polarity_charge_centroid"};
};

struct StreamerHeadDiagnostics {
  double time_s{};
  bool head_valid{};
  std::string head_status{"NO_HEAD"};
  int head_polarity{};
  int head_cell_count{};
  double head_charge_C{std::numeric_limits<double>::quiet_NaN()};
  double head_z_m{std::numeric_limits<double>::quiet_NaN()};
  double head_r_mean_m{std::numeric_limits<double>::quiet_NaN()};
  double head_r_rms_m{std::numeric_limits<double>::quiet_NaN()};
  double head_velocity_z_m_s{std::numeric_limits<double>::quiet_NaN()};
  double head_acceleration_z_m_s2{std::numeric_limits<double>::quiet_NaN()};
  bool velocity_valid{};
  bool acceleration_valid{};
  std::string velocity_status{"INSUFFICIENT_HISTORY_FOR_VELOCITY"};
  std::string acceleration_status{"INSUFFICIENT_HISTORY_FOR_ACCELERATION"};
  double head_peak_rho_C_m3{std::numeric_limits<double>::quiet_NaN()};
  double head_peak_E_V_m{std::numeric_limits<double>::quiet_NaN()};
  double head_peak_EoverN_Td{std::numeric_limits<double>::quiet_NaN()};
  double segmentation_fraction{};
  double segmentation_parameter{};
  std::string position_method;
};

StreamerHeadDiagnostics compute_streamer_head_diagnostics(
    const AxisymmetricGrid& grid, const StreamerState& state,
    const StreamerConfig& config, const StreamerHeadSegmentationConfig& segmentation = {});

class StreamerHeadTracker {
 public:
  explicit StreamerHeadTracker(StreamerHeadSegmentationConfig segmentation = {});
  StreamerHeadDiagnostics sample(const AxisymmetricGrid& grid, const StreamerState& state,
                                 const StreamerConfig& config);
  void reset();

 private:
  struct Sample { double time{}, z{}, charge{}; };
  StreamerHeadSegmentationConfig segmentation_;
  std::array<Sample, 2> previous_{};
  int previous_count_{};
};

}  // namespace streamer_rf::streamer
