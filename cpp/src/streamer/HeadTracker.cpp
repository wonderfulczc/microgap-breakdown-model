#include "streamer_rf/streamer/HeadTracker.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace streamer_rf::streamer {
namespace {
bool is_gas_cell(const AxisymmetricGrid& g, const StreamerConfig& c, int i, int j) {
  return !c.electrode_geometry || c.electrode_geometry->classify(g, i, j) == ElectrodeCellType::Gas;
}

std::size_t index(const AxisymmetricGrid& g, int i, int j) {
  return static_cast<std::size_t>(j) * g.nr() + i;
}

double cell_rho(const StreamerState& s, int i, int j) {
  return s.rho(i, j);
}
}  // namespace

StreamerHeadDiagnostics compute_streamer_head_diagnostics(
    const AxisymmetricGrid& g, const StreamerState& s, const StreamerConfig& c,
    const StreamerHeadSegmentationConfig& seg) {
  StreamerHeadDiagnostics d;
  d.time_s = s.time;
  d.segmentation_parameter = seg.relative_rho_threshold;
  d.position_method = seg.position_method;
  if (!(seg.relative_rho_threshold > 0.0 && seg.relative_rho_threshold <= 1.0) ||
      !std::isfinite(seg.relative_rho_threshold) ||
      (seg.requested_polarity != 0 && seg.requested_polarity != 1 && seg.requested_polarity != -1)) {
    d.head_status = "INVALID_SEGMENTATION_PARAMETER";
    return d;
  }

  int gas_cells = 0, peak_i = -1, peak_j = -1;
  double peak_abs = 0.0, peak_rho = 0.0;
  const int requested = seg.requested_polarity;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (!is_gas_cell(g, c, i, j)) continue;
      ++gas_cells;
      const double rho = cell_rho(s, i, j);
      const double e = s.emag(i, j);
      if (!std::isfinite(rho) || !std::isfinite(e)) {
        d.head_status = "NONFINITE_INPUT";
        return d;
      }
      const double a = std::abs(rho);
      if (requested == 0) {
        if (a > peak_abs) {
          peak_abs = a;
          peak_rho = rho;
          peak_i = i;
          peak_j = j;
        }
      } else if (requested * rho > peak_abs) {
        peak_abs = requested * rho;
        peak_rho = rho;
        peak_i = i;
        peak_j = j;
      }
    }
  }
  if (gas_cells == 0) {
    d.head_status = "NO_HEAD";
    return d;
  }
  if (!(peak_abs > 0.0) || peak_i < 0) {
    d.head_status = "INSUFFICIENT_CHARGE";
    return d;
  }
  const int polarity = peak_rho >= 0.0 ? 1 : -1;
  const double threshold = seg.relative_rho_threshold * peak_abs;
  auto selected = [&](int i, int j) {
    if (!is_gas_cell(g, c, i, j)) return false;
    const double rho = cell_rho(s, i, j);
    return std::isfinite(rho) && polarity * rho >= threshold;
  };
  if (!selected(peak_i, peak_j)) {
    d.head_status = "INSUFFICIENT_CHARGE";
    return d;
  }

  std::vector<unsigned char> visited(g.size(), 0);
  std::vector<std::pair<int, int>> stack{{peak_i, peak_j}};
  visited[index(g, peak_i, peak_j)] = 1;
  double q = 0.0, wsum = 0.0, zsum = 0.0, rsum = 0.0, r2sum = 0.0;
  double max_e = 0.0;
  while (!stack.empty()) {
    auto [i, j] = stack.back();
    stack.pop_back();
    ++d.head_cell_count;
    const double rho = cell_rho(s, i, j);
    const double vol = g.cell_volume(i);
    const double w = std::abs(rho) * vol;
    q += rho * vol;
    wsum += w;
    zsum += w * g.z(j);
    rsum += w * g.r(i);
    r2sum += w * g.r(i) * g.r(i);
    max_e = std::max(max_e, s.emag(i, j));
    const int di[4] = {-1, 1, 0, 0};
    const int dj[4] = {0, 0, -1, 1};
    for (int n = 0; n < 4; ++n) {
      const int ni = i + di[n], nj = j + dj[n];
      if (ni < 0 || ni >= g.nr() || nj < 0 || nj >= g.nz()) continue;
      const std::size_t k = index(g, ni, nj);
      if (visited[k] || !selected(ni, nj)) continue;
      visited[k] = 1;
      stack.push_back({ni, nj});
    }
  }
  if (!(wsum > 0.0) || d.head_cell_count == 0) {
    d.head_status = "INSUFFICIENT_CHARGE";
    return d;
  }
  d.head_valid = true;
  d.head_status = "VALID";
  d.head_polarity = polarity;
  d.head_charge_C = q;
  d.head_z_m = zsum / wsum;
  d.head_r_mean_m = rsum / wsum;
  d.head_r_rms_m = std::sqrt(r2sum / wsum);
  d.head_peak_rho_C_m3 = peak_rho;
  d.head_peak_E_V_m = max_e;
  d.head_peak_EoverN_Td = max_e / c.neutral_density * 1e21;
  d.segmentation_fraction = static_cast<double>(d.head_cell_count) / gas_cells;
  return d;
}

StreamerHeadTracker::StreamerHeadTracker(StreamerHeadSegmentationConfig segmentation)
    : segmentation_(std::move(segmentation)) {}

void StreamerHeadTracker::reset() {
  previous_count_ = 0;
}

StreamerHeadDiagnostics StreamerHeadTracker::sample(const AxisymmetricGrid& g, const StreamerState& s,
                                                    const StreamerConfig& c) {
  auto d = compute_streamer_head_diagnostics(g, s, c, segmentation_);
  if (!d.head_valid) {
    reset();
    return d;
  }
  if (previous_count_ >= 1) {
    const auto& p = previous_[previous_count_ == 1 ? 0 : 1];
    const double dt = d.time_s - p.time;
    if (dt > 0.0 && std::isfinite(dt)) {
      d.head_velocity_z_m_s = (d.head_z_m - p.z) / dt;
      d.velocity_valid = true;
      d.velocity_status = "VALID";
    } else {
      d.velocity_status = "NONFINITE_INPUT";
    }
  }
  if (previous_count_ >= 2 && d.velocity_valid) {
    const auto& p0 = previous_[0];
    const auto& p1 = previous_[1];
    const double dt_prev = p1.time - p0.time;
    const double dt_now = d.time_s - p1.time;
    if (dt_prev > 0.0 && dt_now > 0.0 && std::isfinite(dt_prev) && std::isfinite(dt_now)) {
      const double v_prev = (p1.z - p0.z) / dt_prev;
      const double v_now = (d.head_z_m - p1.z) / dt_now;
      d.head_acceleration_z_m_s2 = 2.0 * (v_now - v_prev) / (dt_now + dt_prev);
      d.acceleration_valid = true;
      d.acceleration_status = "VALID";
    } else {
      d.acceleration_status = "NONFINITE_INPUT";
    }
  }
  const Sample current{d.time_s, d.head_z_m, d.head_charge_C};
  if (previous_count_ == 0) {
    previous_[0] = current;
    previous_count_ = 1;
  } else if (previous_count_ == 1) {
    previous_[1] = current;
    previous_count_ = 2;
  } else {
    previous_[0] = previous_[1];
    previous_[1] = current;
  }
  return d;
}

}  // namespace streamer_rf::streamer
