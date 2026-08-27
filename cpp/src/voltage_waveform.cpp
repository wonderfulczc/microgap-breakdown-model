#include "streamer_rf/voltage_waveform.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace streamer_rf {

SampledVoltage::SampledVoltage(std::vector<double> time_s, std::vector<double> voltage_V,
                               WaveformOutOfRangePolicy policy)
    : time_s_(std::move(time_s)), voltage_V_(std::move(voltage_V)), policy_(policy) {
  if (time_s_.size() != voltage_V_.size() || time_s_.size() < 2) {
    throw std::invalid_argument("sampled voltage requires matching time and voltage arrays with at least two samples");
  }
  for (std::size_t i = 1; i < time_s_.size(); ++i) {
    if (!(time_s_[i] > time_s_[i - 1])) {
      throw std::invalid_argument("sampled voltage time_s must be strictly increasing");
    }
  }
}

SampledVoltage SampledVoltage::from_csv(const std::filesystem::path& path, WaveformOutOfRangePolicy policy) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("failed to open voltage waveform CSV");
  std::string line;
  if (!std::getline(f, line)) throw std::runtime_error("empty voltage waveform CSV");
  std::vector<double> t, v;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string a, b;
    if (!std::getline(ss, a, ',') || !std::getline(ss, b, ',')) {
      throw std::runtime_error("invalid voltage waveform CSV row");
    }
    t.push_back(std::stod(a));
    v.push_back(std::stod(b));
  }
  return SampledVoltage(std::move(t), std::move(v), policy);
}

double SampledVoltage::value(double time_s) const {
  if (time_s < time_s_.front()) {
    if (policy_ == WaveformOutOfRangePolicy::HoldEndpoint) return voltage_V_.front();
    throw std::out_of_range("sampled voltage time is before first sample");
  }
  if (time_s > time_s_.back()) {
    if (policy_ == WaveformOutOfRangePolicy::HoldEndpoint) return voltage_V_.back();
    throw std::out_of_range("sampled voltage time is after last sample");
  }
  if (time_s == time_s_.front()) return voltage_V_.front();
  if (time_s == time_s_.back()) return voltage_V_.back();
  const auto hi = std::upper_bound(time_s_.begin(), time_s_.end(), time_s);
  const std::size_t k = static_cast<std::size_t>(hi - time_s_.begin());
  const double t0 = time_s_[k - 1], t1 = time_s_[k];
  const double v0 = voltage_V_[k - 1], v1 = voltage_V_[k];
  const double w = (time_s - t0) / (t1 - t0);
  return (1.0 - w) * v0 + w * v1;
}

}  // namespace streamer_rf
