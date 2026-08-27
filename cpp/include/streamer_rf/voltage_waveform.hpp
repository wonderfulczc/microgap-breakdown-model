#pragma once
#include <filesystem>
#include <string>
#include <vector>

namespace streamer_rf {

enum class WaveformOutOfRangePolicy { Error, HoldEndpoint };

class VoltageWaveform {
 public:
  virtual ~VoltageWaveform() = default;
  virtual double value(double time_s) const = 0;
};

class ConstantVoltage final : public VoltageWaveform {
 public:
  explicit ConstantVoltage(double voltage_V) : voltage_V_(voltage_V) {}
  double value(double) const override { return voltage_V_; }

 private:
  double voltage_V_{};
};

class SampledVoltage final : public VoltageWaveform {
 public:
  SampledVoltage(std::vector<double> time_s, std::vector<double> voltage_V,
                 WaveformOutOfRangePolicy policy = WaveformOutOfRangePolicy::Error);
  static SampledVoltage from_csv(const std::filesystem::path& path,
                                 WaveformOutOfRangePolicy policy = WaveformOutOfRangePolicy::Error);
  double value(double time_s) const override;
  const std::vector<double>& times_s() const { return time_s_; }
  const std::vector<double>& voltages_V() const { return voltage_V_; }

 private:
  std::vector<double> time_s_;
  std::vector<double> voltage_V_;
  WaveformOutOfRangePolicy policy_{WaveformOutOfRangePolicy::Error};
};

}  // namespace streamer_rf
