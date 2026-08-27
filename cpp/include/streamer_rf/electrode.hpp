#pragma once
#include "streamer_rf/types.hpp"
#include <string>

namespace streamer_rf {

enum class ElectrodeCellType { Gas, HighVoltageElectrode, GroundElectrode };

struct ElectrodeGeometryMetadata {
  std::string geometry_id;
  double gap_m{};
  double tip_radius_m{};
  std::string representation{"structured cell-center mask"};
};

class AxisymmetricNeedlePlaneGeometry {
 public:
  AxisymmetricNeedlePlaneGeometry(std::string geometry_id, double ground_z_m,
                                  double tip_z_m, double tip_radius_m,
                                  double shank_radius_m, double ground_thickness_m);
  ElectrodeCellType classify(double r_m, double z_m) const;
  ElectrodeCellType classify(const AxisymmetricGrid& grid, int i, int j) const {
    return classify(grid.r(i), grid.z(j));
  }
  const ElectrodeGeometryMetadata& metadata() const { return metadata_; }
  double ground_z_m() const { return ground_z_m_; }
  double tip_z_m() const { return tip_z_m_; }
  double tip_radius_m() const { return tip_radius_m_; }
  double shank_radius_m() const { return shank_radius_m_; }
  double ground_thickness_m() const { return ground_thickness_m_; }
  double seed_z_from_tip_offset(double z_offset_from_tip_m) const { return tip_z_m_ + z_offset_from_tip_m; }

 private:
  ElectrodeGeometryMetadata metadata_;
  double ground_z_m_{};
  double tip_z_m_{};
  double tip_radius_m_{};
  double shank_radius_m_{};
  double ground_thickness_m_{};
};

struct ElectrodePoissonDiagnostics {
  double time_s{};
  double applied_voltage_V{};
  double emax_V_m{};
  double phi_hv_residual_V{};
  double phi_ground_residual_V{};
  std::string geometry_id;
  double gap_m{};
  double tip_radius_m{};
  int ksp_iterations{};
};

ElectrodePoissonDiagnostics evaluate_electrode_diagnostics(
    const ScalarField2D& phi, const ScalarField2D& er, const ScalarField2D& ez,
    const AxisymmetricNeedlePlaneGeometry& geometry, double time_s,
    double applied_voltage_V, int ksp_iterations);

}  // namespace streamer_rf
