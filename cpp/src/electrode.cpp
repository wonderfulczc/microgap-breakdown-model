#include "streamer_rf/electrode.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace streamer_rf {

AxisymmetricNeedlePlaneGeometry::AxisymmetricNeedlePlaneGeometry(
    std::string geometry_id, double ground_z_m, double tip_z_m, double tip_radius_m,
    double shank_radius_m, double ground_thickness_m)
    : ground_z_m_(ground_z_m),
      tip_z_m_(tip_z_m),
      tip_radius_m_(tip_radius_m),
      shank_radius_m_(shank_radius_m),
      ground_thickness_m_(ground_thickness_m) {
  if (geometry_id.empty() || tip_z_m <= ground_z_m || tip_radius_m <= 0 ||
      shank_radius_m <= 0 || ground_thickness_m <= 0) {
    throw std::invalid_argument("invalid axisymmetric needle-plane geometry");
  }
  metadata_.geometry_id = std::move(geometry_id);
  metadata_.gap_m = tip_z_m_ - ground_z_m_;
  metadata_.tip_radius_m = tip_radius_m_;
}

ElectrodeCellType AxisymmetricNeedlePlaneGeometry::classify(double r_m, double z_m) const {
  if (z_m <= ground_z_m_ + ground_thickness_m_) {
    return ElectrodeCellType::GroundElectrode;
  }
  const double dr_tip = r_m;
  const double dz_tip = z_m - tip_z_m_;
  if (std::hypot(dr_tip, dz_tip) <= tip_radius_m_) {
    return ElectrodeCellType::HighVoltageElectrode;
  }
  if (r_m <= shank_radius_m_ && z_m >= tip_z_m_) {
    return ElectrodeCellType::HighVoltageElectrode;
  }
  return ElectrodeCellType::Gas;
}

ElectrodePoissonDiagnostics evaluate_electrode_diagnostics(
    const ScalarField2D& phi, const ScalarField2D& er, const ScalarField2D& ez,
    const AxisymmetricNeedlePlaneGeometry& geometry, double time_s,
    double applied_voltage_V, int ksp_iterations) {
  ElectrodePoissonDiagnostics d;
  d.time_s = time_s;
  d.applied_voltage_V = applied_voltage_V;
  d.geometry_id = geometry.metadata().geometry_id;
  d.gap_m = geometry.metadata().gap_m;
  d.tip_radius_m = geometry.metadata().tip_radius_m;
  d.ksp_iterations = ksp_iterations;
  const auto& g = phi.grid();
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      d.emax_V_m = std::max(d.emax_V_m, std::hypot(er(i, j), ez(i, j)));
      const auto c = geometry.classify(g, i, j);
      if (c == ElectrodeCellType::HighVoltageElectrode) {
        d.phi_hv_residual_V = std::max(d.phi_hv_residual_V, std::abs(phi(i, j) - applied_voltage_V));
      } else if (c == ElectrodeCellType::GroundElectrode) {
        d.phi_ground_residual_V = std::max(d.phi_ground_residual_V, std::abs(phi(i, j)));
      }
    }
  }
  return d;
}

}  // namespace streamer_rf
