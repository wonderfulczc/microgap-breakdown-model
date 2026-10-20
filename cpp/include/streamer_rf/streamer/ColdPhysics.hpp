#pragma once
#include "streamer_rf/geometry2d.hpp"
#include <limits>
#include <string>
#include <vector>

namespace streamer_rf::streamer {
struct IonFieldTable {
  std::string species_id,gas_composition,source_doi,database_id,file_sha256,evidence_class;
  std::string mobility_unit{"m2/(V*s)"},diffusion_unit{"m2/s"},field_unit{"Td"};
  std::string interpolation{"LINEAR_EN"},out_of_domain{"REJECT"};
  std::string aggregate_mapping{"EFFECTIVE_SPECIES_APPROXIMATION"};
  double reference_density_m3{NAN}; // Required only for reduced mobility K0=mu*N/N0.
  std::vector<double> EN_Td,mobility,D_parallel,D_perpendicular;
  void validate() const;
};
struct IonTransportCoefficients {double mobility,D_parallel,D_perpendicular;};
enum class ColdPhysicsMode { LegacyElectronOnly, CoupledColdV3 };
enum class PhotonMetalBoundary { Unavailable, OpaqueAbsorbingLarsen };
struct IonTransportModel {
  std::string model_id, source, species_id, approximation, qualification;
  double mobility_m2_Vs{NAN}, diffusion_m2_s{NAN};
  double min_EN_Td{NAN}, max_EN_Td{NAN}, min_temperature_K{NAN}, max_temperature_K{NAN};
  double min_pressure_Pa{NAN}, max_pressure_Pa{NAN};
  IonFieldTable table;
  void validate(double EN_Td, double temperature_K, double pressure_Pa) const;
};
class PositiveIonTransportProvider {
 public:
  explicit PositiveIonTransportProvider(const IonTransportModel& model):model_(model){}
  IonTransportCoefficients evaluate(double EN_Td,double temperature_K,double pressure_Pa,double neutral_density_m3) const;
 private:
  const IonTransportModel& model_;
};
struct SurfaceEmissionModel {
  bool enabled{false};
  std::string model_id, source, qualification;
  double gamma_i{NAN};
  bool source_contract_required{false};
  std::string incident_species_id,gas_composition,cathode_material,surface_state;
  std::string incident_energy_assumption,source_sha256,evidence_class;
  double incident_energy_eV{NAN},min_energy_eV{NAN},max_energy_eV{NAN};
  std::string uncertainty,energy_status{"UNAVAILABLE"};
  void validate_source() const;
  double emitted_flux(double incident_m2s, ElectrodeCellType conductor, double voltage_V) const;
};
struct ColdPhysicsConfig {
  ColdPhysicsMode mode{ColdPhysicsMode::LegacyElectronOnly};
  bool positive_ion_transport{false}, local_charge_port_v2{false};
  IonTransportModel ions;
  SurfaceEmissionModel surface;
  PhotonMetalBoundary photon_boundary{PhotonMetalBoundary::Unavailable};
  std::string excitation_source;
  void validate(bool geometry, bool photoionization) const;
};
struct IonFaceFlux {
  int i,j;
  FaceDirection2D direction;
  GeometryFace2D face;
  double gamma_p_out_m2s, gamma_secondary_in_m2s;
};
struct IonTransportResult {
  ScalarField2D rhs, secondary_rhs;
  std::vector<IonFaceFlux> faces;
  double positivity_dt_s{std::numeric_limits<double>::infinity()};
  double hv_collected_s1{}, return_collected_s1{}, outer_lost_s1{};
  double hv_emitted_s1{}, return_emitted_s1{};
  explicit IonTransportResult(const AxisymmetricGrid& g):rhs(g),secondary_rhs(g){}
};
// First-order upwind drift and centered Fick diffusion, with absorbing half-cell walls.
IonTransportResult positive_ion_transport(const Geometry2DView&, const ScalarField2D& np,
    const ScalarField2D& er, const ScalarField2D& ez, const IonTransportModel&,
    const SurfaceEmissionModel&, double voltage, double neutral_density,
    double temperature, double pressure, bool capture_faces=true);
}
