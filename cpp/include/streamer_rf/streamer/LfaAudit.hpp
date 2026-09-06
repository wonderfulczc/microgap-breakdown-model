#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace streamer_rf::streamer {

struct ElectronRelaxationMetadata {
  std::string source_id;
  std::string source_reference;
  std::string gas_composition;
  double pressure_Pa{};
  double temperature_K{};
  double neutral_density_m3{};
};

struct ElectronRelaxationLookup {
  std::string status{"NO_RELAXATION_TABLE"};
  double reduced_field_Td{};
  bool has_tau_epsilon{};
  bool has_lambda_epsilon{};
  double tau_epsilon_s{std::numeric_limits<double>::quiet_NaN()};
  double lambda_epsilon_m{std::numeric_limits<double>::quiet_NaN()};
};

class ElectronRelaxationTable {
 public:
  ElectronRelaxationTable(ElectronRelaxationMetadata metadata, std::vector<double> reduced_field_Td,
                          std::optional<std::vector<double>> tau_epsilon_s,
                          std::optional<std::vector<double>> lambda_epsilon_m);

  ElectronRelaxationLookup lookup(double reduced_field_Td) const;
  bool has_tau_epsilon() const { return !tau_epsilon_s_.empty(); }
  bool has_lambda_epsilon() const { return !lambda_epsilon_m_.empty(); }
  double min_reduced_field_Td() const { return reduced_field_Td_.front(); }
  double max_reduced_field_Td() const { return reduced_field_Td_.back(); }
  const ElectronRelaxationMetadata& metadata() const { return metadata_; }

 private:
  ElectronRelaxationMetadata metadata_;
  std::vector<double> reduced_field_Td_;
  std::vector<double> tau_epsilon_s_;
  std::vector<double> lambda_epsilon_m_;
};

enum class LfaAuditRegionKind { AllGas, ActiveElectron, HighField, CustomMask };

struct LfaAuditRegionDefinition {
  std::string name{"ALL_GAS"};
  LfaAuditRegionKind kind{LfaAuditRegionKind::AllGas};
  double active_electron_relative_floor{1e-12};
  double high_field_fraction{0.5};
};

std::vector<LfaAuditRegionDefinition> default_lfa_audit_regions();

struct LfaAuditSummary {
  std::string region_name{"ALL_GAS"};
  double eovern_max_Td{};
  double le_min_valid_m{std::numeric_limits<double>::quiet_NaN()};
  double le_p05_m{std::numeric_limits<double>::quiet_NaN()};
  double le_median_m{std::numeric_limits<double>::quiet_NaN()};
  double le_valid_fraction{};
  double tauE_min_valid_s{std::numeric_limits<double>::quiet_NaN()};
  double tauE_p05_s{std::numeric_limits<double>::quiet_NaN()};
  double tauE_median_s{std::numeric_limits<double>::quiet_NaN()};
  double tauE_valid_fraction{};
  int gas_cells{};
  int selected_gas_cells{};
  int le_valid_cells{};
  int tauE_valid_cells{};
  int insufficient_stencil_cells{};
  int near_zero_gradE_cells{};
  int near_zero_dEdt_cells{};
  int nonfinite_input_cells{};
  double chiL_p95{std::numeric_limits<double>::quiet_NaN()};
  double chiL_median{std::numeric_limits<double>::quiet_NaN()};
  double chiT_p95{std::numeric_limits<double>::quiet_NaN()};
  double chiT_median{std::numeric_limits<double>::quiet_NaN()};
  double relaxation_coverage_fraction{};
  double relaxation_table_min_Td{std::numeric_limits<double>::quiet_NaN()};
  double relaxation_table_max_Td{std::numeric_limits<double>::quiet_NaN()};
  double fraction_outside_relaxation_table{};
  int relaxation_covered_cells{};
  int outside_relaxation_table_cells{};
  int chiL_valid_cells{};
  int chiT_valid_cells{};
  std::string temporal_status{"INVALID_INITIAL_SAMPLE"};
  std::string lfa_relaxation_data_status{"NOT_AVAILABLE"};
  std::string lfa_applicability{"UNRESOLVED_RELAXATION_DATA"};
};

LfaAuditSummary compute_lfa_audit(const AxisymmetricGrid& grid, const StreamerState& state,
                                  const StreamerConfig& config, double accepted_dt_s,
                                  const ScalarField2D* previous_emag = nullptr,
                                  const std::vector<unsigned char>* optional_cell_mask = nullptr,
                                  const ElectronRelaxationTable* relaxation_table = nullptr,
                                  LfaAuditRegionDefinition region = {});

class LfaAuditHistory {
 public:
  LfaAuditSummary sample(const AxisymmetricGrid& grid, const StreamerState& state,
                         const StreamerConfig& config, double accepted_dt_s,
                         const std::vector<unsigned char>* optional_cell_mask = nullptr,
                         const ElectronRelaxationTable* relaxation_table = nullptr,
                         LfaAuditRegionDefinition region = {});
  void reset();
  bool has_previous() const { return has_previous_; }

 private:
  std::vector<double> previous_emag_values_;
  bool has_previous_{};
};

}  // namespace streamer_rf::streamer
