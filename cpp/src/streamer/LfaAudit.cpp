#include "streamer_rf/streamer/LfaAudit.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace streamer_rf::streamer {
namespace {

bool gas_cell(const AxisymmetricGrid& g, const StreamerConfig& cfg, int i, int j) {
  return !cfg.electrode_geometry || cfg.electrode_geometry->classify(g, i, j) == ElectrodeCellType::Gas;
}

bool selected_cell(const AxisymmetricGrid& g, int i, int j, const std::vector<unsigned char>* mask) {
  if (!mask) return true;
  if (mask->size() != g.size()) throw std::invalid_argument("LFA audit mask size does not match grid");
  return (*mask)[static_cast<std::size_t>(j) * g.nr() + i] != 0;
}

double interpolate_linear(const std::vector<double>& x, const std::vector<double>& y, double v) {
  auto it = std::upper_bound(x.begin(), x.end(), v);
  if (it == x.begin()) return y.front();
  if (it == x.end()) return y.back();
  const auto k = static_cast<std::size_t>(it - x.begin());
  const double t = (v - x[k - 1]) / (x[k] - x[k - 1]);
  return (1.0 - t) * y[k - 1] + t * y[k];
}

double percentile(std::vector<double> values, double q) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double x = q * static_cast<double>(values.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(x));
  const auto hi = static_cast<std::size_t>(std::ceil(x));
  const double f = x - static_cast<double>(lo);
  return (1.0 - f) * values[lo] + f * values[hi];
}

struct Difference {
  double value{};
  bool valid{};
};

Difference directional_difference(const AxisymmetricGrid& g, const StreamerState& s,
                                  const StreamerConfig& cfg, int i, int j, int di, int dj,
                                  double h) {
  const int im = i - di, jm = j - dj;
  const int ip = i + di, jp = j + dj;
  const bool has_m = im >= 0 && im < g.nr() && jm >= 0 && jm < g.nz() && gas_cell(g, cfg, im, jm) &&
                     std::isfinite(s.emag(im, jm));
  const bool has_p = ip >= 0 && ip < g.nr() && jp >= 0 && jp < g.nz() && gas_cell(g, cfg, ip, jp) &&
                     std::isfinite(s.emag(ip, jp));
  if (has_m && has_p) return {(s.emag(ip, jp) - s.emag(im, jm)) / (2.0 * h), true};
  if (has_p) return {(s.emag(ip, jp) - s.emag(i, j)) / h, true};
  if (has_m) return {(s.emag(i, j) - s.emag(im, jm)) / h, true};
  return {};
}

double max_selected_emag(const AxisymmetricGrid& g, const StreamerState& s, const StreamerConfig& cfg,
                         const std::vector<unsigned char>* mask, const LfaAuditRegionDefinition& region,
                         double active_threshold, double high_field_threshold) {
  double out = 0.0;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (!gas_cell(g, cfg, i, j) || !selected_cell(g, i, j, mask)) continue;
      if (region.kind == LfaAuditRegionKind::ActiveElectron && !(s.ne(i, j) > active_threshold)) continue;
      if (region.kind == LfaAuditRegionKind::HighField && !(std::abs(s.emag(i, j)) >= high_field_threshold)) continue;
      if (std::isfinite(s.emag(i, j))) out = std::max(out, std::abs(s.emag(i, j)));
    }
  }
  return out;
}

double max_gas_ne(const AxisymmetricGrid& g, const StreamerState& s, const StreamerConfig& cfg) {
  double out = 0.0;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (gas_cell(g, cfg, i, j) && std::isfinite(s.ne(i, j))) out = std::max(out, s.ne(i, j));
    }
  }
  return out;
}

double max_gas_emag(const AxisymmetricGrid& g, const StreamerState& s, const StreamerConfig& cfg) {
  double out = 0.0;
  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (gas_cell(g, cfg, i, j) && std::isfinite(s.emag(i, j))) out = std::max(out, std::abs(s.emag(i, j)));
    }
  }
  return out;
}

bool selected_region_cell(const StreamerState& s, int i, int j, const LfaAuditRegionDefinition& region,
                          double active_threshold, double high_field_threshold) {
  switch (region.kind) {
    case LfaAuditRegionKind::AllGas:
      return true;
    case LfaAuditRegionKind::ActiveElectron:
      return std::isfinite(s.ne(i, j)) && s.ne(i, j) > active_threshold;
    case LfaAuditRegionKind::HighField:
      return std::isfinite(s.emag(i, j)) && std::abs(s.emag(i, j)) >= high_field_threshold;
    case LfaAuditRegionKind::CustomMask:
      return true;
  }
  return false;
}

}  // namespace

ElectronRelaxationTable::ElectronRelaxationTable(ElectronRelaxationMetadata metadata,
                                                 std::vector<double> reduced_field_Td,
                                                 std::optional<std::vector<double>> tau_epsilon_s,
                                                 std::optional<std::vector<double>> lambda_epsilon_m)
    : metadata_(std::move(metadata)), reduced_field_Td_(std::move(reduced_field_Td)) {
  if (metadata_.source_id.empty() || metadata_.source_reference.empty() || metadata_.gas_composition.empty() ||
      !std::isfinite(metadata_.pressure_Pa) || metadata_.pressure_Pa <= 0.0 ||
      !std::isfinite(metadata_.temperature_K) || metadata_.temperature_K <= 0.0 ||
      !std::isfinite(metadata_.neutral_density_m3) || metadata_.neutral_density_m3 <= 0.0) {
    throw std::invalid_argument("invalid electron relaxation table metadata");
  }
  if (reduced_field_Td_.size() < 2) throw std::invalid_argument("electron relaxation table requires at least two rows");
  for (std::size_t k = 0; k < reduced_field_Td_.size(); ++k) {
    if (!std::isfinite(reduced_field_Td_[k])) throw std::invalid_argument("nonfinite reduced field in relaxation table");
    if (k > 0 && !(reduced_field_Td_[k] > reduced_field_Td_[k - 1])) {
      throw std::invalid_argument("electron relaxation table E/N must be strictly increasing");
    }
  }
  if (!tau_epsilon_s && !lambda_epsilon_m) {
    throw std::invalid_argument("electron relaxation table must provide tau_epsilon or lambda_epsilon");
  }
  auto validate_quantity = [&](const std::vector<double>& v, const char* name) {
    if (v.size() != reduced_field_Td_.size()) throw std::invalid_argument(std::string(name) + " size mismatch");
    for (double x : v) {
      if (!std::isfinite(x) || x <= 0.0) throw std::invalid_argument(std::string("invalid ") + name);
    }
  };
  if (tau_epsilon_s) {
    validate_quantity(*tau_epsilon_s, "tau_epsilon");
    tau_epsilon_s_ = *tau_epsilon_s;
  }
  if (lambda_epsilon_m) {
    validate_quantity(*lambda_epsilon_m, "lambda_epsilon");
    lambda_epsilon_m_ = *lambda_epsilon_m;
  }
}

ElectronRelaxationLookup ElectronRelaxationTable::lookup(double reduced_field_Td) const {
  ElectronRelaxationLookup out;
  out.reduced_field_Td = reduced_field_Td;
  if (!std::isfinite(reduced_field_Td)) {
    out.status = "NONFINITE_REDUCED_FIELD";
    return out;
  }
  if (reduced_field_Td < reduced_field_Td_.front() || reduced_field_Td > reduced_field_Td_.back()) {
    out.status = "OUTSIDE_RELAXATION_TABLE_RANGE";
    return out;
  }
  out.status = "IN_RANGE";
  if (has_tau_epsilon()) {
    out.has_tau_epsilon = true;
    out.tau_epsilon_s = interpolate_linear(reduced_field_Td_, tau_epsilon_s_, reduced_field_Td);
  }
  if (has_lambda_epsilon()) {
    out.has_lambda_epsilon = true;
    out.lambda_epsilon_m = interpolate_linear(reduced_field_Td_, lambda_epsilon_m_, reduced_field_Td);
  }
  return out;
}

std::vector<LfaAuditRegionDefinition> default_lfa_audit_regions() {
  return {{"ALL_GAS", LfaAuditRegionKind::AllGas, 1e-12, 0.5},
          {"ACTIVE_ELECTRON", LfaAuditRegionKind::ActiveElectron, 1e-12, 0.5},
          {"HIGH_FIELD", LfaAuditRegionKind::HighField, 1e-12, 0.5}};
}

LfaAuditSummary compute_lfa_audit(const AxisymmetricGrid& g, const StreamerState& s,
                                  const StreamerConfig& cfg, double accepted_dt_s,
                                  const ScalarField2D* previous_emag,
                                  const std::vector<unsigned char>* optional_cell_mask,
                                  const ElectronRelaxationTable* relaxation_table,
                                  LfaAuditRegionDefinition region) {
  LfaAuditSummary out;
  out.region_name = region.name.empty() ? "UNNAMED_REGION" : region.name;
  if (optional_cell_mask && optional_cell_mask->size() != g.size()) {
    throw std::invalid_argument("LFA audit mask size does not match grid");
  }
  if (previous_emag && &previous_emag->grid() != &g) {
    throw std::invalid_argument("LFA audit previous field grid mismatch");
  }

  const double ne_max = max_gas_ne(g, s, cfg);
  const double emax = max_gas_emag(g, s, cfg);
  const double active_threshold = region.active_electron_relative_floor * std::max(ne_max, cfg.n_ref);
  const double high_field_threshold = region.high_field_fraction * emax;
  const double e_scale = max_selected_emag(g, s, cfg, optional_cell_mask, region, active_threshold, high_field_threshold);
  const double h_min = std::min(g.dr(), g.dz());
  const double grad_zero_tol = 128.0 * std::numeric_limits<double>::epsilon() *
                               std::max(e_scale / std::max(h_min, std::numeric_limits<double>::min()), 1.0);
  const bool temporal_available = previous_emag && accepted_dt_s > 0.0 && std::isfinite(accepted_dt_s);
  const double dedt_zero_tol = temporal_available
                                   ? 128.0 * std::numeric_limits<double>::epsilon() *
                                         std::max(e_scale / accepted_dt_s, 1.0)
                                   : std::numeric_limits<double>::quiet_NaN();
  out.temporal_status = temporal_available ? "VALID" : "INVALID_INITIAL_SAMPLE";

  std::vector<double> le_values;
  std::vector<double> tau_values;
  std::vector<double> chiL_values;
  std::vector<double> chiT_values;

  if (relaxation_table) {
    out.relaxation_table_min_Td = relaxation_table->min_reduced_field_Td();
    out.relaxation_table_max_Td = relaxation_table->max_reduced_field_Td();
    if (relaxation_table->has_tau_epsilon() && relaxation_table->has_lambda_epsilon()) {
      out.lfa_relaxation_data_status = "AVAILABLE_TAU_AND_LAMBDA";
    } else if (relaxation_table->has_tau_epsilon()) {
      out.lfa_relaxation_data_status = "AVAILABLE_TAU_ONLY";
    } else if (relaxation_table->has_lambda_epsilon()) {
      out.lfa_relaxation_data_status = "AVAILABLE_LAMBDA_ONLY";
    }
    out.lfa_applicability = "UNRESOLVED_NO_APPROVED_CRITERION";
  }

  for (int j = 0; j < g.nz(); ++j) {
    for (int i = 0; i < g.nr(); ++i) {
      if (!gas_cell(g, cfg, i, j)) continue;
      ++out.gas_cells;
      if (!selected_cell(g, i, j, optional_cell_mask)) continue;
      if (!selected_region_cell(s, i, j, region, active_threshold, high_field_threshold)) continue;
      ++out.selected_gas_cells;

      const double e = s.emag(i, j);
      if (!std::isfinite(e) || !std::isfinite(cfg.neutral_density) || cfg.neutral_density <= 0.0) {
        ++out.nonfinite_input_cells;
        continue;
      }
      out.eovern_max_Td = std::max(out.eovern_max_Td, e / cfg.neutral_density * 1e21);

      const double eovern = e / cfg.neutral_density * 1e21;
      ElectronRelaxationLookup relaxation;
      if (relaxation_table) {
        relaxation = relaxation_table->lookup(eovern);
        if (relaxation.status == "IN_RANGE") ++out.relaxation_covered_cells;
        else if (relaxation.status == "OUTSIDE_RELAXATION_TABLE_RANGE") ++out.outside_relaxation_table_cells;
      }

      bool le_valid = false;
      double le_value = std::numeric_limits<double>::quiet_NaN();
      const auto dr = directional_difference(g, s, cfg, i, j, 1, 0, g.dr());
      const auto dz = directional_difference(g, s, cfg, i, j, 0, 1, g.dz());
      if (!dr.valid && !dz.valid) {
        ++out.insufficient_stencil_cells;
      } else {
        const double gr = dr.valid ? dr.value : 0.0;
        const double gz = dz.valid ? dz.value : 0.0;
        const double grad = std::hypot(gr, gz);
        if (!std::isfinite(grad)) {
          ++out.nonfinite_input_cells;
        } else if (grad <= grad_zero_tol) {
          ++out.near_zero_gradE_cells;
        } else {
          le_value = e / grad;
          le_values.push_back(le_value);
          le_valid = true;
        }
      }

      bool tau_valid = false;
      double tau_value = std::numeric_limits<double>::quiet_NaN();
      if (temporal_available) {
        const double ep = (*previous_emag)(i, j);
        if (!std::isfinite(ep)) {
          ++out.nonfinite_input_cells;
        } else {
          const double dedt = (e - ep) / accepted_dt_s;
          if (!std::isfinite(dedt)) {
            ++out.nonfinite_input_cells;
          } else if (std::abs(dedt) <= dedt_zero_tol) {
            ++out.near_zero_dEdt_cells;
          } else {
            tau_value = e / std::abs(dedt);
            tau_values.push_back(tau_value);
            tau_valid = true;
          }
        }
      }
      if (relaxation.status == "IN_RANGE") {
        if (le_valid && relaxation.has_lambda_epsilon) chiL_values.push_back(relaxation.lambda_epsilon_m / le_value);
        if (tau_valid && relaxation.has_tau_epsilon) chiT_values.push_back(relaxation.tau_epsilon_s / tau_value);
      }
    }
  }

  out.le_valid_cells = static_cast<int>(le_values.size());
  out.tauE_valid_cells = static_cast<int>(tau_values.size());
  const double denom = std::max(out.selected_gas_cells, 1);
  out.le_valid_fraction = static_cast<double>(out.le_valid_cells) / denom;
  out.tauE_valid_fraction = static_cast<double>(out.tauE_valid_cells) / denom;
  out.le_min_valid_m = percentile(le_values, 0.0);
  out.le_p05_m = percentile(le_values, 0.05);
  out.le_median_m = percentile(le_values, 0.5);
  out.tauE_min_valid_s = percentile(tau_values, 0.0);
  out.tauE_p05_s = percentile(tau_values, 0.05);
  out.tauE_median_s = percentile(tau_values, 0.5);
  out.chiL_valid_cells = static_cast<int>(chiL_values.size());
  out.chiT_valid_cells = static_cast<int>(chiT_values.size());
  out.chiL_p95 = percentile(chiL_values, 0.95);
  out.chiL_median = percentile(chiL_values, 0.5);
  out.chiT_p95 = percentile(chiT_values, 0.95);
  out.chiT_median = percentile(chiT_values, 0.5);
  out.relaxation_coverage_fraction = static_cast<double>(out.relaxation_covered_cells) / denom;
  out.fraction_outside_relaxation_table = static_cast<double>(out.outside_relaxation_table_cells) / denom;
  return out;
}

LfaAuditSummary LfaAuditHistory::sample(const AxisymmetricGrid& grid, const StreamerState& state,
                                        const StreamerConfig& config, double accepted_dt_s,
                                        const std::vector<unsigned char>* optional_cell_mask,
                                        const ElectronRelaxationTable* relaxation_table,
                                        LfaAuditRegionDefinition region) {
  ScalarField2D previous(grid);
  ScalarField2D* previous_ptr = nullptr;
  if (has_previous_) {
    previous.values() = previous_emag_values_;
    previous_ptr = &previous;
  }
  auto out = compute_lfa_audit(grid, state, config, accepted_dt_s, previous_ptr, optional_cell_mask, relaxation_table,
                               region);
  previous_emag_values_ = state.emag.values();
  has_previous_ = true;
  return out;
}

void LfaAuditHistory::reset() {
  previous_emag_values_.clear();
  has_previous_ = false;
}

}  // namespace streamer_rf::streamer
