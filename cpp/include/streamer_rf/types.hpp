#pragma once
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace streamer_rf {
struct UnitMetadata { std::string quantity, unit, convention; };
struct SolverTolerances { double rtol{1e-10}, atol{1e-14}; int max_iterations{10000}; };
enum class BoundaryKind { Dirichlet, Neumann, BackgroundField, OpenCharge, Robin };
struct BoundaryCondition { BoundaryKind kind{BoundaryKind::Dirichlet}; double value{0.0}; };
struct PressureConfig { double oxygen_torr{150.0}, total_torr{760.0}, quenching_torr{30.0}; };

class AxisymmetricGrid {
 public:
  AxisymmetricGrid(int nr, int nz, double r_max, double z_min, double z_max);
  int nr() const { return nr_; } int nz() const { return nz_; }
  double dr() const { return dr_; } double dz() const { return dz_; }
  double r(int i) const { return (i + 0.5) * dr_; }
  double z(int j) const { return z_min_ + (j + 0.5) * dz_; }
  double radial_face(int i) const { return i * dr_; }
  double cell_volume(int i) const;
  std::size_t size() const { return static_cast<std::size_t>(nr_) * nz_; }
 private: int nr_, nz_; double r_max_, z_min_, z_max_, dr_, dz_;
};

class ScalarField2D {
 public:
  explicit ScalarField2D(const AxisymmetricGrid& grid, UnitMetadata units = {})
    : grid_(&grid), values_(grid.size()), units_(std::move(units)) {}
  double& operator()(int i,int j) { return values_.at(static_cast<std::size_t>(j)*grid_->nr()+i); }
  double operator()(int i,int j) const { return values_.at(static_cast<std::size_t>(j)*grid_->nr()+i); }
  std::vector<double>& values(){return values_;} const std::vector<double>& values()const{return values_;}
  const AxisymmetricGrid& grid()const{return *grid_;}
 private: const AxisymmetricGrid* grid_; std::vector<double> values_; UnitMetadata units_;
};
}

