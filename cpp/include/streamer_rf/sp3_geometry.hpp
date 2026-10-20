#pragma once
#include "streamer_rf/sp3.hpp"
#include "streamer_rf/geometry2d.hpp"
#include <memory>

namespace streamer_rf {
struct GeometrySp3Audit {
  double photo_created_s1{}, effective_source_s1{}, metal_loss_s1{}, outer_loss_s1{};
  double balance_residual_s1{}, max_relative_linear_residual{};
  int iterations{}, metal_faces{};
  std::string boundary_id{"OPAQUE_ABSORBING_METAL_LARSEN_FV_V1"};
};
// Immutable gas/metal topology and matrices; KSP workspace is reusable, not physical history.
// Serial only. Each solve overwrites every RHS and starts from zero.
class GeometrySp3Workspace {
 public:
  GeometrySp3Workspace(const Geometry2DView&,const PressureConfig&,const SolverTolerances&);
  ~GeometrySp3Workspace();
  GeometrySp3Audit solve(const ScalarField2D&,ScalarField2D&,Sp3WarmStart* = nullptr);
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}
