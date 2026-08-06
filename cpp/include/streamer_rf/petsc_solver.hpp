#pragma once
#include "streamer_rf/types.hpp"
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
namespace streamer_rf {
class PetscSolverContext {
 public:
  PetscSolverContext(const AxisymmetricGrid&, SolverTolerances);
  ~PetscSolverContext(); PetscSolverContext(const PetscSolverContext&)=delete;
  DM dm()const{return dm_;} KSP ksp()const{return ksp_;}
 private: DM dm_{nullptr}; KSP ksp_{nullptr};
};
}
