#include "streamer_rf/petsc_solver.hpp"
namespace streamer_rf {
PetscSolverContext::PetscSolverContext(const AxisymmetricGrid& g,SolverTolerances t){
 DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,
   g.nr(),g.nz(),PETSC_DECIDE,PETSC_DECIDE,1,2,nullptr,nullptr,&dm_);
 DMSetUp(dm_); KSPCreate(PETSC_COMM_WORLD,&ksp_); KSPSetType(ksp_,KSPBCGS); KSPSetTolerances(ksp_,t.rtol,t.atol,PETSC_DEFAULT,t.max_iterations);
 PC pc;KSPGetPC(ksp_,&pc);PCSetType(pc,PCGAMG);KSPSetFromOptions(ksp_);
}
PetscSolverContext::~PetscSolverContext(){if(ksp_)KSPDestroy(&ksp_);if(dm_)DMDestroy(&dm_);}
}
