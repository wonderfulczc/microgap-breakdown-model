#include "streamer_rf/sp3_geometry.hpp"
#include <petscksp.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace streamer_rf {
namespace {
using Block=std::array<double,4>;
// Face value = (I + h/2 B)^-1 center; outward normal derivative = -B face.
Block face_operator(double k,double h){
  const double s=std::sqrt(6./5.);
  const Block b{k*5./96*(34+11*s),k*5./96*(2+s),
                k*5./96*(2-s),k*5./96*(34-11*s)};
  const double a=1+h/2*b[0],c=h/2*b[1],d=h/2*b[2],e=1+h/2*b[3];
  const double det=a*e-c*d;
  return {(b[0]*e-b[1]*d)/det,(-b[0]*c+b[1]*a)/det,
          (b[2]*e-b[3]*d)/det,(-b[2]*c+b[3]*a)/det};
}
void check(PetscErrorCode error){if(error)throw std::runtime_error("GEOMETRY_SP3_PETSC_ERROR");}
}
struct GeometrySp3Workspace::Impl {
  const Geometry2DView& v;
  Sp3Constants constants;
  std::vector<int> ids;
  int count=0;
  struct Group {
    Mat matrix=nullptr; Vec rhs=nullptr,x=nullptr,residual=nullptr; KSP ksp=nullptr;
    ~Group(){KSPDestroy(&ksp);VecDestroy(&residual);VecDestroy(&x);VecDestroy(&rhs);MatDestroy(&matrix);}
  };
  std::array<Group,3> groups;
  Impl(const Geometry2DView& view,const PressureConfig& p,const SolverTolerances& t):v(view),constants(sp3_constants(p)),ids(v.grid().size(),-1){
    int size;MPI_Comm_size(PETSC_COMM_WORLD,&size);
    if(size!=1)throw std::invalid_argument("GEOMETRY_SP3_SERIAL_ONLY");
    if(!(p.oxygen_torr>0)||!std::isfinite(p.oxygen_torr))throw std::invalid_argument("INVALID_SP3_OXYGEN_PRESSURE");
    const auto& g=v.grid();
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)
      if(v.classify(i,j)==ElectrodeCellType::Gas)ids[j*g.nr()+i]=count++;
    for(int z=0;z<3;++z){
      auto& a=groups[z];const double k=constants.groups[z].k_si;
      check(MatCreateSeqAIJ(PETSC_COMM_SELF,2*count,2*count,7,nullptr,&a.matrix));
      for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
        const int id=ids[j*g.nr()+i];if(id<0)continue;
        const double vol=g.cell_volume(i);
        for(int c=0;c<2;++c){
          const PetscInt row=2*id+c;
          const double kap=c==0?constants.kappa1_sq:constants.kappa2_sq;
          check(MatSetValue(a.matrix,row,row,k*k/kap,ADD_VALUES));
          for(auto dir:{FaceDirection2D::RadialMinus,FaceDirection2D::RadialPlus,FaceDirection2D::AxialMinus,FaceDirection2D::AxialPlus}){
            const auto f=v.face(i,j,dir);if(f.role==FaceRole2D::Axis)continue;
            const double h=f.normal_r!=0?g.dr():g.dz();
            if(f.role==FaceRole2D::GasInterface){
              const double conductance=f.area_m2/(h*vol);
              check(MatSetValue(a.matrix,row,row,conductance,ADD_VALUES));
              check(MatSetValue(a.matrix,row,2*ids[f.neighbor_j*g.nr()+f.neighbor_i]+c,-conductance,ADD_VALUES));
            }else{
              const auto b=face_operator(k,h);
              for(int cc=0;cc<2;++cc)
                check(MatSetValue(a.matrix,row,2*id+cc,f.area_m2/vol*b[2*c+cc],ADD_VALUES));
            }
          }
        }
      }
      check(MatAssemblyBegin(a.matrix,MAT_FINAL_ASSEMBLY));check(MatAssemblyEnd(a.matrix,MAT_FINAL_ASSEMBLY));
      check(MatCreateVecs(a.matrix,&a.x,&a.rhs));check(VecDuplicate(a.rhs,&a.residual));
      check(KSPCreate(PETSC_COMM_SELF,&a.ksp));check(KSPSetOperators(a.ksp,a.matrix,a.matrix));
      check(KSPSetType(a.ksp,KSPGMRES));PC pc;check(KSPGetPC(a.ksp,&pc));check(PCSetType(pc,PCILU));
      check(KSPSetTolerances(a.ksp,t.rtol,t.atol,PETSC_DEFAULT,t.max_iterations));
      check(KSPSetReusePreconditioner(a.ksp,PETSC_TRUE));
      check(KSPSetInitialGuessNonzero(a.ksp,PETSC_FALSE));
    }
  }
};
GeometrySp3Workspace::GeometrySp3Workspace(const Geometry2DView& v,const PressureConfig& p,const SolverTolerances& t):impl_(std::make_unique<Impl>(v,p,t)){}
GeometrySp3Workspace::~GeometrySp3Workspace()=default;
GeometrySp3Audit GeometrySp3Workspace::solve(const ScalarField2D& source,ScalarField2D& sph,Sp3WarmStart* warm){
  auto& n=*impl_;const auto& g=n.v.grid();
  if(&source.grid()!=&g||&sph.grid()!=&g)throw std::invalid_argument("SP3_GRID_IDENTITY_MISMATCH");
  for(double value:source.values())if(!std::isfinite(value)||value<0)throw std::invalid_argument("INVALID_SP3_SOURCE");
  GeometrySp3Audit audit;std::fill(sph.values().begin(),sph.values().end(),0);
  for(int z=0;z<3;++z){
    auto& a=n.groups[z];const auto& c=n.constants;const double k=c.groups[z].k_si;
    const double weights[2]={c.gamma2/(c.gamma2-c.gamma1),-c.gamma1/(c.gamma2-c.gamma1)};
    PetscScalar* b;check(VecGetArray(a.rhs,&b));
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
      int id=n.ids[j*g.nr()+i];if(id<0)continue;
      b[2*id]=k/c.kappa1_sq*source(i,j);b[2*id+1]=k/c.kappa2_sq*source(i,j);
    }
    check(VecRestoreArray(a.rhs,&b));check(VecSet(a.x,0));check(KSPSolve(a.ksp,a.rhs,a.x));
    KSPConvergedReason reason;check(KSPGetConvergedReason(a.ksp,&reason));
    if(reason<=0)throw std::runtime_error("GEOMETRY_SP3_LINEAR_SOLVE_UNRESOLVED");
    PetscInt iterations;check(KSPGetIterationNumber(a.ksp,&iterations));audit.iterations+=iterations;
    check(MatMult(a.matrix,a.x,a.residual));check(VecAYPX(a.residual,-1,a.rhs));
    PetscReal rn,bn;check(VecNorm(a.residual,NORM_2,&rn));check(VecNorm(a.rhs,NORM_2,&bn));
    audit.max_relative_linear_residual=std::max(audit.max_relative_linear_residual,bn>0?rn/bn:rn);
    const PetscScalar* x;check(VecGetArrayRead(a.x,&x));
    if(warm){warm->components[2*z].assign(g.size(),0);warm->components[2*z+1].assign(g.size(),0);}
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
      int id=n.ids[j*g.nr()+i];if(id<0)continue;
      const double xx[2]={PetscRealPart(x[2*id]),PetscRealPart(x[2*id+1])};
      const double photo=c.groups[z].a_si*(weights[0]*xx[0]+weights[1]*xx[1]);
      sph(i,j)+=photo;audit.photo_created_s1+=photo*g.cell_volume(i);
      audit.effective_source_s1+=c.groups[z].a_si/k*source(i,j)*g.cell_volume(i);
      if(warm)for(int cc=0;cc<2;++cc)warm->components[2*z+cc][j*g.nr()+i]=xx[cc];
      for(auto dir:{FaceDirection2D::RadialMinus,FaceDirection2D::RadialPlus,FaceDirection2D::AxialMinus,FaceDirection2D::AxialPlus}){
        const auto f=n.v.face(i,j,dir);
        if(f.role==FaceRole2D::GasInterface||f.role==FaceRole2D::Axis)continue;
        const auto bface=face_operator(k,f.normal_r!=0?g.dr():g.dz());
        double loss=0;
        for(int cc=0;cc<2;++cc)loss+=weights[cc]*(cc==0?c.kappa1_sq:c.kappa2_sq)/k*
            (bface[2*cc]*xx[0]+bface[2*cc+1]*xx[1])*f.area_m2*c.groups[z].a_si/k;
        if(f.role==FaceRole2D::Electrode){audit.metal_loss_s1+=loss;if(z==0)++audit.metal_faces;}
        else audit.outer_loss_s1+=loss;
      }
    }
    check(VecRestoreArrayRead(a.x,&x));
  }
  audit.balance_residual_s1=audit.effective_source_s1-audit.photo_created_s1-audit.metal_loss_s1-audit.outer_loss_s1;
  for(double value:sph.values())if(!std::isfinite(value)||value<0)throw std::runtime_error("SP3_NONPOSITIVE_RATE_NOT_CLIPPED");
  return audit;
}
}
