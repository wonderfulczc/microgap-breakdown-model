#include "streamer_rf/sp3.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace streamer_rf {
namespace {int g_last_sp3_ksp_iterations=0,g_last_sp3_boundary_iterations=0;}
Sp3Constants sp3_constants(const PressureConfig&p){
 const double s=std::sqrt(6.0/5.0);Sp3Constants c;c.kappa1_sq=3.0/7-2.0*s/7;c.kappa2_sq=3.0/7+2.0*s/7;c.gamma1=5.0/7*(1-3*s);c.gamma2=5.0/7*(1+3*s);const double A[3]={0.0067,0.0346,0.3059},lam[3]={0.0447,0.1121,0.5994};for(int j=0;j<3;++j)c.groups[j]={A[j],lam[j],A[j]*p.oxygen_torr*100.0,lam[j]*p.oxygen_torr*100.0};return c;
}
void solve_photoionization(const ScalarField2D&src,const PressureConfig&p,const PoissonBoundaryConfig&bc,ScalarField2D&sph,const SolverTolerances&t,const Sp3BoundaryOptions&bo,std::array<Sp3GroupDiagnostics,3>*diag,Sp3WarmStart*warm){
 g_last_sp3_ksp_iterations=0;g_last_sp3_boundary_iterations=0;const auto c=sp3_constants(p);std::fill(sph.values().begin(),sph.values().end(),0.0);const double root=std::sqrt(6.0/5.0),alpha1=5.0/96*(34+11*root),alpha2=5.0/96*(34-11*root),beta1=5.0/96*(2-root),beta2=5.0/96*(2+root);int group=0;
 for(const auto&g:c.groups){
  ScalarField2D b1(src.grid()),b2(src.grid()),x1(src.grid()),x2(src.grid()),old1(src.grid()),old2(src.grid());
  if(warm&&warm->components[2*group].size()==src.values().size()){x1.values()=warm->components[2*group];x2.values()=warm->components[2*group+1];}
  for(std::size_t k=0;k<src.values().size();++k){b1.values()[k]=-g.k_si/c.kappa1_sq*src.values()[k];b2.values()[k]=-g.k_si/c.kappa2_sq*src.values()[k];}
  Sp3GroupDiagnostics gd;bool converged=false;
  for(int it=1;it<=bo.max_iterations;++it){
   old1.values()=x1.values();old2.values()=x2.values();
   if(bo.zero_dirichlet){solve_axisymmetric_elliptic(b1,bc,x1,g.k_si*g.k_si/c.kappa1_sq,t);g_last_sp3_ksp_iterations+=last_elliptic_iterations();solve_axisymmetric_elliptic(b2,bc,x2,g.k_si*g.k_si/c.kappa2_sq,t);g_last_sp3_ksp_iterations+=last_elliptic_iterations();}
   else{solve_axisymmetric_robin(b1,x1,g.k_si*g.k_si/c.kappa1_sq,g.k_si*alpha1,g.k_si*beta2,x2,t);g_last_sp3_ksp_iterations+=last_elliptic_iterations();solve_axisymmetric_robin(b2,x2,g.k_si*g.k_si/c.kappa2_sq,g.k_si*alpha2,g.k_si*beta1,x1,t);g_last_sp3_ksp_iterations+=last_elliptic_iterations();}
   double change=0,norm=0;const auto&grid=src.grid();for(int j=0;j<grid.nz();++j)for(int i=0;i<grid.nr();++i)if(i==grid.nr()-1||j==0||j==grid.nz()-1){change=std::max(change,std::max(std::abs(x1(i,j)-old1(i,j)),std::abs(x2(i,j)-old2(i,j))));norm=std::max(norm,std::max(std::abs(x1(i,j)),std::abs(x2(i,j))));}double residual=change/std::max(bo.atol,norm);gd.residual_history.push_back(residual);if(it==1)gd.initial_residual=residual;gd.final_residual=residual;gd.iterations=it;if(change<=bo.atol+bo.rtol*norm){converged=true;break;}for(std::size_t k=0;k<x1.values().size();++k){x1.values()[k]=bo.relaxation*x1.values()[k]+(1-bo.relaxation)*old1.values()[k];x2.values()[k]=bo.relaxation*x2.values()[k]+(1-bo.relaxation)*old2.values()[k];}
  }
  gd.converged=converged;gd.sph_change=gd.final_residual;g_last_sp3_boundary_iterations+=gd.iterations;if(!converged)throw std::runtime_error("SP3 coupled Robin fixed-point did not converge: group="+std::to_string(group+1)+" initial="+std::to_string(gd.initial_residual)+" final="+std::to_string(gd.final_residual));if(diag)(*diag)[group]=gd;if(warm){warm->components[2*group]=x1.values();warm->components[2*group+1]=x2.values();}for(std::size_t k=0;k<src.values().size();++k){double psi=(c.gamma2*x1.values()[k]-c.gamma1*x2.values()[k])/(c.gamma2-c.gamma1);sph.values()[k]+=g.a_si*psi;}++group;
 }
}
int last_sp3_ksp_iterations(){return g_last_sp3_ksp_iterations;}
int last_sp3_boundary_iterations(){return g_last_sp3_boundary_iterations;}
}
