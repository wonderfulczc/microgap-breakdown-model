#include "streamer_rf/poisson.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <cstdint>
namespace streamer_rf {
namespace { constexpr double eps0=8.8541878128e-12;int g_last_elliptic_iterations=0,g_last_poisson_iterations=0;
using Timer=std::chrono::steady_clock;
double elapsed(Timer::time_point a){return std::chrono::duration<double>(Timer::now()-a).count();}
bool profile_poisson=false;
std::vector<PoissonTiming> timings;
void measured_solve(KSP ksp,Mat A,Vec b,Vec u,PoissonTiming& t){
 auto start=Timer::now();KSPSetUp(ksp);t.setup_s=elapsed(start);
 start=Timer::now();KSPSolve(ksp,b,u);t.solve_s=elapsed(start);
 PetscInt its;KSPGetIterationNumber(ksp,&its);t.iterations=its;
 KSPConvergedReason reason;KSPGetConvergedReason(ksp,&reason);t.reason=reason;
 start=Timer::now();Vec residual;VecDuplicate(b,&residual);MatMult(A,u,residual);VecAYPX(residual,-1,b);
 PetscReal rn,bn;VecNorm(residual,NORM_2,&rn);VecNorm(b,NORM_2,&bn);t.true_relative_residual=bn>0?rn/bn:rn;VecDestroy(&residual);t.residual_s=elapsed(start);
}
void assemble(const AxisymmetricGrid& g,DM dm,Mat A,Vec b,const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,double shift,const ScalarField2D* boundary,const AxisymmetricNeedlePlaneGeometry* geometry=nullptr,double applied_voltage=0.0,const Geometry2DView* view=nullptr){
 DMDALocalInfo info;DMDAGetLocalInfo(dm,&info); const double dr=g.dr(),dz=g.dz();
 PetscScalar** barr;DMDAVecGetArray(dm,b,&barr);
 auto st=[](int ii,int jj){MatStencil s{};s.i=ii;s.j=jj;return s;};
 for(int j=info.ys;j<info.ys+info.ym;++j)for(int i=info.xs;i<info.xs+info.xm;++i){
  MatStencil row=st(i,j),col[5];PetscScalar val[5];int n=0;const bool boundary_cell=i==g.nr()-1||j==0||j==g.nz()-1;
  const auto electrode_cell=view?view->classify(i,j):(geometry?geometry->classify(g,i,j):ElectrodeCellType::Gas);
  if(electrode_cell!=ElectrodeCellType::Gas){col[n]=row;val[n++]=1;barr[j][i]=electrode_cell==ElectrodeCellType::HighVoltageElectrode?applied_voltage:0.0;}
  else if(boundary_cell){
   auto dirichlet=[&](double v){col[n]=row;val[n++]=1;barr[j][i]=v;};
   auto neumann=[&](int ni,int nj,double h,double v){col[n]=row;val[n++]=1.0/h;col[n]=st(ni,nj);val[n++]=-1.0/h;barr[j][i]=v;};
   if(boundary)dirichlet((*boundary)(i,j));
   else if(j==0){
    if(bc.z_lower.kind==BoundaryKind::Neumann)neumann(i,j+1,g.dz(),-bc.z_lower.value);
    else dirichlet(bc.z_lower.value);
   }else if(j==g.nz()-1){
    if(bc.z_upper.kind==BoundaryKind::Neumann)neumann(i,j-1,g.dz(),bc.z_upper.value);
    else dirichlet(bc.z_upper.value);
   }else if(i==g.nr()-1){
    if(bc.r_outer.kind==BoundaryKind::Neumann)neumann(i-1,j,g.dr(),bc.r_outer.value);
    else dirichlet(bc.r_outer.value);
   }
  }
  else {const double r=g.r(i),rp=g.radial_face(i+1),rm=g.radial_face(i),ar=rp/(r*dr*dr),al=rm/(r*dr*dr),az=1/(dz*dz);
   col[n]=st(i,j);val[n++]=-(ar+al+2*az)-shift;
   col[n]=st(i+1,j);val[n++]=ar;if(i>0){col[n]=st(i-1,j);val[n++]=al;}
   col[n]=st(i,j+1);val[n++]=az;col[n]=st(i,j-1);val[n++]=az;
   barr[j][i]=rhs(i,j);}
  if(A)MatSetValuesStencil(A,1,&row,n,col,val,INSERT_VALUES);
 }DMDAVecRestoreArray(dm,b,&barr);}
void solve_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,ScalarField2D& x,const SolverTolerances&t,double shift,const ScalarField2D* boundary=nullptr){
 PetscSolverContext ctx(rhs.grid(),t);Mat A;Vec b,u;DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);
 assemble(rhs.grid(),ctx.dm(),A,b,rhs,bc,shift,boundary);MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);VecAssemblyBegin(b);VecAssemblyEnd(b);
 KSPSetOperators(ctx.ksp(),A,A);KSPSolve(ctx.ksp(),b,u);PetscInt its;KSPGetIterationNumber(ctx.ksp(),&its);g_last_elliptic_iterations=static_cast<int>(its);KSPConvergedReason reason;KSPGetConvergedReason(ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("PETSc elliptic solve diverged");
 Vec natural,all;DMDACreateNaturalVector(ctx.dm(),&natural);DMDAGlobalToNaturalBegin(ctx.dm(),u,INSERT_VALUES,natural);DMDAGlobalToNaturalEnd(ctx.dm(),u,INSERT_VALUES,natural);VecScatter scatter;VecScatterCreateToAll(natural,&scatter,&all);VecScatterBegin(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);
 const PetscScalar* arr;VecGetArrayRead(all,&arr);for(std::size_t k=0;k<x.values().size();++k)x.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(all,&arr);
 VecScatterDestroy(&scatter);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);
}
void solve_electrode_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,ScalarField2D& x,const SolverTolerances&t,const AxisymmetricNeedlePlaneGeometry* geometry,const Geometry2DView* view,double applied_voltage,const ScalarField2D* boundary=nullptr){
 PoissonTiming timing;auto clock=Timer::now();
 PetscSolverContext ctx(rhs.grid(),t);Mat A;Vec b,u;DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);
 assemble(rhs.grid(),ctx.dm(),A,b,rhs,bc,0,boundary,geometry,applied_voltage,view);MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);VecAssemblyBegin(b);VecAssemblyEnd(b);
 timing.assembly_s=elapsed(clock);KSPSetOperators(ctx.ksp(),A,A);if(profile_poisson)measured_solve(ctx.ksp(),A,b,u,timing);else KSPSolve(ctx.ksp(),b,u);PetscInt its;KSPGetIterationNumber(ctx.ksp(),&its);g_last_elliptic_iterations=static_cast<int>(its);KSPConvergedReason reason;KSPGetConvergedReason(ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("PETSc electrode elliptic solve diverged");
 clock=Timer::now();
 Vec natural,all;DMDACreateNaturalVector(ctx.dm(),&natural);DMDAGlobalToNaturalBegin(ctx.dm(),u,INSERT_VALUES,natural);DMDAGlobalToNaturalEnd(ctx.dm(),u,INSERT_VALUES,natural);VecScatter scatter;VecScatterCreateToAll(natural,&scatter,&all);VecScatterBegin(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);
 const PetscScalar* arr;VecGetArrayRead(all,&arr);for(std::size_t k=0;k<x.values().size();++k)x.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(all,&arr);
 VecScatterDestroy(&scatter);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);
 timing.transfer_s=elapsed(clock);if(profile_poisson)timings.push_back(timing);
}
void solve_robin_impl(const ScalarField2D& rhs,ScalarField2D& x,double shift,double a,double cross,const ScalarField2D& partner,const SolverTolerances&t){
 PetscSolverContext ctx(rhs.grid(),t);Mat A;Vec b,u;DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);DMDALocalInfo info;DMDAGetLocalInfo(ctx.dm(),&info);const auto&g=rhs.grid();
 PetscScalar** barr;DMDAVecGetArray(ctx.dm(),b,&barr);
 auto st=[](int ii,int jj){MatStencil s{};s.i=ii;s.j=jj;return s;};
 for(int j=info.ys;j<info.ys+info.ym;++j)for(int i=info.xs;i<info.xs+info.xm;++i){MatStencil row=st(i,j),col[5];PetscScalar val[5];int n=0;double q=rhs(i,j);const bool bd=i==g.nr()-1||j==0||j==g.nz()-1;
  if(bd){double h;int i1=i,j1=j,i2=i,j2=j;if(i==g.nr()-1){h=g.dr();i1=i-1;i2=i-2;}else if(j==0){h=g.dz();j1=j+1;j2=j+2;}else{h=g.dz();j1=j-1;j2=j-2;}col[n]=row;val[n++]=3/(2*h*h)+a/h;col[n]=st(i1,j1);val[n++]=-4/(2*h*h);col[n]=st(i2,j2);val[n++]=1/(2*h*h);q=-cross*partner(i,j)/h;}
  else{const double r=g.r(i),rp=g.radial_face(i+1),rm=g.radial_face(i),ar=rp/(r*g.dr()*g.dr()),al=rm/(r*g.dr()*g.dr()),az=1/(g.dz()*g.dz());col[n]=row;val[n++]=ar+al+2*az+shift;col[n]=st(i+1,j);val[n++]=-ar;if(i>0){col[n]=st(i-1,j);val[n++]=-al;}col[n]=st(i,j+1);val[n++]=-az;col[n]=st(i,j-1);val[n++]=-az;q=-q;}
  MatSetValuesStencil(A,1,&row,n,col,val,INSERT_VALUES);barr[j][i]=q;}DMDAVecRestoreArray(ctx.dm(),b,&barr);
 MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);VecAssemblyBegin(b);VecAssemblyEnd(b);KSPSetOperators(ctx.ksp(),A,A);KSPSolve(ctx.ksp(),b,u);PetscInt its;KSPGetIterationNumber(ctx.ksp(),&its);g_last_elliptic_iterations=static_cast<int>(its);KSPConvergedReason reason;KSPGetConvergedReason(ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("PETSc Robin solve diverged");Vec natural,all;DMDACreateNaturalVector(ctx.dm(),&natural);DMDAGlobalToNaturalBegin(ctx.dm(),u,INSERT_VALUES,natural);DMDAGlobalToNaturalEnd(ctx.dm(),u,INSERT_VALUES,natural);VecScatter sc;VecScatterCreateToAll(natural,&sc,&all);VecScatterBegin(sc,natural,all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(sc,natural,all,INSERT_VALUES,SCATTER_FORWARD);const PetscScalar*arr;VecGetArrayRead(all,&arr);for(std::size_t k=0;k<x.values().size();++k)x.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(all,&arr);VecScatterDestroy(&sc);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);
}
}
void set_poisson_profiling(bool value){profile_poisson=value;}
const std::vector<PoissonTiming>& poisson_timing_records(){return timings;}
void clear_poisson_timing_records(){timings.clear();}
PoissonDiagnosticSnapshot snapshot_poisson_diagnostics(){return {g_last_elliptic_iterations,g_last_poisson_iterations,profile_poisson,timings};}
void restore_poisson_diagnostics(const PoissonDiagnosticSnapshot& s){g_last_elliptic_iterations=s.elliptic_iterations;g_last_poisson_iterations=s.poisson_iterations;profile_poisson=s.profiling;timings=s.records;}
std::string electrode_operator_identity(const Geometry2DView& v,const PoissonBoundaryConfig& bc,const SolverTolerances& t){
 const auto& g=v.grid();std::uint64_t hash=1469598103934665603ULL;
 for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){hash^=static_cast<unsigned>(v.classify(i,j));hash*=1099511628211ULL;}
 std::ostringstream s;s<<std::setprecision(17)<<"STAGE_C_PETSC_FV_V1:epsilon="<<eps0<<":"<<g.nr()<<":"<<g.nz()<<":"<<g.dr()<<":"<<g.dz()<<":"<<g.z(0)<<":"<<hash<<":"<<v.geometry_id()<<":"<<v.revision()<<":"<<v.reference_plane_id()<<":basis_HV1_RETURN0";
 for(auto b:{bc.r_outer,bc.z_lower,bc.z_upper})s<<":"<<int(b.kind)<<":"<<b.value;
 s<<":"<<t.rtol<<":"<<t.atol<<":"<<t.max_iterations;
 char* options=nullptr;PetscOptionsGetAll(nullptr,&options);if(options){s<<":"<<options;PetscFree(options);}return s.str();
}
struct FixedElectrodePoisson::Impl {
 const Geometry2DView& view;PoissonBoundaryConfig bc;SolverTolerances tolerance;std::string key;
 PetscSolverContext ctx;Mat A{};Vec b{},u{},natural{},all{};VecScatter scatter{};ScalarField2D rhs;bool first=true;double creation_s=0;
 Impl(const Geometry2DView& v,const PoissonBoundaryConfig& c,const SolverTolerances& t):view(v),bc(c),tolerance(t),key(electrode_operator_identity(v,c,t)),ctx(v.grid(),t),rhs(v.grid()){
  auto clock=Timer::now();for(auto bd:{bc.r_outer,bc.z_lower,bc.z_upper})if(bd.kind==BoundaryKind::OpenCharge)throw std::invalid_argument("fixed workspace requires explicit Dirichlet/Neumann; OpenCharge not enabled");
  DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);
  assemble(v.grid(),ctx.dm(),A,b,rhs,bc,0,nullptr,nullptr,0,&v);MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  KSPSetOperators(ctx.ksp(),A,A);KSPSetInitialGuessNonzero(ctx.ksp(),PETSC_FALSE);
  DMDACreateNaturalVector(ctx.dm(),&natural);VecScatterCreateToAll(natural,&scatter,&all);creation_s=elapsed(clock);
 }
 ~Impl(){VecScatterDestroy(&scatter);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);}
};
FixedElectrodePoisson::FixedElectrodePoisson(const Geometry2DView& v,const PoissonBoundaryConfig& b,const SolverTolerances&t):impl_(std::make_unique<Impl>(v,b,t)){}
FixedElectrodePoisson::~FixedElectrodePoisson()=default;
void FixedElectrodePoisson::solve(const ScalarField2D& rho,double voltage,ScalarField2D& phi){
 auto& p=*impl_;if(&rho.grid()!=&p.view.grid()||&phi.grid()!=&p.view.grid()||electrode_operator_identity(p.view,p.bc,p.tolerance)!=p.key)throw std::invalid_argument("fixed Poisson cache identity changed; construct a new workspace");
 PoissonTiming timing;timing.reused=!p.first;auto clock=Timer::now();
 for(std::size_t k=0;k<rho.values().size();++k)p.rhs.values()[k]=-rho.values()[k]/eps0;
 assemble(rho.grid(),p.ctx.dm(),nullptr,p.b,p.rhs,p.bc,0,nullptr,nullptr,voltage,&p.view);
 VecAssemblyBegin(p.b);VecAssemblyEnd(p.b);VecSet(p.u,0);timing.assembly_s=elapsed(clock)+(p.first?p.creation_s:0);
 if(profile_poisson)measured_solve(p.ctx.ksp(),p.A,p.b,p.u,timing);else KSPSolve(p.ctx.ksp(),p.b,p.u);
 PetscInt its;KSPGetIterationNumber(p.ctx.ksp(),&its);g_last_elliptic_iterations=g_last_poisson_iterations=its;
 KSPConvergedReason reason;KSPGetConvergedReason(p.ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("cached native Poisson diverged");
 clock=Timer::now();DMDAGlobalToNaturalBegin(p.ctx.dm(),p.u,INSERT_VALUES,p.natural);DMDAGlobalToNaturalEnd(p.ctx.dm(),p.u,INSERT_VALUES,p.natural);
 VecScatterBegin(p.scatter,p.natural,p.all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(p.scatter,p.natural,p.all,INSERT_VALUES,SCATTER_FORWARD);
 const PetscScalar* arr;VecGetArrayRead(p.all,&arr);for(std::size_t k=0;k<phi.values().size();++k)phi.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(p.all,&arr);
 timing.transfer_s=elapsed(clock);if(profile_poisson)timings.push_back(timing);p.first=false;
}
bool FixedElectrodePoisson::matches_fresh_matrix(const ScalarField2D& rho,double voltage)const{
 auto& p=*impl_;Mat fresh;Vec b;DMCreateMatrix(p.ctx.dm(),&fresh);VecDuplicate(p.b,&b);
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 assemble(rho.grid(),p.ctx.dm(),fresh,b,rhs,p.bc,0,nullptr,nullptr,voltage,&p.view);MatAssemblyBegin(fresh,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(fresh,MAT_FINAL_ASSEMBLY);
 PetscBool equal;MatEqual(p.A,fresh,&equal);MatDestroy(&fresh);VecDestroy(&b);return equal;
}
void solve_potential(const ScalarField2D& rho,double background,const PoissonBoundaryConfig& bc,ScalarField2D& phi,const SolverTolerances&t,const OpenBoundaryOptions&o){
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 const bool open=bc.r_outer.kind==BoundaryKind::OpenCharge||bc.z_lower.kind==BoundaryKind::OpenCharge||bc.z_upper.kind==BoundaryKind::OpenCharge;ScalarField2D bd(rho.grid());if(open)bd=open_charge_boundary(rho,o);
 solve_elliptic(rhs,bc,phi,t,0,open?&bd:nullptr);g_last_poisson_iterations=g_last_elliptic_iterations;for(int j=0;j<rho.grid().nz();++j)for(int i=0;i<rho.grid().nr();++i)phi(i,j)+=-background*rho.grid().z(j);
}
void solve_potential_with_electrodes(const ScalarField2D& rho,const AxisymmetricNeedlePlaneGeometry& geometry,double applied_voltage,const PoissonBoundaryConfig& bc,ScalarField2D& phi,const SolverTolerances&t,const OpenBoundaryOptions&o){
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 const bool open=bc.r_outer.kind==BoundaryKind::OpenCharge||bc.z_lower.kind==BoundaryKind::OpenCharge||bc.z_upper.kind==BoundaryKind::OpenCharge;ScalarField2D bd(rho.grid());if(open)bd=open_charge_boundary(rho,o);
 solve_electrode_elliptic(rhs,bc,phi,t,&geometry,nullptr,applied_voltage,open?&bd:nullptr);g_last_poisson_iterations=g_last_elliptic_iterations;
}
void solve_potential_with_electrodes(const ScalarField2D& rho,const Geometry2DView& geometry,double applied_voltage,const PoissonBoundaryConfig& bc,ScalarField2D& phi,const SolverTolerances&t,const OpenBoundaryOptions&o){
 if (&rho.grid()!=&geometry.grid()) throw std::invalid_argument("geometry/grid identity mismatch");
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 const bool open=bc.r_outer.kind==BoundaryKind::OpenCharge||bc.z_lower.kind==BoundaryKind::OpenCharge||bc.z_upper.kind==BoundaryKind::OpenCharge;ScalarField2D bd(rho.grid());if(open)bd=open_charge_boundary(rho,o);
 solve_electrode_elliptic(rhs,bc,phi,t,nullptr,&geometry,applied_voltage,open?&bd:nullptr);g_last_poisson_iterations=g_last_elliptic_iterations;
}
void solve_axisymmetric_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,ScalarField2D& solution,double shift,const SolverTolerances&t){solve_elliptic(rhs,bc,solution,t,shift);}
void solve_axisymmetric_dirichlet(const ScalarField2D&rhs,const ScalarField2D&boundary,ScalarField2D&solution,double shift,const SolverTolerances&t){PoissonBoundaryConfig bc;solve_elliptic(rhs,bc,solution,t,shift,&boundary);}
double ring_potential_on_axis(double q,double a,double z){return q/(4*3.14159265358979323846*eps0*std::sqrt(a*a+z*z));}
double complete_elliptic_k(double m){if(m<0||m>=1)return std::numeric_limits<double>::quiet_NaN();return std::comp_ellint_1(std::sqrt(m));}
double ring_potential(double q,double a,double R,double z){if(R==0)return ring_potential_on_axis(q,a,z);const double den=(R+a)*(R+a)+z*z,m=std::min(1-1e-15,4*R*a/den);return q*complete_elliptic_k(m)/(2*3.14159265358979323846*3.14159265358979323846*eps0*std::sqrt(den));}
ScalarField2D open_charge_boundary(const ScalarField2D&rho,const OpenBoundaryOptions&o){
 const auto&g=rho.grid();ScalarField2D out(g);double maxrho=0;for(double v:rho.values())maxrho=std::max(maxrho,std::abs(v));int rank,size;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&size);constexpr double xi=0.5773502691896257645;
 struct Source{double r,z,q;};std::vector<Source> sources;sources.reserve(g.size()/std::max(size,1));
 for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){const int id=j*g.nr()+i;if(id%size!=rank||rho(i,j)==0||(maxrho>0&&std::abs(rho(i,j))<o.threshold*maxrho))continue;if(o.mode==RingQuadratureMode::CellCenterRing)sources.push_back({g.r(i),g.z(j),rho(i,j)*g.cell_volume(i)});else for(int qr:{-1,1})for(int qz:{-1,1}){double r=g.r(i)+qr*xi*g.dr()/2;sources.push_back({r,g.z(j)+qz*xi*g.dz()/2,rho(i,j)*2*3.14159265358979323846*r*g.dr()*g.dz()/4});}}
 for(int J=0;J<g.nz();++J)for(int I=0;I<g.nr();++I)if(I==g.nr()-1||J==0||J==g.nz()-1){double local=0;for(const auto&s:sources)local+=ring_potential(s.q,s.r,g.r(I),g.z(J)-s.z);double global;MPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);out(I,J)=global;}return out;
}
void solve_axisymmetric_robin(const ScalarField2D&r,ScalarField2D&x,double shift,double self,double cross,const ScalarField2D&p,const SolverTolerances&t){solve_robin_impl(r,x,shift,self,cross,p,t);}
int last_elliptic_iterations(){return g_last_elliptic_iterations;}
int last_poisson_iterations(){return g_last_poisson_iterations;}
}
