#include "streamer_rf/poisson.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace streamer_rf {
namespace { constexpr double eps0=8.8541878128e-12;int g_last_elliptic_iterations=0,g_last_poisson_iterations=0;
void assemble(const AxisymmetricGrid& g,DM dm,Mat A,Vec b,const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,double shift,const ScalarField2D* boundary,const AxisymmetricNeedlePlaneGeometry* geometry=nullptr,double applied_voltage=0.0){
 DMDALocalInfo info;DMDAGetLocalInfo(dm,&info); const double dr=g.dr(),dz=g.dz();
 PetscScalar** barr;DMDAVecGetArray(dm,b,&barr);
 auto st=[](int ii,int jj){MatStencil s{};s.i=ii;s.j=jj;return s;};
 for(int j=info.ys;j<info.ys+info.ym;++j)for(int i=info.xs;i<info.xs+info.xm;++i){
  MatStencil row=st(i,j),col[5];PetscScalar val[5];int n=0;const bool boundary_cell=i==g.nr()-1||j==0||j==g.nz()-1;
  const auto electrode_cell=geometry?geometry->classify(g,i,j):ElectrodeCellType::Gas;
  if(electrode_cell!=ElectrodeCellType::Gas){col[n]=row;val[n++]=1;barr[j][i]=electrode_cell==ElectrodeCellType::HighVoltageElectrode?applied_voltage:0.0;}
  else if(boundary_cell){col[n]=row;val[n++]=1;barr[j][i]=boundary?(*boundary)(i,j):bc.r_outer.value;}
  else {const double r=g.r(i),rp=g.radial_face(i+1),rm=g.radial_face(i),ar=rp/(r*dr*dr),al=rm/(r*dr*dr),az=1/(dz*dz);
   col[n]=st(i,j);val[n++]=-(ar+al+2*az)-shift;
   col[n]=st(i+1,j);val[n++]=ar;if(i>0){col[n]=st(i-1,j);val[n++]=al;}
   col[n]=st(i,j+1);val[n++]=az;col[n]=st(i,j-1);val[n++]=az;
   barr[j][i]=rhs(i,j);}
  MatSetValuesStencil(A,1,&row,n,col,val,INSERT_VALUES);
 }DMDAVecRestoreArray(dm,b,&barr);}
void solve_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,ScalarField2D& x,const SolverTolerances&t,double shift,const ScalarField2D* boundary=nullptr){
 PetscSolverContext ctx(rhs.grid(),t);Mat A;Vec b,u;DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);
 assemble(rhs.grid(),ctx.dm(),A,b,rhs,bc,shift,boundary);MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);VecAssemblyBegin(b);VecAssemblyEnd(b);
 KSPSetOperators(ctx.ksp(),A,A);KSPSolve(ctx.ksp(),b,u);PetscInt its;KSPGetIterationNumber(ctx.ksp(),&its);g_last_elliptic_iterations=static_cast<int>(its);KSPConvergedReason reason;KSPGetConvergedReason(ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("PETSc elliptic solve diverged");
 Vec natural,all;DMDACreateNaturalVector(ctx.dm(),&natural);DMDAGlobalToNaturalBegin(ctx.dm(),u,INSERT_VALUES,natural);DMDAGlobalToNaturalEnd(ctx.dm(),u,INSERT_VALUES,natural);VecScatter scatter;VecScatterCreateToAll(natural,&scatter,&all);VecScatterBegin(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);
 const PetscScalar* arr;VecGetArrayRead(all,&arr);for(std::size_t k=0;k<x.values().size();++k)x.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(all,&arr);
 VecScatterDestroy(&scatter);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);
}
void solve_electrode_elliptic(const ScalarField2D& rhs,const PoissonBoundaryConfig& bc,ScalarField2D& x,const SolverTolerances&t,const AxisymmetricNeedlePlaneGeometry& geometry,double applied_voltage,const ScalarField2D* boundary=nullptr){
 PetscSolverContext ctx(rhs.grid(),t);Mat A;Vec b,u;DMCreateMatrix(ctx.dm(),&A);DMCreateGlobalVector(ctx.dm(),&b);VecDuplicate(b,&u);
 assemble(rhs.grid(),ctx.dm(),A,b,rhs,bc,0,boundary,&geometry,applied_voltage);MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);VecAssemblyBegin(b);VecAssemblyEnd(b);
 KSPSetOperators(ctx.ksp(),A,A);KSPSolve(ctx.ksp(),b,u);PetscInt its;KSPGetIterationNumber(ctx.ksp(),&its);g_last_elliptic_iterations=static_cast<int>(its);KSPConvergedReason reason;KSPGetConvergedReason(ctx.ksp(),&reason);if(reason<0)throw std::runtime_error("PETSc electrode elliptic solve diverged");
 Vec natural,all;DMDACreateNaturalVector(ctx.dm(),&natural);DMDAGlobalToNaturalBegin(ctx.dm(),u,INSERT_VALUES,natural);DMDAGlobalToNaturalEnd(ctx.dm(),u,INSERT_VALUES,natural);VecScatter scatter;VecScatterCreateToAll(natural,&scatter,&all);VecScatterBegin(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);VecScatterEnd(scatter,natural,all,INSERT_VALUES,SCATTER_FORWARD);
 const PetscScalar* arr;VecGetArrayRead(all,&arr);for(std::size_t k=0;k<x.values().size();++k)x.values()[k]=PetscRealPart(arr[k]);VecRestoreArrayRead(all,&arr);
 VecScatterDestroy(&scatter);VecDestroy(&all);VecDestroy(&natural);VecDestroy(&u);VecDestroy(&b);MatDestroy(&A);
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
void solve_potential(const ScalarField2D& rho,double background,const PoissonBoundaryConfig& bc,ScalarField2D& phi,const SolverTolerances&t,const OpenBoundaryOptions&o){
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 const bool open=bc.r_outer.kind==BoundaryKind::OpenCharge||bc.z_lower.kind==BoundaryKind::OpenCharge||bc.z_upper.kind==BoundaryKind::OpenCharge;ScalarField2D bd(rho.grid());if(open)bd=open_charge_boundary(rho,o);
 solve_elliptic(rhs,bc,phi,t,0,open?&bd:nullptr);g_last_poisson_iterations=g_last_elliptic_iterations;for(int j=0;j<rho.grid().nz();++j)for(int i=0;i<rho.grid().nr();++i)phi(i,j)+=-background*rho.grid().z(j);
}
void solve_potential_with_electrodes(const ScalarField2D& rho,const AxisymmetricNeedlePlaneGeometry& geometry,double applied_voltage,const PoissonBoundaryConfig& bc,ScalarField2D& phi,const SolverTolerances&t,const OpenBoundaryOptions&o){
 ScalarField2D rhs(rho.grid());for(std::size_t k=0;k<rho.values().size();++k)rhs.values()[k]=-rho.values()[k]/eps0;
 const bool open=bc.r_outer.kind==BoundaryKind::OpenCharge||bc.z_lower.kind==BoundaryKind::OpenCharge||bc.z_upper.kind==BoundaryKind::OpenCharge;ScalarField2D bd(rho.grid());if(open)bd=open_charge_boundary(rho,o);
 solve_electrode_elliptic(rhs,bc,phi,t,geometry,applied_voltage,open?&bd:nullptr);g_last_poisson_iterations=g_last_elliptic_iterations;
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
