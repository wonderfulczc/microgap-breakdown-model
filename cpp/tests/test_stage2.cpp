#include "streamer_rf/transport.hpp"
#include "streamer_rf/sp3.hpp"
#include "streamer_rf/transport_table.hpp"
#include <petscsys.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace streamer_rf;
namespace {int count=0;void check(bool q,const char*n){++count;if(!q)throw std::runtime_error(n);std::cout<<"ok "<<count<<" - "<<n<<'\n';}}
int main(int argc,char**argv){PetscInitialize(&argc,&argv,nullptr,nullptr);try{
 check(std::abs(bernoulli(0)-1)<1e-15,"Bernoulli zero");check(std::abs(bernoulli(1e-7)-(1-5e-8))<1e-14,"Bernoulli series");
 check(std::abs(sg_flux(2,1,0,3,.5)-6)<1e-13,"SG pure diffusion");check(std::abs(sg_flux(2,2,0,3,.5))<1e-14,"SG uniform zero field");
 check(sg_flux(2,1,100,1,.1)>0,"SG upwind direction");check(std::abs(sg_flux(2,1,3,1,.1)+sg_flux(1,2,-3,1,.1))<1e-12,"SG antisymmetry");
 IsgDiagnostics d;double f=isg0_flux(2,1,0,100,1,.1,.01,1,&d);check(std::isfinite(f)&&d.branch=="isg0","ISG-0 regular branch");
 check(std::abs(isg0_flux(2,1,3,3,1,.1,.01,1,&d)-sg_flux(2,1,3,1,.1))<1e-14,"ISG fallback");
 double z=isg0_zero_width_flux(2,1,3,1,.1,1);check(std::isfinite(z),"ISG zero-width");
 double fs=isg0_flux(2e9,1e9,0,100,1,.1,.01,1e9,nullptr);check(std::abs(fs/f-1e9)<1e-5,"n_ref scaling");
 AxisymmetricGrid g(4,8,.04,-.04,.04);check(std::abs(g.cell_volume(0)-3.14159265358979323846*.01*.01*.01)<1e-15,"exact axis cell volume");
 check(g.r(0)==.005&&g.radial_face(0)==0,"axisymmetric coordinates");
 auto c=sp3_constants({});check(std::abs(c.groups[0].A_original-.0067)<1e-15,"SP3 Table 3 coefficients");check(c.kappa1_sq>0&&c.kappa2_sq>c.kappa1_sq,"SP3 kappa formula");
 check(std::isfinite(complete_elliptic_k(.5)),"open boundary elliptic kernel uses parameter m");check(std::abs(ring_potential_on_axis(1e-12,.01,0)-.89875517923)<1e-8,"ring charge analytic potential");
 TransportTable tab({1,2,3},{2,4,6});check(tab.interpolate(1.5)==3,"transport table interpolation");
 PoissonBoundaryConfig bc;ScalarField2D rho(g),phi(g),src(g),sph(g);solve_potential(rho,0,bc,phi);check(std::isfinite(phi(1,1)),"Poisson matrix and solve");
 rho(1,4)=1e-9;auto ob1=open_charge_boundary(rho,{RingQuadratureMode::CellCenterRing,0});double before=ob1(g.nr()-1,4);rho(1,4)*=2;auto ob2=open_charge_boundary(rho,{RingQuadratureMode::CellCenterRing,0});check(std::abs(ob2(g.nr()-1,4)/before-2)<1e-13,"OpenCharge updates with rho");
 const double q=rho(1,4)*g.cell_volume(1),ref=ring_potential(q,g.r(1),g.r(g.nr()-1),g.z(4)-g.z(4));check(std::abs(ob2(g.nr()-1,4)-ref)/std::abs(ref)<1e-13,"OpenCharge cell-center ring kernel");
 ScalarField2D rhs(g),partner0(g),partner1(g),robin0(g),robin1(g);for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)if(i==g.nr()-1||j==0||j==g.nz()-1)partner1(i,j)=1;solve_axisymmetric_robin(rhs,robin0,100,20,10,partner0);solve_axisymmetric_robin(rhs,robin1,100,20,10,partner1);double cross_response=0;for(std::size_t k=0;k<g.size();++k)cross_response=std::max(cross_response,std::abs(robin1.values()[k]-robin0.values()[k]));check(cross_response>1e-8,"Robin cross-field contribution changes solution");
 src(1,4)=1e15;std::array<Sp3GroupDiagnostics,3> sd;solve_photoionization(src,{},bc,sph,{},Sp3BoundaryOptions{},&sd);check(std::isfinite(sph(1,4)),"six SP3 Helmholtz solves");check(sd[0].converged&&sd[0].final_residual<1e-8,"coupled Robin fixed-point convergence");check(sd[0].residual_history.size()>1&&sd[0].residual_history.back()<sd[0].residual_history.front(),"Robin residual decreases");
 ScalarField2D exact(g),mrhs(g),mnum(g);for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){exact(i,j)=g.r(i)*g.r(i)+g.z(j)*g.z(j);mrhs(i,j)=6;}solve_axisymmetric_dirichlet(mrhs,exact,mnum,0);double merr=0;for(std::size_t k=0;k<g.size();++k)merr=std::max(merr,std::abs(mnum.values()[k]-exact.values()[k]));check(merr<1e-10,"axisymmetric manufactured quadratic solution");
 std::cout<<"passed "<<count<<" C++ checks\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';PetscFinalize();return 1;}PetscFinalize();return 0;}
