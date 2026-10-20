#include "streamer_rf/poisson.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sys/resource.h>

using namespace streamer_rf;
void require(bool condition,const char* message) {
  if(!condition) throw std::runtime_error(message);
}
void self_test() {
  const double pi=std::acos(-1.0);
  for(double scale:{0.5,1.0,2.0}) {
    RoundedConeProfile2D p{35e-6*scale,200e-6*scale,pi/12,1e-3*scale,.03*scale,1};
    p.validate();
    const double s=p.join_s_m(),h=1e-10*scale;
    require(std::abs(p.radius_at_s(s-h)-p.radius_at_s(s+h))<1e-9*scale,"arc-cone continuity");
    const auto n=p.side_normal(s),m=p.side_normal(s+h);
    require(std::hypot(n[0]-m[0],n[1]-m[1])<1e-6,"arc-cone normal continuity");
    require(std::abs(p.radius_at_s(p.cone_end_s_m())-p.shank_radius_m)<1e-15,"cone-cylinder continuity");
    require(!p.contains(0,p.tip_z_m-h)&&p.contains(0,p.tip_z_m+h),"tip orientation");
    require(!p.contains(0,p.tip_z_m+p.total_length_m+h),"back cap");
    RoundedConeProfile2D lower=p;lower.tip_z_m=-p.tip_z_m;lower.extension_direction=-1;
    AxisymmetricGrid g(24,80,.6e-3*scale,-.6e-3*scale,.6e-3*scale);
    auto v=Geometry2DView::facing_profiles(g,lower,p,"test","1","HV-return");
    require(std::abs(v.net_gap_m()-70e-6*scale)<1e-18,"gap");
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i) {
      const auto role=v.classify(i,j);
      if(role==ElectrodeCellType::Gas)for(auto dir:{FaceDirection2D::RadialMinus,FaceDirection2D::RadialPlus,FaceDirection2D::AxialMinus,FaceDirection2D::AxialPlus}) {
        const auto f=v.face(i,j,dir);require(f.area_m2>=0,"face area");
        require(std::abs(std::hypot(f.normal_r,f.normal_z)-1)<1e-14,"face normal");
      }
    }
    bool rejected=false;try{(void)Geometry2DView::facing_profiles(g,lower,p,"test","1","HV-return",false);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"3D must be rejected");
    rejected=false;try{(void)Geometry2DView::facing_profiles(g,lower,p,"test","1","HV-return",true,true);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"dielectric must be rejected");
    lower.tip_radius_m=150e-6*scale;lower.cone_half_angle_rad=pi/9;
    lower.shank_radius_m=1.1e-3*scale;lower.total_length_m=.025*scale;
    const auto asymmetric=Geometry2DView::facing_profiles(g,lower,p,"asymmetric","1","HV-return");
    require(asymmetric.classify(0,0)==ElectrodeCellType::GroundElectrode,"independent return profile");
    require(asymmetric.classify(0,g.nz()-1)==ElectrodeCellType::HighVoltageElectrode,"independent HV profile");
    p.total_length_m=p.cone_end_s_m()/2;
    rejected=false;try{p.validate();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"negative cylinder length accepted");
  }
  std::cout<<"PROFILE_TEST PASS scale, connections, normals, gap, masks, rejection\n";
}
int main(int argc,char** argv) {
  if(argc==2&&std::string(argv[1])=="--self-test") {self_test();return 0;}
  PetscInitialize(&argc,&argv,nullptr,nullptr);
  try {
    if(argc==2 && std::string(argv[1])=="--boundary-test") {
      AxisymmetricGrid g(12,32,1e-3,-1e-3,1e-3);
      ScalarField2D zero(g),phi(g);
      PoissonBoundaryConfig bc{{BoundaryKind::Neumann,0},{BoundaryKind::Dirichlet,0},{BoundaryKind::Dirichlet,1}};
      solve_potential(zero,0,bc,phi,{1e-12,1e-16,3000});
      double error=0;
      for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)
        error=std::max(error,std::abs(phi(i,j)-static_cast<double>(j)/(g.nz()-1)));
      require(error<1e-9,"native Neumann/Dirichlet plate analytic error");
      std::cout<<"BOUNDARY_TEST PASS max_psi_error="<<error<<'\n';
      PetscFinalize();return 0;
    }
    require(argc>=4,"usage: input.txt output.bin output.json [PETSc options]");
    std::ifstream in(argv[1]); int nr,nz;double rmax,zmin,zmax;std::string kind;
    in>>nr>>nz>>rmax>>zmin>>zmax>>kind;
    RoundedConeProfile2D lo,hi;
    for(auto* p:{&lo,&hi})in>>p->tip_z_m>>p->tip_radius_m>>p->cone_half_angle_rad>>p->shank_radius_m>>p->total_length_m>>p->extension_direction;
    require(bool(in),"invalid static input");
    AxisymmetricGrid g(nr,nz,rmax,zmin,zmax);
    const auto view=Geometry2DView::facing_profiles(g,lo,hi,"B0-profile-from-exports","1","HV-return");
    require(kind=="NEUMANN_ZERO"||kind=="DIRICHLET_ZERO","unsupported boundary");
    const auto bk=kind=="NEUMANN_ZERO"?BoundaryKind::Neumann:BoundaryKind::Dirichlet;
    PoissonBoundaryConfig bc{{bk,0},{bk,0},{bk,0}};
    ScalarField2D rho(g),phi(g);
    const auto start=std::chrono::steady_clock::now();
    solve_potential_with_electrodes(rho,view,1.,bc,phi,{1e-11,1e-15,3000});
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream bin(argv[2],std::ios::binary);
    for(int j=0;j<nz;++j)for(int i=0;i<nr;++i) {
      // Same cell-centered derivative convention as the native solver.
      const double er=i==0?0:-(phi(std::min(i+1,nr-1),j)-phi(i-1,j))/((i==nr-1?1:2)*g.dr());
      const double ez=-(phi(i,std::min(j+1,nz-1))-phi(i,std::max(j-1,0)))/((j==0||j==nz-1?1:2)*g.dz());
      const double row[]{phi(i,j),er,ez,static_cast<double>(view.classify(i,j))};
      bin.write(reinterpret_cast<const char*>(row),sizeof(row));
    }
    struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::ofstream meta(argv[3]);meta.precision(17);
    meta<<"{\"solver_seconds\":"<<seconds<<",\"peak_rss_kib\":"<<usage.ru_maxrss
        <<",\"iterations\":"<<last_poisson_iterations()<<",\"cells\":"<<g.size()
        <<",\"binary_layout\":\"native float64 rows: phi,Er,Ez,material; radial index fastest\"}\n";
    std::cout<<"STATIC_SOLVE "<<nr<<'x'<<nz<<" seconds="<<seconds<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';PetscFinalize();return 1;}
  PetscFinalize();return 0;
}
