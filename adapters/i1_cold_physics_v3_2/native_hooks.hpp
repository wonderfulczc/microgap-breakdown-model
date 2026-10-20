#pragma once
#include "native_full_stable.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace i1v32 {
using namespace streamer_rf;
using namespace streamer_rf::streamer;
inline void configure(StreamerConfig& cfg,const std::map<std::string,std::string>& input){
  auto& m=cfg.cold_physics.ions;auto& t=m.table;
  m.model_id="SOURCE_BOUND_EN_TABLE";m.approximation="EFFECTIVE_SPECIES_APPROXIMATION";
  t.species_id=m.species_id;t.gas_composition=input.at("table_gas");
  t.source_doi=input.at("table_doi");t.database_id=input.at("table_database");
  t.file_sha256=input.at("table_sha256");t.evidence_class=input.at("table_class");
  t.mobility_unit=input.at("table_mu_unit");t.reference_density_m3=std::stod(input.at("table_Nref"));
  std::ifstream f(input.at("table_file"));std::string line;std::getline(f,line);
  if(line!="EN_Td,mobility,D_parallel_m2_s,D_perpendicular_m2_s")throw std::runtime_error("TABLE_HEADER_MISMATCH");
  while(std::getline(f,line)){
    std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double x,mu,dp,dt;std::string extra;
    if(!(row>>x>>mu>>dp>>dt)||(row>>extra))throw std::runtime_error("TABLE_ROW_INVALID");
    t.EN_Td.push_back(x);t.mobility.push_back(mu);t.D_parallel.push_back(dp);t.D_perpendicular.push_back(dt);
  }
  auto& s=cfg.cold_physics.surface;s.source_contract_required=true;
  s.incident_species_id=m.species_id;s.gas_composition=t.gas_composition;
  s.cathode_material=input.at("surface_material");s.surface_state=input.at("surface_state");
  s.incident_energy_assumption=input.at("surface_energy_assumption");s.source_sha256=input.at("surface_sha256");
  s.evidence_class=input.at("surface_class");s.uncertainty=input.at("surface_uncertainty");
  if(s.evidence_class!="C")throw std::runtime_error("V32_SMOKE_REQUIRES_SYNTHETIC_MODE");
  m.validate(0,300,101325);s.validate_source();
}
class Diagnostics {
  const fp3full::ActiveNetwork& net;const StreamerConfig& cfg;std::filesystem::path out;
  int calls=0;fp3full::StableNetwork stable;
 public:
  Diagnostics(const fp3full::ActiveNetwork& n,const StreamerConfig& c,const std::filesystem::path& o):net(n),cfg(c),out(o),stable(n,c,calls){}
  void accepted(const StreamerState& base,const AcceptedStepRecord& r,std::uint64_t step){
    auto m=stable.budget(base,r,net.budget(base,r,cfg));
    m["step"]=step;m["electron_H_collected_N"]=r.hv_absorbed_electrons;m["electron_Return_collected_N"]=r.return_absorbed_electrons;
    m["ion_H_collected_N"]=r.hv_collected_ions;m["ion_Return_collected_N"]=r.return_collected_ions;
    m["secondary_H_N"]=r.hv_secondary_electrons;m["secondary_Return_N"]=r.return_secondary_electrons;
    double lo=INFINITY,hi=0;for(std::size_t k=0;k<net.g.size();++k)if(net.gas(k%net.g.nr(),k/net.g.nr())){
      double EN=r.pre.emag.values()[k]/cfg.neutral_density*1e21;lo=std::min(lo,EN);hi=std::max(hi,EN);
    }
    m["pre_EN_min_Td"]=lo;m["pre_EN_max_Td"]=hi;
    std::ofstream metrics(out/("weak_step_"+std::to_string(step)+".csv"));metrics<<std::setprecision(17)<<"quantity,value\n";
    for(auto [k,v]:m)metrics<<k<<','<<v<<'\n';
    if(step!=3)return;
    export_source(r);
  }
  void export_source(const AcceptedStepRecord& r){
    const auto& g=net.g;const double qe=1.602176634e-19;
    std::vector<std::array<double,4>> je(g.size()),jp(g.size());
    std::ofstream faces(out/"species_resolved_current_faces.csv");faces<<std::setprecision(17);
    faces<<"i,j,ni,nj,direction,role,conductor,area_m2,gamma_e_net_m2s,gamma_p_m2s,gamma_secondary_m2s,Je_out_A_m2,Jp_out_A_m2,Jtotal_out_A_m2,time_pre_s,time_post_s\n";
    std::map<std::tuple<int,int,int>,const IonFaceFlux*> ions;
    for(const auto& f:r.positive_ion_faces)ions[{f.i,f.j,int(f.direction)}]=&f;
    for(const auto& f:r.electron_faces){
      int d=int(f.direction),k=f.j*g.nr()+f.i;double gp=0,ge=0;
      auto found=ions.find({f.i,f.j,d});
      if(found!=ions.end()){gp=found->second->gamma_p_out_m2s;ge=found->second->gamma_secondary_in_m2s;}
      else if(f.role!=FaceRole2D::Axis)throw std::runtime_error("SPECIES_FACE_IDENTITY_MISMATCH");
      const double e=-qe*f.gamma_e_m2s,p=qe*gp;
      je[k][d]=e;jp[k][d]=p;
      faces<<f.i<<','<<f.j<<','<<f.neighbor_i<<','<<f.neighbor_j<<','<<d<<','<<int(f.role)<<','<<int(f.neighbor_material)<<','<<f.area_m2<<','<<f.gamma_e_m2s<<','<<gp<<','<<ge<<','<<e<<','<<p<<','<<e+p<<','<<r.pre_time_s<<','<<r.post_time_s<<'\n';
    }
    std::ofstream cells(out/"species_center_source.csv");cells<<std::setprecision(17);
    cells<<"cell_id,r_m,z_m,dr_m,dz_m,volume_m3,rho_C_m3,Je_r_A_m2,Je_z_A_m2,Jp_r_A_m2,Jp_z_A_m2,time_pre_s,time_post_s\n";
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)if(net.gas(i,j)){
      int k=j*g.nr()+i;
      // Normal components are outward; center vector is the mean of the two oriented face values.
      cells<<k<<','<<g.r(i)<<','<<g.z(j)<<','<<g.dr()<<','<<g.dz()<<','<<g.cell_volume(i)<<','<<r.pre.rho(i,j)<<','
        <<.5*(je[k][1]-je[k][0])<<','<<.5*(je[k][3]-je[k][2])<<','<<.5*(jp[k][1]-jp[k][0])<<','<<.5*(jp[k][3]-jp[k][2])<<','<<r.pre_time_s<<','<<r.post_time_s<<'\n';
    }
    if(!faces||!cells)throw std::runtime_error("STAGE_F_EXPORT_FAILED");
  }
};
}
