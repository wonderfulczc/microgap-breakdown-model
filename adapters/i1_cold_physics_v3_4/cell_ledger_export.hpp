#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>

namespace i1v34 {
using namespace streamer_rf;
using namespace streamer_rf::streamer;
inline void export_geometry(const Geometry2DView& geometry,const RoundedConeProfile2D& lo,
                            const RoundedConeProfile2D& hi,const std::filesystem::path& out){
 const auto& g=geometry.grid();
 std::ofstream meta(out/"geometry_native.json");meta<<std::setprecision(17);
 meta<<"{\"geometry_id\":\""<<geometry.geometry_id()<<"\",\"revision\":\""<<geometry.revision()
     <<"\",\"reference_plane_id\":\""<<geometry.reference_plane_id()<<"\",\"nr\":"<<g.nr()<<",\"nz\":"<<g.nz()
     <<",\"dr_m\":"<<g.dr()<<",\"dz_m\":"<<g.dz()<<",\"z0_center_m\":"<<g.z(0)<<",\"profiles\":[";
 bool first=true;for(const auto* p:{&lo,&hi}){
  if(!first)meta<<',';first=false;
  meta<<"{\"tip_z_m\":"<<p->tip_z_m<<",\"tip_radius_m\":"<<p->tip_radius_m
      <<",\"cone_half_angle_rad\":"<<p->cone_half_angle_rad<<",\"shank_radius_m\":"<<p->shank_radius_m
      <<",\"total_length_m\":"<<p->total_length_m<<",\"extension_direction\":"<<p->extension_direction<<'}';
 }meta<<"]}\n";
 std::ofstream mask(out/"geometry_mask.bin",std::ios::binary);
 std::ofstream faces(out/"geometry_metal_faces.csv");faces<<std::setprecision(17)
   <<"i,j,ni,nj,direction,conductor,r_m,z_m,normal_r,normal_z,area_m2\n";
 for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
  const auto material=geometry.classify(i,j);unsigned char m=static_cast<unsigned char>(material);mask.write(reinterpret_cast<char*>(&m),1);
  if(material!=ElectrodeCellType::Gas)continue;
  for(int d=0;d<4;++d){auto f=geometry.face(i,j,static_cast<FaceDirection2D>(d));if(f.role!=FaceRole2D::Electrode)continue;
   faces<<i<<','<<j<<','<<f.neighbor_i<<','<<f.neighbor_j<<','<<d<<','<<int(f.neighbor_material)<<','
        <<g.r(i)+.5*g.dr()*f.normal_r<<','<<g.z(j)+.5*g.dz()*f.normal_z<<','<<f.normal_r<<','<<f.normal_z<<','<<f.area_m2<<'\n';
  }
 }
 if(!meta||!faces||!mask)throw std::runtime_error("GEOMETRY_EXPORT_FAILED");
}
inline void export_accepted(const Geometry2DView& geometry,const AcceptedStepRecord& r,
                            std::uint64_t step,const std::filesystem::path& out){
 const auto& g=geometry.grid();auto stem="cell_step_"+std::to_string(step);
 std::ofstream cells(out/(stem+".csv"));cells<<std::setprecision(17)
  <<"cell_id,volume_m3,pre_time_s,post_time_s,dt_s,raw_ne_m3,raw_np_m3,raw_nn_m3,pre_ne_m3,pre_np_m3,pre_nn_m3,post_ne_m3,post_np_m3,post_nn_m3,impact_m3s,photo_m3s,attach2_m3s,attach3_m3s,ep_recomb_m3s,pn_recomb_m3s,source_e_m3s,source_p_m3s,source_n_m3s,reaction_defect_m3s,pre_floor_e_m3,pre_floor_p_m3,pre_floor_n_m3,clip_e_m3,clip_p_m3,clip_n_m3,post_floor_e_m3,post_floor_p_m3,post_floor_n_m3\n";
 std::size_t count=0;
 for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)if(geometry.classify(i,j)==ElectrodeCellType::Gas){
  auto k=std::size_t(j)*g.nr()+i;++count;
  cells<<k<<','<<g.cell_volume(i)<<','<<r.pre_time_s<<','<<r.post_time_s<<','<<r.dt_s;
  for(double v:{r.raw_pre_ne_m3[k],r.raw_pre_np_m3[k],r.raw_pre_nn_m3[k],r.pre.ne.values()[k],r.pre.np.values()[k],r.pre.nn.values()[k],
   r.post.ne.values()[k],r.post.np.values()[k],r.post.nn.values()[k],r.impact_m3s[k],r.photo_m3s[k],r.attachment2_m3s[k],r.attachment3_m3s[k],r.recombination_m3s[k],r.ion_recombination_m3s[k],
   r.reaction_e_m3s[k],r.reaction_p_m3s[k],r.reaction_n_m3s[k],r.reaction_charge_defect_m3s[k],
   r.pre_floor_ne_delta_m3[k],r.pre_floor_np_delta_m3[k],r.pre_floor_nn_delta_m3[k],r.clip_ne_delta_m3[k],r.clip_np_delta_m3[k],r.clip_nn_delta_m3[k],
   r.post_floor_ne_delta_m3[k],r.post_floor_np_delta_m3[k],r.post_floor_nn_delta_m3[k]})cells<<','<<v;
  cells<<'\n';
 }
 std::map<std::tuple<int,int,int>,const IonFaceFlux*> ions;
 for(const auto& f:r.positive_ion_faces)ions[{f.i,f.j,int(f.direction)}]=&f;
 std::ofstream faces(out/(stem+"_faces.csv"));faces<<std::setprecision(17)
  <<"i,j,ni,nj,direction,role,conductor,area_m2,gamma_e_net_m2s,gamma_p_m2s,gamma_secondary_m2s\n";
 for(const auto& f:r.electron_faces){
  auto it=ions.find({f.i,f.j,int(f.direction)});
  if(it==ions.end())throw std::runtime_error("CELL_LEDGER_ION_FACE_MISSING");
  faces<<f.i<<','<<f.j<<','<<f.neighbor_i<<','<<f.neighbor_j<<','<<int(f.direction)<<','<<int(f.role)<<','<<int(f.neighbor_material)<<','<<f.area_m2<<','<<f.gamma_e_m2s<<','<<it->second->gamma_p_out_m2s<<','<<it->second->gamma_secondary_in_m2s<<'\n';
 }
 std::ofstream index(out/(stem+"_index.json"));index<<std::setprecision(17)
  <<"{\"version\":\"I1_NATIVE_CELL_CHARGE_LEDGER_V1\",\"accepted_step_id\":"<<step
  <<",\"pre_time_s\":"<<r.pre_time_s<<",\"post_time_s\":"<<r.post_time_s
  <<",\"geometry_id\":\""<<r.geometry_id<<"\",\"mesh_id\":\""<<r.mesh_id
  <<"\",\"operator_id\":\""<<r.operator_id<<"\",\"reference_plane_id\":\""<<r.reference_plane_id
  <<"\",\"coverage\":\"ALL_GAS_CELLS_ALL_RECORDED_FACES\",\"gas_cells\":"<<count
  <<",\"faces\":"<<r.electron_faces.size()<<",\"time_layer\":\"RAW_PRE_AND_ACCEPTED_POST_WITH_PRE_RHS\",\"negative_ion_transport\":\"IMMOBILE_BY_MODEL\",\"secondary_in_net_face\":true}\n";
 if(!cells||!faces||!index)throw std::runtime_error("CELL_LEDGER_EXPORT_FAILED");
}
}
