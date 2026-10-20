#include "streamer_rf/streamer/ColdFullProvider.hpp"
#include "native_full_network.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#ifdef I1_V32_HOOKS
#include "native_hooks.hpp"
#endif
#ifdef I1_V34_HOOKS
#include "cell_ledger_export.hpp"
#endif
using namespace streamer_rf;
using namespace streamer_rf::streamer;
using Clock=std::chrono::steady_clock;
void dump_state(const StreamerState& s,const std::filesystem::path& p){
 std::ofstream f(p,std::ios::binary);for(auto* a:{&s.ne,&s.np,&s.nn,&s.rho,&s.phi,&s.er,&s.ez,&s.emag,&s.sph})f.write(reinterpret_cast<const char*>(a->values().data()),a->values().size()*sizeof(double));
 if(!f)throw std::runtime_error("STATE_WRITE_FAILED");
}
void dump_vector(const std::vector<double>& v,const std::filesystem::path& p){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));if(!f)throw std::runtime_error("VECTOR_WRITE_FAILED");}
double seed(const AxisymmetricGrid& g,int i,int j){
 const double s=12e-6,a=g.radial_face(i),b=g.radial_face(i+1),z0=g.z(j)-g.dz()/2,z1=z0+g.dz();
 const double zd=z0>=0?std::erfc(z0/s)-std::erfc(z1/s):z1<=0?std::erfc(-z1/s)-std::erfc(-z0/s):std::erf(z1/s)-std::erf(z0/s);
 return 1e16*std::exp(-a*a/(s*s))*(-std::expm1(-(b*b-a*a)/(s*s)))*s*s/(b*b-a*a)*s*std::sqrt(std::acos(-1.))*zd/(2*g.dz());
}
int main(int argc,char** argv){
 PetscInitialize(&argc,&argv,nullptr,nullptr);int code=0;auto start=Clock::now();int calls=0;
 try{
  if(argc!=6)throw std::invalid_argument("mode identity native_input output payload");
  std::string mode=argv[1],identity=argv[2];const std::filesystem::path out=argv[4],payload=argv[5];
  if(mode!="reference"&&mode!="resume")throw std::invalid_argument("INVALID_MODE");
  std::map<std::string,std::string> input;std::ifstream infile(argv[3]);std::string line;
  while(std::getline(infile,line)){auto k=line.find('=');if(k==std::string::npos)throw std::runtime_error("INPUT_FORMAT");input.emplace(line.substr(0,k),line.substr(k+1));}
  auto num=[&](const std::string& key){auto& s=input.at(key);std::size_t used;double x=std::stod(s,&used);if(used!=s.size()||!std::isfinite(x))throw std::runtime_error("INPUT_NONFINITE");return x;};
  AxisymmetricGrid g(200,400,.0005,-.0005,.0005);
  RoundedConeProfile2D lo{-35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,-1},hi{35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,1};
  auto geometry=Geometry2DView::facing_profiles(g,lo,hi,"B0_R02_A30_g70_Nominal","export-profile-1","LOCAL_HV_RETURN");
#ifdef I1_V34_HOOKS
  i1v34::export_geometry(geometry,lo,hi,out);
#endif
  ConstantVoltage voltage(500);StreamerConfig cfg;cfg.geometry_view=&geometry;cfg.voltage_waveform=&voltage;
  cfg.poisson_boundary={{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0}};
  cfg.elliptic={1e-12,1e-20,5000};cfg.photoionization=true;
  auto& c=cfg.cold_physics;c.mode=ColdPhysicsMode::CoupledColdV3;c.positive_ion_transport=true;c.local_charge_port_v2=true;
  c.ions={"CONSTANT_EFFECTIVE_POSITIVE_ION",input.at("ion_source"),"EFFECTIVE_SINGLE_POSITIVE_ION","CONSTANT_MU_D",input.at("qualification"),
    num("mobility_m2_Vs"),num("diffusion_m2_s"),num("min_EN_Td"),num("max_EN_Td"),num("min_temperature_K"),num("max_temperature_K"),num("min_pressure_Pa"),num("max_pressure_Pa")};
  c.surface={true,"ION_IMPACT_EFFECTIVE_YIELD",input.at("surface_source"),input.at("qualification"),num("gamma_i")};
  c.photon_boundary=PhotonMetalBoundary::OpaqueAbsorbingLarsen;c.excitation_source=input.at("excitation_source");
#ifdef I1_V32_HOOKS
  i1v32::configure(cfg,input);
#endif
  StreamerSolver solver(g,cfg,{true,true,false});ColdFullProvider provider(solver);AcceptedCheckpointLedger ledger;
  if(mode=="resume")ledger=provider.load_accepted_payload(payload,identity);
  else {ledger.capacitance_mode=input.at("capacitance_mode");for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)solver.state().ne(i,j)=solver.state().np(i,j)=seed(g,i,j);StreamerTransactionAccess::initial_floor(solver);}
  if(ledger.capacitance_mode!=input.at("capacitance_mode"))throw std::runtime_error("CHECKPOINT_CAPACITANCE_IDENTITY_MISMATCH");
  ScalarField2D port_weight(g);if(!ledger.port_weight.empty()){if(ledger.port_weight.size()!=g.size())throw std::runtime_error("PORT_GRID_MISMATCH");port_weight.values()=ledger.port_weight;}
  fp3full::ActiveNetwork net(geometry,cfg,ledger.port_weight.empty()?nullptr:&port_weight);
  ledger.port_weight=net.psi.values();
#ifdef I1_V32_HOOKS
  i1v32::Diagnostics v32diagnostics(net,cfg,out);
#endif
  const auto tracepath=out/"trace.csv";
  if(mode=="resume"&&std::filesystem::file_size(tracepath)!=ledger.output_offsets.at("trace.csv"))throw std::runtime_error("OUTPUT_OFFSET_MISMATCH");
  std::ofstream trace(tracepath,std::ios::app);trace<<std::setprecision(17);
  if(mode=="reference")trace<<"step,time_s,dt_s,Ne,Np,Nn,gas_charge_C,ion_H_N,ion_R_N,emit_H_N,emit_R_N,e_abs_H_N,e_abs_R_N,photo_N,photon_boundary_s1,photon_balance_s1,stable_charge_C,Wport_J,Wcond_J,DeltaU_J,Wouter_J,Wcorrection_J,Wreaction_J,energy_residual_J\n";
  const int steps=mode=="reference"?3:2;
  for(int k=0;k<steps;++k){
   if(std::chrono::duration<double>(Clock::now()-start).count()>450)throw std::runtime_error("COMPUTATIONAL_BUDGET_REACHED");
   const auto base=StreamerTransactionAccess::capture(solver);
   const double dt=std::min(2e-14,solver.timestep_limits().selected);
   if(dt<=0||solver.state().time+dt>1e-11)throw std::runtime_error("STEP_DOMAIN");
   ++calls;++ledger.native_calls;auto t=provider.evaluate_trial(dt);const auto& r=*t->record;auto m=net.budget(base.state,r,cfg);
   provider.commit(*t);ledger.accepted_step_id=provider.version();ledger.dt_s=dt;ledger.previous_time_s=base.state.time;
#ifdef I1_V34_HOOKS
   i1v34::export_accepted(geometry,r,provider.version(),out);
#endif
#ifdef I1_V32_HOOKS
   v32diagnostics.accepted(base.state,r,provider.version());
#endif
   long double ne=0,np=0,nn=0,ph=0,dqrho=0;
   for(std::size_t a=0;a<g.size();++a)if(geometry.classify(a%g.nr(),a/g.nr())==ElectrodeCellType::Gas){
    const double vol=g.cell_volume(a%g.nr());ne+=r.post.ne.values()[a]*vol;np+=r.post.np.values()[a]*vol;nn+=r.post.nn.values()[a]*vol;ph+=r.photo_m3s[a]*dt*vol;
    if(net.active[a])dqrho-=(long double)net.psi.values()[a]*(r.post.rho.values()[a]-base.state.rho.values()[a])*vol;
    for(auto item:{std::make_pair("floor_clip_e_N",r.post_floor_ne_delta_m3[a]+r.pre_floor_ne_delta_m3[a]+r.clip_ne_delta_m3[a]),
                   std::make_pair("floor_clip_p_N",r.post_floor_np_delta_m3[a]+r.pre_floor_np_delta_m3[a]+r.clip_np_delta_m3[a]),
                   std::make_pair("floor_clip_n_N",r.post_floor_nn_delta_m3[a]+r.pre_floor_nn_delta_m3[a]+r.clip_nn_delta_m3[a])})ledger.cumulative[item.first]+=item.second*vol;
   }
   const double dq=double(dqrho)+m.at("Deltaq_Hg_C");
   for(auto key:{"Wport_native_J","Wcond_face_J","DeltaU_E_J","Wouter_J","Wcorrection_J","Wreaction_J","correction_charge_C","reaction_charge_C","Deltaq_Hg_C"})ledger.cumulative[key]+=m.at(key);
   ledger.cumulative["stable_port_charge_C"]+=dq;ledger.cumulative["photo_N"]+=double(ph);
   ledger.cumulative["ion_H_N"]+=r.hv_collected_ions;ledger.cumulative["ion_R_N"]+=r.return_collected_ions;
   ledger.cumulative["emit_H_N"]+=r.hv_secondary_electrons;ledger.cumulative["emit_R_N"]+=r.return_secondary_electrons;
   ledger.last_electron_faces.clear();ledger.last_ion_faces.clear();
   for(auto f:r.electron_faces)for(double x:{double(f.i),double(f.j),double(f.neighbor_i),double(f.neighbor_j),double(f.direction),double(f.role),double(f.neighbor_material),f.gamma_e_m2s,f.area_m2})ledger.last_electron_faces.push_back(x);
   for(auto f:r.positive_ion_faces)for(double x:{double(f.i),double(f.j),double(f.face.neighbor_i),double(f.face.neighbor_j),double(f.direction),double(f.face.role),double(f.face.neighbor_material),f.gamma_p_out_m2s,f.gamma_secondary_in_m2s,f.face.area_m2})ledger.last_ion_faces.push_back(x);
   trace<<provider.version()<<','<<r.post_time_s<<','<<dt<<','<<ne<<','<<np<<','<<nn<<','<<1.602176634e-19L*(np-ne-nn)<<','<<r.hv_collected_ions<<','<<r.return_collected_ions<<','<<r.hv_secondary_electrons<<','<<r.return_secondary_electrons<<','<<r.hv_absorbed_electrons<<','<<r.return_absorbed_electrons<<','<<ph<<','<<r.photon_boundary_audit.metal_loss_s1+r.photon_boundary_audit.outer_loss_s1<<','<<r.photon_boundary_audit.balance_residual_s1<<','<<dq<<','<<m.at("Wport_native_J")<<','<<m.at("Wcond_face_J")<<','<<m.at("DeltaU_E_J")<<','<<m.at("Wouter_J")<<','<<m.at("Wcorrection_J")<<','<<m.at("Wreaction_J")<<','<<m.at("energy_residual_J")<<'\n';trace.flush();
   if(!trace)throw std::runtime_error("TRACE_WRITE_FAILED");ledger.output_offsets["trace.csv"]=std::filesystem::file_size(tracepath);
   if(mode=="reference"&&k==0){
    provider.save_accepted_payload(payload,identity,ledger);
    std::ofstream meta(out/"checkpoint_metadata.json");meta<<std::setprecision(17)<<"{\"accepted_step_id\":"<<provider.version()<<",\"time_s\":"<<r.post_time_s<<",\"dt_s\":"<<dt<<",\"previous_time_s\":"<<base.state.time<<",\"native_calls\":"<<ledger.native_calls<<",\"capacitance_mode\":\""<<ledger.capacitance_mode<<"\",\"output_offsets\":{\"trace.csv\":"<<ledger.output_offsets["trace.csv"]<<"}}\n";
   }else{
    auto stem="step_"+std::to_string(provider.version());dump_state(r.post,out/(stem+"_state.bin"));
    const auto saved=StreamerTransactionAccess::capture(solver);
    dump_vector({saved.head_position,saved.head_time,saved.q_hv,saved.q_return,saved.cap_native,saved.cap_conjugate,double(saved.has_history),double(saved.has_weight)},out/(stem+"_history.bin"));
    std::vector<double> warm;for(const auto& component:saved.warm.components)warm.insert(warm.end(),component.begin(),component.end());dump_vector(warm,out/(stem+"_warm.bin"));
    std::vector<double> weight;for(const auto* f:{&saved.psi,&saved.wr,&saved.wz})weight.insert(weight.end(),f->values().begin(),f->values().end());dump_vector(weight,out/(stem+"_weight.bin"));
    dump_vector(ledger.last_electron_faces,out/(stem+"_electron_faces.bin"));dump_vector(ledger.last_ion_faces,out/(stem+"_ion_faces.bin"));
    provider.save_accepted_payload(out/(stem+"_accepted.bin"),identity,ledger);
    std::ofstream budget(out/(stem+"_cumulative.csv"));budget<<std::setprecision(17)<<"key,value\n";for(auto [key,value]:ledger.cumulative)budget<<key<<','<<value<<'\n';
   }
  }
  std::ofstream summary(out/"execution_native.json");summary<<std::setprecision(17)<<"{\"native_step_calls_this_process\":"<<calls<<",\"accepted_step_id\":"<<provider.version()<<",\"time_s\":"<<solver.state().time<<",\"wall_s\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<",\"normal_exit\":true}\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';code=1;}
 PetscFinalize();return code;
}
