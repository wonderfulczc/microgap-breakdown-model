#include "streamer_rf/streamer/ColdFullProvider.hpp"
#include "native_full_stable.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

using namespace streamer_rf;
using namespace streamer_rf::streamer;
using Clock=std::chrono::steady_clock;
constexpr double echarge=1.602176634e-19;
double cell_average(const AxisymmetricGrid& g,int i,int j){
  const double s=12e-6,a=g.radial_face(i),b=g.radial_face(i+1),z0=g.z(j)-g.dz()/2,z1=z0+g.dz();
  const double zd=z0>=0?std::erfc(z0/s)-std::erfc(z1/s):z1<=0?std::erfc(-z1/s)-std::erfc(-z0/s):std::erf(z1/s)-std::erf(z0/s);
  return 1e16*std::exp(-a*a/(s*s))*(-std::expm1(-(b*b-a*a)/(s*s)))*s*s/(b*b-a*a)*s*std::sqrt(std::acos(-1.))*zd/(2*g.dz());
}
void state_file(const StreamerState& s,const std::filesystem::path& p){
  std::ofstream f(p,std::ios::binary);
  for(const auto* v:{&s.ne,&s.np,&s.nn,&s.rho,&s.phi,&s.er,&s.ez,&s.emag,&s.sph})
    f.write(reinterpret_cast<const char*>(v->values().data()),v->values().size()*sizeof(double));
}
int main(int argc,char** argv){
  PetscInitialize(&argc,&argv,nullptr,nullptr);int code=0;
  try{
    if(argc!=3)throw std::invalid_argument("usage: i1_cold_scenario validated_native_input fresh_output_directory");
    std::map<std::string,std::string> values;std::ifstream input(argv[1]);std::string line;
    while(std::getline(input,line)){const auto pos=line.find('=');if(pos==std::string::npos)throw std::invalid_argument("INVALID_SCENARIO_LINE");values.emplace(line.substr(0,pos),line.substr(pos+1));}
    auto num=[&](const std::string& key){std::size_t used;const auto& s=values.at(key);double v=std::stod(s,&used);if(used!=s.size()||!std::isfinite(v))throw std::invalid_argument("NONFINITE_SCENARIO_PARAMETER");return v;};
    const double end=num("end_time_s"),cap=num("max_step_s"),budget=num("wall_budget_s");
    const int max_calls=int(num("max_native_steps"));
    if(end<=0||end>1e-11||cap<=0||cap>2e-14||budget<=0||budget>2400||max_calls<1||max_calls>1500)throw std::invalid_argument("UNAUTHORIZED_RESOURCE_LIMIT");
    const std::filesystem::path out=argv[2];if(std::filesystem::exists(out/"native_species_charge_conservation.csv"))throw std::runtime_error("REFUSE_OVERWRITE");
    auto start=Clock::now();int calls=0,accepted=0;
    AxisymmetricGrid g(200,400,.0005,-.0005,.0005);
    RoundedConeProfile2D lo{-35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,-1},hi{35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,1};
    auto geometry=Geometry2DView::facing_profiles(g,lo,hi,"B0_R02_A30_g70_Nominal","export-profile-1","LOCAL_HV_RETURN");
    ConstantVoltage voltage(500);StreamerConfig cfg;cfg.geometry_view=&geometry;cfg.voltage_waveform=&voltage;
    cfg.poisson_boundary={{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0}};
    cfg.elliptic={1e-12,1e-20,5000};cfg.photoionization=true;
    auto& c=cfg.cold_physics;c.mode=ColdPhysicsMode::CoupledColdV3;c.positive_ion_transport=true;c.local_charge_port_v2=true;
    c.ions={"CONSTANT_EFFECTIVE_POSITIVE_ION",values.at("ion_source"),"EFFECTIVE_SINGLE_POSITIVE_ION","CONSTANT_MU_D",values.at("qualification"),
      num("mobility_m2_Vs"),num("diffusion_m2_s"),num("min_EN_Td"),num("max_EN_Td"),num("min_temperature_K"),num("max_temperature_K"),num("min_pressure_Pa"),num("max_pressure_Pa")};
    c.surface={true,"ION_IMPACT_EFFECTIVE_YIELD",values.at("surface_source"),values.at("qualification"),num("gamma_i")};
    c.photon_boundary=PhotonMetalBoundary::OpaqueAbsorbingLarsen;c.excitation_source=values.at("excitation_source");
    StreamerSolver solver(g,cfg,{true,true,false});
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)solver.state().ne(i,j)=solver.state().np(i,j)=cell_average(g,i,j);
    StreamerTransactionAccess::initial_floor(solver);
    state_file(solver.state(),out/"state_initial.bin");
    fp3full::ActiveNetwork net(geometry,cfg);
    ColdFullProvider provider(solver);
    std::ofstream ledger(out/"native_species_charge_conservation.csv"),wall(out/"secondary_emission_boundary_audit.csv"),
      ions(out/"ion_transport_and_wall_flux_budget.csv"),photo(out/"b0_sp3_internal_boundary_audit.csv"),port(out/"local_charge_and_work_budget.csv"),transaction(out/"transaction_audit.csv");
    for(auto* f:{&ledger,&wall,&ions,&photo,&port,&transaction})*f<<std::setprecision(17);
    ledger<<"pre_s,post_s,dt_s,Ne_pre,Ne_post,Np_pre,Np_post,Nn_pre,Nn_post,impact_N,photo_N,attach2_N,attach3_N,recomb_ep_N,recomb_pn_N,correction_e_N,correction_p_N,correction_n_N,residual_e_N,residual_p_N,residual_n_N,residual_charge_C,GR_pre_S,GR_post_S\n";
    wall<<"pre_s,post_s,i,j,direction,conductor,area_m2,ion_incident_m2s,electron_absorbed_m2s,secondary_emitted_m2s,net_electron_out_m2s,ion_collected_C,emission_electrode_C,net_gas_charge_out_C\n";
    ions<<"pre_s,post_s,ion_H_N,ion_Return_N,ion_outer_N,secondary_H_N,secondary_Return_N,I_Hg_A,I_Return_g_A,ion_positivity_dt_s\n";
    photo<<"pre_s,post_s,boundary_id,metal_faces,source_s1,photo_created_s1,metal_loss_s1,outer_loss_s1,balance_residual_s1,linear_relative_residual,iterations\n";
    port<<"pre_s,post_s,C_local_F,q_A_C,q_B_native_links_C,q_B_active_V1_C,q_B_stable_V2_C,Deltaq_Hg_C,stable_Deltaq_port_C,stable_charge_residual_C,transport_abs_C,Wcond_native_faces_J,Wport_J,DeltaU_field_J,native_energy_residual_J,TOTAL_minus_SPLIT_C\n";
    transaction<<"time_s,check,result\n";
    std::string reason="END_TIME_REACHED";double max_species=0,max_charge=0,max_photo=0,max_port=0;
    auto allowed=[&](){return calls<max_calls&&std::chrono::duration<double>(Clock::now()-start).count()<budget;};
    while(solver.state().time<end){
      if(!allowed()){reason="COMPUTATIONAL_BUDGET_REACHED";break;}
      double dt=std::min({cap,solver.timestep_limits().selected,end-solver.state().time});
      if(dt<=128*std::numeric_limits<double>::epsilon()*end){reason="END_TIME_WITHIN_FLOAT_ALIGNMENT";break;}
      const auto base=StreamerTransactionAccess::capture(solver);
      if(accepted==3&&calls+4<max_calls){
        ++calls;auto a=provider.evaluate_trial(dt);
        ++calls;auto b=provider.evaluate_trial(dt*.5);provider.reject(*b);
        ++calls;auto again=provider.evaluate_trial(dt);
        bool same=a->record->post.ne.values()==again->record->post.ne.values()&&a->record->post.np.values()==again->record->post.np.values();
        provider.reject(*a);provider.reject(*again);
        bool fail=false;try{++calls;provider.evaluate_trial(dt,true);}catch(const std::runtime_error& error){fail=std::string(error.what())=="CONTROLLED_FAILURE_AFTER_NATIVE_FIELDS";}
        same&=native_hash(StreamerTransactionAccess::capture(solver))==native_hash(base)&&native_hash(StreamerTransactionAccess::capture(solver),true)==native_hash(base,true);
        transaction<<base.state.time<<",repeat_reject_controlled_failure,"<<(same&&fail?"PASS":"FAIL")<<'\n';
        transaction.flush();
        if(!same||!fail)throw std::runtime_error("V3_NATIVE_TRANSACTION_MISMATCH");
      }
      ++calls;auto trial=provider.evaluate_trial(dt);const auto& r=*trial->record;
      auto raw=net.budget(base.state,r,cfg);
      provider.commit(*trial);++accepted;
      long double ne0=0,ne1=0,np0=0,np1=0,nn0=0,nn1=0,impact=0,ph=0,a2=0,a3=0,ep=0,pn=0,ce=0,cp=0,cn=0,dqrho=0;
      for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
        if(geometry.classify(i,j)!=ElectrodeCellType::Gas)continue;
        auto k=j*g.nr()+i;const double v=g.cell_volume(i);
        ne0+=r.pre.ne(i,j)*v;ne1+=r.post.ne(i,j)*v;np0+=r.pre.np(i,j)*v;np1+=r.post.np(i,j)*v;nn0+=r.pre.nn(i,j)*v;nn1+=r.post.nn(i,j)*v;
        impact+=dt*r.impact_m3s[k]*v;ph+=dt*r.photo_m3s[k]*v;a2+=dt*r.attachment2_m3s[k]*v;a3+=dt*r.attachment3_m3s[k]*v;ep+=dt*r.recombination_m3s[k]*v;
        pn+=dt*beta_np(cfg.temperature)*r.pre.nn(i,j)*r.pre.np(i,j)*v;
        ce+=(r.clip_ne_delta_m3[k]+r.post_floor_ne_delta_m3[k])*v;cp+=(r.clip_np_delta_m3[k]+r.post_floor_np_delta_m3[k])*v;cn+=(r.clip_nn_delta_m3[k]+r.post_floor_nn_delta_m3[k])*v;
        dqrho-=(long double)net.psi(i,j)*((long double)r.post.rho(i,j)-base.state.rho(i,j))*v*(net.active[k]?1:0);
      }
      const long double re=ne1-ne0-impact-ph+a2+a3+ep+r.hv_absorbed_electrons+r.return_absorbed_electrons+r.outer_lost_electrons-r.hv_secondary_electrons-r.return_secondary_electrons-ce;
      const long double rp=np1-np0-impact-ph+ep+pn+r.hv_collected_ions+r.return_collected_ions+r.outer_lost_ions-cp;
      const long double rn=nn1-nn0-a2-a3+pn-cn;
      const double rq=double(echarge*(rp-re-rn));max_species=std::max({max_species,double(std::abs(re)),double(std::abs(rp)),double(std::abs(rn))});max_charge=std::max(max_charge,std::abs(rq));
      ledger<<r.pre_time_s<<','<<r.post_time_s<<','<<dt<<','<<ne0<<','<<ne1<<','<<np0<<','<<np1<<','<<nn0<<','<<nn1<<','<<impact<<','<<ph<<','<<a2<<','<<a3<<','<<ep<<','<<pn<<','<<ce<<','<<cp<<','<<cn<<','<<re<<','<<rp<<','<<rn<<','<<rq<<','<<r.G_R_pre_S<<','<<r.G_R_post_S<<'\n';
      std::vector<double> net_e(4*g.size());for(auto f:r.electron_faces)net_e[4*(f.j*g.nr()+f.i)+int(f.direction)]=f.gamma_e_m2s;
      for(auto f:r.positive_ion_faces)if(f.face.role==FaceRole2D::Electrode){
        double ge=net_e[4*(f.j*g.nr()+f.i)+int(f.direction)];
        wall<<r.pre_time_s<<','<<r.post_time_s<<','<<f.i<<','<<f.j<<','<<int(f.direction)<<','<<int(f.face.neighbor_material)<<','<<f.face.area_m2<<','<<f.gamma_p_out_m2s<<','<<ge+f.gamma_secondary_in_m2s<<','<<f.gamma_secondary_in_m2s<<','<<ge<<','<<echarge*dt*f.gamma_p_out_m2s*f.face.area_m2<<','<<echarge*dt*f.gamma_secondary_in_m2s*f.face.area_m2<<','<<echarge*dt*(f.gamma_p_out_m2s-ge)*f.face.area_m2<<'\n';
      }
      ions<<r.pre_time_s<<','<<r.post_time_s<<','<<r.hv_collected_ions<<','<<r.return_collected_ions<<','<<r.outer_lost_ions<<','<<r.hv_secondary_electrons<<','<<r.return_secondary_electrons<<','<<r.net_I_Hg_A<<','<<r.net_I_Return_g_A<<','<<r.ion_positivity_dt_s<<'\n';
      const auto& pa=r.photon_boundary_audit;
      max_photo=std::max(max_photo,pa.effective_source_s1>0?std::abs(pa.balance_residual_s1)/pa.effective_source_s1:0);
      photo<<r.pre_time_s<<','<<r.post_time_s<<','<<pa.boundary_id<<','<<pa.metal_faces<<','<<pa.effective_source_s1<<','<<pa.photo_created_s1<<','<<pa.metal_loss_s1<<','<<pa.outer_loss_s1<<','<<pa.balance_residual_s1<<','<<pa.max_relative_linear_residual<<','<<pa.iterations<<'\n';
      const long double qpl=dqrho+raw.at("Deltaq_Hg_C"),qcap=(long double)net.capacitance*(r.post_voltage_V-r.pre_voltage_V),dq=qcap+qpl;
      const double stable_res=double(qpl-raw.at("Deltaq_w_C")-raw.at("Deltaq_outer_C")+raw.at("DeltaQ_psi_correction_C")+raw.at("DeltaQ_psi_reaction_C"));
      max_port=std::max(max_port,std::abs(stable_res));
      long double qp=0;for(std::size_t k=0;k<g.size();++k)if(net.active[k])qp-=(long double)net.psi.values()[k]*r.post.rho.values()[k]*g.cell_volume(k%g.nr());
      port<<r.pre_time_s<<','<<r.post_time_s<<','<<net.capacitance<<','<<r.q_native_hv_C<<','<<r.q_conjugate_hv_C<<','<<raw.at("q_B_active_end_C")<<','<<net.capacitance*r.post_voltage_V+qp<<','<<raw.at("Deltaq_Hg_C")<<','<<dq<<','<<stable_res<<','<<raw.at("unsigned_native_transport_C")<<','<<raw.at("Wcond_face_J")<<','<<raw.at("Wport_native_J")<<','<<raw.at("DeltaU_E_J")<<','<<raw.at("energy_residual_J")<<','<<double(dq-(qcap+qpl))<<'\n';
      for(auto* f:{&ledger,&wall,&ions,&photo,&port,&transaction}){f->flush();if(!*f)throw std::runtime_error("ACCEPTED_OUTPUT_WRITE_FAILED");}
      // Publish progress only after every interval ledger has been flushed.
      {std::ofstream progress(out/"accepted_progress.json.tmp");progress<<std::setprecision(17)
        <<"{\"accepted_steps\":"<<accepted<<",\"native_step_calls\":"<<calls
        <<",\"accepted_time_s\":"<<solver.state().time<<",\"restartable_full_snapshot\":false}\n";
        progress.close();if(!progress)throw std::runtime_error("PROGRESS_WRITE_FAILED");
        std::filesystem::rename(out/"accepted_progress.json.tmp",out/"accepted_progress.json");}
      if(accepted%25==0){state_file(solver.state(),out/"state_latest.bin.tmp");std::filesystem::rename(out/"state_latest.bin.tmp",out/"state_latest.bin");std::cout<<"accepted "<<accepted<<" t_ps "<<solver.state().time/1e-12<<std::endl;}
    }
    state_file(solver.state(),out/"state_final.bin");
    std::ofstream summary(out/"b0_cold_full_smoke_summary.json");summary<<std::setprecision(17)
      <<"{\n\"scenario\":\"COLD_FULL_B0\",\n\"physics_qualification\":\"SYNTHETIC_OR_CONDITIONAL_PHYSICS\",\n\"stop_reason\":\""<<reason<<"\",\n\"native_step_calls\":"<<calls<<",\n\"accepted_steps\":"<<accepted<<",\n\"final_time_s\":"<<solver.state().time
      <<",\n\"wall_s\":"<<std::chrono::duration<double>(Clock::now()-start).count()<<",\n\"max_species_residual_particles\":"<<max_species<<",\n\"max_charge_residual_C\":"<<max_charge<<",\n\"max_photon_balance_relative\":"<<max_photo<<",\n\"max_stable_port_residual_C\":"<<max_port
      <<",\n\"surface_emission_energy\":null,\n\"surface_emission_energy_status\":\"UNAVAILABLE\",\n\"gas_heat\":null,\n\"global_physical_port\":false,\n\"RF_validated\":false\n}\n";
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';code=1;}
  PetscFinalize();return code;
}
