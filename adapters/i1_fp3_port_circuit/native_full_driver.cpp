#include "native_full_network.hpp"
#include "native_full_stable.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>

using namespace fp3full;
using Clock=std::chrono::steady_clock;

struct LinearInput final:VoltageWaveform {
  double t0=0,t1=0,v0=500,v1=500;
  double value(double t)const override{return t<=t0?v0:t>=t1?v1:v0+(v1-v0)*(t-t0)/(t1-t0);}
};
struct Buffer final:AcceptedStepObserver {
  mutable std::unique_ptr<AcceptedStepRecord> last;
  void on_accepted_step(const AcceptedStepRecord&r)const override{last=std::make_unique<AcceptedStepRecord>(r);}
};
struct Trial {
  int id;
  std::string base;
  NativeFullSnapshot frame;
  LinearInput wave;
  Metrics m;
  std::vector<Metrics> rows;
  std::vector<std::shared_ptr<AcceptedStepRecord>> records;
  explicit Trial(const AxisymmetricGrid&g):frame(g){}
};
void emit_metrics(const Metrics&m){bool first=true;for(auto [k,v]:m){if(!first)std::cout<<',';first=false;std::cout<<'"'<<k<<"\":";if(std::isfinite(v))std::cout<<v;else std::cout<<"null";}}
double average(const AxisymmetricGrid&g,int i,int j){
  const double s=12e-6,a=g.radial_face(i),b=g.radial_face(i+1),z0=g.z(j)-g.dz()/2,z1=g.z(j)+g.dz()/2;
  double zdiff=z0>=0?std::erfc(z0/s)-std::erfc(z1/s):z1<=0?std::erfc(-z1/s)-std::erfc(-z0/s):std::erf(z1/s)-std::erf(z0/s);
  double radial=std::exp(-a*a/(s*s))*(-std::expm1(-(b*b-a*a)/(s*s)))*s*s/(b*b-a*a);
  return 1e16*radial*s*std::sqrt(std::acos(-1.))*zdiff/(2*g.dz());
}
void dump_record(const AcceptedStepRecord&r,const std::filesystem::path&p){
  std::ofstream f(p.string()+"_faces.bin",std::ios::binary);
  for(const auto&x:r.electron_faces){int fields[]={x.i,x.j,x.neighbor_i,x.neighbor_j,int(x.direction),int(x.role),int(x.neighbor_material)};f.write(reinterpret_cast<char*>(fields),sizeof fields);for(double v:{x.gamma_e_m2s,x.area_m2})f.write(reinterpret_cast<char*>(&v),sizeof v);}
  std::ofstream c(p.string()+"_corrections.bin",std::ios::binary);
  for(const auto*v:{&r.pre_floor_ne_delta_m3,&r.pre_floor_np_delta_m3,&r.pre_floor_nn_delta_m3,&r.post_floor_ne_delta_m3,&r.post_floor_np_delta_m3,&r.post_floor_nn_delta_m3,&r.clip_ne_delta_m3,&r.clip_np_delta_m3,&r.clip_nn_delta_m3})c.write(reinterpret_cast<const char*>(v->data()),v->size()*sizeof(double));
}
void dump_state(const StreamerState&s,const std::filesystem::path&p){
  std::ofstream f(p,std::ios::binary);for(const auto*v:{&s.ne,&s.np,&s.nn,&s.rho,&s.phi,&s.er,&s.ez,&s.emag,&s.sph})f.write(reinterpret_cast<const char*>(v->values().data()),v->values().size()*sizeof(double));
}

int main(int argc,char**argv){
  PetscInitialize(&argc,&argv,nullptr,nullptr);
  try {
    if(argc!=2)throw std::invalid_argument("usage: native_full_driver fresh_output_dir");
    const std::filesystem::path out=argv[1];if(!std::filesystem::create_directory(out))throw std::runtime_error("REFUSE_OVERWRITE");
    std::filesystem::create_directory(out/"accepted_records");
    const auto begun=Clock::now();int calls=0,next_id=0,version=0,branch=0,commits=0;
    auto guard=[&](){if(calls>=120)throw std::runtime_error("NATIVE_STEP_BUDGET_EXHAUSTED");if(std::chrono::duration<double>(Clock::now()-begun).count()>1800)throw std::runtime_error("NATIVE_WALL_BUDGET_EXHAUSTED");};
    AxisymmetricGrid g(200,400,.0005,-.0005,.0005);
    RoundedConeProfile2D lo{-35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,-1},hi{35e-6,200e-6,std::acos(-1.)/12,1e-3,.03,1};
    const auto view=Geometry2DView::facing_profiles(g,lo,hi,"B0_R02_A30_g70_Nominal","export-profile-1","LOCAL_HV_RETURN");
    LinearInput wave;StreamerConfig cfg;cfg.geometry_view=&view;cfg.voltage_waveform=&wave;cfg.photoionization=false;
    const bool stable_mode=std::getenv("FP3_STABLE_PORT")&&std::string(std::getenv("FP3_STABLE_PORT"))=="1";
    cfg.elliptic={stable_mode?1e-14:1e-11,stable_mode?1e-25:1e-15,3000};cfg.poisson_boundary={{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0},{BoundaryKind::Neumann,0}};
    ScalarField2D cached_weight(g);
    const char* cache_path=stable_mode?std::getenv("FP3_CACHED_WEIGHT"):nullptr;
    if(cache_path){std::ifstream f(cache_path,std::ios::binary);f.read(reinterpret_cast<char*>(cached_weight.values().data()),g.size()*sizeof(double));if(!f||f.peek()!=EOF)throw std::runtime_error("INVALID_CACHED_WEIGHT_SIZE");}
    set_poisson_profiling(false);ActiveNetwork net(view,cfg,cache_path?&cached_weight:nullptr);
    std::ofstream identity(out/"operator_identity.txt");identity<<electrode_operator_identity(view,cfg.poisson_boundary,cfg.elliptic)<<'\n';identity.close();
    ScalarField2D zero(g),vacuum(g);if(!cache_path)solve_potential_with_electrodes(zero,view,500.,cfg.poisson_boundary,vacuum,cfg.elliptic);
    int static_rhs_calls=2;
    if(cache_path){const char* prior=std::getenv("FP3_PRIOR_STATIC_CALLS");if(!prior)throw std::runtime_error("MISSING_CACHE_RESOURCE_PROVENANCE");static_rhs_calls=std::stoi(prior);if(static_rhs_calls<0||static_rhs_calls>24)throw std::runtime_error("INVALID_PRIOR_STATIC_COUNT");}
    std::unique_ptr<StableNetwork> stable;
    if(stable_mode)stable=std::make_unique<StableNetwork>(net,cfg,static_rhs_calls);
    Metrics vacuum_audit{{"C_local_F",net.capacitance},{"Q_vacuum_500_C",double(net.charge(vacuum,500))},
      {"Q_minus_CV_C",double(net.charge(vacuum,500)-net.capacitance*500.)},
      {"unit_Q_minus_2U_C",double(net.charge(net.psi,1)-2*net.energy(net.psi,1))},
      {"U_vacuum_J",double(net.energy(vacuum,500))},{"active_network_edges",double(net.edges.size())}};
    if(cache_path){vacuum_audit={{"C_local_F",net.capacitance},{"active_network_edges",double(net.edges.size())},{"weight_cache_reused",1}};}
    StreamerSolver solver(g,cfg,{true,true,false});Buffer buffer;solver.set_accepted_step_observer(&buffer);
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)if(view.classify(i,j)==ElectrodeCellType::Gas)solver.state().ne(i,j)=solver.state().np(i,j)=average(g,i,j);
    solver.refresh_electrostatic_fields();solver.reset_electrode_history();
    const auto raw=StreamerTransactionAccess::capture(solver);
    StreamerTransactionAccess::initial_floor(solver);
    auto accepted=StreamerTransactionAccess::capture(solver);const auto initial=accepted;const LinearInput initial_wave=wave;
    if(stable){dump_state(accepted.state,out/"initial_state.bin");std::ofstream pf(out/"psi.bin",std::ios::binary);pf.write(reinterpret_cast<const char*>(net.psi.values().data()),g.size()*sizeof(double));}
    double initial_floor_N=0,initial_floor_Q=0;for(std::size_t k=0;k<g.size();++k){double vol=g.cell_volume(k%g.nr());initial_floor_N+=(accepted.state.ne.values()[k]-raw.state.ne.values()[k])*vol;initial_floor_Q+=(accepted.state.rho.values()[k]-raw.state.rho.values()[k])*vol;}
    const std::string instance="FULL_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    auto token=[&](){return instance+":"+std::to_string(branch)+":"+std::to_string(version);};
    std::map<int,std::unique_ptr<Trial>> candidates;
    std::map<int,std::unique_ptr<NativeFullSnapshot>> references;int next_reference=0;
    std::ofstream charge(out/"accepted_charge_budget.csv"),energy(out/"accepted_energy_budget.csv");
    charge<<std::setprecision(17);energy<<std::setprecision(17);bool csv_header=false;
    auto current_reply=[&](){auto actual=StreamerTransactionAccess::capture(solver);std::cout<<"{\"ok\":true,\"token\":\""<<token()<<"\",\"version\":"<<version<<",\"time_s\":"<<actual.state.time<<",\"voltage_V\":"<<wave.value(actual.state.time)
      <<",\"physical_hash\":\""<<native_hash(actual)<<"\",\"history_hash\":\""<<native_hash(actual,true)<<"\",\"step_calls\":"<<calls<<",\"commits\":"<<commits<<",\"static_rhs_calls\":"<<static_rhs_calls<<",\"metrics\":{";emit_metrics(vacuum_audit);std::cout<<"},\"initial_floor_N\":"<<initial_floor_N<<",\"initial_floor_Q_C\":"<<initial_floor_Q<<"}"<<std::endl;};
    std::cout<<std::setprecision(17);current_reply();
    std::string line;
    while(std::getline(std::cin,line)){
      std::istringstream in(line);std::string cmd;in>>cmd;
      try {
        if(cmd=="QUIT")break;
        if(cmd=="INFO"){current_reply();continue;}
        if(cmd=="STATIC_REFERENCE"){
          if(!stable)throw std::runtime_error("STABLE_MODE_NOT_ENABLED");
          std::string path;double voltage,rtol;in>>path>>voltage>>rtol;
          if(!in||(rtol!=1e-11&&rtol!=1e-13&&rtol!=1e-14)||static_rhs_calls+2>24)throw std::runtime_error("INVALID_STATIC_REFERENCE_OR_BUDGET");
          StreamerState s=accepted.state;
          if(path!="ACCEPTED"){
            std::ifstream file(path,std::ios::binary);if(!file)throw std::runtime_error("MISSING_STATIC_ENDPOINT");
            for(auto*p:{&s.ne,&s.np,&s.nn,&s.rho,&s.phi,&s.er,&s.ez,&s.emag,&s.sph})file.read(reinterpret_cast<char*>(p->values().data()),p->values().size()*sizeof(double));
            if(!file)throw std::runtime_error("INCOMPLETE_STATIC_ENDPOINT");
          }
          auto refcfg=cfg;refcfg.elliptic={rtol,rtol==1e-11?1e-15:1e-25,3000};
          ++static_rhs_calls;ActiveNetwork refnet(view,refcfg);StableNetwork ref(refnet,refcfg,static_rhs_calls);
          auto m=ref.endpoint(s,voltage);m["rtol"]=rtol;m["C_local_F"]=refnet.capacitance;m["explicit_static_rhs_calls"]=static_rhs_calls;
          std::cout<<"{\"ok\":true,\"step_calls\":"<<calls<<",\"metrics\":{";emit_metrics(m);std::cout<<"}}"<<std::endl;continue;
        }
        if(cmd=="COMPARE"){
          int ti,ri;in>>ti>>ri;auto a=candidates.find(ti);auto b=references.find(ri);
          if(a==candidates.end()||b==references.end())throw std::runtime_error("UNKNOWN_COMPARISON_ID");
          const auto&x=a->second->frame;const auto&y=*b->second;double mx=0,mr=0;int different=0;
          const ScalarField2D* xx[]={&x.state.ne,&x.state.np,&x.state.nn,&x.state.rho,&x.state.phi,&x.state.er,&x.state.ez,&x.state.emag,&x.state.sph};
          const ScalarField2D* yy[]={&y.state.ne,&y.state.np,&y.state.nn,&y.state.rho,&y.state.phi,&y.state.er,&y.state.ez,&y.state.emag,&y.state.sph};
          for(int n=0;n<9;++n){double err=0,scale=0;for(std::size_t k=0;k<g.size();++k){err=std::max(err,std::abs(xx[n]->values()[k]-yy[n]->values()[k]));scale=std::max(scale,std::abs(yy[n]->values()[k]));if(xx[n]->values()[k]!=yy[n]->values()[k])++different;}mx=std::max(mx,err);if(scale>0)mr=std::max(mr,err/scale);else if(err>0)mr=std::numeric_limits<double>::infinity();}
          std::cout<<"{\"ok\":true,\"different_values\":"<<different<<",\"maximum_absolute_field_difference\":"<<mx<<",\"maximum_field_relative_difference\":"<<mr<<",\"history_equal\":"<<(native_hash(x,true)==native_hash(y,true)?"true":"false")<<",\"step_calls\":"<<calls<<"}"<<std::endl;continue;
        }
        if(cmd=="RESET_INITIAL"){
          ++branch;version=0;candidates.clear();wave=initial_wave;StreamerTransactionAccess::restore(solver,initial);accepted=StreamerTransactionAccess::capture(solver);buffer.last.reset();current_reply();continue;
        }
        if(cmd=="REFERENCE"){
          int count,observe;double dt;in>>count>>dt>>observe;if(!in||count<1||count>2||!(dt>0&&dt<=20e-15))throw std::runtime_error("INVALID_REFERENCE");
          LinearInput rw=wave;rw.t0=accepted.state.time;rw.t1=rw.t0+count*dt;rw.v0=rw.v1=wave.value(accepted.state.time);
          auto rcfg=cfg;rcfg.voltage_waveform=&rw;StreamerSolver direct(g,rcfg,{true,true,false});
          auto fork=accepted;fork.owner=&direct;fork.config=rcfg;fork.observer=observe?&buffer:nullptr;StreamerTransactionAccess::restore(direct,fork);
          for(int n=0;n<count;++n){guard();++calls;StreamerDiagnostics d;if(!direct.step(dt,d))throw std::runtime_error("DIRECT_STEP_REJECTED");}
          auto ref=StreamerTransactionAccess::capture(direct);StreamerTransactionAccess::restore(solver,accepted);buffer.last.reset();
          int rid=++next_reference;references.emplace(rid,std::make_unique<NativeFullSnapshot>(ref));
          std::cout<<"{\"ok\":true,\"reference_id\":"<<rid<<",\"physical_hash\":\""<<native_hash(ref)<<"\",\"history_hash\":\""<<native_hash(ref,true)<<"\",\"step_calls\":"<<calls<<"}"<<std::endl;continue;
        }
        if(cmd=="COMMIT"||cmd=="REJECT"){
          std::string base;int id;in>>base>>id;
          auto found=candidates.find(id);if(!in||base!=token()||found==candidates.end()||found->second->base!=token())throw std::runtime_error("STALE_FOREIGN_OR_DUPLICATE_TOKEN");
          if(cmd=="REJECT"){candidates.erase(found);current_reply();continue;}
          auto& t=*found->second;
          wave=t.wave;StreamerTransactionAccess::restore(solver,t.frame);accepted=StreamerTransactionAccess::capture(solver);++version;++commits;
          for(std::size_t n=0;n<t.rows.size();++n){
            if(!csv_header){charge<<"branch,version,substep";for(auto[k,v]:t.rows[n])charge<<','<<k;charge<<'\n';energy<<"branch,version,substep";for(auto[k,v]:t.rows[n])energy<<','<<k;energy<<'\n';csv_header=true;}
            for(auto*f:{&charge,&energy}){*f<<branch<<','<<version<<','<<n;for(auto[k,v]:t.rows[n])*f<<','<<v;*f<<'\n';f->flush();}
            dump_record(*t.records[n],out/"accepted_records"/("b"+std::to_string(branch)+"_v"+std::to_string(version)+"_s"+std::to_string(n)));
            if(stable){auto prefix=out/"accepted_records"/("b"+std::to_string(branch)+"_v"+std::to_string(version)+"_s"+std::to_string(n));dump_state(t.records[n]->pre,prefix.string()+"_pre.bin");dump_state(t.records[n]->post,prefix.string()+"_post.bin");}
          }
          buffer.last.reset();candidates.clear();current_reply();continue;
        }
        if(cmd!="TRIAL")throw std::runtime_error("UNKNOWN_COMMAND");
        std::string base;double tn,dt,v0,v1;int inject;in>>base>>tn>>dt>>v0>>v1>>inject;
        if(!in||base!=token())throw std::runtime_error("STALE_OR_FOREIGN_TOKEN");
        if(!std::isfinite(dt)||!std::isfinite(v0)||!std::isfinite(v1)||!(dt>0&&dt<=20e-15)||tn!=accepted.state.time||v0!=wave.value(tn)||v1<499||v1>501)throw std::runtime_error("INVALID_TRIAL_TIME_OR_VOLTAGE");
        guard();const auto saved_wave=wave;wave.t0=tn;wave.t1=tn+dt;wave.v0=v0;wave.v1=v1;
        auto trial=std::make_unique<Trial>(g);trial->base=base;trial->id=++next_id;
        try {
          StreamerTransactionAccess::restore(solver,accepted);solver.set_accepted_step_observer(&buffer);
          while(solver.state().time<wave.t1){
            guard();const auto start=solver.state();const auto limits=solver.timestep_limits();
            double subdt=std::min({20e-15,limits.selected,wave.t1-start.time});
            if(!std::isfinite(subdt)||subdt<=0)throw std::runtime_error("INVALID_NATIVE_SUBSTEP");
            ++calls;StreamerDiagnostics diagnostic;buffer.last.reset();if(!solver.step(subdt,diagnostic))throw std::runtime_error("NATIVE_STEP_REJECTED");
            if(inject)throw std::runtime_error("CONTROLLED_FAILURE_AFTER_NATIVE_STATE_CHANGE");
            auto m=net.budget(start,*buffer.last,cfg);
            if(stable)m=stable->budget(start,*buffer.last,std::move(m));
            m["native_dt_limit_s"]=limits.selected;m["native_poisson_iterations"]=diagnostic.poisson_iterations;
            trial->rows.push_back(m);trial->records.emplace_back(buffer.last.release());
          }
          trial->frame=StreamerTransactionAccess::capture(solver);trial->wave=wave;
          auto&m=trial->m;m["time_n_s"]=tn;m["time_n1_s"]=tn+dt;m["dt_s"]=dt;m["Vp_n_V"]=v0;m["Vp_n1_V"]=v1;m["substeps"]=trial->rows.size();m["C_local_F"]=net.capacitance;
          for(const auto&r:trial->rows)for(auto key:{"Deltaq_port_C","Deltaq_Hg_C","Deltaq_w_C","Deltaq_outer_C","DeltaQ_o_C","DeltaQ_psi_correction_C","DeltaQ_psi_reaction_C","charge_weighted_residual_C","unsigned_native_transport_C","correction_charge_C","DeltaU_E_J","Wport_native_J","Wcond_face_J","Wouter_J","Wcorrection_J","Wreaction_J","energy_residual_J","excluded_boundary_rho_change_TV_C","nonfinite_active_mu_cells"})m[key]+=r.at(key);
          m["U_E_start_J"]=trial->rows.front().at("U_E_start_J");m["U_E_end_J"]=trial->rows.back().at("U_E_end_J");
          m["q_B_active_start_C"]=trial->rows.front().at("q_B_active_start_C");m["q_B_active_end_C"]=trial->rows.back().at("q_B_active_end_C");
          m["q_A_end_C"]=trial->rows.back().at("q_A_end_C");m["q_B_native_links_end_C"]=trial->rows.back().at("q_B_native_links_end_C");
          m["Wport_circuit_J"]=.5*(v0+v1)*m["Deltaq_port_C"];m["macro_minus_native_work_J"]=m["Wport_circuit_J"]-m["Wport_native_J"];
          if(stable){
            for(const auto&r:trial->rows)for(const auto&kv:r){
              const auto& key=kv.first;
              if(key.rfind("stable_",0)==0&&key.find("remaining_start")==std::string::npos&&key.find("remaining_end")==std::string::npos)m[key]+=kv.second;
            }
            m["stable_U_remaining_start_J"]=trial->rows.front().at("stable_U_remaining_start_J");
            m["stable_U_remaining_end_J"]=trial->rows.back().at("stable_U_remaining_end_J");
          }
          auto ph=native_hash(trial->frame),hh=native_hash(trial->frame,true);int id=trial->id;candidates.emplace(id,std::move(trial));
          wave=saved_wave;StreamerTransactionAccess::restore(solver,accepted);buffer.last.reset();
          std::cout<<"{\"ok\":true,\"trial_id\":"<<id<<",\"physical_hash\":\""<<ph<<"\",\"history_hash\":\""<<hh<<"\",\"accepted_hash\":\""<<native_hash(StreamerTransactionAccess::capture(solver))<<"\",\"step_calls\":"<<calls<<",\"metrics\":{";emit_metrics(candidates.at(id)->m);std::cout<<"}}"<<std::endl;
        }catch(...){wave=saved_wave;StreamerTransactionAccess::restore(solver,accepted);buffer.last.reset();throw;}
      }catch(const std::exception&error){std::cout<<"{\"ok\":false,\"error\":\""<<error.what()<<"\",\"step_calls\":"<<calls<<",\"accepted_hash\":\""<<native_hash(accepted)<<"\"}"<<std::endl;}
    }
    std::ofstream final(out/"native_final_state.bin",std::ios::binary);for(const auto*p:{&accepted.state.ne,&accepted.state.np,&accepted.state.nn,&accepted.state.rho,&accepted.state.phi,&accepted.state.er,&accepted.state.ez,&accepted.state.emag,&accepted.state.sph})final.write(reinterpret_cast<const char*>(p->values().data()),p->values().size()*sizeof(double));
    std::ofstream status(out/"native_resource.json");status<<std::setprecision(17)<<"{\"native_step_calls\":"<<calls<<",\"commits\":"<<commits<<",\"explicit_static_rhs_calls\":"<<static_rhs_calls<<",\"wall_s\":"<<std::chrono::duration<double>(Clock::now()-begun).count()<<"}\n";
  }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;PetscFinalize();return 1;}
  PetscFinalize();return 0;
}
