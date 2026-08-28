#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <stdexcept>
namespace streamer_rf::streamer {
namespace {constexpr double qe=1.602176634e-19,eps0=8.8541878128e-12;}
StreamerSolver::StreamerSolver(const AxisymmetricGrid&g,StreamerConfig c):g_(g),c_(c),state_(g){}
double axisymmetric_radial_face_area(const AxisymmetricGrid& g,int face_i){return 2*3.14159265358979323846*g.radial_face(face_i)*g.dz();}
double axisymmetric_axial_face_area(const AxisymmetricGrid& g,int cell_i){return 3.14159265358979323846*(std::pow(g.radial_face(cell_i+1),2)-std::pow(g.radial_face(cell_i),2));}
double absorbing_electrode_flux(double ne,double velocity_positive,double diffusion,double h,bool positive_face){
 if(ne<0||diffusion<0||h<=0)throw std::invalid_argument("invalid absorbing electrode flux input");
 const double diffusive=2.0*diffusion*ne/h;
 return positive_face?std::max(velocity_positive,0.0)*ne+diffusive:std::min(velocity_positive,0.0)*ne-diffusive;
}
double gas_gas_electron_flux(double nl,double nr,double vl,double vr,double d,double h,double eps,double ref){
 if(d<=0&&std::abs(vl)+std::abs(vr)<1e-300)return 0.0;
 return isg0_flux(nl,nr,vl,vr,d,h,eps,ref);
}
ElectrodeCellType StreamerSolver::cell_type(int i,int j)const{return electrode_mode()?c_.electrode_geometry->classify(g_,i,j):ElectrodeCellType::Gas;}
void StreamerSolver::enforce_plasma_mask(){
 if(!electrode_mode())return;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(!is_gas(i,j)){state_.ne(i,j)=0;state_.np(i,j)=0;state_.nn(i,j)=0;state_.sph(i,j)=0;state_.rho(i,j)=0;}
}
void StreamerSolver::save_checkpoint(const std::filesystem::path&p)const{std::ofstream f(p,std::ios::binary);const char magic[8]={'S','R','F','S','3','C','P','1'};f.write(magic,8);int nr=g_.nr(),nz=g_.nz();f.write(reinterpret_cast<const char*>(&nr),sizeof nr);f.write(reinterpret_cast<const char*>(&nz),sizeof nz);f.write(reinterpret_cast<const char*>(&state_.time),sizeof state_.time);for(const auto*field:{&state_.ne,&state_.np,&state_.nn})f.write(reinterpret_cast<const char*>(field->values().data()),field->values().size()*sizeof(double));if(!f)throw std::runtime_error("checkpoint write failed");}
void StreamerSolver::load_checkpoint(const std::filesystem::path&p){std::ifstream f(p,std::ios::binary);char magic[8];int nr,nz;f.read(magic,8);f.read(reinterpret_cast<char*>(&nr),sizeof nr);f.read(reinterpret_cast<char*>(&nz),sizeof nz);if(std::string(magic,8)!="SRFS3CP1"||nr!=g_.nr()||nz!=g_.nz())throw std::runtime_error("checkpoint mismatch");f.read(reinterpret_cast<char*>(&state_.time),sizeof state_.time);for(auto*field:{&state_.ne,&state_.np,&state_.nn})f.read(reinterpret_cast<char*>(field->values().data()),field->values().size()*sizeof(double));if(!f)throw std::runtime_error("checkpoint read failed");last_head_position_=last_head_time_=std::numeric_limits<double>::quiet_NaN();fields();reset_electrode_history();}
void StreamerSolver::initialize_gaussian(double n0,double sigma,double z0){
 initialize_gaussians({GaussianSeed{n0,sigma,z0}});
}
void StreamerSolver::initialize_gaussian_at_tip_offset(double n0,double sigma,double z_offset_from_tip){
 if(!c_.electrode_geometry)throw std::runtime_error("tip-relative seed requires electrode geometry");
 initialize_gaussian(n0,sigma,c_.electrode_geometry->seed_z_from_tip_offset(z_offset_from_tip));
}
void StreamerSolver::initialize_gaussians(const std::vector<GaussianSeed>& seeds){
 for(int j=0;j<g_.nz();++j){
  for(int i=0;i<g_.nr();++i){
   double n=0.0;
   for(const auto& seed:seeds)n+=seed.n0*std::exp(-(g_.r(i)*g_.r(i)+std::pow(g_.z(j)-seed.z0,2))/(seed.sigma*seed.sigma));
   state_.ne(i,j)=state_.np(i,j)=n;
   state_.nn(i,j)=0;
  }
 }
 state_.time=0.0;
 last_head_position_=last_head_time_=std::numeric_limits<double>::quiet_NaN();
 enforce_plasma_mask();
 fields();
 reset_electrode_history();
}
void StreamerSolver::fields(){
 enforce_plasma_mask();
 for(std::size_t k=0;k<state_.rho.values().size();++k)state_.rho.values()[k]=qe*(state_.np.values()[k]-state_.ne.values()[k]-state_.nn.values()[k]);
 enforce_plasma_mask();
 PoissonBoundaryConfig b=c_.poisson_boundary;if(c_.electrode_geometry&&c_.voltage_waveform)solve_potential_with_electrodes(state_.rho,*c_.electrode_geometry,c_.voltage_waveform->value(state_.time),b,state_.phi,c_.elliptic,c_.open_boundary);else solve_potential(state_.rho,c_.background_field,b,state_.phi,c_.elliptic,c_.open_boundary);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){auto dr=[&](int a,int b){return(state_.phi(b,j)-state_.phi(a,j))/((b-a)*g_.dr());};auto dz=[&](int a,int b){return(state_.phi(i,b)-state_.phi(i,a))/((b-a)*g_.dz());};state_.er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));state_.ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));state_.emag(i,j)=std::hypot(state_.er(i,j),state_.ez(i,j));}
}
void StreamerSolver::fill_electrode_diagnostics(StreamerDiagnostics&diag)const{
 if(!c_.electrode_geometry||!c_.voltage_waveform)return;
 auto e=evaluate_electrode_diagnostics(state_.phi,state_.er,state_.ez,*c_.electrode_geometry,state_.time,c_.voltage_waveform->value(state_.time),last_poisson_iterations());
 diag.applied_voltage=e.applied_voltage_V;diag.phi_hv_residual=e.phi_hv_residual_V;diag.phi_ground_residual=e.phi_ground_residual_V;diag.geometry_id=e.geometry_id;diag.gap=e.gap_m;diag.tip_radius=e.tip_radius_m;
}
double StreamerSolver::head_position()const{
 double best_ne=-1.0,best_z=0.0;
 bool threshold_hit=false;
 for(int j=0;j<g_.nz();++j)if(is_gas(0,j)){
  const double n=state_.ne(0,j);
  if(n>=c_.head_ne_threshold){if(!threshold_hit){best_z=g_.z(j);threshold_hit=true;}else best_z=std::min(best_z,g_.z(j));}
  if(!threshold_hit&&n>best_ne){best_ne=n;best_z=g_.z(j);}
 }
 return best_z;
}
bool StreamerSolver::bridge_flag()const{
 if(!electrode_mode())return false;
 bool in_gap=false;
 for(int j=0;j<g_.nz();++j){
  if(!is_gas(0,j))continue;
  const double z=g_.z(j);
  if(z<=c_.electrode_geometry->ground_z_m()+c_.electrode_geometry->ground_thickness_m()||z>=c_.electrode_geometry->tip_z_m())continue;
  in_gap=true;
  if(state_.ne(0,j)<c_.bridge_ne_threshold)return false;
 }
 return in_gap;
}
ElectrodeSurfaceDiagnostics StreamerSolver::electrode_surface_diagnostics()const{
 ElectrodeSurfaceDiagnostics d;if(!electrode_mode())return d;
 auto add=[&](ElectrodeCellType t,double en,double area){if(t==ElectrodeCellType::HighVoltageElectrode){d.q_hv+=eps0*en*area;d.area_hv+=area;}else if(t==ElectrodeCellType::GroundElectrode){d.q_ground+=eps0*en*area;d.area_ground+=area;}};
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j))continue;
  if(i<g_.nr()-1&&!is_gas(i+1,j))add(cell_type(i+1,j),-state_.er(i,j),axisymmetric_radial_face_area(g_,i+1));
  if(i>0&&!is_gas(i-1,j))add(cell_type(i-1,j),state_.er(i,j),axisymmetric_radial_face_area(g_,i));
  if(j<g_.nz()-1&&!is_gas(i,j+1))add(cell_type(i,j+1),-state_.ez(i,j),axisymmetric_axial_face_area(g_,i));
  if(j>0&&!is_gas(i,j-1))add(cell_type(i,j-1),state_.ez(i,j),axisymmetric_axial_face_area(g_,i));
 }
 return d;
}
ConductanceDiagnostics StreamerSolver::conductance_diagnostics(double voltage)const{
 ConductanceDiagnostics d;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);const double sigma=qe*q.mobility*state_.ne(i,j);d.p_cond+=sigma*state_.emag(i,j)*state_.emag(i,j)*g_.cell_volume(i);}
 const double vtol=std::max(c_.voltage_tolerance,1024*std::numeric_limits<double>::epsilon());
 d.valid=std::isfinite(voltage)&&std::abs(voltage)>vtol;
 if(d.valid){d.gb=d.p_cond/(voltage*voltage);d.rb=d.gb>0?1.0/d.gb:std::numeric_limits<double>::infinity();}
 else{d.gb=std::numeric_limits<double>::quiet_NaN();d.rb=std::numeric_limits<double>::quiet_NaN();}
 return d;
}
double StreamerSolver::vacuum_gap_capacitance()const{
 if(!electrode_mode())return 0.0;
 if(std::isfinite(c_gap_vacuum_cache_))return c_gap_vacuum_cache_;
 ScalarField2D rho(g_),phi(g_),er(g_),ez(g_);PoissonBoundaryConfig b;b.r_outer.kind=b.z_lower.kind=b.z_upper.kind=BoundaryKind::OpenCharge;solve_potential_with_electrodes(rho,*c_.electrode_geometry,1.0,b,phi,c_.elliptic,c_.open_boundary);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){auto dr=[&](int a,int b){return(phi(b,j)-phi(a,j))/((b-a)*g_.dr());};auto dz=[&](int a,int b){return(phi(i,b)-phi(i,a))/((b-a)*g_.dz());};er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));}
 double q_hv=0.0;auto add=[&](ElectrodeCellType t,double en,double area){if(t==ElectrodeCellType::HighVoltageElectrode)q_hv+=eps0*en*area;};
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(c_.electrode_geometry->classify(g_,i,j)!=ElectrodeCellType::Gas)continue;if(i<g_.nr()-1&&c_.electrode_geometry->classify(g_,i+1,j)!=ElectrodeCellType::Gas)add(c_.electrode_geometry->classify(g_,i+1,j),-er(i,j),axisymmetric_radial_face_area(g_,i+1));if(i>0&&c_.electrode_geometry->classify(g_,i-1,j)!=ElectrodeCellType::Gas)add(c_.electrode_geometry->classify(g_,i-1,j),er(i,j),axisymmetric_radial_face_area(g_,i));if(j<g_.nz()-1&&c_.electrode_geometry->classify(g_,i,j+1)!=ElectrodeCellType::Gas)add(c_.electrode_geometry->classify(g_,i,j+1),-ez(i,j),axisymmetric_axial_face_area(g_,i));if(j>0&&c_.electrode_geometry->classify(g_,i,j-1)!=ElectrodeCellType::Gas)add(c_.electrode_geometry->classify(g_,i,j-1),ez(i,j),axisymmetric_axial_face_area(g_,i));}
 c_gap_vacuum_cache_=std::abs(q_hv);return c_gap_vacuum_cache_;
}
void StreamerSolver::reset_electrode_history(){if(!electrode_mode()){has_electrode_history_=false;return;}auto q=electrode_surface_diagnostics();last_q_hv_=q.q_hv;last_q_ground_=q.q_ground;has_electrode_history_=true;}
void StreamerSolver::refresh_electrostatic_fields(){fields();}
double StreamerSolver::electrode_displacement_current(double q1,double q0,double dt)const{return dt>0?(q1-q0)/dt:std::numeric_limits<double>::quiet_NaN();}
StreamerDiagnostics StreamerSolver::sample_terminal_diagnostics_from_history(double dt){
 StreamerDiagnostics diag{};diag.time=state_.time;diag.dt=dt;diag.poisson_iterations=last_poisson_iterations();diag.emax=*std::max_element(state_.emag.values().begin(),state_.emag.values().end());diag.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());diag.np_max=*std::max_element(state_.np.values().begin(),state_.np.values().end());diag.nn_max=*std::max_element(state_.nn.values().begin(),state_.nn.values().end());fill_electrode_diagnostics(diag);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;double vol=g_.cell_volume(i);diag.total_electrons+=state_.ne(i,j)*vol;diag.total_charge+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);diag.sigma_max=std::max(diag.sigma_max,qe*q.mobility*state_.ne(i,j));}
 if(electrode_mode()){auto qsurf=electrode_surface_diagnostics();diag.q_hv=qsurf.q_hv;diag.q_ground=qsurf.q_ground;diag.displacement_current_valid=has_electrode_history_&&dt>0;diag.c_gap_vacuum=vacuum_gap_capacitance();if(diag.displacement_current_valid){diag.i_disp_hv=electrode_displacement_current(qsurf.q_hv,last_q_hv_,dt);diag.i_disp_ground=electrode_displacement_current(qsurf.q_ground,last_q_ground_,dt);diag.i_total_hv=diag.i_cond_hv+diag.i_disp_hv;diag.i_total_ground=diag.i_cond_ground+diag.i_disp_ground;}last_q_hv_=qsurf.q_hv;last_q_ground_=qsurf.q_ground;has_electrode_history_=true;}
 auto gd=conductance_diagnostics(diag.applied_voltage);diag.p_cond=gd.p_cond;diag.gb=gd.gb;diag.rb=gd.rb;diag.rb_valid=gd.valid;diag.head_position=head_position();diag.bridge_flag=bridge_flag();return diag;
}
double StreamerSolver::numerical_density_tolerance()const{
 double ne_max=0.0;
 if(electrode_mode()){for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(is_gas(i,j))ne_max=std::max(ne_max,state_.ne(i,j));}
 else ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());
 return 1e-12*std::max(ne_max,c_.n_ref);
}
double StreamerSolver::stage_c_activity_tolerance()const{
 double ne_max=0.0;for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(is_gas(i,j))ne_max=std::max(ne_max,state_.ne(i,j));
 double tol=1000.0*numerical_density_tolerance();
 if(ne_max>1000.0*c_.n_ref)tol=std::max(tol,c_.isg_epsilon*c_.n_ref);
 return tol;
}
void StreamerSolver::enforce_stage_c_activity_floor(){
 if(!electrode_mode())return;
 const double tol=stage_c_activity_tolerance();
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;if(state_.ne(i,j)<=tol)state_.ne(i,j)=0;if(state_.np(i,j)<=tol)state_.np(i,j)=0;if(state_.nn(i,j)<=tol)state_.nn(i,j)=0;}
}
ReactionTimestepDiagnostic StreamerSolver::reaction_timestep_diagnostic()const{
 ReactionTimestepDiagnostic d;d.numerical_density_tolerance=numerical_density_tolerance();d.reaction_dt=std::numeric_limits<double>::infinity();
 if(electrode_mode()){for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(is_gas(i,j))d.ne_max=std::max(d.ne_max,state_.ne(i,j));}
 else d.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){const auto t=cell_type(i,j);if(t!=ElectrodeCellType::Gas)continue;const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);const double recomb=beta_ep(q.mobility,q.diffusion)*state_.np.values()[k];const double loss=q.attachment_two_body_frequency+q.attachment_three_body_frequency+recomb;const double rdt=.2/std::max(loss,1e-300);if(rdt<d.reaction_dt){d.i=i;d.j=j;d.r=g_.r(i);d.z=g_.z(j);d.E=state_.emag.values()[k];d.E_over_N_Td=d.E/c_.neutral_density*1e21;d.ne=state_.ne.values()[k];d.np=state_.np.values()[k];d.nn=state_.nn.values()[k];d.ionization_frequency=q.ionization_frequency;d.attachment_two_body_frequency=q.attachment_two_body_frequency;d.attachment_three_body_frequency=q.attachment_three_body_frequency;d.recombination_frequency=recomb;d.reaction_dt=rdt;d.cell_classification="Gas";}}
 return d;
}
TimeStepLimits StreamerSolver::timestep_limits()const{
 if(!electrode_mode()){
 double vr=0,vz=0,dmax=0,numax=0,cond=0,loss=0;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);vr=std::max(vr,q.mobility*std::abs(state_.er.values()[k]));vz=std::max(vz,q.mobility*std::abs(state_.ez.values()[k]));dmax=std::max(dmax,q.diffusion);numax=std::max(numax,q.ionization_frequency);cond=std::max(cond,qe*q.mobility*state_.ne.values()[k]);loss=std::max(loss,q.attachment_two_body_frequency+q.attachment_three_body_frequency+beta_ep(q.mobility,q.diffusion)*state_.np.values()[k]);}
 TimeStepLimits x;x.drift=.5*std::min(g_.dr()/std::max(vr,1e-300),g_.dz()/std::max(vz,1e-300));x.diffusion=.5/(2*std::max(dmax,1e-300)*(1/(g_.dr()*g_.dr())+1/(g_.dz()*g_.dz())));x.ionization=.05/std::max(numax,1e-300);x.dielectric=.2*eps0/std::max(cond,1e-300);x.reaction=.2/std::max(loss,1e-300);x.selected=x.drift;x.controller="drift";for(auto p:{std::pair{x.diffusion,"diffusion"},std::pair{x.ionization,"ionization"},std::pair{x.dielectric,"dielectric"},std::pair{x.reaction,"reaction"}})if(p.first<x.selected){x.selected=p.first;x.controller=p.second;}return x;
 }
 double vr=0,vz=0,dmax=0,numax=0,cond=0,reaction=std::numeric_limits<double>::infinity();
 const double tol=stage_c_activity_tolerance();
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j))continue;
  const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;
  const double ne=state_.ne.values()[k],np=state_.np.values()[k],nn=state_.nn.values()[k];
  const bool electron_active=ne>tol,positive_active=np>tol,negative_active=nn>tol;
  if(!electron_active&&!positive_active&&!negative_active)continue;
  auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);
  if(electron_active){
   vr=std::max(vr,q.mobility*std::abs(state_.er.values()[k]));vz=std::max(vz,q.mobility*std::abs(state_.ez.values()[k]));dmax=std::max(dmax,q.diffusion);numax=std::max(numax,q.ionization_frequency);cond=std::max(cond,qe*q.mobility*ne);
  }
  auto s=evaluate_reactions(ne,np,nn,state_.sph.values()[k],q,c_.temperature);
  if(electron_active&&s.electron<0)reaction=std::min(reaction,.2*ne/(-s.electron));
  if(positive_active&&s.positive_ion<0)reaction=std::min(reaction,.2*np/(-s.positive_ion));
  if(negative_active&&s.negative_ion<0)reaction=std::min(reaction,.2*nn/(-s.negative_ion));
 }
 if(!std::isfinite(reaction))reaction=std::numeric_limits<double>::max();
 TimeStepLimits x;x.drift=.5*std::min(g_.dr()/std::max(vr,1e-300),g_.dz()/std::max(vz,1e-300));x.diffusion=.5/(2*std::max(dmax,1e-300)*(1/(g_.dr()*g_.dr())+1/(g_.dz()*g_.dz())));x.ionization=.05/std::max(numax,1e-300);x.dielectric=.2*eps0/std::max(cond,1e-300);x.reaction=reaction;x.selected=x.drift;x.controller="drift";for(auto p:{std::pair{x.diffusion,"diffusion"},std::pair{x.ionization,"ionization"},std::pair{x.dielectric,"dielectric"},std::pair{x.reaction,"reaction"}})if(p.first<x.selected){x.selected=p.first;x.controller=p.second;}return x;
}
bool StreamerSolver::step(double dt,StreamerDiagnostics&diag){
 enforce_plasma_mask();
 enforce_stage_c_activity_floor();
 const bool had_electrode_history=has_electrode_history_;
 const double prev_q_hv=last_q_hv_,prev_q_ground=last_q_ground_;
 ScalarField2D emission(g_);for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;if(!is_gas(i,j)){emission.values()[k]=0;continue;}auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);emission.values()[k]=c_.excitation_ratio*q.ionization_frequency*state_.ne.values()[k];}
 if(c_.photoionization){PoissonBoundaryConfig b;solve_photoionization(emission,{},b,state_.sph,c_.elliptic,c_.sp3_boundary,nullptr,&sp3_warm_);if(electrode_mode())for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(!is_gas(i,j))state_.sph(i,j)=0;}else std::fill(state_.sph.values().begin(),state_.sph.values().end(),0);
 ScalarField2D dne(g_),dnp(g_),dnn(g_);double old_ne=0,old_q=0,charge_inventory=0,re_int=0,boundary_out=0,hv_abs=0,ground_abs=0;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j)){dne(i,j)=dnp(i,j)=dnn(i,j)=0;continue;}
  auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);auto s=evaluate_reactions(state_.ne(i,j),state_.np(i,j),state_.nn(i,j),state_.sph(i,j),q,c_.temperature);double div=0;
  if(i<g_.nr()-1&&is_gas(i+1,j)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i+1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i+1,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i,j),vr=-r.mobility*state_.er(i+1,j);div+=g_.radial_face(i+1)*gas_gas_electron_flux(state_.ne(i,j),state_.ne(i+1,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref)/(g_.r(i)*g_.dr());}
  else if(i<g_.nr()-1){double v=-q.mobility*state_.er(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dr(),true);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());double loss=flux*2*3.14159265358979323846*g_.radial_face(i+1)*g_.dz();if(cell_type(i+1,j)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.er(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());boundary_out+=flux*2*3.14159265358979323846*g_.radial_face(i+1)*g_.dz();}
  if(i>0&&is_gas(i-1,j)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i-1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i-1,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i-1,j),vr=-r.mobility*state_.er(i,j);div-=g_.radial_face(i)*gas_gas_electron_flux(state_.ne(i-1,j),state_.ne(i,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref)/(g_.r(i)*g_.dr());}
  else if(i>0){double v=-q.mobility*state_.er(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dr(),false);div-=g_.radial_face(i)*flux/(g_.r(i)*g_.dr());double loss=-flux*2*3.14159265358979323846*g_.radial_face(i)*g_.dz();if(cell_type(i-1,j)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  if(j<g_.nz()-1&&is_gas(i,j+1)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j+1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j+1),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j),vr=-r.mobility*state_.ez(i,j+1);div+=gas_gas_electron_flux(state_.ne(i,j),state_.ne(i,j+1),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref)/g_.dz();}
  else if(j<g_.nz()-1){double v=-q.mobility*state_.ez(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dz(),true);div+=flux/g_.dz();double loss=flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));if(cell_type(i,j+1)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=flux/g_.dz();boundary_out+=flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  if(j>0&&is_gas(i,j-1)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j-1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j-1),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j-1),vr=-r.mobility*state_.ez(i,j);div-=gas_gas_electron_flux(state_.ne(i,j-1),state_.ne(i,j),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref)/g_.dz();}
  else if(j>0){double v=-q.mobility*state_.ez(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dz(),false);div-=flux/g_.dz();double loss=-flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));if(cell_type(i,j-1)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::min(v,0.0)*state_.ne(i,j);div-=flux/g_.dz();boundary_out+=-flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  dne(i,j)=-div+s.electron;dnp(i,j)=s.positive_ion;dnn(i,j)=s.negative_ion;double vol=g_.cell_volume(i);old_ne+=state_.ne(i,j)*vol;old_q+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;charge_inventory+=qe*(state_.np(i,j)+state_.ne(i,j)+state_.nn(i,j))*vol;re_int+=s.electron*vol;
 }
 double tol=electrode_mode()?stage_c_activity_tolerance():numerical_density_tolerance();
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;if(!is_gas(i,j)){state_.ne.values()[k]=state_.np.values()[k]=state_.nn.values()[k]=0;continue;}double ne=state_.ne.values()[k]+dt*dne.values()[k],np=state_.np.values()[k]+dt*dnp.values()[k],nn=state_.nn.values()[k]+dt*dnn.values()[k];if(ne<-tol||np<-tol||nn<-tol)return false;state_.ne.values()[k]=std::max(0.,ne);state_.np.values()[k]=std::max(0.,np);state_.nn.values()[k]=std::max(0.,nn);}
 enforce_stage_c_activity_floor();
 state_.time+=dt;fields();diag={};diag.time=state_.time;diag.dt=dt;diag.poisson_iterations=last_poisson_iterations();diag.sp3_ksp_iterations=c_.photoionization?last_sp3_ksp_iterations():0;diag.sp3_boundary_iterations=c_.photoionization?last_sp3_boundary_iterations():0;diag.emax=*std::max_element(state_.emag.values().begin(),state_.emag.values().end());diag.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());diag.np_max=*std::max_element(state_.np.values().begin(),state_.np.values().end());diag.nn_max=*std::max_element(state_.nn.values().begin(),state_.nn.values().end());diag.controller=timestep_limits().controller;fill_electrode_diagnostics(diag);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;double vol=g_.cell_volume(i);diag.total_electrons+=state_.ne(i,j)*vol;diag.total_charge+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);diag.sigma_max=std::max(diag.sigma_max,qe*q.mobility*state_.ne(i,j));}
 diag.absorbed_electron_hv=dt*hv_abs;diag.absorbed_electron_ground=dt*ground_abs;diag.i_cond_hv=qe*hv_abs;diag.i_cond_ground=qe*ground_abs;diag.outer_boundary_current=qe*boundary_out;
 if(electrode_mode()){auto qsurf=electrode_surface_diagnostics();diag.q_hv=qsurf.q_hv;diag.q_ground=qsurf.q_ground;diag.displacement_current_valid=had_electrode_history&&dt>0;diag.c_gap_vacuum=vacuum_gap_capacitance();if(diag.displacement_current_valid){diag.i_disp_hv=electrode_displacement_current(qsurf.q_hv,prev_q_hv,dt);diag.i_disp_ground=electrode_displacement_current(qsurf.q_ground,prev_q_ground,dt);diag.i_total_hv=diag.i_cond_hv+diag.i_disp_hv;diag.i_total_ground=diag.i_cond_ground+diag.i_disp_ground;}last_q_hv_=qsurf.q_hv;last_q_ground_=qsurf.q_ground;has_electrode_history_=true;}
 auto gd=conductance_diagnostics(diag.applied_voltage);diag.p_cond=gd.p_cond;diag.gb=gd.gb;diag.rb=gd.rb;diag.rb_valid=gd.valid;
 diag.head_position=head_position();diag.bridge_flag=bridge_flag();if(std::isfinite(last_head_position_)&&diag.time>last_head_time_)diag.head_velocity=(diag.head_position-last_head_position_)/(diag.time-last_head_time_);last_head_position_=diag.head_position;last_head_time_=diag.time;
 const double total_loss=boundary_out+hv_abs+ground_abs;double ne_expected=old_ne+dt*(re_int-total_loss),q_expected=old_q+qe*dt*total_loss;double re=std::abs(diag.total_electrons-ne_expected)/std::max({std::abs(diag.total_electrons),std::abs(ne_expected),1.0}),rq=std::abs(diag.total_charge-q_expected)/std::max(charge_inventory,1e-30);diag.conservation_residual=std::max(re,rq);diag.plasma_charge_derivative=(diag.total_charge-old_q)/dt;const double expected_charge_rate=diag.i_cond_hv+diag.i_cond_ground+diag.outer_boundary_current;diag.current_continuity_residual=std::abs(diag.plasma_charge_derivative-expected_charge_rate)/std::max({std::abs(diag.plasma_charge_derivative),std::abs(expected_charge_rate),1e-30});return true;
}
}
