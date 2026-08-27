#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
namespace streamer_rf::streamer {
namespace {constexpr double qe=1.602176634e-19,eps0=8.8541878128e-12;}
StreamerSolver::StreamerSolver(const AxisymmetricGrid&g,StreamerConfig c):g_(g),c_(c),state_(g){}
void StreamerSolver::save_checkpoint(const std::filesystem::path&p)const{std::ofstream f(p,std::ios::binary);const char magic[8]={'S','R','F','S','3','C','P','1'};f.write(magic,8);int nr=g_.nr(),nz=g_.nz();f.write(reinterpret_cast<const char*>(&nr),sizeof nr);f.write(reinterpret_cast<const char*>(&nz),sizeof nz);f.write(reinterpret_cast<const char*>(&state_.time),sizeof state_.time);for(const auto*field:{&state_.ne,&state_.np,&state_.nn})f.write(reinterpret_cast<const char*>(field->values().data()),field->values().size()*sizeof(double));if(!f)throw std::runtime_error("checkpoint write failed");}
void StreamerSolver::load_checkpoint(const std::filesystem::path&p){std::ifstream f(p,std::ios::binary);char magic[8];int nr,nz;f.read(magic,8);f.read(reinterpret_cast<char*>(&nr),sizeof nr);f.read(reinterpret_cast<char*>(&nz),sizeof nz);if(std::string(magic,8)!="SRFS3CP1"||nr!=g_.nr()||nz!=g_.nz())throw std::runtime_error("checkpoint mismatch");f.read(reinterpret_cast<char*>(&state_.time),sizeof state_.time);for(auto*field:{&state_.ne,&state_.np,&state_.nn})f.read(reinterpret_cast<char*>(field->values().data()),field->values().size()*sizeof(double));if(!f)throw std::runtime_error("checkpoint read failed");fields();}
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
 fields();
}
void StreamerSolver::fields(){
 for(std::size_t k=0;k<state_.rho.values().size();++k)state_.rho.values()[k]=qe*(state_.np.values()[k]-state_.ne.values()[k]-state_.nn.values()[k]);
 PoissonBoundaryConfig b;b.r_outer.kind=b.z_lower.kind=b.z_upper.kind=BoundaryKind::OpenCharge;if(c_.electrode_geometry&&c_.voltage_waveform)solve_potential_with_electrodes(state_.rho,*c_.electrode_geometry,c_.voltage_waveform->value(state_.time),b,state_.phi,c_.elliptic,c_.open_boundary);else solve_potential(state_.rho,c_.background_field,b,state_.phi,c_.elliptic,c_.open_boundary);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){auto dr=[&](int a,int b){return(state_.phi(b,j)-state_.phi(a,j))/((b-a)*g_.dr());};auto dz=[&](int a,int b){return(state_.phi(i,b)-state_.phi(i,a))/((b-a)*g_.dz());};state_.er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));state_.ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));state_.emag(i,j)=std::hypot(state_.er(i,j),state_.ez(i,j));}
}
void StreamerSolver::fill_electrode_diagnostics(StreamerDiagnostics&diag)const{
 if(!c_.electrode_geometry||!c_.voltage_waveform)return;
 auto e=evaluate_electrode_diagnostics(state_.phi,state_.er,state_.ez,*c_.electrode_geometry,state_.time,c_.voltage_waveform->value(state_.time),last_poisson_iterations());
 diag.applied_voltage=e.applied_voltage_V;diag.phi_hv_residual=e.phi_hv_residual_V;diag.phi_ground_residual=e.phi_ground_residual_V;diag.geometry_id=e.geometry_id;diag.gap=e.gap_m;diag.tip_radius=e.tip_radius_m;
}
TimeStepLimits StreamerSolver::timestep_limits()const{
 double vr=0,vz=0,dmax=0,numax=0,cond=0,loss=0;
 for(std::size_t k=0;k<state_.ne.values().size();++k){auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);vr=std::max(vr,q.mobility*std::abs(state_.er.values()[k]));vz=std::max(vz,q.mobility*std::abs(state_.ez.values()[k]));dmax=std::max(dmax,q.diffusion);numax=std::max(numax,q.ionization_frequency);cond=std::max(cond,qe*q.mobility*state_.ne.values()[k]);loss=std::max(loss,q.attachment_two_body_frequency+q.attachment_three_body_frequency+beta_ep(q.mobility,q.diffusion)*state_.np.values()[k]);}
 TimeStepLimits x;x.drift=.5*std::min(g_.dr()/std::max(vr,1e-300),g_.dz()/std::max(vz,1e-300));x.diffusion=.5/(2*std::max(dmax,1e-300)*(1/(g_.dr()*g_.dr())+1/(g_.dz()*g_.dz())));x.ionization=.05/std::max(numax,1e-300);x.dielectric=.2*eps0/std::max(cond,1e-300);x.reaction=.2/std::max(loss,1e-300);x.selected=x.drift;x.controller="drift";for(auto p:{std::pair{x.diffusion,"diffusion"},std::pair{x.ionization,"ionization"},std::pair{x.dielectric,"dielectric"},std::pair{x.reaction,"reaction"}})if(p.first<x.selected){x.selected=p.first;x.controller=p.second;}return x;
}
bool StreamerSolver::step(double dt,StreamerDiagnostics&diag){
 ScalarField2D emission(g_);for(std::size_t k=0;k<emission.values().size();++k){auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);emission.values()[k]=c_.excitation_ratio*q.ionization_frequency*state_.ne.values()[k];}
 if(c_.photoionization){PoissonBoundaryConfig b;solve_photoionization(emission,{},b,state_.sph,c_.elliptic,c_.sp3_boundary,nullptr,&sp3_warm_);}else std::fill(state_.sph.values().begin(),state_.sph.values().end(),0);
 ScalarField2D dne(g_),dnp(g_),dnn(g_);double old_ne=0,old_q=0,charge_inventory=0,re_int=0,boundary_out=0;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);auto s=evaluate_reactions(state_.ne(i,j),state_.np(i,j),state_.nn(i,j),state_.sph(i,j),q,c_.temperature);double div=0;
  if(i<g_.nr()-1){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i+1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i+1,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i,j),vr=-r.mobility*state_.er(i+1,j);div+=g_.radial_face(i+1)*isg0_flux(state_.ne(i,j),state_.ne(i+1,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref)/(g_.r(i)*g_.dr());}
  else{double v=-q.mobility*state_.er(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());boundary_out+=flux*2*3.14159265358979323846*g_.radial_face(i+1)*g_.dz();}
  if(i>0){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i-1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i-1,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i-1,j),vr=-r.mobility*state_.er(i,j);div-=g_.radial_face(i)*isg0_flux(state_.ne(i-1,j),state_.ne(i,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref)/(g_.r(i)*g_.dr());}
  if(j<g_.nz()-1){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j+1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j+1),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j),vr=-r.mobility*state_.ez(i,j+1);div+=isg0_flux(state_.ne(i,j),state_.ne(i,j+1),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref)/g_.dz();}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=flux/g_.dz();boundary_out+=flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  if(j>0){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j-1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j-1),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j-1),vr=-r.mobility*state_.ez(i,j);div-=isg0_flux(state_.ne(i,j-1),state_.ne(i,j),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref)/g_.dz();}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::min(v,0.0)*state_.ne(i,j);div-=flux/g_.dz();boundary_out+=-flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  dne(i,j)=-div+s.electron;dnp(i,j)=s.positive_ion;dnn(i,j)=s.negative_ion;double vol=g_.cell_volume(i);old_ne+=state_.ne(i,j)*vol;old_q+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;charge_inventory+=qe*(state_.np(i,j)+state_.ne(i,j)+state_.nn(i,j))*vol;re_int+=s.electron*vol;
 }
 double tol=1e-12*std::max(*std::max_element(state_.ne.values().begin(),state_.ne.values().end()),c_.n_ref);
 for(std::size_t k=0;k<state_.ne.values().size();++k){double ne=state_.ne.values()[k]+dt*dne.values()[k],np=state_.np.values()[k]+dt*dnp.values()[k],nn=state_.nn.values()[k]+dt*dnn.values()[k];if(ne<-tol||np<-tol||nn<-tol)return false;state_.ne.values()[k]=std::max(0.,ne);state_.np.values()[k]=std::max(0.,np);state_.nn.values()[k]=std::max(0.,nn);}
 state_.time+=dt;fields();diag={};diag.time=state_.time;diag.dt=dt;diag.poisson_iterations=last_poisson_iterations();diag.sp3_ksp_iterations=c_.photoionization?last_sp3_ksp_iterations():0;diag.sp3_boundary_iterations=c_.photoionization?last_sp3_boundary_iterations():0;diag.emax=*std::max_element(state_.emag.values().begin(),state_.emag.values().end());diag.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());diag.np_max=*std::max_element(state_.np.values().begin(),state_.np.values().end());diag.nn_max=*std::max_element(state_.nn.values().begin(),state_.nn.values().end());diag.controller=timestep_limits().controller;fill_electrode_diagnostics(diag);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){double vol=g_.cell_volume(i);diag.total_electrons+=state_.ne(i,j)*vol;diag.total_charge+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;}
 double ne_expected=old_ne+dt*(re_int-boundary_out),q_expected=old_q+qe*dt*boundary_out;double re=std::abs(diag.total_electrons-ne_expected)/std::max({std::abs(diag.total_electrons),std::abs(ne_expected),1.0}),rq=std::abs(diag.total_charge-q_expected)/std::max(charge_inventory,1e-30);diag.conservation_residual=std::max(re,rq);return true;
}
}
