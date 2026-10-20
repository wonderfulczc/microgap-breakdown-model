#include "streamer_rf/streamer/StreamerSolver.hpp"
#include "streamer_rf/transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <memory>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <chrono>
namespace streamer_rf::streamer {
namespace {constexpr double qe=1.602176634e-19,eps0=8.8541878128e-12;}
StreamerSolver::StreamerSolver(const AxisymmetricGrid&g,StreamerConfig c,StreamerPerformanceOptions performance):performance_(performance),g_(g),c_(c),state_(g){
 c_.cold_physics.validate(c_.geometry_view!=nullptr,c_.photoionization);
 if(c_.geometry_view&&c_.photoionization&&c_.cold_physics.photon_boundary==PhotonMetalBoundary::OpaqueAbsorbingLarsen&&c_.sp3_boundary.zero_dirichlet)
  throw std::invalid_argument("GEOMETRY_SP3_BOUNDARY_CONFLICT");
 if((performance_.fixed_weight_cache||performance_.fixed_poisson_workspace)&&!c_.geometry_view)throw std::invalid_argument("fixed caches require immutable generic geometry");
 if(performance_.fixed_weight_cache&&(c_.poisson_boundary.r_outer.kind==BoundaryKind::OpenCharge||c_.poisson_boundary.z_lower.kind==BoundaryKind::OpenCharge||c_.poisson_boundary.z_upper.kind==BoundaryKind::OpenCharge))throw std::invalid_argument("weight cache requires explicit fixed boundary conditions");
 if(c_.electrode_geometry&&c_.geometry_view)throw std::invalid_argument("ambiguous legacy and generic geometry entries");
 if(c_.geometry_view){
  if(&c_.geometry_view->grid()!=&g)throw std::invalid_argument("geometry/grid identity mismatch");
  if(!c_.voltage_waveform)throw std::invalid_argument("generic geometry requires voltage waveform");
  if(c_.photoionization&&c_.cold_physics.photon_boundary!=PhotonMetalBoundary::OpaqueAbsorbingLarsen)throw std::invalid_argument("dual-geometry SP3 photon boundary not validated; use SP3 OFF");
 }
}
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
ElectrodeCellType StreamerSolver::cell_type(int i,int j)const{return c_.geometry_view?c_.geometry_view->classify(i,j):(electrode_mode()?c_.electrode_geometry->classify(g_,i,j):ElectrodeCellType::Gas);}
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
 if(!c_.electrode_geometry&&!c_.geometry_view)throw std::runtime_error("tip-relative seed requires electrode geometry");
 initialize_gaussian(n0,sigma,c_.geometry_view?c_.geometry_view->seed_reference_z_m()+z_offset_from_tip:c_.electrode_geometry->seed_z_from_tip_offset(z_offset_from_tip));
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
 PoissonBoundaryConfig b=c_.poisson_boundary;if(performance_.fixed_poisson_workspace){if(!poisson_workspace_)poisson_workspace_=std::make_unique<FixedElectrodePoisson>(*c_.geometry_view,b,c_.elliptic);poisson_workspace_->solve(state_.rho,c_.voltage_waveform->value(state_.time),state_.phi);}else if(c_.geometry_view)solve_potential_with_electrodes(state_.rho,*c_.geometry_view,c_.voltage_waveform->value(state_.time),b,state_.phi,c_.elliptic,c_.open_boundary);else if(c_.electrode_geometry&&c_.voltage_waveform)solve_potential_with_electrodes(state_.rho,*c_.electrode_geometry,c_.voltage_waveform->value(state_.time),b,state_.phi,c_.elliptic,c_.open_boundary);else solve_potential(state_.rho,c_.background_field,b,state_.phi,c_.elliptic,c_.open_boundary);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){auto dr=[&](int a,int b){return(state_.phi(b,j)-state_.phi(a,j))/((b-a)*g_.dr());};auto dz=[&](int a,int b){return(state_.phi(i,b)-state_.phi(i,a))/((b-a)*g_.dz());};state_.er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));state_.ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));state_.emag(i,j)=std::hypot(state_.er(i,j),state_.ez(i,j));}
}
void StreamerSolver::fill_electrode_diagnostics(StreamerDiagnostics&diag)const{
 if(c_.geometry_view){
  diag.applied_voltage=c_.voltage_waveform->value(state_.time);
  diag.geometry_id=c_.geometry_view->geometry_id();diag.gap=c_.geometry_view->net_gap_m();diag.tip_radius=c_.geometry_view->tip_radius_m();
  for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
   const auto t=cell_type(i,j);
   if(t==ElectrodeCellType::HighVoltageElectrode)diag.phi_hv_residual=std::max(diag.phi_hv_residual,std::abs(state_.phi(i,j)-diag.applied_voltage));
   else if(t==ElectrodeCellType::GroundElectrode)diag.phi_ground_residual=std::max(diag.phi_ground_residual,std::abs(state_.phi(i,j)));
  }
  if(c_.geometry_view->family()!="NEEDLE_PLANE"){
   diag.head_status="NOT_APPLICABLE";diag.bridge_status="NOT_APPLICABLE";
  }
  return;
 }
 if(!c_.electrode_geometry||!c_.voltage_waveform)return;
 auto e=evaluate_electrode_diagnostics(state_.phi,state_.er,state_.ez,*c_.electrode_geometry,state_.time,c_.voltage_waveform->value(state_.time),last_poisson_iterations());
 diag.applied_voltage=e.applied_voltage_V;diag.phi_hv_residual=e.phi_hv_residual_V;diag.phi_ground_residual=e.phi_ground_residual_V;diag.geometry_id=e.geometry_id;diag.gap=e.gap_m;diag.tip_radius=e.tip_radius_m;
}
double StreamerSolver::head_position()const{
 if(c_.geometry_view&&c_.geometry_view->family()!="NEEDLE_PLANE")return std::numeric_limits<double>::quiet_NaN();
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
 if(c_.geometry_view&&c_.geometry_view->family()!="NEEDLE_PLANE")return false;
 bool in_gap=false;
 for(int j=0;j<g_.nz();++j){
  if(!is_gas(0,j))continue;
  const double z=g_.z(j);
  if(c_.geometry_view){if(z<=c_.geometry_view->return_surface_z_m()||z>=c_.geometry_view->seed_reference_z_m())continue;}
  else if(z<=c_.electrode_geometry->ground_z_m()+c_.electrode_geometry->ground_thickness_m()||z>=c_.electrode_geometry->tip_z_m())continue;
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
ReactionSourceDecompositionDiagnostics StreamerSolver::reaction_source_decomposition_diagnostics()const{
 ReactionSourceDecompositionDiagnostics d;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j))continue;
  ++d.gas_cells;
  const auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);
  const auto s=electron_reaction_source_components(state_.ne(i,j),state_.np(i,j),state_.nn(i,j),state_.sph(i,j),q,c_.temperature);
  const double vol=g_.cell_volume(i);
  d.impact_rate_s1+=s.impact_source_m3s*vol;
  d.photo_rate_s1+=s.photo_source_m3s*vol;
  d.attach2_rate_s1+=s.attachment2_loss_m3s*vol;
  d.attach3_rate_s1+=s.attachment3_loss_m3s*vol;
  d.recomb_e_rate_s1+=s.electron_recombination_loss_m3s*vol;
  d.net_electron_reaction_rate_s1+=s.net_electron_reaction_source_m3s*vol;
  d.local_source_closure_max_abs_m3s=std::max(d.local_source_closure_max_abs_m3s,std::abs(s.algebraic_closure_residual_m3s));
 }
 const double component_net=d.impact_rate_s1+d.photo_rate_s1-d.attach2_rate_s1-d.attach3_rate_s1-d.recomb_e_rate_s1;
 d.source_closure_abs_s1=std::abs(component_net-d.net_electron_reaction_rate_s1);
 const double denom=std::abs(d.impact_rate_s1)+std::abs(d.photo_rate_s1)+std::abs(d.attach2_rate_s1)+std::abs(d.attach3_rate_s1)+std::abs(d.recomb_e_rate_s1)+std::abs(d.net_electron_reaction_rate_s1);
 d.source_closure_rel=denom>0.0?d.source_closure_abs_s1/denom:(d.source_closure_abs_s1==0.0?0.0:std::numeric_limits<double>::infinity());
 return d;
}
ElectronTransportCurrentSource StreamerSolver::electron_transport_current_source()const{
 ElectronTransportCurrentSource out(g_);
 auto radial_face_flux=[&](int left_i,int right_i,int j)->double{
  const bool left_gas=left_i>=0&&left_i<g_.nr()&&is_gas(left_i,j);
  const bool right_gas=right_i>=0&&right_i<g_.nr()&&is_gas(right_i,j);
  if(left_gas&&right_gas){
   auto f=evaluate_morrow_lowke(.5*(state_.emag(left_i,j)+state_.emag(right_i,j)),c_.neutral_density,c_.pressure,c_.temperature);
   auto l=evaluate_morrow_lowke(state_.emag(left_i,j),c_.neutral_density,c_.pressure,c_.temperature);
   auto r=evaluate_morrow_lowke(state_.emag(right_i,j),c_.neutral_density,c_.pressure,c_.temperature);
   return gas_gas_electron_flux(state_.ne(left_i,j),state_.ne(right_i,j),-l.mobility*state_.er(left_i,j),-r.mobility*state_.er(right_i,j),f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref);
  }
  if(left_gas&&!right_gas){
   auto q=evaluate_morrow_lowke(state_.emag(left_i,j),c_.neutral_density,c_.pressure,c_.temperature);
   const double v=-q.mobility*state_.er(left_i,j);
   return right_i<g_.nr()?absorbing_electrode_flux(state_.ne(left_i,j),v,q.diffusion,g_.dr(),true):std::max(v,0.0)*state_.ne(left_i,j);
  }
  if(!left_gas&&right_gas){
   auto q=evaluate_morrow_lowke(state_.emag(right_i,j),c_.neutral_density,c_.pressure,c_.temperature);
   const double v=-q.mobility*state_.er(right_i,j);
   return left_i>=0?absorbing_electrode_flux(state_.ne(right_i,j),v,q.diffusion,g_.dr(),false):0.0;
  }
  return 0.0;
 };
 auto axial_face_flux=[&](int i,int lower_j,int upper_j)->double{
  const bool lower_gas=lower_j>=0&&lower_j<g_.nz()&&is_gas(i,lower_j);
  const bool upper_gas=upper_j>=0&&upper_j<g_.nz()&&is_gas(i,upper_j);
  if(lower_gas&&upper_gas){
   auto f=evaluate_morrow_lowke(.5*(state_.emag(i,lower_j)+state_.emag(i,upper_j)),c_.neutral_density,c_.pressure,c_.temperature);
   auto l=evaluate_morrow_lowke(state_.emag(i,lower_j),c_.neutral_density,c_.pressure,c_.temperature);
   auto r=evaluate_morrow_lowke(state_.emag(i,upper_j),c_.neutral_density,c_.pressure,c_.temperature);
   return gas_gas_electron_flux(state_.ne(i,lower_j),state_.ne(i,upper_j),-l.mobility*state_.ez(i,lower_j),-r.mobility*state_.ez(i,upper_j),f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref);
  }
  if(lower_gas&&!upper_gas){
   auto q=evaluate_morrow_lowke(state_.emag(i,lower_j),c_.neutral_density,c_.pressure,c_.temperature);
   const double v=-q.mobility*state_.ez(i,lower_j);
   return upper_j<g_.nz()?absorbing_electrode_flux(state_.ne(i,lower_j),v,q.diffusion,g_.dz(),true):std::max(v,0.0)*state_.ne(i,lower_j);
  }
  if(!lower_gas&&upper_gas){
   auto q=evaluate_morrow_lowke(state_.emag(i,upper_j),c_.neutral_density,c_.pressure,c_.temperature);
   const double v=-q.mobility*state_.ez(i,upper_j);
   return lower_j>=0?absorbing_electrode_flux(state_.ne(i,upper_j),v,q.diffusion,g_.dz(),false):std::min(v,0.0)*state_.ne(i,upper_j);
  }
  return 0.0;
 };
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j)){out.jr(i,j)=0.0;out.jz(i,j)=0.0;continue;}
  const double gamma_r_minus=radial_face_flux(i-1,i,j);
  const double gamma_r_plus=radial_face_flux(i,i+1,j);
  const double gamma_z_minus=axial_face_flux(i,j-1,j);
  const double gamma_z_plus=axial_face_flux(i,j,j+1);
  out.jr(i,j)=-0.5*qe*(gamma_r_minus+gamma_r_plus);
  out.jz(i,j)=-0.5*qe*(gamma_z_minus+gamma_z_plus);
  const double vol=g_.cell_volume(i);
  out.current_moment_r+=out.jr(i,j)*vol;
  out.current_moment_z+=out.jz(i,j)*vol;
  out.integral_abs_jz+=std::abs(out.jz(i,j))*vol;
  out.max_abs_jz=std::max(out.max_abs_jz,std::abs(out.jz(i,j)));
  out.max_abs_j=std::max({out.max_abs_j,std::abs(out.jr(i,j)),std::abs(out.jz(i,j))});
 }
 return out;
}
double StreamerSolver::vacuum_gap_capacitance()const{
 if(!electrode_mode())return 0.0;
 if(std::isfinite(c_gap_vacuum_cache_))return c_gap_vacuum_cache_;
 ScalarField2D rho(g_),phi(g_),er(g_),ez(g_);PoissonBoundaryConfig b;b.r_outer.kind=b.z_lower.kind=b.z_upper.kind=BoundaryKind::OpenCharge;if(c_.geometry_view)solve_potential_with_electrodes(rho,*c_.geometry_view,1.0,b,phi,c_.elliptic,c_.open_boundary);else solve_potential_with_electrodes(rho,*c_.electrode_geometry,1.0,b,phi,c_.elliptic,c_.open_boundary);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){auto dr=[&](int a,int b){return(phi(b,j)-phi(a,j))/((b-a)*g_.dr());};auto dz=[&](int a,int b){return(phi(i,b)-phi(i,a))/((b-a)*g_.dz());};er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));}
 double q_hv=0.0;auto add=[&](ElectrodeCellType t,double en,double area){if(t==ElectrodeCellType::HighVoltageElectrode)q_hv+=eps0*en*area;};
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(cell_type(i,j)!=ElectrodeCellType::Gas)continue;if(i<g_.nr()-1&&cell_type(i+1,j)!=ElectrodeCellType::Gas)add(cell_type(i+1,j),-er(i,j),axisymmetric_radial_face_area(g_,i+1));if(i>0&&cell_type(i-1,j)!=ElectrodeCellType::Gas)add(cell_type(i-1,j),er(i,j),axisymmetric_radial_face_area(g_,i));if(j<g_.nz()-1&&cell_type(i,j+1)!=ElectrodeCellType::Gas)add(cell_type(i,j+1),-ez(i,j),axisymmetric_axial_face_area(g_,i));if(j>0&&cell_type(i,j-1)!=ElectrodeCellType::Gas)add(cell_type(i,j-1),ez(i,j),axisymmetric_axial_face_area(g_,i));}
 c_gap_vacuum_cache_=std::abs(q_hv);return c_gap_vacuum_cache_;
}
double StreamerSolver::conjugate_hv_charge()const{
 if(!electrode_mode())return std::numeric_limits<double>::quiet_NaN();
 return conjugate_hv_charge_for(state_.phi,c_.voltage_waveform->value(state_.time));
}
double StreamerSolver::conjugate_hv_charge_for(const ScalarField2D& phi,double v)const{
 constexpr double pi=3.14159265358979323846;
 double q=0;
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j))continue;
  if(i+1<g_.nr()&&cell_type(i+1,j)==ElectrodeCellType::HighVoltageElectrode)q+=eps0*(v-phi(i,j))*2*pi*g_.radial_face(i+1)*g_.dz()/g_.dr();
  if(i>0&&cell_type(i-1,j)==ElectrodeCellType::HighVoltageElectrode)q+=eps0*(v-phi(i,j))*2*pi*g_.radial_face(i)*g_.dz()/g_.dr();
  const double area=pi*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));
  if(j+1<g_.nz()&&cell_type(i,j+1)==ElectrodeCellType::HighVoltageElectrode)q+=eps0*(v-phi(i,j))*area/g_.dz();
  if(j>0&&cell_type(i,j-1)==ElectrodeCellType::HighVoltageElectrode)q+=eps0*(v-phi(i,j))*area/g_.dz();
 }
 return q;
}
double StreamerSolver::conjugate_vacuum_capacitance()const{
 if(!electrode_mode())return std::numeric_limits<double>::quiet_NaN();
 if(std::isfinite(c_gap_conjugate_cache_))return c_gap_conjugate_cache_;
 ScalarField2D rho(g_),phi(g_);
 if(c_.geometry_view)solve_potential_with_electrodes(rho,*c_.geometry_view,1.0,c_.poisson_boundary,phi,c_.elliptic,c_.open_boundary);
 else solve_potential_with_electrodes(rho,*c_.electrode_geometry,1.0,c_.poisson_boundary,phi,c_.elliptic,c_.open_boundary);
 c_gap_conjugate_cache_=conjugate_hv_charge_for(phi,1.0);return c_gap_conjugate_cache_;
}
void StreamerSolver::reset_electrode_history(){if(!electrode_mode()){has_electrode_history_=false;return;}auto q=electrode_surface_diagnostics();last_q_hv_=q.q_hv;last_q_ground_=q.q_ground;has_electrode_history_=true;}
void StreamerSolver::refresh_electrostatic_fields(){fields();}
double StreamerSolver::electrode_displacement_current(double q1,double q0,double dt)const{return dt>0?(q1-q0)/dt:std::numeric_limits<double>::quiet_NaN();}
StreamerDiagnostics StreamerSolver::sample_terminal_diagnostics_from_history(double dt){
 StreamerDiagnostics diag{};diag.time=state_.time;diag.dt=dt;diag.poisson_iterations=last_poisson_iterations();diag.emax=*std::max_element(state_.emag.values().begin(),state_.emag.values().end());diag.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());diag.np_max=*std::max_element(state_.np.values().begin(),state_.np.values().end());diag.nn_max=*std::max_element(state_.nn.values().begin(),state_.nn.values().end());fill_electrode_diagnostics(diag);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;double vol=g_.cell_volume(i);diag.total_electrons+=state_.ne(i,j)*vol;diag.total_charge+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);diag.sigma_max=std::max(diag.sigma_max,qe*q.mobility*state_.ne(i,j));}
 if(electrode_mode()){auto qsurf=electrode_surface_diagnostics();diag.q_hv=qsurf.q_hv;diag.q_ground=qsurf.q_ground;if(c_.geometry_view){diag.q_conjugate_hv=conjugate_hv_charge();diag.c_gap_conjugate=conjugate_vacuum_capacitance();}diag.displacement_current_valid=has_electrode_history_&&dt>0;diag.c_gap_vacuum=vacuum_gap_capacitance();if(diag.displacement_current_valid){diag.i_disp_hv=electrode_displacement_current(qsurf.q_hv,last_q_hv_,dt);diag.i_disp_ground=electrode_displacement_current(qsurf.q_ground,last_q_ground_,dt);diag.i_total_hv=diag.i_cond_hv+diag.i_disp_hv;diag.i_total_ground=diag.i_cond_ground+diag.i_disp_ground;}last_q_hv_=qsurf.q_hv;last_q_ground_=qsurf.q_ground;has_electrode_history_=true;}
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
 TimeStepLimits x;x.drift=.5*std::min(g_.dr()/std::max(vr,1e-300),g_.dz()/std::max(vz,1e-300));x.diffusion=.5/(2*std::max(dmax,1e-300)*(1/(g_.dr()*g_.dr())+1/(g_.dz()*g_.dz())));x.ionization=.05/std::max(numax,1e-300);x.dielectric=.2*eps0/std::max(cond,1e-300);x.reaction=reaction;x.selected=x.drift;x.controller="drift";for(auto p:{std::pair{x.diffusion,"diffusion"},std::pair{x.ionization,"ionization"},std::pair{x.dielectric,"dielectric"},std::pair{x.reaction,"reaction"}})if(p.first<x.selected){x.selected=p.first;x.controller=p.second;}
 if(c_.cold_physics.positive_ion_transport){
  const auto ions=positive_ion_transport(*c_.geometry_view,state_.np,state_.er,state_.ez,
      c_.cold_physics.ions,c_.cold_physics.surface,c_.voltage_waveform->value(state_.time),
      c_.neutral_density,c_.temperature,c_.pressure,false);
  if(ions.positivity_dt_s<x.selected){x.selected=ions.positivity_dt_s;x.controller="positive_ion_transport";}
  // Reserve a separate bound for reactions; their losses add to transport losses.
  x.selected=std::min(x.selected,.5*ions.positivity_dt_s);
 }
 return x;
}
bool StreamerSolver::step(double dt,StreamerDiagnostics&diag){
 if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("INVALID_NATIVE_STEP_DT");
 using ProfileClock=std::chrono::steady_clock;
 auto profile_mark=performance_.profile?ProfileClock::now():ProfileClock::time_point{};
 StreamerStepTiming profile;
 auto mark=[&](){if(!performance_.profile)return 0.0;const auto now=ProfileClock::now();double seconds=std::chrono::duration<double>(now-profile_mark).count();profile_mark=now;return seconds;};
 std::unique_ptr<AcceptedStepRecord> record;
 std::vector<double> raw_ne,raw_np,raw_nn;
 if(observer_){record=std::make_unique<AcceptedStepRecord>(g_);raw_ne=state_.ne.values();raw_np=state_.np.values();raw_nn=state_.nn.values();record->pre_time_s=state_.time;record->dt_s=dt;record->geometry_id=c_.geometry_view?c_.geometry_view->geometry_id():(c_.electrode_geometry?c_.electrode_geometry->metadata().geometry_id:"NONE");record->geometry_revision=c_.geometry_view?c_.geometry_view->revision():"LEGACY";record->reference_plane_id=c_.geometry_view?c_.geometry_view->reference_plane_id():"LEGACY_HV_TO_GROUND";std::ostringstream mesh;mesh<<g_.nr()<<'x'<<g_.nz()<<':'<<std::setprecision(17)<<g_.dr()<<':'<<g_.dz()<<':'<<g_.z(0);record->mesh_id=mesh.str();record->pre_voltage_V=electrode_mode()?c_.voltage_waveform->value(state_.time):0;}
 enforce_plasma_mask();
 enforce_stage_c_activity_floor();
 if(record){record->pre=state_;const auto n=g_.size();record->pre_floor_ne_delta_m3.resize(n);record->pre_floor_np_delta_m3.resize(n);record->pre_floor_nn_delta_m3.resize(n);record->impact_m3s.resize(n);record->photo_m3s.resize(n);record->attachment2_m3s.resize(n);record->attachment3_m3s.resize(n);record->recombination_m3s.resize(n);record->clip_ne_delta_m3.assign(n,0);record->clip_np_delta_m3.assign(n,0);record->clip_nn_delta_m3.assign(n,0);for(std::size_t k=0;k<n;++k){record->pre_floor_ne_delta_m3[k]=state_.ne.values()[k]-raw_ne[k];record->pre_floor_np_delta_m3[k]=state_.np.values()[k]-raw_np[k];record->pre_floor_nn_delta_m3[k]=state_.nn.values()[k]-raw_nn[k];record->pre.rho.values()[k]=qe*(record->pre.np.values()[k]-record->pre.ne.values()[k]-record->pre.nn.values()[k]);}}
 if(record){
  record->raw_pre_ne_m3=raw_ne;record->raw_pre_np_m3=raw_np;record->raw_pre_nn_m3=raw_nn;
  record->ion_recombination_m3s.resize(g_.size());record->reaction_e_m3s.resize(g_.size());
  record->reaction_p_m3s.resize(g_.size());record->reaction_n_m3s.resize(g_.size());
  record->reaction_charge_defect_m3s.resize(g_.size());
 }
 auto capture_face=[&](int i,int j,FaceDirection2D dir,double gamma_positive){
  if(!record)return;
  int ni=i,nj=j;double sign=1,area=0;constexpr double pi=3.14159265358979323846;
  switch(dir){
   case FaceDirection2D::RadialMinus:--ni;sign=-1;area=2*pi*g_.radial_face(i)*g_.dz();break;
   case FaceDirection2D::RadialPlus:++ni;area=2*pi*g_.radial_face(i+1)*g_.dz();break;
   case FaceDirection2D::AxialMinus:--nj;sign=-1;area=axisymmetric_axial_face_area(g_,i);break;
   case FaceDirection2D::AxialPlus:++nj;area=axisymmetric_axial_face_area(g_,i);break;
  }
  const bool inside=ni>=0&&ni<g_.nr()&&nj>=0&&nj<g_.nz();
  const auto neighbor=inside?cell_type(ni,nj):ElectrodeCellType::Gas;
  const auto role=!inside?(ni<0?FaceRole2D::Axis:FaceRole2D::OuterBoundary):
     (neighbor==ElectrodeCellType::Gas?FaceRole2D::GasInterface:FaceRole2D::Electrode);
  record->electron_faces.push_back({i,j,ni,nj,dir,role,neighbor,sign*gamma_positive,area});
 };
 const bool had_electrode_history=has_electrode_history_;
 const double prev_q_hv=last_q_hv_,prev_q_ground=last_q_ground_;
 ScalarField2D emission(g_);for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;if(!is_gas(i,j)){emission.values()[k]=0;continue;}auto q=evaluate_morrow_lowke(state_.emag.values()[k],c_.neutral_density,c_.pressure,c_.temperature);emission.values()[k]=c_.excitation_ratio*q.ionization_frequency*state_.ne.values()[k];}
 GeometrySp3Audit photon_audit;
 if(c_.photoionization){
  if(c_.geometry_view&&c_.cold_physics.photon_boundary==PhotonMetalBoundary::OpaqueAbsorbingLarsen){
   if(!photon_workspace_)photon_workspace_=std::make_shared<GeometrySp3Workspace>(*c_.geometry_view,PressureConfig{},c_.elliptic);
   photon_audit=photon_workspace_->solve(emission,state_.sph,&sp3_warm_);
  }else{PoissonBoundaryConfig b;solve_photoionization(emission,{},b,state_.sph,c_.elliptic,c_.sp3_boundary,nullptr,&sp3_warm_);if(electrode_mode())for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i)if(!is_gas(i,j))state_.sph(i,j)=0;}
 }else std::fill(state_.sph.values().begin(),state_.sph.values().end(),0);
 std::unique_ptr<IonTransportResult> ions;
 if(c_.cold_physics.positive_ion_transport){
  ions=std::make_unique<IonTransportResult>(positive_ion_transport(*c_.geometry_view,state_.np,state_.er,state_.ez,
      c_.cold_physics.ions,c_.cold_physics.surface,c_.voltage_waveform->value(state_.time),
      c_.neutral_density,c_.temperature,c_.pressure));
  if(dt>ions->positivity_dt_s)throw std::invalid_argument("ION_POSITIVITY_TIMESTEP_EXCEEDED");
 }
 if(record)record->pre.sph=state_.sph;
 ScalarField2D dne(g_),dnp(g_),dnn(g_);double old_ne=0,old_q=0,charge_inventory=0,re_int=0,boundary_out=0,hv_abs=0,ground_abs=0;
 profile.prelude_s=mark();
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  if(!is_gas(i,j)){dne(i,j)=dnp(i,j)=dnn(i,j)=0;continue;}
  auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);auto s=evaluate_reactions(state_.ne(i,j),state_.np(i,j),state_.nn(i,j),state_.sph(i,j),q,c_.temperature);double div=0;
  if(record){
   const auto k=static_cast<std::size_t>(j)*g_.nr()+i;
   record->ion_recombination_m3s[k]=s.ion_recombination;
   record->reaction_e_m3s[k]=s.electron;record->reaction_p_m3s[k]=s.positive_ion;
   record->reaction_n_m3s[k]=s.negative_ion;record->reaction_charge_defect_m3s[k]=s.charge_balance;
  }
  if(record){const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;const auto parts=electron_reaction_source_components(state_.ne(i,j),state_.np(i,j),state_.nn(i,j),state_.sph(i,j),q,c_.temperature);record->impact_m3s[k]=parts.impact_source_m3s;record->photo_m3s[k]=parts.photo_source_m3s;record->attachment2_m3s[k]=parts.attachment2_loss_m3s;record->attachment3_m3s[k]=parts.attachment3_loss_m3s;record->recombination_m3s[k]=parts.electron_recombination_loss_m3s;}
  if(i<g_.nr()-1&&is_gas(i+1,j)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i+1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i+1,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i,j),vr=-r.mobility*state_.er(i+1,j);const double flux=gas_gas_electron_flux(state_.ne(i,j),state_.ne(i+1,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());capture_face(i,j,FaceDirection2D::RadialPlus,flux);}
  else if(i<g_.nr()-1){double v=-q.mobility*state_.er(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dr(),true);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());capture_face(i,j,FaceDirection2D::RadialPlus,flux);double loss=flux*2*3.14159265358979323846*g_.radial_face(i+1)*g_.dz();if(cell_type(i+1,j)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.er(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=g_.radial_face(i+1)*flux/(g_.r(i)*g_.dr());capture_face(i,j,FaceDirection2D::RadialPlus,flux);boundary_out+=flux*2*3.14159265358979323846*g_.radial_face(i+1)*g_.dz();}
  if(i>0&&is_gas(i-1,j)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i-1,j)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i-1,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.er(i-1,j),vr=-r.mobility*state_.er(i,j);const double flux=gas_gas_electron_flux(state_.ne(i-1,j),state_.ne(i,j),vl,vr,f.diffusion,g_.dr(),c_.isg_epsilon,c_.n_ref);div-=g_.radial_face(i)*flux/(g_.r(i)*g_.dr());capture_face(i,j,FaceDirection2D::RadialMinus,flux);}
  else if(i>0){double v=-q.mobility*state_.er(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dr(),false);div-=g_.radial_face(i)*flux/(g_.r(i)*g_.dr());capture_face(i,j,FaceDirection2D::RadialMinus,flux);double loss=-flux*2*3.14159265358979323846*g_.radial_face(i)*g_.dz();if(cell_type(i-1,j)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  if(j<g_.nz()-1&&is_gas(i,j+1)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j+1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j+1),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j),vr=-r.mobility*state_.ez(i,j+1);const double flux=gas_gas_electron_flux(state_.ne(i,j),state_.ne(i,j+1),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref);div+=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialPlus,flux);}
  else if(j<g_.nz()-1){double v=-q.mobility*state_.ez(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dz(),true);div+=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialPlus,flux);double loss=flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));if(cell_type(i,j+1)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::max(v,0.0)*state_.ne(i,j);div+=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialPlus,flux);boundary_out+=flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  if(j>0&&is_gas(i,j-1)){auto f=evaluate_morrow_lowke(.5*(state_.emag(i,j)+state_.emag(i,j-1)),c_.neutral_density,c_.pressure,c_.temperature),l=evaluate_morrow_lowke(state_.emag(i,j-1),c_.neutral_density,c_.pressure,c_.temperature),r=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);double vl=-l.mobility*state_.ez(i,j-1),vr=-r.mobility*state_.ez(i,j);const double flux=gas_gas_electron_flux(state_.ne(i,j-1),state_.ne(i,j),vl,vr,f.diffusion,g_.dz(),c_.isg_epsilon,c_.n_ref);div-=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialMinus,flux);}
  else if(j>0){double v=-q.mobility*state_.ez(i,j),flux=absorbing_electrode_flux(state_.ne(i,j),v,q.diffusion,g_.dz(),false);div-=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialMinus,flux);double loss=-flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));if(cell_type(i,j-1)==ElectrodeCellType::HighVoltageElectrode)hv_abs+=loss;else ground_abs+=loss;}
  else{double v=-q.mobility*state_.ez(i,j),flux=std::min(v,0.0)*state_.ne(i,j);div-=flux/g_.dz();capture_face(i,j,FaceDirection2D::AxialMinus,flux);boundary_out+=-flux*3.14159265358979323846*(std::pow(g_.radial_face(i+1),2)-std::pow(g_.radial_face(i),2));}
  dne(i,j)=-div+s.electron;dnp(i,j)=s.positive_ion;dnn(i,j)=s.negative_ion;double vol=g_.cell_volume(i);old_ne+=state_.ne(i,j)*vol;old_q+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;charge_inventory+=qe*(state_.np(i,j)+state_.ne(i,j)+state_.nn(i,j))*vol;re_int+=s.electron*vol;
 }
 if(ions){
  for(std::size_t k=0;k<g_.size();++k){dnp.values()[k]+=ions->rhs.values()[k];dne.values()[k]+=ions->secondary_rhs.values()[k];}
  // Emission is an inward boundary flux, never a second reaction volume source.
  if(record){
   std::vector<double> secondary_by_face(4*g_.size(),0);
   for(const auto& face:ions->faces)secondary_by_face[4*(face.j*g_.nr()+face.i)+int(face.direction)]=face.gamma_secondary_in_m2s;
   for(auto& face:record->electron_faces)
    face.gamma_e_m2s-=secondary_by_face[4*(face.j*g_.nr()+face.i)+int(face.direction)];
  }
 }
 double tol=electrode_mode()?stage_c_activity_tolerance():numerical_density_tolerance();
 profile.transport_and_record_s=mark();
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
  const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;
  if(!is_gas(i,j)){state_.ne.values()[k]=state_.np.values()[k]=state_.nn.values()[k]=0;continue;}
  double ne=state_.ne.values()[k]+dt*dne.values()[k],np=state_.np.values()[k]+dt*dnp.values()[k],nn=state_.nn.values()[k]+dt*dnn.values()[k];
  if(c_.cold_physics.mode==ColdPhysicsMode::CoupledColdV3&&(!std::isfinite(ne)||!std::isfinite(np)||!std::isfinite(nn)))throw std::runtime_error("NONFINITE_NATIVE_SPECIES");
  if(ne<-tol||np<-tol||nn<-tol){
   if(c_.cold_physics.mode==ColdPhysicsMode::CoupledColdV3){
    std::ostringstream why;why<<std::setprecision(17)<<"NATIVE_POSITIVITY_REJECTED i="<<i<<" j="<<j<<" ne="<<ne<<" np="<<np<<" nn="<<nn<<" dt="<<dt<<" tolerance="<<tol;
    throw std::runtime_error(why.str());
   }
   return false;
  }
  if(record){record->clip_ne_delta_m3[k]=std::max(0.,ne)-ne;record->clip_np_delta_m3[k]=std::max(0.,np)-np;record->clip_nn_delta_m3[k]=std::max(0.,nn)-nn;}
  state_.ne.values()[k]=std::max(0.,ne);state_.np.values()[k]=std::max(0.,np);state_.nn.values()[k]=std::max(0.,nn);
 }
 std::vector<double> pre_post_floor_ne,pre_post_floor_np,pre_post_floor_nn;
 if(record){pre_post_floor_ne=state_.ne.values();pre_post_floor_np=state_.np.values();pre_post_floor_nn=state_.nn.values();}
 enforce_stage_c_activity_floor();
 if(record){const auto n=g_.size();record->post_floor_ne_delta_m3.resize(n);record->post_floor_np_delta_m3.resize(n);record->post_floor_nn_delta_m3.resize(n);for(std::size_t k=0;k<n;++k){record->post_floor_ne_delta_m3[k]=state_.ne.values()[k]-pre_post_floor_ne[k];record->post_floor_np_delta_m3[k]=state_.np.values()[k]-pre_post_floor_np[k];record->post_floor_nn_delta_m3[k]=state_.nn.values()[k]-pre_post_floor_nn[k];}}
 profile.update_s=mark();state_.time+=dt;fields();profile.fields_s=mark();diag={};diag.time=state_.time;diag.dt=dt;diag.poisson_iterations=last_poisson_iterations();diag.sp3_ksp_iterations=c_.photoionization?last_sp3_ksp_iterations():0;diag.sp3_boundary_iterations=c_.photoionization?last_sp3_boundary_iterations():0;diag.emax=*std::max_element(state_.emag.values().begin(),state_.emag.values().end());diag.ne_max=*std::max_element(state_.ne.values().begin(),state_.ne.values().end());diag.np_max=*std::max_element(state_.np.values().begin(),state_.np.values().end());diag.nn_max=*std::max_element(state_.nn.values().begin(),state_.nn.values().end());diag.controller=timestep_limits().controller;fill_electrode_diagnostics(diag);
 for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){if(!is_gas(i,j))continue;double vol=g_.cell_volume(i);diag.total_electrons+=state_.ne(i,j)*vol;diag.total_charge+=qe*(state_.np(i,j)-state_.ne(i,j)-state_.nn(i,j))*vol;auto q=evaluate_morrow_lowke(state_.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature);diag.sigma_max=std::max(diag.sigma_max,qe*q.mobility*state_.ne(i,j));}
 if(c_.geometry_view&&c_.photoionization&&c_.cold_physics.photon_boundary==PhotonMetalBoundary::OpaqueAbsorbingLarsen){
  diag.sp3_ksp_iterations=photon_audit.iterations;
  diag.sp3_boundary_iterations=0; // Direct coupled boundary block, no outer boundary iteration.
 }
 diag.absorbed_electron_hv=dt*hv_abs;diag.absorbed_electron_ground=dt*ground_abs;diag.i_cond_hv=qe*hv_abs;diag.i_cond_ground=qe*ground_abs;diag.outer_boundary_current=qe*boundary_out;
 if(ions){
  diag.i_cond_hv-=qe*(ions->hv_collected_s1+ions->hv_emitted_s1);
  diag.i_cond_ground-=qe*(ions->return_collected_s1+ions->return_emitted_s1);
  diag.outer_boundary_current-=qe*ions->outer_lost_s1;
 }
 if(electrode_mode()){auto qsurf=electrode_surface_diagnostics();diag.q_hv=qsurf.q_hv;diag.q_ground=qsurf.q_ground;if(c_.geometry_view){diag.q_conjugate_hv=conjugate_hv_charge();diag.c_gap_conjugate=conjugate_vacuum_capacitance();}diag.displacement_current_valid=had_electrode_history&&dt>0;diag.c_gap_vacuum=vacuum_gap_capacitance();if(diag.displacement_current_valid){diag.i_disp_hv=electrode_displacement_current(qsurf.q_hv,prev_q_hv,dt);diag.i_disp_ground=electrode_displacement_current(qsurf.q_ground,prev_q_ground,dt);diag.i_total_hv=diag.i_cond_hv+diag.i_disp_hv;diag.i_total_ground=diag.i_cond_ground+diag.i_disp_ground;}last_q_hv_=qsurf.q_hv;last_q_ground_=qsurf.q_ground;has_electrode_history_=true;}
 auto gd=conductance_diagnostics(diag.applied_voltage);diag.p_cond=gd.p_cond;diag.gb=gd.gb;diag.rb=gd.rb;diag.rb_valid=gd.valid;
 diag.head_position=head_position();diag.bridge_flag=bridge_flag();if(std::isfinite(last_head_position_)&&diag.time>last_head_time_)diag.head_velocity=(diag.head_position-last_head_position_)/(diag.time-last_head_time_);last_head_position_=diag.head_position;last_head_time_=diag.time;
 const double total_loss=boundary_out+hv_abs+ground_abs;double ne_expected=old_ne+dt*(re_int-total_loss),q_expected=old_q+qe*dt*total_loss;double re=std::abs(diag.total_electrons-ne_expected)/std::max({std::abs(diag.total_electrons),std::abs(ne_expected),1.0}),rq=std::abs(diag.total_charge-q_expected)/std::max(charge_inventory,1e-30);diag.conservation_residual=std::max(re,rq);diag.plasma_charge_derivative=(diag.total_charge-old_q)/dt;const double expected_charge_rate=diag.i_cond_hv+diag.i_cond_ground+diag.outer_boundary_current;diag.current_continuity_residual=std::abs(diag.plasma_charge_derivative-expected_charge_rate)/std::max({std::abs(diag.plasma_charge_derivative),std::abs(expected_charge_rate),1e-30});
 profile.terminal_s=mark();
 if(ions){
  const double total_exchange=qe*dt*(hv_abs+ground_abs+boundary_out-
      ions->hv_collected_s1-ions->return_collected_s1-ions->outer_lost_s1-
      ions->hv_emitted_s1-ions->return_emitted_s1);
  const double electron_expected=old_ne+dt*(re_int-boundary_out-hv_abs-ground_abs+
      ions->hv_emitted_s1+ions->return_emitted_s1);
  diag.conservation_residual=std::max(std::abs(diag.total_charge-old_q-total_exchange)/std::max(charge_inventory,1e-30),
      std::abs(diag.total_electrons-electron_expected)/std::max({std::abs(diag.total_electrons),std::abs(electron_expected),1.}));
 }
 if(record){
  record->physics_mode=c_.cold_physics.mode;record->photon_boundary_audit=photon_audit;
  record->net_I_Hg_A=diag.i_cond_hv;record->net_I_Return_g_A=diag.i_cond_ground;
  if(ions){
   record->positive_ion_faces=std::move(ions->faces);
   record->ion_positivity_dt_s=ions->positivity_dt_s;
   record->hv_collected_ions=dt*ions->hv_collected_s1;record->return_collected_ions=dt*ions->return_collected_s1;
   record->outer_lost_ions=dt*ions->outer_lost_s1;
   record->hv_secondary_electrons=dt*ions->hv_emitted_s1;record->return_secondary_electrons=dt*ions->return_emitted_s1;
   record->species_current_id="ELECTRON_NET_PLUS_POSITIVE_ION_V3_NEGATIVE_ION_IMMOBILE";
  }
  record->post=state_;record->post_time_s=state_.time;record->post_voltage_V=electrode_mode()?c_.voltage_waveform->value(state_.time):0;
  record->hv_absorbed_electrons=dt*hv_abs;record->return_absorbed_electrons=dt*ground_abs;record->outer_lost_electrons=dt*boundary_out;
  record->q_native_hv_C=diag.q_hv;record->q_conjugate_hv_C=conjugate_hv_charge();record->native_capacitance_F=diag.c_gap_vacuum;record->conjugate_capacitance_F=conjugate_vacuum_capacitance();
  if(electrode_mode()){
   const auto weight_start=performance_.profile?ProfileClock::now():ProfileClock::time_point{};
   const std::string weight_key=performance_.fixed_weight_cache?electrode_operator_identity(*c_.geometry_view,c_.poisson_boundary,c_.elliptic):std::string{};
   if(weight_cache_){
    if(weight_cache_->key!=weight_key)throw std::runtime_error("fixed weight cache invalidated; recreate solver");
    record->weight_psi=weight_cache_->psi;record->weight_er=weight_cache_->er;record->weight_ez=weight_cache_->ez;
   }else{
   ScalarField2D zero(g_);
   if(c_.geometry_view)solve_potential_with_electrodes(zero,*c_.geometry_view,1.0,c_.poisson_boundary,record->weight_psi,c_.elliptic,c_.open_boundary);
   else solve_potential_with_electrodes(zero,*c_.electrode_geometry,1.0,c_.poisson_boundary,record->weight_psi,c_.elliptic,c_.open_boundary);
   for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
    auto dr=[&](int a,int b){return(record->weight_psi(b,j)-record->weight_psi(a,j))/((b-a)*g_.dr());};
    auto dz=[&](int a,int b){return(record->weight_psi(i,b)-record->weight_psi(i,a))/((b-a)*g_.dz());};
    record->weight_er(i,j)=i==0?0:-(i==g_.nr()-1?dr(i-1,i):dr(i-1,i+1));
    record->weight_ez(i,j)=-(j==0?dz(0,1):j==g_.nz()-1?dz(j-1,j):dz(j-1,j+1));
   }
   if(performance_.fixed_weight_cache){weight_cache_=std::make_unique<WeightCache>(g_);weight_cache_->key=weight_key;weight_cache_->psi=record->weight_psi;weight_cache_->er=record->weight_er;weight_cache_->ez=record->weight_ez;}
   }
   if(performance_.profile)profile.weight_s=std::chrono::duration<double>(ProfileClock::now()-weight_start).count();
   for(int j=0;j<g_.nz();++j)for(int i=0;i<g_.nr();++i){
    if(!is_gas(i,j))continue;
    const std::size_t k=static_cast<std::size_t>(j)*g_.nr()+i;const double vol=g_.cell_volume(i);
    const double ew2=std::pow(record->weight_er(i,j),2)+std::pow(record->weight_ez(i,j),2);
    const double mu_pre=evaluate_morrow_lowke(record->pre.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature).mobility;
    const double mu_post=evaluate_morrow_lowke(record->post.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature).mobility;
    const double wpre=qe*mu_pre*ew2,wpost=qe*mu_post*ew2;
    record->G_R_pre_S+=wpre*record->pre.ne(i,j)*vol;record->G_R_post_S+=wpost*record->post.ne(i,j)*vol;
    record->W_rho_pre_C+=record->weight_psi(i,j)*record->pre.rho(i,j)*vol;
    record->W_rho_post_C+=record->weight_psi(i,j)*record->post.rho(i,j)*vol;
    const double source=record->impact_m3s[k]+record->photo_m3s[k]-record->attachment2_m3s[k]-record->attachment3_m3s[k]-record->recombination_m3s[k];
    record->delta_G_reaction_S+=dt*wpre*source*vol;
    record->delta_G_mobility_S+=(wpost-wpre)*record->post.ne(i,j)*vol;
    record->delta_G_correction_S+=wpre*(record->clip_ne_delta_m3[k]+record->post_floor_ne_delta_m3[k])*vol;
   }
   for(const auto& face:record->electron_faces){
    const int i=face.i,j=face.j;
    const double ew2=std::pow(record->weight_er(i,j),2)+std::pow(record->weight_ez(i,j),2);
    const double mu=evaluate_morrow_lowke(record->pre.emag(i,j),c_.neutral_density,c_.pressure,c_.temperature).mobility;
    const double term=-dt*qe*mu*ew2*face.gamma_e_m2s*face.area_m2;
    if(face.role==FaceRole2D::GasInterface)record->delta_G_transport_S+=term;
    else record->delta_G_boundary_S+=term;
   }
   record->delta_G_closure_residual_S=(record->G_R_post_S-record->G_R_pre_S)-
     (record->delta_G_reaction_S+record->delta_G_transport_S+record->delta_G_boundary_S+
      record->delta_G_mobility_S+record->delta_G_correction_S);
  }
  profile.observer_assembly_s=mark();
  observer_->on_accepted_step(*record);
  profile.callback_s=mark();
 }
 if(performance_.profile)step_timings_.push_back(profile);
 return true;
}
}
