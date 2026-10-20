#pragma once
#include "native_full_snapshot.hpp"
#include "streamer_rf/transport.hpp"
#include <map>

namespace fp3full {
using namespace streamer_rf;
using namespace streamer_rf::streamer;
constexpr double echarge=1.602176634e-19,epsilon=8.8541878128e-12;
using Metrics=std::map<std::string,double>;

struct Edge {int a,b;double c;};
class ActiveNetwork {
 public:
  const Geometry2DView& view;
  const AxisymmetricGrid& g;
  ScalarField2D psi;
  std::vector<Edge> edges;
  std::vector<bool> active;
  double capacitance;
  explicit ActiveNetwork(const Geometry2DView& v,const StreamerConfig& cfg,const ScalarField2D* cached_psi=nullptr):view(v),g(v.grid()),psi(g),active(g.size(),false){
    for(auto b:{cfg.poisson_boundary.r_outer,cfg.poisson_boundary.z_lower,cfg.poisson_boundary.z_upper})
      if(b.kind!=BoundaryKind::Neumann||b.value!=0)throw std::invalid_argument("ACTIVE_NETWORK_REQUIRES_ZERO_NEUMANN");
    if(cached_psi)psi=*cached_psi;
    else {ScalarField2D zero(g);solve_potential_with_electrodes(zero,v,1.,cfg.poisson_boundary,psi,cfg.elliptic);}
    for(int j=1;j<g.nz()-1;++j)for(int i=0;i<g.nr()-1;++i)active[idx(i,j)]=gas(i,j);
    for(int j=1;j<g.nz()-1;++j)for(int i=0;i<g.nr()-1;++i)if(active[idx(i,j)]){
      for(auto dir:{FaceDirection2D::RadialMinus,FaceDirection2D::RadialPlus,FaceDirection2D::AxialMinus,FaceDirection2D::AxialPlus}){
        auto f=v.face(i,j,dir);if(f.role==FaceRole2D::Axis)continue;
        int a=idx(i,j),b=idx(f.neighbor_i,f.neighbor_j);
        if(f.neighbor_material==ElectrodeCellType::Gas){if(!active[b]||b<a)continue;}
        double h=(dir==FaceDirection2D::RadialMinus||dir==FaceDirection2D::RadialPlus)?g.dr():g.dz();
        edges.push_back({a,b,epsilon*f.area_m2/h});
      }
    }
    capacitance=double(charge(psi,1));
  }
  int idx(int i,int j)const{return j*g.nr()+i;}
  bool gas(int i,int j)const{return view.classify(i,j)==ElectrodeCellType::Gas;}
  double potential(const ScalarField2D& phi,int k,double voltage)const{
    auto t=view.classify(k%g.nr(),k/g.nr());
    return t==ElectrodeCellType::HighVoltageElectrode?voltage:t==ElectrodeCellType::GroundElectrode?0.:phi.values()[k];
  }
  long double charge(const ScalarField2D& phi,double v)const{
    long double q=0;for(auto edge:edges)if(view.classify(edge.b%g.nr(),edge.b/g.nr())==ElectrodeCellType::HighVoltageElectrode)
      q+=(long double)edge.c*(v-phi.values()[edge.a]);return q;
  }
  long double energy(const ScalarField2D& phi,double v)const{
    long double u=0;for(auto edge:edges){long double d=potential(phi,edge.a,v)-potential(phi,edge.b,v);u+=.5L*edge.c*d*d;}return u;
  }
  long double delta_energy(const ScalarField2D& a,double va,const ScalarField2D& b,double vb)const{
    long double u=0;for(auto edge:edges){long double x=potential(a,edge.a,va)-potential(a,edge.b,va),y=potential(b,edge.a,vb)-potential(b,edge.b,vb);u+=.5L*edge.c*(y+x)*(y-x);}return u;
  }
  double gauss(const StreamerState& s,double v)const{
    std::vector<long double> r(g.size());
    for(auto edge:edges){long double q=edge.c*(potential(s.phi,edge.a,v)-potential(s.phi,edge.b,v));r[edge.a]+=q;if(active[edge.b])r[edge.b]-=q;}
    long double result=0;for(std::size_t k=0;k<g.size();++k)if(active[k])result=std::max(result,std::abs(r[k]-s.rho.values()[k]*g.cell_volume(k%g.nr())));return double(result);
  }
  Metrics budget(const StreamerState& start,const AcceptedStepRecord& r,const StreamerConfig& cfg)const{
    const double dt=r.dt_s,v0=r.pre_voltage_V,v1=r.post_voltage_V,vm=.5*(v0+v1);
    long double qw=0,fo=0,qhg=0,wcond=0,wouter=0,qcorr=0,wcorr=0,qreact=0,wreact=0,inventory=0,transport_tv=0;
    long double dqplasma=0,rho_source=0,corrcharge=0,boundarycharge=0,excluded_delta=0;
    double mu_bad=0,active_low=0,active_nonzero=0;
    for(std::size_t k=0;k<g.size();++k){int i=k%g.nr(),j=k/g.nr();if(!gas(i,j))continue;
      const double vol=g.cell_volume(i);
      if(!active[k]){excluded_delta+=std::abs((r.post.rho.values()[k]-start.rho.values()[k])*vol);continue;}
      const double pbar=.5*(start.phi.values()[k]+r.post.phi.values()[k]);
      const double corr=echarge*(r.pre_floor_np_delta_m3[k]+r.post_floor_np_delta_m3[k]+r.clip_np_delta_m3[k]
                        -r.pre_floor_ne_delta_m3[k]-r.post_floor_ne_delta_m3[k]-r.clip_ne_delta_m3[k]
                        -r.pre_floor_nn_delta_m3[k]-r.post_floor_nn_delta_m3[k]-r.clip_nn_delta_m3[k])*vol;
      auto tr=evaluate_morrow_lowke(r.pre.emag.values()[k],cfg.neutral_density,cfg.pressure,cfg.temperature);
      auto source=evaluate_reactions(r.pre.ne.values()[k],r.pre.np.values()[k],r.pre.nn.values()[k],r.pre.sph.values()[k],tr,cfg.temperature);
      const double qsrc=dt*echarge*(source.positive_ion-source.electron-source.negative_ion)*vol;
      qcorr+=psi.values()[k]*corr;wcorr+=pbar*corr;corrcharge+=corr;
      qreact+=psi.values()[k]*qsrc;wreact+=pbar*qsrc;rho_source+=qsrc;
      inventory+=echarge*(start.ne.values()[k]+start.np.values()[k]+start.nn.values()[k])*vol;
      dqplasma+=(r.post.rho.values()[k]-start.rho.values()[k])*vol;
      if(r.pre.ne.values()[k]>0){++active_nonzero;if(r.pre.emag.values()[k]/cfg.neutral_density*1e4<2.6e-17)++active_low;if(!std::isfinite(tr.mobility))++mu_bad;}
    }
    auto accumulate_face=[&](int i,int j,int ni,int nj,ElectrodeCellType metal,double gamma,double area,double species_charge){
      int a=idx(i,j);if(!active[a]||area==0)return;
      const bool inside=ni>=0&&ni<g.nr()&&nj>=0&&nj<g.nz();
      int b=inside?idx(ni,nj):-1;
      const long double q=species_charge*dt*gamma*area;
      if(inside&&active[b]&&b<a)return;
      transport_tv+=std::abs(q);
      const double pa=.5*(potential(start.phi,a,v0)+potential(r.post.phi,a,v1));
      const double wa=psi.values()[a];
      if(inside&&active[b]){
        const double pb=.5*(potential(start.phi,b,v0)+potential(r.post.phi,b,v1));
        qw+=q*(wa-psi.values()[b]);wcond+=q*(pa-pb);
      }else if(inside&&metal!=ElectrodeCellType::Gas){
        const bool hv=metal==ElectrodeCellType::HighVoltageElectrode;
        qw+=q*(wa-(hv?1.:0.));wcond+=q*(pa-(hv?vm:0.));
        if(hv)qhg-=q;boundarycharge-=q;
      }else {fo+=wa*q;wouter-=pa*q;boundarycharge-=q;}
    };
    for(const auto& f:r.electron_faces)accumulate_face(f.i,f.j,f.neighbor_i,f.neighbor_j,f.neighbor_material,f.gamma_e_m2s,f.area_m2,-echarge);
    for(const auto& f:r.positive_ion_faces)accumulate_face(f.i,f.j,f.face.neighbor_i,f.face.neighbor_j,f.face.neighbor_material,f.gamma_p_out_m2s,f.face.area_m2,echarge);
    const long double q0=charge(start.phi,v0),q1=charge(r.post.phi,v1),dq=q1-q0+qhg;
    const long double du=delta_energy(start.phi,v0,r.post.phi,v1),wp=vm*dq;
    Metrics m{{"time_n_s",start.time},{"time_n1_s",r.post.time},{"dt_s",dt},{"Vp_n_V",v0},{"Vp_n1_V",v1},
      {"q_B_active_start_C",double(q0)},{"q_B_active_end_C",double(q1)},{"q_A_end_C",r.q_native_hv_C},{"q_B_native_links_end_C",r.q_conjugate_hv_C},
      {"C_local_F",capacitance},{"Deltaq_port_C",double(dq)},{"Deltaq_Hg_C",double(qhg)},
      {"Deltaq_Hg_all_native_C",echarge*(r.hv_absorbed_electrons-r.hv_secondary_electrons-r.hv_collected_ions)},
      {"Deltaq_w_C",double(qw)},{"Deltaq_outer_C",double(fo)},{"DeltaQ_o_C",0.},
      {"DeltaQ_psi_correction_C",double(qcorr)},{"DeltaQ_psi_reaction_C",double(qreact)},
      {"charge_weighted_residual_C",double(dq-capacitance*(v1-v0)-qw-fo+qcorr+qreact)},
      {"DeltaQ_plasma_C",double(dqplasma)},{"transport_boundary_charge_C",double(boundarycharge)},
      {"reaction_charge_C",double(rho_source)},{"correction_charge_C",double(corrcharge)},
      {"plasma_charge_residual_C",double(dqplasma-boundarycharge-rho_source-corrcharge)},
      {"unsigned_native_transport_C",double(transport_tv)},{"initial_charge_inventory_C",double(inventory)},
      {"excluded_boundary_rho_change_TV_C",double(excluded_delta)},
      {"U_E_start_J",double(energy(start.phi,v0))},{"U_E_end_J",double(energy(r.post.phi,v1))},{"DeltaU_E_J",double(du)},
      {"Wport_native_J",double(wp)},{"Wcond_face_J",double(wcond)},{"Wouter_J",double(wouter)},
      {"Wcorrection_J",double(wcorr)},{"Wreaction_J",double(wreact)},
      {"energy_residual_J",double(wp+wouter+wcorr+wreact-du-wcond)},
      {"gauss_residual_max_C",gauss(r.post,v1)},{"active_lowfield_cells",active_low},{"active_electron_cells",active_nonzero},{"nonfinite_active_mu_cells",mu_bad},
      {"G_R_pre_S",r.G_R_pre_S},{"G_R_post_S",r.G_R_post_S},{"W_rho_FULL_GAS_post_C",r.W_rho_post_C}};
    return m;
  }
};
}
