#pragma once
#include "native_full_network.hpp"
#include <memory>

namespace fp3full {
// Same active links and homogeneous conductor response; never changes the native field.
class StableNetwork {
 public:
  const ActiveNetwork& n;
  const StreamerConfig& cfg;
  int& rhs_calls;
  std::vector<long double> rpsi;
  long double cp_link=0,ce_minus_cp=0;
  struct Response { ScalarField2D rho,u; explicit Response(const AxisymmetricGrid&g):rho(g),u(g){} };
  std::vector<std::shared_ptr<Response>> cache;
  StableNetwork(const ActiveNetwork& net,const StreamerConfig& c,int& calls):n(net),cfg(c),rhs_calls(calls){
    rpsi=residual(n.psi,1,nullptr);
    cp_link=n.charge(n.psi,1);
    long double dot=0;for(std::size_t k=0;k<n.g.size();++k)if(n.active[k])dot+=(long double)n.psi.values()[k]*rpsi[k];
    ce_minus_cp=-dot+(cp_link-(long double)n.capacitance);
  }
  std::vector<long double> residual(const ScalarField2D& phi,double v,const ScalarField2D* rho)const{
    std::vector<long double> r(n.g.size());
    for(std::size_t k=0;k<r.size();++k)if(n.active[k]&&rho)r[k]=(long double)rho->values()[k]*n.g.cell_volume(k%n.g.nr());
    for(auto e:n.edges){long double q=(long double)e.c*((long double)n.potential(phi,e.a,v)-n.potential(phi,e.b,v));r[e.a]-=q;if(n.active[e.b])r[e.b]+=q;}
    return r;
  }
  const ScalarField2D& response(const StreamerState& s){
    for(const auto& p:cache)if(p->rho.values()==s.rho.values())return p->u;
    if(rhs_calls>=24)throw std::runtime_error("EXPLICIT_STATIC_RHS_BUDGET_EXHAUSTED");
    auto p=std::make_shared<Response>(n.g);p->rho=s.rho;++rhs_calls;
    solve_potential_with_electrodes(s.rho,n.view,0.,cfg.poisson_boundary,p->u,cfg.elliptic);
    cache.push_back(p);return p->u;
  }
  long double projection(const ScalarField2D& rho)const{
    long double q=0;for(std::size_t k=0;k<n.g.size();++k)if(n.active[k])q-=(long double)n.psi.values()[k]*rho.values()[k]*n.g.cell_volume(k%n.g.nr());return q;
  }
  long double cross(const ScalarField2D& u)const{
    long double x=0;for(std::size_t k=0;k<n.g.size();++k)if(n.active[k])x-=(long double)u.values()[k]*rpsi[k];return x;
  }
  long double rem(const ScalarField2D& u,double v)const{return n.energy(u,0)+v*cross(u)+.5L*ce_minus_cp*v*v;}
  ScalarField2D reconstruct(const ScalarField2D& u,double v)const{
    ScalarField2D p(n.g);for(std::size_t k=0;k<n.g.size();++k)p.values()[k]=v*n.psi.values()[k]+u.values()[k];return p;
  }
  Metrics endpoint(const StreamerState& s,double v){
    const auto& u=response(s);auto rp=residual(s.phi,v,&s.rho);
    long double eq=0,dot=0,maxr=0,maxw=0,diff=0,scale=0;
    for(std::size_t k=0;k<n.g.size();++k)if(n.active[k]){
      eq+=(long double)n.psi.values()[k]*rp[k]-rpsi[k]*s.phi.values()[k];
      dot+=(long double)n.psi.values()[k]*rpsi[k];maxr=std::max(maxr,std::abs(rp[k]));maxw=std::max(maxw,std::abs(rpsi[k]));
      diff=std::max(diff,std::abs((long double)s.phi.values()[k]-v*n.psi.values()[k]-u.values()[k]));scale=std::max(scale,std::abs((long double)s.phi.values()[k]));
    }
    long double qp=projection(s.rho),ql=n.charge(u,0),q=n.charge(s.phi,v);
    return {{"Q_rho_projection_C",double(qp)},{"Q_rho_zero_voltage_links_C",double(ql)},
      {"response_link_minus_projection_C",double(ql-qp)},
      {"Q_link_C",double(q)},{"Cp_V_C",n.capacitance*v},{"e_Q_residual_contraction_C",double(eq)},
      {"Cp_link_minus_stored_F",double(cp_link-n.capacitance)},
      {"Q_identity_residual_C",double(q-(long double)n.capacitance*v-qp-eq-(cp_link-n.capacitance)*v)},
      {"r_phi_max_C",double(maxr)},{"r_psi_max_F",double(maxw)},
      {"C_E_minus_Cp_F",double(ce_minus_cp)},
      {"C_E_minus_Cp_link_check_F",double(2*n.energy(n.psi,1)-n.capacitance-ce_minus_cp)},
      {"X_cross_C",double(cross(u))},{"U_sc_J",double(n.energy(u,0))},
      {"U_remaining_stable_J",double(rem(u,v))},
      {"native_minus_response_potential_max_V",double(diff)},
      {"native_minus_response_potential_relative",double(diff/scale)}};
  }
  Metrics budget(const StreamerState& start,const AcceptedStepRecord& rec,Metrics old){
    const double v0=rec.pre_voltage_V,v1=rec.post_voltage_V,vm=.5*(v0+v1);
    const auto& a=response(start);const auto& b=response(rec.post);
    long double dq_rho=0;
    for(std::size_t k=0;k<n.g.size();++k)if(n.active[k])dq_rho-=(long double)n.psi.values()[k]*((long double)rec.post.rho.values()[k]-start.rho.values()[k])*n.g.cell_volume(k%n.g.nr());
    const long double qpl=dq_rho+old.at("Deltaq_Hg_C"),qcap=(long double)n.capacitance*(v1-v0),dq=qcap+qpl;
    auto p0=reconstruct(a,v0),p1=reconstruct(b,v1);
    StreamerState s=start;s.phi=p0;AcceptedStepRecord r=rec;r.post.phi=p1;
    auto paired=n.budget(s,r,cfg);
    long double wcond=0,wabs=0,wouter=0,wcorr=0,wreact=0;
    // Work is evaluated from original native face charge, not a budget residual.
    auto accumulate_work=[&](int i,int j,int ni,int nj,ElectrodeCellType metal,double gamma,double area,long double species_charge){
      int x=n.idx(i,j);if(!n.active[x]||area==0)return;
      bool in=ni>=0&&ni<n.g.nr()&&nj>=0&&nj<n.g.nz();
      int y=in?n.idx(ni,nj):-1;if(in&&n.active[y]&&y<x)return;
      long double q=species_charge*rec.dt_s*gamma*area;
      auto mean=[&](int k){return .5L*((long double)v0*n.potential(n.psi,k,1)+n.potential(a,k,0)+(long double)v1*n.potential(n.psi,k,1)+n.potential(b,k,0));};
      if(in&&(n.active[y]||metal!=ElectrodeCellType::Gas)){
        long double w=q*(mean(x)-mean(y));wcond+=w;wabs+=std::abs(w);
      }else wouter-=mean(x)*q;
    };
    for(const auto& f:rec.electron_faces)accumulate_work(f.i,f.j,f.neighbor_i,f.neighbor_j,f.neighbor_material,f.gamma_e_m2s,f.area_m2,-(long double)echarge);
    for(const auto& f:rec.positive_ion_faces)accumulate_work(f.i,f.j,f.face.neighbor_i,f.face.neighbor_j,f.face.neighbor_material,f.gamma_p_out_m2s,f.face.area_m2,(long double)echarge);
    // Native correction/source arrays and original formulas, paired with response potentials.
    wcorr=paired.at("Wcorrection_J");wreact=paired.at("Wreaction_J");
    long double du_sc=n.delta_energy(a,0,b,0),x0=cross(a),x1=cross(b);
    long double du_rem=du_sc+(long double)v1*x1-(long double)v0*x0+ce_minus_cp*vm*(v1-v0);
    long double du_cap=(long double)n.capacitance*vm*(v1-v0),wp=(long double)vm*dq;
    long double re=(long double)vm*qpl+wouter+wcorr+wreact-du_rem-wcond;
    long double rq=qpl-old.at("Deltaq_w_C")-old.at("Deltaq_outer_C")+old.at("DeltaQ_psi_correction_C")+old.at("DeltaQ_psi_reaction_C");
    auto ea=endpoint(start,v0),eb=endpoint(rec.post,v1);
    for(auto [key,val]:ea)old["start_"+key]=val;
    for(auto [key,val]:eb)old["end_"+key]=val;
    auto r0=residual(start.phi,v0,&start.rho),r1=residual(rec.post.phi,v1,&rec.post.rho);
    long double deq=0,der=0;
    for(std::size_t k=0;k<n.g.size();++k)if(n.active[k]){
      deq+=(long double)n.psi.values()[k]*(r1[k]-r0[k])-rpsi[k]*((long double)rec.post.phi.values()[k]-start.phi.values()[k]);
      der+=.5L*((long double)start.phi.values()[k]+rec.post.phi.values()[k])*(r1[k]-r0[k]);
    }
    old.insert({{"stable_DeltaQ_rho_C",double(dq_rho)},{"stable_q_pl_C",double(qpl)},
      {"stable_Deltaq_port_C",double(dq)},{"stable_charge_residual_C",double(rq)},
      {"stable_Wport_J",double(wp)},{"stable_DeltaU_remaining_J",double(du_rem)},
      {"stable_DeltaU_E_J",double(du_rem+du_cap)},{"stable_Wcond_J",double(wcond)},
      {"stable_Wcond_abs_J",double(wabs)},{"stable_Wouter_J",double(wouter)},
      {"stable_Wcorrection_J",double(wcorr)},{"stable_Wreaction_J",double(wreact)},
      {"stable_energy_residual_J",double(re)},
      {"stable_U_remaining_start_J",double(rem(a,v0))},{"stable_U_remaining_end_J",double(rem(b,v1))},
      {"Delta_e_Q_C",double(deq)},
      {"old_RQ_minus_Delta_e_Q_C",old.at("charge_weighted_residual_C")-double(deq)-(double(cp_link)-n.capacitance)*(v1-v0)},
      {"ubar_Delta_r_phi_J",double(der)},
      {"old_RE_minus_ubar_Delta_r_phi_J",old.at("energy_residual_J")-double(der)},
      {"explicit_static_rhs_calls",double(rhs_calls)}});
    return old;
  }
};
}
