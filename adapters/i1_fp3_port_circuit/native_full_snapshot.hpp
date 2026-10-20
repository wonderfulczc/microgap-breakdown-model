#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace streamer_rf::streamer {

// No PETSc objects or ownership-bearing raw pointers are copied by this frame.
struct NativeFullSnapshot {
  const StreamerSolver* owner;
  StreamerState state;
  StreamerConfig config;
  StreamerPerformanceOptions performance;
  Sp3WarmStart warm;
  std::shared_ptr<GeometrySp3Workspace> photon_workspace;
  double head_position,head_time,q_hv,q_return,cap_native,cap_conjugate;
  bool has_history,has_weight;
  ScalarField2D psi,wr,wz;
  std::string weight_key;
  const AcceptedStepObserver* observer;
  std::vector<StreamerStepTiming> timings;
  PoissonDiagnosticSnapshot poisson_diagnostics;
  explicit NativeFullSnapshot(const AxisymmetricGrid& g):owner(nullptr),state(g),psi(g),wr(g),wz(g){}
};

class StreamerTransactionAccess {
 public:
  static NativeFullSnapshot capture(const StreamerSolver& s){
    NativeFullSnapshot f(s.g_);f.owner=&s;f.state=s.state_;f.config=s.c_;f.performance=s.performance_;
    f.warm=s.sp3_warm_;f.head_position=s.last_head_position_;f.head_time=s.last_head_time_;
    f.photon_workspace=s.photon_workspace_;
    f.q_hv=s.last_q_hv_;f.q_return=s.last_q_ground_;f.cap_native=s.c_gap_vacuum_cache_;
    f.cap_conjugate=s.c_gap_conjugate_cache_;f.has_history=s.has_electrode_history_;
    f.observer=s.observer_;f.timings=s.step_timings_;f.poisson_diagnostics=snapshot_poisson_diagnostics();
    f.has_weight=bool(s.weight_cache_);
    if(f.has_weight){f.psi=s.weight_cache_->psi;f.wr=s.weight_cache_->er;f.wz=s.weight_cache_->ez;f.weight_key=s.weight_cache_->key;}
    return f;
  }
  static void restore(StreamerSolver& s,const NativeFullSnapshot& f){
    if(f.owner!=&s||&f.state.ne.grid()!=&s.g_)throw std::invalid_argument("FOREIGN_NATIVE_SNAPSHOT");
    s.state_=f.state;s.c_=f.config;s.performance_=f.performance;s.sp3_warm_=f.warm;
    s.photon_workspace_=f.photon_workspace;
    s.last_head_position_=f.head_position;s.last_head_time_=f.head_time;
    s.last_q_hv_=f.q_hv;s.last_q_ground_=f.q_return;s.c_gap_vacuum_cache_=f.cap_native;
    s.c_gap_conjugate_cache_=f.cap_conjugate;s.has_electrode_history_=f.has_history;
    s.observer_=f.observer;s.step_timings_=f.timings;restore_poisson_diagnostics(f.poisson_diagnostics);
    if(f.has_weight){
      if(!s.weight_cache_)s.weight_cache_=std::make_unique<StreamerSolver::WeightCache>(s.g_);
      s.weight_cache_->psi=f.psi;s.weight_cache_->er=f.wr;s.weight_cache_->ez=f.wz;s.weight_cache_->key=f.weight_key;
    }else s.weight_cache_.reset();
    // Fixed workspace is not copied: same immutable matrix, fresh RHS/u per solve.
    // Its preconditioner is reusable; recreation counters are not physical state.
  }
  static void initial_floor(StreamerSolver& s){s.enforce_plasma_mask();s.enforce_stage_c_activity_floor();s.fields();s.reset_electrode_history();}
};

inline std::string native_hash(const NativeFullSnapshot& f,bool history=false){
  std::uint64_t h=1469598103934665603ULL;
  auto bytes=[&](const void*p,std::size_t n){const auto*c=static_cast<const unsigned char*>(p);for(std::size_t k=0;k<n;++k){h^=c[k];h*=1099511628211ULL;}};
  if(!history){for(const auto*p:{&f.state.ne,&f.state.np,&f.state.nn,&f.state.rho,&f.state.phi,&f.state.er,&f.state.ez,&f.state.emag,&f.state.sph})bytes(p->values().data(),p->values().size()*sizeof(double));bytes(&f.state.time,sizeof(double));}
  else {for(double x:{f.head_position,f.head_time,f.q_hv,f.q_return,f.cap_native,f.cap_conjugate})bytes(&x,sizeof x);bytes(&f.has_history,sizeof(bool));for(const auto&v:f.warm.components)bytes(v.data(),v.size()*sizeof(double));}
  return std::to_string(h);
}
}
