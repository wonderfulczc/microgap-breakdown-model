#include "streamer_rf/streamer/ColdFullProvider.hpp"
#include "native_full_snapshot.hpp"
#include "native_persistent_checkpoint.hpp"
#include <stdexcept>

namespace streamer_rf::streamer {
struct ColdFullProvider::Trial::Candidate {
  NativeFullSnapshot snapshot;
  explicit Candidate(const NativeFullSnapshot& s):snapshot(s){}
};
struct ColdFullProvider::Impl {
  StreamerSolver& solver;
  std::uint64_t version=0;
  struct Buffer:AcceptedStepObserver {
    mutable std::shared_ptr<const AcceptedStepRecord> record;
    void on_accepted_step(const AcceptedStepRecord& r)const override{record=std::make_shared<AcceptedStepRecord>(r);}
  };
  explicit Impl(StreamerSolver& s):solver(s){}
};
ColdFullProvider::ColdFullProvider(StreamerSolver& s):impl_(std::make_unique<Impl>(s)){}
ColdFullProvider::~ColdFullProvider()=default;
std::uint64_t ColdFullProvider::version()const{return impl_->version;}
std::unique_ptr<ColdFullProvider::Trial> ColdFullProvider::evaluate_trial(double dt,bool fail){
  auto& s=impl_->solver;const auto accepted=StreamerTransactionAccess::capture(s);
  Impl::Buffer buffer;s.set_accepted_step_observer(&buffer);
  auto trial=std::make_unique<Trial>();trial->owner=this;trial->version=impl_->version;
  trial->base_physical_hash=native_hash(accepted);trial->base_history_hash=native_hash(accepted,true);
  try{
    if(!s.step(dt,trial->diagnostics))throw std::runtime_error("NATIVE_STEP_REJECTED");
    if(fail)throw std::runtime_error("CONTROLLED_FAILURE_AFTER_NATIVE_FIELDS");
    auto candidate=StreamerTransactionAccess::capture(s);
    candidate.observer=accepted.observer;
    trial->candidate=std::make_shared<Trial::Candidate>(candidate);trial->record=buffer.record;
    StreamerTransactionAccess::restore(s,accepted);
    return trial;
  }catch(...){StreamerTransactionAccess::restore(s,accepted);throw;}
}
void ColdFullProvider::commit(Trial& trial){
  if(trial.owner!=this)throw std::invalid_argument("FOREIGN_TRIAL");
  if(trial.consumed)throw std::invalid_argument("TRIAL_ALREADY_CONSUMED");
  if(trial.version!=impl_->version)throw std::invalid_argument("STALE_TRIAL");
  if(!trial.candidate||!trial.record)throw std::invalid_argument("MISSING_NATIVE_CANDIDATE");
  const auto current=StreamerTransactionAccess::capture(impl_->solver);
  if(native_hash(current)!=trial.base_physical_hash||native_hash(current,true)!=trial.base_history_hash)
    throw std::invalid_argument("ACCEPTED_STATE_CHANGED_OUTSIDE_PROVIDER");
  StreamerTransactionAccess::restore(impl_->solver,trial.candidate->snapshot);
  trial.consumed=true;++impl_->version;
  // Accepted output is published once, after committing the actual candidate state.
  if(trial.candidate->snapshot.observer)trial.candidate->snapshot.observer->on_accepted_step(*trial.record);
}
void ColdFullProvider::reject(Trial& trial)const{
  if(trial.owner!=this||trial.consumed)throw std::invalid_argument("INVALID_REJECT");
  trial.consumed=true;
}
void ColdFullProvider::save_accepted_payload(const std::filesystem::path& path,const std::string& identity,const AcceptedCheckpointLedger& ledger)const{
  auto frame=StreamerTransactionAccess::capture(impl_->solver);
  if(frame.observer)throw std::invalid_argument("PERSISTENT_OBSERVER_REBIND_REQUIRED");
  persistent_v31::save(path,frame,ledger,impl_->version,identity);
}
AcceptedCheckpointLedger ColdFullProvider::load_accepted_payload(const std::filesystem::path& path,const std::string& identity){
  if(impl_->version!=0||impl_->solver.state().time!=0)throw std::invalid_argument("RESTORE_REQUIRES_FRESH_PROVIDER");
  auto frame=StreamerTransactionAccess::capture(impl_->solver);std::uint64_t version=0;
  auto ledger=persistent_v31::load(path,frame,version,identity);
  StreamerTransactionAccess::restore(impl_->solver,frame);impl_->version=version;return ledger;
}
}
