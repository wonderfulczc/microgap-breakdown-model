#pragma once
#include "streamer_rf/streamer/StreamerSolver.hpp"
#include <memory>
#include <cstdint>
#include <map>

namespace streamer_rf::streamer {
struct AcceptedCheckpointLedger {
  double dt_s{},previous_time_s{};
  std::uint64_t accepted_step_id{},native_calls{};
  std::string capacitance_mode{"TOTAL_CHARGE"};
  std::map<std::string,double> cumulative;
  std::map<std::string,std::uint64_t> output_offsets;
  std::vector<double> last_electron_faces,last_ion_faces,port_weight;
};
// Thin transaction boundary over the existing FP3 native snapshot implementation.
// The immutable configuration, geometry and voltage owner must outlive the provider.
class ColdFullProvider {
 public:
  struct Trial;
  explicit ColdFullProvider(StreamerSolver&);
  ~ColdFullProvider();
  std::unique_ptr<Trial> evaluate_trial(double dt_s,bool controlled_failure=false);
  void commit(Trial&);
  void reject(Trial&) const;
  std::uint64_t version()const;
  void save_accepted_payload(const std::filesystem::path&,const std::string& identity,
                             const AcceptedCheckpointLedger&)const;
  AcceptedCheckpointLedger load_accepted_payload(const std::filesystem::path&,const std::string& identity);
  struct Trial {
    const ColdFullProvider* owner{};
    std::uint64_t version{};
    bool consumed{false};
    std::string base_physical_hash,base_history_hash;
    StreamerDiagnostics diagnostics;
    std::shared_ptr<const AcceptedStepRecord> record;
    struct Candidate;
    std::shared_ptr<Candidate> candidate;
  };
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}
