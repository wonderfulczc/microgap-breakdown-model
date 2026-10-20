#pragma once
#include "native_full_snapshot.hpp"
#include <fstream>
#include <sstream>
#include <type_traits>

// Same-build, serial payload. The scenario envelope supplies SHA256, immutable
// config/source identities and atomic completion. No PETSc pointer is serialized.
namespace streamer_rf::streamer::persistent_v31 {
struct Stream {
  std::stringstream s;
  bool reading=false;
  template<class T> void scalar(T& x){
    static_assert(std::is_arithmetic_v<T>);
    if(reading)s.read(reinterpret_cast<char*>(&x),sizeof x);
    else s.write(reinterpret_cast<const char*>(&x),sizeof x);
    if(!s)throw std::runtime_error("CHECKPOINT_TRUNCATED");
  }
  void text(std::string& x){
    std::uint64_t n=x.size();scalar(n);
    if(n>1048576)throw std::runtime_error("CHECKPOINT_STRING_SIZE");
    if(reading){x.resize(n);s.read(x.data(),n);}else s.write(x.data(),n);
    if(!s)throw std::runtime_error("CHECKPOINT_TRUNCATED");
  }
  void vec(std::vector<double>& x){
    std::uint64_t n=x.size();scalar(n);
    if(n>20000000)throw std::runtime_error("CHECKPOINT_VECTOR_SIZE");
    if(reading)x.resize(n);
    for(auto& v:x)scalar(v);
  }
  void field(ScalarField2D& f){
    auto n=f.values().size();vec(f.values());
    if(f.values().size()!=n)throw std::runtime_error("CHECKPOINT_GRID_MISMATCH");
  }
  template<class T>void map(std::map<std::string,T>& x){
    std::uint64_t n=x.size();scalar(n);if(n>10000)throw std::runtime_error("CHECKPOINT_MAP_SIZE");
    if(reading){x.clear();for(std::uint64_t i=0;i<n;++i){std::string k;T v{};text(k);scalar(v);if(!x.emplace(k,v).second)throw std::runtime_error("CHECKPOINT_DUPLICATE_KEY");}}
    else for(const auto& item:x){std::string k=item.first;T v=item.second;text(k);scalar(v);}
  }
};
inline std::uint64_t checksum(const std::string& b){std::uint64_t h=1469598103934665603ULL;for(unsigned char c:b){h^=c;h*=1099511628211ULL;}return h;}
inline void transfer(Stream& b,NativeFullSnapshot& f,AcceptedCheckpointLedger& l,std::uint64_t& version,std::string& identity){
  std::string magic="I1_ACCEPTED_CHECKPOINT_V3.1";b.text(magic);
  if(magic!="I1_ACCEPTED_CHECKPOINT_V3.1")throw std::runtime_error("CHECKPOINT_VERSION_OR_KIND_MISMATCH");
  std::uint64_t endian=0x0102030405060708ULL;b.scalar(endian);
  if(endian!=0x0102030405060708ULL)throw std::runtime_error("CHECKPOINT_ABI_MISMATCH");
  b.text(identity);b.scalar(version);
  for(auto* p:{&f.state.ne,&f.state.np,&f.state.nn,&f.state.rho,&f.state.phi,&f.state.er,&f.state.ez,&f.state.emag,&f.state.sph})b.field(*p);
  b.scalar(f.state.time);
  for(auto& v:f.warm.components)b.vec(v);
  for(auto* p:{&f.head_position,&f.head_time,&f.q_hv,&f.q_return,&f.cap_native,&f.cap_conjugate})b.scalar(*p);
  b.scalar(f.has_history);b.scalar(f.has_weight);b.field(f.psi);b.field(f.wr);b.field(f.wz);b.text(f.weight_key);
  b.scalar(f.performance.fixed_weight_cache);b.scalar(f.performance.fixed_poisson_workspace);b.scalar(f.performance.profile);
  std::uint64_t nt=f.timings.size();b.scalar(nt);if(nt>1500)throw std::runtime_error("CHECKPOINT_TIMING_SIZE");
  if(b.reading)f.timings.resize(nt);
  for(auto& t:f.timings)for(auto* p:{&t.prelude_s,&t.transport_and_record_s,&t.update_s,&t.fields_s,&t.terminal_s,&t.observer_assembly_s,&t.callback_s,&t.weight_s})b.scalar(*p);
  auto& d=f.poisson_diagnostics;b.scalar(d.elliptic_iterations);b.scalar(d.poisson_iterations);b.scalar(d.profiling);
  nt=d.records.size();b.scalar(nt);if(nt>10000)throw std::runtime_error("CHECKPOINT_TIMING_SIZE");
  if(b.reading)d.records.resize(nt);
  for(auto& t:d.records){for(auto* p:{&t.assembly_s,&t.setup_s,&t.solve_s,&t.transfer_s,&t.residual_s,&t.true_relative_residual})b.scalar(*p);b.scalar(t.iterations);b.scalar(t.reason);b.scalar(t.reused);}
  b.scalar(l.dt_s);b.scalar(l.previous_time_s);b.scalar(l.accepted_step_id);b.scalar(l.native_calls);
  b.text(l.capacitance_mode);b.map(l.cumulative);b.map(l.output_offsets);b.vec(l.last_electron_faces);b.vec(l.last_ion_faces);b.vec(l.port_weight);
}
inline void validate(const NativeFullSnapshot& f,const AcceptedCheckpointLedger& l,std::uint64_t v){
  if(v!=l.accepted_step_id||l.native_calls<v||!std::isfinite(f.state.time)||f.state.time<=0||
     !std::isfinite(l.dt_s)||l.dt_s<=0||l.previous_time_s<0||
     std::abs(f.state.time-l.previous_time_s-l.dt_s)>1e-14*l.dt_s)
    throw std::runtime_error("CHECKPOINT_ACCEPTED_INTERVAL_MISMATCH");
  if(l.capacitance_mode!="TOTAL_CHARGE"&&l.capacitance_mode!="SPLIT_GEOMETRIC_CAPACITANCE")throw std::runtime_error("CHECKPOINT_OWNERSHIP_MISMATCH");
  for(const auto* p:{&f.state.ne,&f.state.np,&f.state.nn,&f.state.rho,&f.state.phi,&f.state.er,&f.state.ez,&f.state.emag,&f.state.sph})
    for(double x:p->values())if(!std::isfinite(x))throw std::runtime_error("CHECKPOINT_NONFINITE_FIELD");
  for(auto [k,x]:l.cumulative)if(!std::isfinite(x))throw std::runtime_error("CHECKPOINT_NONFINITE_LEDGER");
}
inline void save(const std::filesystem::path& path,NativeFullSnapshot f,AcceptedCheckpointLedger l,std::uint64_t version,std::string identity){
  validate(f,l,version);Stream b;transfer(b,f,l,version,identity);auto bytes=b.s.str();auto hash=checksum(bytes);
  if(std::filesystem::exists(path))throw std::runtime_error("CHECKPOINT_REFUSE_OVERWRITE");
  auto tmp=path;tmp+=".incomplete";std::ofstream out(tmp,std::ios::binary|std::ios::trunc);
  out.write(reinterpret_cast<const char*>(&hash),sizeof hash);out.write(bytes.data(),bytes.size());out.close();
  if(!out)throw std::runtime_error("CHECKPOINT_WRITE_FAILED");
  std::filesystem::rename(tmp,path);
}
inline AcceptedCheckpointLedger load(const std::filesystem::path& path,NativeFullSnapshot& f,std::uint64_t& version,const std::string& expected){
  if(std::filesystem::file_size(path)>512ULL*1024*1024)throw std::runtime_error("CHECKPOINT_SIZE");
  std::ifstream in(path,std::ios::binary);std::uint64_t hash=0;in.read(reinterpret_cast<char*>(&hash),sizeof hash);
  std::string bytes((std::istreambuf_iterator<char>(in)),{});
  if(!in||checksum(bytes)!=hash)throw std::runtime_error("CHECKPOINT_CHECKSUM_MISMATCH");
  Stream b;b.reading=true;b.s.str(bytes);std::string identity;AcceptedCheckpointLedger l;
  transfer(b,f,l,version,identity);
  if(identity!=expected)throw std::runtime_error("CHECKPOINT_IDENTITY_MISMATCH");
  if(b.s.peek()!=std::char_traits<char>::eof())throw std::runtime_error("CHECKPOINT_TRAILING_DATA");
  validate(f,l,version);f.photon_workspace.reset();return l;
}
}
