#pragma once
#include "streamer_rf/poisson.hpp"
#include "streamer_rf/sp3.hpp"
#include "streamer_rf/sp3_geometry.hpp"
#include "streamer_rf/streamer/ColdPhysics.hpp"
#include "streamer_rf/streamer/ReactionModel.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <string>
#include <filesystem>
#include <limits>
#include <vector>
namespace streamer_rf::streamer {
struct StreamerPerformanceOptions {bool fixed_weight_cache{false},fixed_poisson_workspace{false},profile{false};};
struct StreamerStepTiming {double prelude_s{},transport_and_record_s{},update_s{},fields_s{},terminal_s{},observer_assembly_s{},callback_s{},weight_s{};};
struct StreamerConfig {double pressure{101325},temperature{300},neutral_density{101325/(1.380649e-23*300)},background_field{4.8e6},n_ref{1e8},isg_epsilon{.01};bool photoionization{true};double excitation_ratio{.1};OpenBoundaryOptions open_boundary{};PoissonBoundaryConfig poisson_boundary{{BoundaryKind::OpenCharge,0.0},{BoundaryKind::OpenCharge,0.0},{BoundaryKind::OpenCharge,0.0}};SolverTolerances elliptic{};Sp3BoundaryOptions sp3_boundary{1e-8,1e-14,1.0,500,false};const AxisymmetricNeedlePlaneGeometry* electrode_geometry{};const VoltageWaveform* voltage_waveform{};double head_ne_threshold{1e16},bridge_ne_threshold{1e18},voltage_tolerance{1024*std::numeric_limits<double>::epsilon()};const Geometry2DView* geometry_view{};ColdPhysicsConfig cold_physics{};};
struct StreamerState {ScalarField2D ne,np,nn,rho,phi,er,ez,emag,sph;double time{};explicit StreamerState(const AxisymmetricGrid&g):ne(g),np(g),nn(g),rho(g),phi(g),er(g),ez(g),emag(g),sph(g){}};
struct ReactionTimestepDiagnostic {int i{-1},j{-1};double r{},z{},E{},E_over_N_Td{},ne{},np{},nn{},ionization_frequency{},attachment_two_body_frequency{},attachment_three_body_frequency{},recombination_frequency{},reaction_dt{},ne_max{},numerical_density_tolerance{};std::string cell_classification;};
struct TimeStepLimits {double drift{},diffusion{},ionization{},dielectric{},reaction{},selected{};std::string controller;};
struct ElectrodeSurfaceDiagnostics {double q_hv{},q_ground{},area_hv{},area_ground{};};
struct ConductanceDiagnostics {double p_cond{},gb{},rb{std::numeric_limits<double>::infinity()};bool valid{};};
struct ReactionSourceDecompositionDiagnostics {
  double impact_rate_s1{},photo_rate_s1{},attach2_rate_s1{},attach3_rate_s1{},recomb_e_rate_s1{};
  double net_electron_reaction_rate_s1{},source_closure_abs_s1{},source_closure_rel{},local_source_closure_max_abs_m3s{};
  int gas_cells{};
  std::string source_stage{"PRE_LIMITER"};
};
struct ElectronTransportCurrentSource {
  ScalarField2D jr, jz;
  double current_moment_r{}, current_moment_z{}, integral_abs_jz{}, max_abs_j{}, max_abs_jz{};
  explicit ElectronTransportCurrentSource(const AxisymmetricGrid& g) : jr(g), jz(g) {}
};
struct StreamerDiagnostics {double time{},dt{},emax{},ne_max{},np_max{},nn_max{},total_electrons{},total_charge{},conservation_residual{},applied_voltage{},phi_hv_residual{},phi_ground_residual{},gap{},tip_radius{},sigma_max{},head_position{},head_velocity{},absorbed_electron_hv{},absorbed_electron_ground{},i_cond_hv{},i_disp_hv{std::numeric_limits<double>::quiet_NaN()},i_total_hv{std::numeric_limits<double>::quiet_NaN()},i_cond_ground{},i_disp_ground{std::numeric_limits<double>::quiet_NaN()},i_total_ground{std::numeric_limits<double>::quiet_NaN()},q_hv{},q_ground{},q_conjugate_hv{std::numeric_limits<double>::quiet_NaN()},c_gap_vacuum{},c_gap_conjugate{std::numeric_limits<double>::quiet_NaN()},p_cond{},gb{},rb{std::numeric_limits<double>::infinity()},outer_boundary_current{},plasma_charge_derivative{},current_continuity_residual{};bool bridge_flag{},displacement_current_valid{},rb_valid{};int poisson_iterations{},sp3_ksp_iterations{},sp3_boundary_iterations{};std::string controller,geometry_id,head_status{"LEGACY"},bridge_status{"LEGACY"};};
struct GaussianSeed {double n0{},sigma{},z0{};};
struct NativeFaceFlux {int i{},j{},neighbor_i{},neighbor_j{};FaceDirection2D direction{};FaceRole2D role{};ElectrodeCellType neighbor_material{};double gamma_e_m2s{},area_m2{};};
struct AcceptedStepRecord {
  ColdPhysicsMode physics_mode{ColdPhysicsMode::LegacyElectronOnly};
  std::vector<IonFaceFlux> positive_ion_faces;
  GeometrySp3Audit photon_boundary_audit;
  double hv_collected_ions{},return_collected_ions{},outer_lost_ions{};
  double hv_secondary_electrons{},return_secondary_electrons{};
  double net_I_Hg_A{},net_I_Return_g_A{};
  double ion_positivity_dt_s{std::numeric_limits<double>::infinity()};
  std::string species_current_id{"LEGACY_ELECTRON_ONLY"};
  std::string geometry_id,geometry_revision,reference_plane_id,mesh_id,operator_id{"STAGE_C_PETSC_FV_V1"};
  double pre_time_s{},post_time_s{},dt_s{},pre_voltage_V{},post_voltage_V{};
  StreamerState pre,post;
  ScalarField2D weight_psi,weight_er,weight_ez;
  std::vector<NativeFaceFlux> electron_faces;
  std::vector<double> impact_m3s,photo_m3s,attachment2_m3s,attachment3_m3s,recombination_m3s;
  // I1_NATIVE_CELL_CHARGE_LEDGER_V1: actual reaction RHS, not a balance-derived source.
  std::vector<double> raw_pre_ne_m3,raw_pre_np_m3,raw_pre_nn_m3;
  std::vector<double> ion_recombination_m3s,reaction_e_m3s,reaction_p_m3s,reaction_n_m3s,reaction_charge_defect_m3s;
  std::vector<double> pre_floor_ne_delta_m3,pre_floor_np_delta_m3,pre_floor_nn_delta_m3;
  std::vector<double> post_floor_ne_delta_m3,post_floor_np_delta_m3,post_floor_nn_delta_m3;
  std::vector<double> clip_ne_delta_m3,clip_np_delta_m3,clip_nn_delta_m3;
  double hv_absorbed_electrons{},return_absorbed_electrons{},outer_lost_electrons{};
  double q_native_hv_C{},q_conjugate_hv_C{},native_capacitance_F{},conjugate_capacitance_F{};
  double G_R_pre_S{},G_R_post_S{},W_rho_pre_C{},W_rho_post_C{};
  double delta_G_reaction_S{},delta_G_transport_S{},delta_G_boundary_S{},delta_G_mobility_S{},delta_G_correction_S{},delta_G_closure_residual_S{};
  std::string outer_electrostatic_correction_status{"NOT_RESOLVED_FOR_OPEN_BOUNDARY"};
  std::string physical_port_status{"PENDING"};
  explicit AcceptedStepRecord(const AxisymmetricGrid& g):pre(g),post(g),weight_psi(g),weight_er(g),weight_ez(g){}
};
class AcceptedStepObserver {public:virtual ~AcceptedStepObserver()=default;virtual void on_accepted_step(const AcceptedStepRecord&) const=0;};
double absorbing_electrode_flux(double ne,double velocity_positive,double diffusion,double h,bool positive_face);
double axisymmetric_radial_face_area(const AxisymmetricGrid& g,int face_i);
double axisymmetric_axial_face_area(const AxisymmetricGrid& g,int cell_i);
class StreamerSolver {
 friend class StreamerTransactionAccess;
 public:
 const std::vector<StreamerStepTiming>& performance_timings()const{return step_timings_;}
 StreamerSolver(const AxisymmetricGrid&,StreamerConfig,StreamerPerformanceOptions = {});
 void initialize_gaussian(double n0,double sigma,double z0);
 void initialize_gaussian_at_tip_offset(double n0,double sigma,double z_offset_from_tip);
 void initialize_gaussians(const std::vector<GaussianSeed>& seeds);
 TimeStepLimits timestep_limits()const;
 ReactionTimestepDiagnostic reaction_timestep_diagnostic()const;
 double numerical_density_tolerance()const;
 ElectrodeSurfaceDiagnostics electrode_surface_diagnostics()const;
 ConductanceDiagnostics conductance_diagnostics(double voltage)const;
 ReactionSourceDecompositionDiagnostics reaction_source_decomposition_diagnostics()const;
 ElectronTransportCurrentSource electron_transport_current_source()const;
 double vacuum_gap_capacitance()const;
 double conjugate_hv_charge()const;
 double conjugate_vacuum_capacitance()const;
 void set_accepted_step_observer(const AcceptedStepObserver* observer){observer_=observer;}
 void reset_electrode_history();
 void refresh_electrostatic_fields();
 StreamerDiagnostics sample_terminal_diagnostics_from_history(double dt);
 bool step(double dt,StreamerDiagnostics&);
 void save_checkpoint(const std::filesystem::path&)const;
 void load_checkpoint(const std::filesystem::path&);
 const StreamerState&state()const{return state_;}
 StreamerState&state(){return state_;}
 private:
 std::vector<StreamerStepTiming> step_timings_;
 struct WeightCache {ScalarField2D psi,er,ez;std::string key;explicit WeightCache(const AxisymmetricGrid&g):psi(g),er(g),ez(g){}};
 StreamerPerformanceOptions performance_{};
 std::unique_ptr<FixedElectrodePoisson> poisson_workspace_;
 std::unique_ptr<WeightCache> weight_cache_;
 std::shared_ptr<GeometrySp3Workspace> photon_workspace_;
 const AxisymmetricGrid&g_;
 StreamerConfig c_;
 StreamerState state_;
 Sp3WarmStart sp3_warm_;
 double last_head_position_{std::numeric_limits<double>::quiet_NaN()},last_head_time_{std::numeric_limits<double>::quiet_NaN()},last_q_hv_{},last_q_ground_{};
 mutable double c_gap_vacuum_cache_{std::numeric_limits<double>::quiet_NaN()},c_gap_conjugate_cache_{std::numeric_limits<double>::quiet_NaN()};
 bool has_electrode_history_{};
 const AcceptedStepObserver* observer_{};
 bool electrode_mode()const{return (c_.electrode_geometry||c_.geometry_view)&&c_.voltage_waveform;}
 ElectrodeCellType cell_type(int i,int j)const;
 bool is_gas(int i,int j)const{return cell_type(i,j)==ElectrodeCellType::Gas;}
 double conjugate_hv_charge_for(const ScalarField2D&,double)const;
 double stage_c_activity_tolerance()const;
 void enforce_plasma_mask();void enforce_stage_c_activity_floor();void fields();
 void fill_electrode_diagnostics(StreamerDiagnostics&)const;
 double electrode_displacement_current(double q1,double q0,double dt)const;
 double head_position()const;bool bridge_flag()const;
};
}
