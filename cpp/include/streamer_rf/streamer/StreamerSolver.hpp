#pragma once
#include "streamer_rf/poisson.hpp"
#include "streamer_rf/sp3.hpp"
#include "streamer_rf/streamer/ReactionModel.hpp"
#include "streamer_rf/voltage_waveform.hpp"
#include <string>
#include <filesystem>
#include <limits>
#include <vector>
namespace streamer_rf::streamer {
struct StreamerConfig {double pressure{101325},temperature{300},neutral_density{101325/(1.380649e-23*300)},background_field{4.8e6},n_ref{1e8},isg_epsilon{.01};bool photoionization{true};double excitation_ratio{.1};OpenBoundaryOptions open_boundary{};PoissonBoundaryConfig poisson_boundary{{BoundaryKind::OpenCharge,0.0},{BoundaryKind::OpenCharge,0.0},{BoundaryKind::OpenCharge,0.0}};SolverTolerances elliptic{};Sp3BoundaryOptions sp3_boundary{1e-8,1e-14,1.0,500,false};const AxisymmetricNeedlePlaneGeometry* electrode_geometry{};const VoltageWaveform* voltage_waveform{};double head_ne_threshold{1e16},bridge_ne_threshold{1e18},voltage_tolerance{1024*std::numeric_limits<double>::epsilon()};};
struct StreamerState {ScalarField2D ne,np,nn,rho,phi,er,ez,emag,sph;double time{};explicit StreamerState(const AxisymmetricGrid&g):ne(g),np(g),nn(g),rho(g),phi(g),er(g),ez(g),emag(g),sph(g){}};
struct ReactionTimestepDiagnostic {int i{-1},j{-1};double r{},z{},E{},E_over_N_Td{},ne{},np{},nn{},ionization_frequency{},attachment_two_body_frequency{},attachment_three_body_frequency{},recombination_frequency{},reaction_dt{},ne_max{},numerical_density_tolerance{};std::string cell_classification;};
struct TimeStepLimits {double drift{},diffusion{},ionization{},dielectric{},reaction{},selected{};std::string controller;};
struct ElectrodeSurfaceDiagnostics {double q_hv{},q_ground{},area_hv{},area_ground{};};
struct ConductanceDiagnostics {double p_cond{},gb{},rb{std::numeric_limits<double>::infinity()};bool valid{};};
struct ElectronTransportCurrentSource {
  ScalarField2D jr, jz;
  double current_moment_r{}, current_moment_z{}, integral_abs_jz{}, max_abs_j{}, max_abs_jz{};
  explicit ElectronTransportCurrentSource(const AxisymmetricGrid& g) : jr(g), jz(g) {}
};
struct StreamerDiagnostics {double time{},dt{},emax{},ne_max{},np_max{},nn_max{},total_electrons{},total_charge{},conservation_residual{},applied_voltage{},phi_hv_residual{},phi_ground_residual{},gap{},tip_radius{},sigma_max{},head_position{},head_velocity{},absorbed_electron_hv{},absorbed_electron_ground{},i_cond_hv{},i_disp_hv{std::numeric_limits<double>::quiet_NaN()},i_total_hv{std::numeric_limits<double>::quiet_NaN()},i_cond_ground{},i_disp_ground{std::numeric_limits<double>::quiet_NaN()},i_total_ground{std::numeric_limits<double>::quiet_NaN()},q_hv{},q_ground{},c_gap_vacuum{},p_cond{},gb{},rb{std::numeric_limits<double>::infinity()},outer_boundary_current{},plasma_charge_derivative{},current_continuity_residual{};bool bridge_flag{},displacement_current_valid{},rb_valid{};int poisson_iterations{},sp3_ksp_iterations{},sp3_boundary_iterations{};std::string controller,geometry_id;};
struct GaussianSeed {double n0{},sigma{},z0{};};
double absorbing_electrode_flux(double ne,double velocity_positive,double diffusion,double h,bool positive_face);
double axisymmetric_radial_face_area(const AxisymmetricGrid& g,int face_i);
double axisymmetric_axial_face_area(const AxisymmetricGrid& g,int cell_i);
class StreamerSolver {public: StreamerSolver(const AxisymmetricGrid&,StreamerConfig);void initialize_gaussian(double n0,double sigma,double z0);void initialize_gaussian_at_tip_offset(double n0,double sigma,double z_offset_from_tip);void initialize_gaussians(const std::vector<GaussianSeed>& seeds);TimeStepLimits timestep_limits()const;ReactionTimestepDiagnostic reaction_timestep_diagnostic()const;double numerical_density_tolerance()const;ElectrodeSurfaceDiagnostics electrode_surface_diagnostics()const;ConductanceDiagnostics conductance_diagnostics(double voltage)const;ElectronTransportCurrentSource electron_transport_current_source()const;double vacuum_gap_capacitance()const;void reset_electrode_history();void refresh_electrostatic_fields();StreamerDiagnostics sample_terminal_diagnostics_from_history(double dt);bool step(double dt,StreamerDiagnostics&);void save_checkpoint(const std::filesystem::path&)const;void load_checkpoint(const std::filesystem::path&);const StreamerState&state()const{return state_;}StreamerState&state(){return state_;}private:const AxisymmetricGrid&g_;StreamerConfig c_;StreamerState state_;Sp3WarmStart sp3_warm_;double last_head_position_{std::numeric_limits<double>::quiet_NaN()},last_head_time_{std::numeric_limits<double>::quiet_NaN()},last_q_hv_{},last_q_ground_{};mutable double c_gap_vacuum_cache_{std::numeric_limits<double>::quiet_NaN()};bool has_electrode_history_{};bool electrode_mode()const{return c_.electrode_geometry&&c_.voltage_waveform;}ElectrodeCellType cell_type(int i,int j)const;bool is_gas(int i,int j)const{return cell_type(i,j)==ElectrodeCellType::Gas;}double stage_c_activity_tolerance()const;void enforce_plasma_mask();void enforce_stage_c_activity_floor();void fields();void fill_electrode_diagnostics(StreamerDiagnostics&)const;double electrode_displacement_current(double q1,double q0,double dt)const;double head_position()const;bool bridge_flag()const;};
}
