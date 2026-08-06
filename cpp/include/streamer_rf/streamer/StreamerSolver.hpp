#pragma once
#include "streamer_rf/poisson.hpp"
#include "streamer_rf/sp3.hpp"
#include "streamer_rf/streamer/ReactionModel.hpp"
#include <string>
#include <filesystem>
#include <vector>
namespace streamer_rf::streamer {
struct StreamerConfig {double pressure{101325},temperature{300},neutral_density{101325/(1.380649e-23*300)},background_field{4.8e6},n_ref{1e8},isg_epsilon{.01};bool photoionization{true};double excitation_ratio{.1};OpenBoundaryOptions open_boundary{};SolverTolerances elliptic{};Sp3BoundaryOptions sp3_boundary{1e-8,1e-14,1.0,500,false};};
struct StreamerState {ScalarField2D ne,np,nn,rho,phi,er,ez,emag,sph;double time{};explicit StreamerState(const AxisymmetricGrid&g):ne(g),np(g),nn(g),rho(g),phi(g),er(g),ez(g),emag(g),sph(g){}};
struct TimeStepLimits {double drift{},diffusion{},ionization{},dielectric{},reaction{},selected{};std::string controller;};
struct StreamerDiagnostics {double time{},dt{},emax{},ne_max{},np_max{},nn_max{},total_electrons{},total_charge{},conservation_residual{};int poisson_iterations{},sp3_ksp_iterations{},sp3_boundary_iterations{};std::string controller;};
struct GaussianSeed {double n0{},sigma{},z0{};};
class StreamerSolver {public: StreamerSolver(const AxisymmetricGrid&,StreamerConfig);void initialize_gaussian(double n0,double sigma,double z0);void initialize_gaussians(const std::vector<GaussianSeed>& seeds);TimeStepLimits timestep_limits()const;bool step(double dt,StreamerDiagnostics&);void save_checkpoint(const std::filesystem::path&)const;void load_checkpoint(const std::filesystem::path&);const StreamerState&state()const{return state_;}StreamerState&state(){return state_;}private:const AxisymmetricGrid&g_;StreamerConfig c_;StreamerState state_;Sp3WarmStart sp3_warm_;void fields();};
}
