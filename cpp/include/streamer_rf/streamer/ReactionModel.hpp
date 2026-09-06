#pragma once
#include "streamer_rf/streamer/MorrowLowke.hpp"
namespace streamer_rf::streamer {
struct ReactionSources { double ionization{},attachment_two_body{},attachment_three_body{},electron_positive_recombination{},ion_recombination{},electron{},positive_ion{},negative_ion{},charge_balance{}; };
struct ElectronReactionSourceComponents {
  double impact_source_m3s{},photo_source_m3s{},attachment2_loss_m3s{},attachment3_loss_m3s{};
  double electron_recombination_loss_m3s{},net_electron_reaction_source_m3s{},algebraic_closure_residual_m3s{};
};
double electron_temperature(double mobility,double diffusion);
double beta_ep(double mobility,double diffusion);
double beta_np(double gas_temperature);
ReactionSources evaluate_reactions(double ne,double np,double nn,double sph,const TransportCoefficients&,double gas_temperature);
ElectronReactionSourceComponents electron_reaction_source_components(double ne,double np,double nn,double sph,const TransportCoefficients&,double gas_temperature);
}
