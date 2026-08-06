#pragma once
#include "streamer_rf/streamer/MorrowLowke.hpp"
namespace streamer_rf::streamer {
struct ReactionSources { double ionization{},attachment_two_body{},attachment_three_body{},electron_positive_recombination{},ion_recombination{},electron{},positive_ion{},negative_ion{},charge_balance{}; };
double electron_temperature(double mobility,double diffusion);
double beta_ep(double mobility,double diffusion);
double beta_np(double gas_temperature);
ReactionSources evaluate_reactions(double ne,double np,double nn,double sph,const TransportCoefficients&,double gas_temperature);
}
