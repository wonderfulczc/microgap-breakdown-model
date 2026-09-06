#include "streamer_rf/streamer/ReactionModel.hpp"
#include <cmath>
#include <stdexcept>
namespace streamer_rf::streamer {
namespace {constexpr double e=1.602176634e-19,kb=1.380649e-23;}
double electron_temperature(double mu,double D){if(mu<=0||D<0)throw std::invalid_argument("invalid electron transport");return e*D/(mu*kb);}
double beta_ep(double mu,double D){return 1.138e-11*std::pow(std::max(electron_temperature(mu,D),1.0),-.7);}
double beta_np(double T){if(T<=0)throw std::invalid_argument("invalid temperature");return 2e-13*std::sqrt(300/T);}
ReactionSources evaluate_reactions(double ne,double np,double nn,double sph,const TransportCoefficients& c,double T){
 if(ne<0||np<0||nn<0||sph<0)throw std::invalid_argument("negative state/source");ReactionSources s;
 s.ionization=c.ionization_frequency*ne;s.attachment_two_body=c.attachment_two_body_frequency*ne;s.attachment_three_body=c.attachment_three_body_frequency*ne;
 s.electron_positive_recombination=beta_ep(c.mobility,c.diffusion)*ne*np;s.ion_recombination=beta_np(T)*nn*np;
 s.electron=s.ionization-s.attachment_two_body-s.attachment_three_body-s.electron_positive_recombination+sph;
 s.positive_ion=s.ionization-s.electron_positive_recombination-s.ion_recombination+sph;
 s.negative_ion=s.attachment_two_body+s.attachment_three_body-s.ion_recombination;
 s.charge_balance=s.positive_ion-s.electron-s.negative_ion;return s;
}
ElectronReactionSourceComponents electron_reaction_source_components(double ne,double np,double nn,double sph,const TransportCoefficients& c,double T){
 const auto s=evaluate_reactions(ne,np,nn,sph,c,T);
 ElectronReactionSourceComponents out;
 out.impact_source_m3s=s.ionization;
 out.photo_source_m3s=sph;
 out.attachment2_loss_m3s=s.attachment_two_body;
 out.attachment3_loss_m3s=s.attachment_three_body;
 out.electron_recombination_loss_m3s=s.electron_positive_recombination;
 out.net_electron_reaction_source_m3s=s.electron;
 out.algebraic_closure_residual_m3s=out.impact_source_m3s+out.photo_source_m3s-out.attachment2_loss_m3s-out.attachment3_loss_m3s-out.electron_recombination_loss_m3s-out.net_electron_reaction_source_m3s;
 return out;
}
}
