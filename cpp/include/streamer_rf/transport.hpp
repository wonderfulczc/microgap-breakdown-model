#pragma once
#include <string>
namespace streamer_rf {
double bernoulli(double x);
double sg_flux(double n_left,double n_right,double velocity,double diffusion,double h);
struct IsgDiagnostics { bool ordinary_sg{false}, zero_width{false}; double h_virtual{0}, ratio{0}; std::string branch; };
double isg0_zero_width_flux(double n_left,double n_right,double velocity_mid,double diffusion,double h,double n_ref);
double isg0_flux(double n_left,double n_right,double velocity_left,double velocity_right,
                 double diffusion,double h,double epsilon,double n_ref,IsgDiagnostics* diagnostics=nullptr);
}

