#include "streamer_rf/streamer/MorrowLowke.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace streamer_rf::streamer {
TransportCoefficients evaluate_morrow_lowke(double E,double N,double pressure,double temperature){
 if(!std::isfinite(E)||!std::isfinite(N)||N<=0||pressure<=0||temperature<=0)throw std::invalid_argument("invalid Morrow-Lowke state");
 E=std::abs(E);TransportCoefficients c{};if(E==0){c.mobility=6.87e24/N;c.diffusion=0;return c;}
 const double x=E/N*1e4,Nc=N/1e6;
 const double alphaN=x>1.5e-15?2e-16*std::exp(-7.248e-15/x):6.619e-17*std::exp(-5.593e-15/x);
 const double eta2N=std::max(0.0,x>1.05e-15?8.889e-5*x+2.567e-19:6.089e-4*x-2.893e-19);
 const double eta3N2=4.7778e-59*std::pow(x,-1.2749);
 double wcm;
 if(x>2e-15)wcm=7.4e21*x+7.1e6;else if(x>=1e-16)wcm=1.03e22*x+1.3e6;else if(x>=2.6e-17){const double x0=2.6e-17,x1=1e-16,t=(x-x0)/(x1-x0),y0=6.87e22*x0+3.38e4,y1=1.03e22*x1+1.3e6,m0=6.87e22*(x1-x0),m1=1.03e22*(x1-x0);wcm=(2*t*t*t-3*t*t+1)*y0+(t*t*t-2*t*t+t)*m0+(-2*t*t*t+3*t*t)*y1+(t*t*t-t*t)*m1;}else wcm=6.87e22*x+3.38e4;
 c.drift_speed=wcm/100;c.mobility=c.drift_speed/E;c.diffusion=0.3341e9*std::pow(x,.54069)*wcm/E*1e-2;
 c.ionization_townsend=alphaN*Nc*100;c.attachment_two_body_townsend=eta2N*Nc*100;c.attachment_three_body_townsend=eta3N2*Nc*Nc*100;
 c.ionization_frequency=c.ionization_townsend*c.drift_speed;c.attachment_two_body_frequency=c.attachment_two_body_townsend*c.drift_speed;c.attachment_three_body_frequency=c.attachment_three_body_townsend*c.drift_speed;return c;
}
double morrow_lowke_breakdown_field(double N,double p,double T){double lo=1e5,hi=2e7;for(int k=0;k<160;++k){double m=.5*(lo+hi);auto c=evaluate_morrow_lowke(m,N,p,T);if(c.ionization_frequency-c.attachment_two_body_frequency-c.attachment_three_body_frequency>0)hi=m;else lo=m;}return .5*(lo+hi);}
}
