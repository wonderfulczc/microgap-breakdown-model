#include "streamer_rf/transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace streamer_rf {
double bernoulli(double x){
 const double ax=std::abs(x);
 if(ax<1e-3){const double x2=x*x;return 1-x/2+x2/12-x2*x2/720+x2*x2*x2/30240;}
 if(x>50)return x*std::exp(-x); if(x<-50)return -x; return x/std::expm1(x);
}
double sg_flux(double nl,double nr,double w,double d,double h){
 if(d<=0||h<=0)throw std::invalid_argument("D and h must be positive");
 const double p=w*h/d; return d/h*(bernoulli(-p)*nl-bernoulli(p)*nr);
}
double isg0_zero_width_flux(double nl,double nr,double w,double d,double h,double ref){
 if(ref<=0||nl<0||nr<0)throw std::invalid_argument("densities and n_ref invalid");
 const double ul=nl/ref,ur=nr/ref;
 const double nmid=ref*(std::sqrt((1+ul)*(1+ur))-1);
 const double a=(std::log1p(ur)-std::log1p(ul))/h;
 return nmid*(w-d*a);
}
double isg0_flux(double nl,double nr,double wl,double wr,double d,double h,double eps,double ref,IsgDiagnostics* out){
 if(ref<=0||eps<=0||d<=0||h<=0||nl<0||nr<0)throw std::invalid_argument("invalid ISG-0 input");
 IsgDiagnostics q; const double dw=wr-wl,scale=std::max({1.0,std::abs(wl),std::abs(wr)});
 if(std::abs(dw)<64*std::numeric_limits<double>::epsilon()*scale){q.ordinary_sg=true;q.branch="small_velocity_gradient";if(out)*out=q;return sg_flux(nl,nr,.5*(wl+wr),d,h);}
 q.h_virtual=std::sqrt(2*eps*d*h/std::abs(dw));q.ratio=q.h_virtual/h;
 if(q.h_virtual>=h){q.ordinary_sg=true;q.branch="virtual_width_ge_cell";if(out)*out=q;return sg_flux(nl,nr,.5*(wl+wr),d,h);}
 if(q.ratio<1e-6){q.zero_width=true;q.branch="zero_width";if(out)*out=q;return isg0_zero_width_flux(nl,nr,.5*(wl+wr),d,h,ref);}
 const double ul=nl/ref,ur=nr/ref,a=(std::log1p(ur)-std::log1p(ul))/h;
 const double x=.5*q.h_virtual;
 const double nvl=ref*((1+ul)*std::exp(a*(.5*h-x))-1);
 const double nvr=ref*((1+ul)*std::exp(a*(.5*h+x))-1);
 q.branch="isg0";if(out)*out=q;
 return sg_flux(nvl,nvr,.5*(wl+wr),d,q.h_virtual);
}
}

