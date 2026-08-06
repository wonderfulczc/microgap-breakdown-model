#include "streamer_rf/types.hpp"
#include <cmath>
namespace streamer_rf {
AxisymmetricGrid::AxisymmetricGrid(int nr,int nz,double rmax,double zmin,double zmax)
 :nr_(nr),nz_(nz),r_max_(rmax),z_min_(zmin),z_max_(zmax),dr_(rmax/nr),dz_((zmax-zmin)/nz){
 if(nr<2||nz<2||rmax<=0||zmax<=zmin) throw std::invalid_argument("invalid axisymmetric grid");
}
double AxisymmetricGrid::cell_volume(int i)const{
 constexpr double pi=3.141592653589793238462643383279502884;
 const double rl=radial_face(i),rr=radial_face(i+1); return pi*(rr*rr-rl*rl)*dz_;
}
}

