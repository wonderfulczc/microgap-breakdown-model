#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace streamer_rf {

// Axial distance s increases from the facing tip into the conductor.
// total_length_m is tip-to-back length, not cylindrical-section length.
struct RoundedConeProfile2D {
  double tip_z_m{}, tip_radius_m{}, cone_half_angle_rad{}, shank_radius_m{}, total_length_m{};
  int extension_direction{1};

  void validate() const {
    const double pi=3.14159265358979323846;
    if (!std::isfinite(tip_z_m) || !std::isfinite(tip_radius_m) ||
        !std::isfinite(cone_half_angle_rad) || !std::isfinite(shank_radius_m) ||
        !std::isfinite(total_length_m) || tip_radius_m<=0 ||
        cone_half_angle_rad<=0 || cone_half_angle_rad>=pi/2 ||
        shank_radius_m<join_radius_m() || total_length_m<=cone_end_s_m() ||
        (extension_direction!=1 && extension_direction!=-1))
      throw std::invalid_argument("invalid rounded-cone electrode profile");
  }
  double join_s_m() const {return tip_radius_m*(1-std::sin(cone_half_angle_rad));}
  double join_radius_m() const {return tip_radius_m*std::cos(cone_half_angle_rad);}
  double cone_end_s_m() const {
    return join_s_m()+(shank_radius_m-join_radius_m())/std::tan(cone_half_angle_rad);
  }
  double cylinder_length_m() const {return total_length_m-cone_end_s_m();}
  double radius_at_s(double s) const {
    if (s<0 || s>total_length_m) return -1;
    if (s<=join_s_m()) return std::sqrt(std::max(0.0,s*(2*tip_radius_m-s)));
    if (s<=cone_end_s_m()) return join_radius_m()+(s-join_s_m())*std::tan(cone_half_angle_rad);
    return shank_radius_m;
  }
  bool contains(double r,double z) const {
    const double s=extension_direction*(z-tip_z_m);
    return r>=0 && r<=radius_at_s(s);
  }
  // Side-surface normal out of metal into gas. Back cap is a separate face.
  std::array<double,2> side_normal(double s) const {
    if(s<0 || s>=total_length_m) throw std::invalid_argument("side normal outside profile");
    if(s<=join_s_m()) return {radius_at_s(s)/tip_radius_m,
                              extension_direction*(s-tip_radius_m)/tip_radius_m};
    if(s<cone_end_s_m()) return {std::cos(cone_half_angle_rad),
                                -extension_direction*std::sin(cone_half_angle_rad)};
    return {1,0};
  }
};
}
