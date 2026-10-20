#pragma once
#include "streamer_rf/electrode.hpp"
#include "streamer_rf/electrode_profile.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace streamer_rf {

enum class FaceRole2D { GasInterface, Electrode, OuterBoundary, Axis };
enum class FaceDirection2D { RadialMinus, RadialPlus, AxialMinus, AxialPlus };

struct GeometryFace2D {
  FaceRole2D role;
  ElectrodeCellType neighbor_material;
  int neighbor_i, neighbor_j;
  double normal_r, normal_z, area_m2;
};

struct FacingNeedles2D {
  std::string geometry_id, revision{"1"}, reference_plane_id{"HV-to-return"};
  double return_tip_z_m{}, hv_tip_z_m{};
  double return_tip_radius_m{}, hv_tip_radius_m{};
  double return_shank_radius_m{}, hv_shank_radius_m{};
  double return_shank_length_m{}, hv_shank_length_m{};
  bool axisymmetric{true}, dielectric{false};
};

// Immutable once constructed; the grid must outlive this view and the solver.
class Geometry2DView final {
 public:
  static Geometry2DView needle_plane(const AxisymmetricGrid& grid,
                                     const AxisymmetricNeedlePlaneGeometry& old,
                                     std::string revision = "legacy-mask-v1") {
    std::vector<ElectrodeCellType> cells;
    cells.reserve(grid.size());
    for (int j=0;j<grid.nz();++j) for (int i=0;i<grid.nr();++i)
      cells.push_back(old.classify(grid,i,j));
    return Geometry2DView(grid,std::move(cells),old.metadata().geometry_id,
                          std::move(revision),"HV-to-ground",
                          old.ground_z_m()+old.ground_thickness_m(),old.tip_z_m()-old.tip_radius_m(),
                          old.tip_z_m(),old.metadata().tip_radius_m,"NEEDLE_PLANE");
  }
  static Geometry2DView facing_needles(const AxisymmetricGrid& grid,const FacingNeedles2D& p) {
    if (!p.axisymmetric) throw std::invalid_argument("NONAXISYMMETRIC_3D_REQUIRED");
    if (p.dielectric) throw std::invalid_argument("UNSUPPORTED_DIELECTRIC");
    if (p.geometry_id.empty() || p.revision.empty() || p.reference_plane_id.empty() ||
        !(p.hv_tip_z_m>p.return_tip_z_m) ||
        !(p.return_tip_radius_m>0 && p.hv_tip_radius_m>0 &&
          p.return_shank_radius_m>0 && p.hv_shank_radius_m>0 &&
          p.return_shank_length_m>0 && p.hv_shank_length_m>0) ||
        p.return_shank_radius_m>p.return_tip_radius_m ||
        p.hv_shank_radius_m>p.hv_tip_radius_m)
      throw std::invalid_argument("invalid facing-needle geometry");
    std::vector<ElectrodeCellType> cells;
    cells.reserve(grid.size());
    const double lower_center=p.return_tip_z_m-p.return_tip_radius_m;
    const double upper_center=p.hv_tip_z_m+p.hv_tip_radius_m;
    for (int j=0;j<grid.nz();++j) for (int i=0;i<grid.nr();++i) {
      const double r=grid.r(i),z=grid.z(j);
      const bool lower=std::hypot(r,z-lower_center)<=p.return_tip_radius_m ||
        (r<=p.return_shank_radius_m && z<=lower_center && z>=lower_center-p.return_shank_length_m);
      const bool upper=std::hypot(r,z-upper_center)<=p.hv_tip_radius_m ||
        (r<=p.hv_shank_radius_m && z>=upper_center && z<=upper_center+p.hv_shank_length_m);
      if (lower && upper) throw std::invalid_argument("overlapping electrode masks");
      cells.push_back(lower?ElectrodeCellType::GroundElectrode:
                      upper?ElectrodeCellType::HighVoltageElectrode:ElectrodeCellType::Gas);
    }
    return Geometry2DView(grid,std::move(cells),p.geometry_id,p.revision,p.reference_plane_id,
                          p.return_tip_z_m,p.hv_tip_z_m,p.hv_tip_z_m,p.hv_tip_radius_m,"FACING_NEEDLES");
  }
  static Geometry2DView facing_profiles(const AxisymmetricGrid& grid,
      const RoundedConeProfile2D& lower,const RoundedConeProfile2D& upper,
      std::string id,std::string revision,std::string reference,
      bool axisymmetric=true,bool dielectric=false) {
    if(!axisymmetric) throw std::invalid_argument("NONAXISYMMETRIC_3D_REQUIRED");
    if(dielectric) throw std::invalid_argument("UNSUPPORTED_DIELECTRIC");
    lower.validate(); upper.validate();
    if(lower.extension_direction!=-1 || upper.extension_direction!=1 ||
       upper.tip_z_m<=lower.tip_z_m || id.empty() || revision.empty() || reference.empty())
      throw std::invalid_argument("invalid facing profile identity or orientation");
    std::vector<ElectrodeCellType> cells;
    cells.reserve(grid.size());
    for(int j=0;j<grid.nz();++j) for(int i=0;i<grid.nr();++i) {
      const bool lo=lower.contains(grid.r(i),grid.z(j));
      const bool hi=upper.contains(grid.r(i),grid.z(j));
      if(lo&&hi) throw std::invalid_argument("overlapping profile electrodes");
      cells.push_back(lo?ElectrodeCellType::GroundElectrode:
                      hi?ElectrodeCellType::HighVoltageElectrode:ElectrodeCellType::Gas);
    }
    return Geometry2DView(grid,std::move(cells),std::move(id),std::move(revision),
        std::move(reference),lower.tip_z_m,upper.tip_z_m,upper.tip_z_m,
        upper.tip_radius_m,"FACING_ROUNDED_CONE_PROFILES");
  }
  ElectrodeCellType classify(int i,int j) const {
    if (i<0 || i>=grid_->nr() || j<0 || j>=grid_->nz()) throw std::out_of_range("geometry cell");
    return cells_[static_cast<std::size_t>(j)*grid_->nr()+i];
  }
  const AxisymmetricGrid& grid() const {return *grid_;}
  const std::string& geometry_id() const {return geometry_id_;}
  const std::string& revision() const {return revision_;}
  const std::string& reference_plane_id() const {return reference_plane_id_;}
  const std::string& family() const {return family_;}
  double net_gap_m() const {return hv_surface_z_m_-return_surface_z_m_;}
  double return_surface_z_m() const {return return_surface_z_m_;}
  double hv_surface_z_m() const {return hv_surface_z_m_;}
  double seed_reference_z_m() const {return seed_reference_z_m_;}
  double tip_radius_m() const {return tip_radius_m_;}
  GeometryFace2D face(int i,int j,FaceDirection2D direction) const {
    if (classify(i,j)!=ElectrodeCellType::Gas) throw std::invalid_argument("face requires gas cell");
    int ni=i,nj=j; double nr=0,nz=0,area=0;
    constexpr double pi=3.14159265358979323846;
    switch(direction) {
      case FaceDirection2D::RadialMinus: --ni;nr=-1;area=2*pi*grid_->radial_face(i)*grid_->dz();break;
      case FaceDirection2D::RadialPlus: ++ni;nr=1;area=2*pi*grid_->radial_face(i+1)*grid_->dz();break;
      case FaceDirection2D::AxialMinus: --nj;nz=-1;area=pi*(std::pow(grid_->radial_face(i+1),2)-std::pow(grid_->radial_face(i),2));break;
      case FaceDirection2D::AxialPlus: ++nj;nz=1;area=pi*(std::pow(grid_->radial_face(i+1),2)-std::pow(grid_->radial_face(i),2));break;
    }
    if (ni<0) return {FaceRole2D::Axis,ElectrodeCellType::Gas,ni,nj,nr,nz,area};
    if (ni>=grid_->nr()||nj<0||nj>=grid_->nz())
      return {FaceRole2D::OuterBoundary,ElectrodeCellType::Gas,ni,nj,nr,nz,area};
    const auto material=classify(ni,nj);
    return {material==ElectrodeCellType::Gas?FaceRole2D::GasInterface:FaceRole2D::Electrode,
            material,ni,nj,nr,nz,area};
  }
 private:
  Geometry2DView(const AxisymmetricGrid& g,std::vector<ElectrodeCellType> cells,
                 std::string id,std::string rev,std::string reference,double return_surface_z,
                 double hv_surface_z,double seed_reference_z,double radius,std::string family)
      :grid_(&g),cells_(std::move(cells)),geometry_id_(std::move(id)),revision_(std::move(rev)),
       reference_plane_id_(std::move(reference)),family_(std::move(family)),
       return_surface_z_m_(return_surface_z),hv_surface_z_m_(hv_surface_z),
       seed_reference_z_m_(seed_reference_z),tip_radius_m_(radius) {
    if (!(hv_surface_z_m_>return_surface_z_m_)) throw std::invalid_argument("nonpositive net gas gap");
    bool hv=false,ret=false,gas=false;
    for (const auto role:cells_) {hv|=role==ElectrodeCellType::HighVoltageElectrode;
      ret|=role==ElectrodeCellType::GroundElectrode;gas|=role==ElectrodeCellType::Gas;}
    if (!hv||!ret||!gas) throw std::invalid_argument("geometry lacks gas/HV/return cells on target grid");
  }
  const AxisymmetricGrid* grid_;
  const std::vector<ElectrodeCellType> cells_;
  const std::string geometry_id_,revision_,reference_plane_id_,family_;
  const double return_surface_z_m_,hv_surface_z_m_,seed_reference_z_m_,tip_radius_m_;
};
}  // namespace streamer_rf
