#include "streamer_rf/streamer/ColdPhysics.hpp"
#include <algorithm>
#include <stdexcept>

namespace streamer_rf::streamer {
namespace {
void require(bool ok,const char* reason){if(!ok)throw std::invalid_argument(reason);}
bool present(const std::string& s){return !s.empty();}
}
void IonFieldTable::validate()const{
  require(present(species_id)&&present(gas_composition)&&present(source_doi)&&present(database_id)&&
          file_sha256.size()==64&&(evidence_class=="A"||evidence_class=="B"||evidence_class=="C"),"ION_TABLE_PROVENANCE_MISSING");
  require(field_unit=="Td"&&diffusion_unit=="m2/s"&&
          (mobility_unit=="m2/(V*s)"||mobility_unit=="K0_m2/(V*s)"),"ION_TABLE_UNIT_MISMATCH");
  require(interpolation=="LINEAR_EN"&&out_of_domain=="REJECT"&&aggregate_mapping=="EFFECTIVE_SPECIES_APPROXIMATION","ION_TABLE_METHOD_UNSUPPORTED");
  require(EN_Td.size()>=2&&mobility.size()==EN_Td.size(),"ION_TABLE_SHAPE");
  require(D_parallel.size()==EN_Td.size()&&D_perpendicular.size()==EN_Td.size(),"ION_DIFFUSION_UNAVAILABLE");
  if(mobility_unit=="K0_m2/(V*s)")require(std::isfinite(reference_density_m3)&&reference_density_m3>0,"REDUCED_MOBILITY_REFERENCE_REQUIRED");
  for(std::size_t k=0;k<EN_Td.size();++k){
    require(std::isfinite(EN_Td[k])&&EN_Td[k]>=0&&(k==0||EN_Td[k]>EN_Td[k-1]),"ION_TABLE_NONMONOTONE");
    require(std::isfinite(mobility[k])&&mobility[k]>0&&std::isfinite(D_parallel[k])&&D_parallel[k]>=0&&
            std::isfinite(D_perpendicular[k])&&D_perpendicular[k]>=0,"ION_TABLE_INVALID_VALUE");
  }
}
void IonTransportModel::validate(double EN,double T,double p)const{
  require(present(model_id)&&present(source)&&present(species_id)&&present(approximation)&&
          present(qualification),"ION_PARAMETER_PROVENANCE_MISSING");
  if(model_id=="SOURCE_BOUND_EN_TABLE"){
    table.validate();
    require(approximation=="EFFECTIVE_SPECIES_APPROXIMATION"&&species_id==table.species_id,"ION_TABLE_SPECIES_MISMATCH");
    require(EN>=table.EN_Td.front()&&EN<=table.EN_Td.back(),"ION_TABLE_OUT_OF_DOMAIN");
  }else{
    require(model_id=="CONSTANT_EFFECTIVE_POSITIVE_ION"&&approximation=="CONSTANT_MU_D","ION_TRANSPORT_MODEL_NOT_IMPLEMENTED");
    require(std::isfinite(mobility_m2_Vs)&&mobility_m2_Vs>0&&std::isfinite(diffusion_m2_s)&&diffusion_m2_s>=0,"INVALID_ION_TRANSPORT");
  }
  require(std::isfinite(min_EN_Td)&&std::isfinite(max_EN_Td)&&min_EN_Td>=0&&max_EN_Td>=min_EN_Td&&
          std::isfinite(min_temperature_K)&&std::isfinite(max_temperature_K)&&min_temperature_K>0&&
          max_temperature_K>=min_temperature_K&&std::isfinite(min_pressure_Pa)&&
          std::isfinite(max_pressure_Pa)&&min_pressure_Pa>0&&max_pressure_Pa>=min_pressure_Pa,
          "ION_VALIDITY_DOMAIN_MISSING");
  require(std::isfinite(EN)&&EN>=min_EN_Td&&EN<=max_EN_Td&&std::isfinite(T)&&
          T>=min_temperature_K&&T<=max_temperature_K&&std::isfinite(p)&&
          p>=min_pressure_Pa&&p<=max_pressure_Pa,"ION_TRANSPORT_OUT_OF_MODEL_DOMAIN");
}
IonTransportCoefficients PositiveIonTransportProvider::evaluate(double EN,double T,double p,double N)const{
  model_.validate(EN,T,p);require(std::isfinite(N)&&N>0,"INVALID_NEUTRAL_DENSITY");
  if(model_.model_id=="CONSTANT_EFFECTIVE_POSITIVE_ION")return {model_.mobility_m2_Vs,model_.diffusion_m2_s,model_.diffusion_m2_s};
  const auto& t=model_.table;
  auto hi=std::upper_bound(t.EN_Td.begin(),t.EN_Td.end(),EN);
  std::size_t k=hi==t.EN_Td.end()?t.EN_Td.size()-2:std::size_t(hi-t.EN_Td.begin()-1);
  const double a=(EN-t.EN_Td[k])/(t.EN_Td[k+1]-t.EN_Td[k]);
  auto interp=[&](const std::vector<double>& v){return v[k]+a*(v[k+1]-v[k]);};
  double mu=interp(t.mobility);if(t.mobility_unit=="K0_m2/(V*s)")mu*=t.reference_density_m3/N;
  return {mu,interp(t.D_parallel),interp(t.D_perpendicular)};
}
void SurfaceEmissionModel::validate_source()const{
  if(!source_contract_required)return;
  require(present(incident_species_id)&&present(gas_composition)&&present(cathode_material)&&present(surface_state)&&
          present(incident_energy_assumption)&&source_sha256.size()==64&&present(uncertainty),"SURFACE_SOURCE_CONTRACT_MISSING");
  require(evidence_class=="A"||evidence_class=="B"||evidence_class=="C","SURFACE_PARAMETERS_UNKNOWN");
  require(energy_status=="UNAVAILABLE","SURFACE_ENERGY_NOT_IMPLEMENTED");
  if(evidence_class!="C"){
    require(cathode_material!="UNKNOWN"&&surface_state!="UNKNOWN"&&std::isfinite(incident_energy_eV)&&
            std::isfinite(min_energy_eV)&&std::isfinite(max_energy_eV)&&min_energy_eV>=0&&
            incident_energy_eV>=min_energy_eV&&incident_energy_eV<=max_energy_eV,"SURFACE_ENERGY_DOMAIN_UNAVAILABLE");
  }
}
double SurfaceEmissionModel::emitted_flux(double incident,ElectrodeCellType metal,double voltage)const{
  require(std::isfinite(incident)&&incident>=0&&std::isfinite(voltage),"INVALID_WALL_INPUT");
  if(!enabled)return 0;
  validate_source();
  require(model_id=="ION_IMPACT_EFFECTIVE_YIELD"&&present(source)&&present(qualification)&&
          std::isfinite(gamma_i)&&gamma_i>=0,"SURFACE_PARAMETER_PROVENANCE_MISSING");
  require(metal!=ElectrodeCellType::Gas,"EMISSION_REQUIRES_METAL_FACE");
  const bool cathode=(voltage>0&&metal==ElectrodeCellType::GroundElectrode)||
                     (voltage<0&&metal==ElectrodeCellType::HighVoltageElectrode);
  return cathode?gamma_i*incident:0;
}
void ColdPhysicsConfig::validate(bool geometry,bool photo)const{
  if(mode==ColdPhysicsMode::LegacyElectronOnly){
    require(!positive_ion_transport&&!surface.enabled&&!local_charge_port_v2&&
            photon_boundary==PhotonMetalBoundary::Unavailable,"LEGACY_CAPABILITY_CONFLICT");
    return;
  }
  require(geometry,"V3_REQUIRES_GEOMETRY_VIEW");
  require(!surface.enabled||positive_ion_transport,"SECONDARY_REQUIRES_REAL_ION_FLUX");
  require(!local_charge_port_v2||positive_ion_transport,"PORT_V2_REQUIRES_ENABLED_SPECIES");
  if(positive_ion_transport)ions.validate(ions.min_EN_Td,ions.min_temperature_K,ions.min_pressure_Pa);
  if(surface.enabled)(void)surface.emitted_flux(0,ElectrodeCellType::GroundElectrode,1);
  if(surface.enabled&&surface.source_contract_required)
    require(surface.incident_species_id==ions.species_id,"SURFACE_ION_IDENTITY_MISMATCH");
  require(!photo||(photon_boundary==PhotonMetalBoundary::OpaqueAbsorbingLarsen&&
                   !excitation_source.empty()),"GEOMETRY_SP3_BOUNDARY_OR_SOURCE_UNAVAILABLE");
}
IonTransportResult positive_ion_transport(const Geometry2DView& v,const ScalarField2D& np,
    const ScalarField2D& er,const ScalarField2D& ez,const IonTransportModel& m,
    const SurfaceEmissionModel& emission,double voltage,double N,double T,double p,bool capture){
  const auto& g=v.grid();
  require(&np.grid()==&g&&&er.grid()==&g&&&ez.grid()==&g,"ION_GRID_IDENTITY_MISMATCH");
  require(std::isfinite(N)&&N>0,"INVALID_NEUTRAL_DENSITY");
  IonTransportResult out(g);
  const bool table_mode=m.model_id=="SOURCE_BOUND_EN_TABLE";
  std::vector<IonTransportCoefficients> coeff(g.size());
  if(table_mode){
    PositiveIonTransportProvider provider(m);
    for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)if(v.classify(i,j)==ElectrodeCellType::Gas){
      auto a=provider.evaluate(std::hypot(er(i,j),ez(i,j))/N*1e21,T,p,N);
      require(a.D_parallel==a.D_perpendicular,"ANISOTROPIC_DIFFUSION_FV_NOT_IMPLEMENTED");
      coeff[j*g.nr()+i]=a;
    }
  }
  for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i){
    if(v.classify(i,j)!=ElectrodeCellType::Gas)continue;
    require(std::isfinite(np(i,j))&&np(i,j)>=0,"INVALID_ION_DENSITY");
    if(!table_mode)m.validate(std::hypot(er(i,j),ez(i,j))/N*1e21,T,p);
    const double mu=table_mode?coeff[j*g.nr()+i].mobility:m.mobility_m2_Vs;
    const double D=table_mode?coeff[j*g.nr()+i].D_parallel:m.diffusion_m2_s;
    double outgoing=0;
    for(auto dir:{FaceDirection2D::RadialMinus,FaceDirection2D::RadialPlus,
                 FaceDirection2D::AxialMinus,FaceDirection2D::AxialPlus}){
      const auto f=v.face(i,j,dir);
      if(f.role==FaceRole2D::Axis)continue;
      const double h=f.normal_r!=0?g.dr():g.dz();
      double vn=mu*(er(i,j)*f.normal_r+ez(i,j)*f.normal_z);
      double coefficient=0,flux=0,secondary=0;
      if(f.role==FaceRole2D::GasInterface){
        vn=.5*m.mobility_m2_Vs*((er(i,j)+er(f.neighbor_i,f.neighbor_j))*f.normal_r+
                              (ez(i,j)+ez(f.neighbor_i,f.neighbor_j))*f.normal_z);
        double Dface=D;
        if(table_mode){
          const auto b=coeff[f.neighbor_j*g.nr()+f.neighbor_i];
          vn=.5*(mu*(er(i,j)*f.normal_r+ez(i,j)*f.normal_z)+b.mobility*(er(f.neighbor_i,f.neighbor_j)*f.normal_r+ez(f.neighbor_i,f.neighbor_j)*f.normal_z));
          Dface=.5*(D+b.D_parallel);
        }
        coefficient=std::max(vn,0.)+Dface/h;
        flux=coefficient*np(i,j)+(std::min(vn,0.)-Dface/h)*np(f.neighbor_i,f.neighbor_j);
      }else{
        coefficient=std::max(vn,0.)+(f.role==FaceRole2D::Electrode?2*D/h:0);
        flux=coefficient*np(i,j);
        if(f.role==FaceRole2D::Electrode){
          secondary=emission.emitted_flux(flux,f.neighbor_material,voltage);
          if(f.neighbor_material==ElectrodeCellType::HighVoltageElectrode){
            out.hv_collected_s1+=flux*f.area_m2;out.hv_emitted_s1+=secondary*f.area_m2;
          }else{out.return_collected_s1+=flux*f.area_m2;out.return_emitted_s1+=secondary*f.area_m2;}
        }else out.outer_lost_s1+=flux*f.area_m2;
      }
      out.rhs(i,j)-=flux*f.area_m2/g.cell_volume(i);
      out.secondary_rhs(i,j)+=secondary*f.area_m2/g.cell_volume(i);
      outgoing+=coefficient*f.area_m2/g.cell_volume(i);
      if(capture)out.faces.push_back({i,j,dir,f,flux,secondary});
    }
    if(outgoing>0)out.positivity_dt_s=std::min(out.positivity_dt_s,.45/outgoing);
  }
  return out;
}
}
