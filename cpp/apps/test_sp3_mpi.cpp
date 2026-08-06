#include "streamer_rf/sp3.hpp"
#include <petscsys.h>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace streamer_rf;
int main(int argc,char**argv){
 PetscInitialize(&argc,&argv,nullptr,nullptr);int rank,size;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&size);
 AxisymmetricGrid g(12,20,.01,-.01,.01);ScalarField2D src(g),sph(g);for(int j=0;j<g.nz();++j)for(int i=0;i<g.nr();++i)src(i,j)=1e20*std::exp(-(g.r(i)*g.r(i)+g.z(j)*g.z(j))/1e-6);
 PoissonBoundaryConfig bc;solve_photoionization(src,{},bc,sph,{1e-12,1e-16,20000});double checksum=0;for(std::size_t k=0;k<sph.values().size();++k)checksum+=(k+1)*sph.values()[k];double lo,hi;MPI_Allreduce(&checksum,&lo,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);MPI_Allreduce(&checksum,&hi,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
 if(rank==0)std::cout<<std::setprecision(17)<<"sp3 ranks="<<size<<" checksum="<<checksum<<" relative_difference="<<std::abs(hi-lo)/std::max(1.,std::abs(hi))<<'\n';PetscFinalize();return std::abs(hi-lo)/std::max(1.,std::abs(hi))<1e-11?0:1;
}
