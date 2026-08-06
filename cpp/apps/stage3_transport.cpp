#include "streamer_rf/streamer/MorrowLowke.hpp"
#include <petscsys.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace streamer_rf::streamer;
int main(int argc,char**argv){
 PetscInitialize(&argc,&argv,nullptr,nullptr);if(argc!=2){std::cerr<<"output csv required\n";return 2;}
 double N=101325/(1.380649e-23*300);std::ofstream f(argv[1]);f<<"electric_field_V_m,E_over_N_Td,mobility_m2_V_s,diffusion_m2_s,nu_i_s_1,nu_a2_s_1,nu_a3_s_1\n"<<std::setprecision(17);
 for(int i=0;i<=200;++i){double E=1e4*std::pow(4000.,i/200.);auto c=evaluate_morrow_lowke(E,N,101325,300);f<<E<<','<<E/N/1e-21<<','<<c.mobility<<','<<c.diffusion<<','<<c.ionization_frequency<<','<<c.attachment_two_body_frequency<<','<<c.attachment_three_body_frequency<<'\n';}
 double ek=morrow_lowke_breakdown_field(N,101325,300);std::cout<<std::setprecision(17)<<"breakdown_field_V_m="<<ek<<'\n';PetscFinalize();return 0;
}
