#include "streamer_rf/streamer/MorrowLowke.hpp"
#include "streamer_rf/streamer/ReactionModel.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
constexpr double kB = 1.380649e-23;

struct Case {
  const char* id;
  double td;
  double ne;
  double np;
  double nn;
  double dt;
};
}  // namespace

int main(int argc, char** argv) {
  namespace fs = std::filesystem;
  const fs::path out = argc > 1 ? fs::path(argv[1])
                                : fs::path("solver3d/afivo_reference/common_benchmark/d3/results/d3_reaction_reference.csv");
  fs::create_directories(out.parent_path());
  constexpr double pressure = 101325.0;
  constexpr double temperature = 300.0;
  constexpr double neutral_density = pressure / (kB * temperature);
  constexpr double volume = 4.0e-5 * 4.0e-5 * 4.0e-5;
  const Case cases[] = {{"ionization", 250.0, 1.0e16, 1.0e16, 0.0, 1.0e-14},
                        {"attachment", 20.0, 1.0e16, 1.0e16, 0.0, 1.0e-13},
                        {"recombination", 90.0, 1.0e18, 2.0e18, 1.0e18, 1.0e-13}};
  std::ofstream f(out);
  f << std::setprecision(17)
    << "case_id,E_over_N_Td,E_Vpm,ne_m3,np_m3,nn_m3,dt_s,volume_m3,"
       "dne_dt_m3s,dnp_dt_m3s,dnn_dt_m3s,initial_ne_amount,initial_np_amount,"
       "initial_nn_amount,expected_ne_amount,expected_np_amount,expected_nn_amount\n";
  for (const auto& c : cases) {
    const double efield = c.td * 1e-21 * neutral_density;
    const auto tr = streamer_rf::streamer::evaluate_morrow_lowke(efield, neutral_density, pressure, temperature);
    const auto s = streamer_rf::streamer::evaluate_reactions(c.ne, c.np, c.nn, 0.0, tr, temperature);
    f << c.id << ',' << c.td << ',' << efield << ',' << c.ne << ',' << c.np << ',' << c.nn << ','
      << c.dt << ',' << volume << ',' << s.electron << ',' << s.positive_ion << ',' << s.negative_ion << ','
      << c.ne * volume << ',' << c.np * volume << ',' << c.nn * volume << ','
      << (c.ne + c.dt * s.electron) * volume << ','
      << (c.np + c.dt * s.positive_ion) * volume << ','
      << (c.nn + c.dt * s.negative_ion) * volume << '\n';
  }
  std::cout << "stage_d3_reaction_reference status=PASS out=" << out << '\n';
  return 0;
}
