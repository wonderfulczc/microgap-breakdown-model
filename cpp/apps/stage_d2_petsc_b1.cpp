#include "streamer_rf/electrode.hpp"
#include "streamer_rf/types.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <petscksp.h>

namespace {
struct Options {
  std::string resolution{"coarse"};
  std::filesystem::path out{"solver3d/afivo_reference/common_benchmark/petsc/generated/d2_b1_petsc_coarse.csv"};
};

Options parse(int argc, char** argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--resolution" && i + 1 < argc) {
      o.resolution = argv[++i];
    } else if (a == "--out" && i + 1 < argc) {
      o.out = argv[++i];
    } else {
      throw std::invalid_argument("usage: stage_d2_petsc_b1 [--resolution coarse|medium] [--out file.csv]");
    }
  }
  return o;
}

int idx(const streamer_rf::AxisymmetricGrid& grid, int i, int j) {
  return j * grid.nr() + i;
}

void solve_b1_homogeneous(const streamer_rf::AxisymmetricGrid& grid,
                          const streamer_rf::AxisymmetricNeedlePlaneGeometry& geom,
                          double voltage, streamer_rf::ScalarField2D& phi,
                          int& iterations) {
  const int n = grid.nr() * grid.nz();
  Mat a;
  Vec b, x;
  KSP ksp;
  MatCreate(PETSC_COMM_WORLD, &a);
  MatSetSizes(a, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(a);
  MatMPIAIJSetPreallocation(a, 5, nullptr, 5, nullptr);
  MatSeqAIJSetPreallocation(a, 5, nullptr);
  MatSetUp(a);
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetSizes(b, PETSC_DECIDE, n);
  VecSetFromOptions(b);
  VecDuplicate(b, &x);

  const double dr = grid.dr();
  const double dz = grid.dz();
  for (int j = 0; j < grid.nz(); ++j) {
    for (int i = 0; i < grid.nr(); ++i) {
      const PetscInt row = idx(grid, i, j);
      const auto c = geom.classify(grid, i, j);
      const bool z_low = j == 0;
      const bool z_high = j == grid.nz() - 1;
      const bool r_outer = i == grid.nr() - 1;
      if (c != streamer_rf::ElectrodeCellType::Gas || z_low || z_high) {
        const double value = c == streamer_rf::ElectrodeCellType::HighVoltageElectrode || z_high ? voltage : 0.0;
        const PetscScalar one = 1.0;
        MatSetValue(a, row, row, one, INSERT_VALUES);
        VecSetValue(b, row, value, INSERT_VALUES);
      } else if (r_outer) {
        const PetscInt cols[2] = {row, idx(grid, i - 1, j)};
        const PetscScalar vals[2] = {1.0, -1.0};
        MatSetValues(a, 1, &row, 2, cols, vals, INSERT_VALUES);
        VecSetValue(b, row, 0.0, INSERT_VALUES);
      } else {
        const double r = grid.r(i);
        const double rp = grid.radial_face(i + 1);
        const double rm = grid.radial_face(i);
        const double ar = rp / (r * dr * dr);
        const double al = rm / (r * dr * dr);
        const double az = 1.0 / (dz * dz);
        PetscInt cols[5];
        PetscScalar vals[5];
        int m = 0;
        cols[m] = row;
        vals[m++] = -(ar + al + 2.0 * az);
        cols[m] = idx(grid, i + 1, j);
        vals[m++] = ar;
        if (i > 0) {
          cols[m] = idx(grid, i - 1, j);
          vals[m++] = al;
        }
        cols[m] = idx(grid, i, j + 1);
        vals[m++] = az;
        cols[m] = idx(grid, i, j - 1);
        vals[m++] = az;
        MatSetValues(a, 1, &row, m, cols, vals, INSERT_VALUES);
        VecSetValue(b, row, 0.0, INSERT_VALUES);
      }
    }
  }
  MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, a, a);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, b, x);
  PetscInt its = 0;
  KSPGetIterationNumber(ksp, &its);
  iterations = static_cast<int>(its);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  if (reason < 0) throw std::runtime_error("D2 PETSc B1 solve diverged");
  const PetscScalar* values = nullptr;
  VecGetArrayRead(x, &values);
  for (int j = 0; j < grid.nz(); ++j) {
    for (int i = 0; i < grid.nr(); ++i) {
      phi(i, j) = PetscRealPart(values[idx(grid, i, j)]);
    }
  }
  VecRestoreArrayRead(x, &values);
  KSPDestroy(&ksp);
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&a);
}
}  // namespace

int main(int argc, char** argv) {
  bool petsc_initialized = false;
  try {
    const auto opt = parse(argc, argv);
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    petsc_initialized = true;
    int nr = 32;
    int nz = 74;
    if (opt.resolution == "medium") {
      nr = 64;
      nz = 148;
    } else if (opt.resolution != "coarse") {
      throw std::invalid_argument("resolution must be coarse or medium");
    }

    constexpr double voltage = 500.0;
    constexpr double r_max = 80e-6;
    constexpr double z_min = -2.5e-6;
    constexpr double z_max = 90e-6;
    constexpr double ground_z = -2.5e-6;
    constexpr double ground_thickness = 2.5e-6;
    constexpr double tip_center_z = 75e-6;
    constexpr double tip_radius = 5e-6;
    constexpr double shank_radius = 5e-6;
    constexpr double gas_gap = 70e-6;

    streamer_rf::AxisymmetricGrid grid(nr, nz, r_max, z_min, z_max);
    streamer_rf::AxisymmetricNeedlePlaneGeometry geom(
        "D2-B1-axisymmetric-electrostatic", ground_z, tip_center_z,
        tip_radius, shank_radius, ground_thickness);
    streamer_rf::ScalarField2D phi(grid), er(grid), ez(grid), emag(grid);
    int poisson_iterations = 0;
    solve_b1_homogeneous(grid, geom, voltage, phi, poisson_iterations);

    for (int j = 0; j < grid.nz(); ++j) {
      for (int i = 0; i < grid.nr(); ++i) {
        auto dr = [&](int a, int b) { return (phi(b, j) - phi(a, j)) / ((b - a) * grid.dr()); };
        auto dz = [&](int a, int b) { return (phi(i, b) - phi(i, a)) / ((b - a) * grid.dz()); };
        er(i, j) = i == 0 ? 0 : -(i == grid.nr() - 1 ? dr(i - 1, i) : dr(i - 1, i + 1));
        ez(i, j) = -(j == 0 ? dz(0, 1) : j == grid.nz() - 1 ? dz(j - 1, j) : dz(j - 1, j + 1));
        emag(i, j) = std::hypot(er(i, j), ez(i, j));
      }
    }

    double emax = 0.0;
    for (int j = 0; j < grid.nz(); ++j) {
      for (int i = 0; i < grid.nr(); ++i) {
        if (geom.classify(grid, i, j) == streamer_rf::ElectrodeCellType::Gas) {
          emax = std::max(emax, emag(i, j));
        }
      }
    }

    int rank = 0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
      std::filesystem::create_directories(opt.out.parent_path());
      std::ofstream f(opt.out);
      f << std::setprecision(17);
      f << "# stage_d2_petsc_b1 resolution=" << opt.resolution << " nr=" << nr << " nz=" << nz
        << " voltage_V=" << voltage << " emax_gas_Vpm=" << emax
        << " poisson_iterations=" << poisson_iterations
        << " boundary=homogeneous_top_voltage_side_neumann_bottom_ground\n";
      f << "z_m,z_over_gap,phi_V,Ez_Vpm,Eabs_Vpm,cell_class\n";
      for (int j = 0; j < grid.nz(); ++j) {
        const int i = 0;
        const auto c = geom.classify(grid, i, j);
        const char* cls = c == streamer_rf::ElectrodeCellType::Gas
                              ? "gas"
                              : (c == streamer_rf::ElectrodeCellType::HighVoltageElectrode ? "hv" : "ground");
        f << grid.z(j) << ',' << grid.z(j) / gas_gap << ',' << phi(i, j) << ','
          << ez(i, j) << ',' << emag(i, j) << ',' << cls << '\n';
      }
      std::cout << "stage_d2_petsc_b1 status=PASS resolution=" << opt.resolution
                << " out=" << opt.out << " emax_gas_Vpm=" << emax
                << " poisson_iterations=" << poisson_iterations << '\n';
    }
  } catch (const std::exception& e) {
    std::cerr << "stage_d2_petsc_b1 status=FAIL error=" << e.what() << '\n';
    if (petsc_initialized) PetscFinalize();
    return 2;
  }
  PetscFinalize();
  return 0;
}
