#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"
#include "Solver.hpp"

#include <iostream>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nx = 48, ny = 48, nz = 1;
    mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
    mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

    int n_comp = 3;
    std::vector<std::vector<double>> W(n_comp, std::vector<double>(n_comp, 1.0));
    for (int i = 0; i < n_comp; ++i) W[i][i] = 0.0;
    auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);

    std::vector<std::vector<double>> kappa(n_comp, std::vector<double>(n_comp, 0.0));
    for (int i = 0; i < n_comp; ++i) kappa[i][i] = 1.0;

    // Test Case A: Constant mobility
    {
        auto mob_const = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);
        mcch::Solver solver(ctx, grid, n_comp, fe, mob_const, kappa);
        solver.initialize_random({0.3, 0.3, 0.4}, 0.05, 12345);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.01;

        for (int step = 0; step < 100; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);

        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            if (rank == 0) {
                std::cout << "[Constant Mobility] Component " << m 
                          << " mass drift: " << drift << "\n";
            }
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);
    }

    // Test Case B: Degenerate mobility
    {
        auto mob_deg = std::make_shared<mcch::DegenerateMobility>(n_comp, 1.0);
        mcch::Solver solver(ctx, grid, n_comp, fe, mob_deg, kappa);
        solver.initialize_random({0.35, 0.35, 0.3}, 0.05, 54321);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.01;

        for (int step = 0; step < 100; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);

        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            if (rank == 0) {
                std::cout << "[Degenerate Mobility] Component " << m 
                          << " mass drift: " << drift << "\n";
            }
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);
    }

    if (rank == 0) {
        std::cout << "[PASS] Mass conservation and sum-to-one tests passed across " << size << " ranks.\n";
    }

    MPI_Finalize();
    return 0;
}
