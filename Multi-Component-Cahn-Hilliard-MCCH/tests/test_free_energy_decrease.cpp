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
    std::vector<std::vector<double>> W = {
        {0.0, 1.5, 1.5},
        {1.5, 0.0, 1.5},
        {1.5, 1.5, 0.0}
    };
    auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
    auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

    std::vector<std::vector<double>> kappa = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
    solver.initialize_random({0.33, 0.33, 0.34}, 0.08, 999);

    double dt = 0.02; // well within CFL
    double current_energy = solver.compute_total_free_energy();
    double initial_energy = current_energy;

    int increases = 0;
    for (int step = 1; step <= 150; ++step) {
        solver.step(dt, mcch::TimeIntegratorType::RK2);
        double new_energy = solver.compute_total_free_energy();
        if (new_energy > current_energy + 1e-9) {
            increases++;
        }
        current_energy = new_energy;
    }

    if (rank == 0) {
        std::cout << "Initial Free Energy: " << initial_energy << "\n";
        std::cout << "Final Free Energy:   " << current_energy << "\n";
        std::cout << "Energy Reduction:    " << (initial_energy - current_energy) << "\n";
        std::cout << "Increases detected:  " << increases << "\n";
    }

    assert(increases == 0);
    assert(current_energy < initial_energy);

    if (rank == 0) {
        std::cout << "[PASS] Thermodynamic consistency: Monotonic free energy dissipation confirmed.\n";
    }

    MPI_Finalize();
    return 0;
}
