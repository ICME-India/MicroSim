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

    // 1. Quaternary 2D Test (N = 4)
    {
        int n_comp = 4;
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        std::vector<std::vector<double>> W(n_comp, std::vector<double>(n_comp, 1.0));
        for (int i = 0; i < n_comp; ++i) W[i][i] = 0.0;
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa(n_comp, std::vector<double>(n_comp, 0.0));
        for (int i = 0; i < n_comp; ++i) kappa[i][i] = 1.0;

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.25, 0.25, 0.25, 0.25}, 0.02, 777);

        for (int s = 0; s < 50; ++s) {
            solver.step(0.02, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics diag = solver.compute_diagnostics(0.02);
        assert(diag.max_unity_deviation < 1e-12);
        if (rank == 0) {
            std::cout << "[PASS] 4-component 2D solver executed successfully.\n";
        }
    }

    // 2. 3D Simulation Test
    {
        int n_comp = 3;
        int nx = 24, ny = 24, nz = 24;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        mcch::Grid grid(3, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        std::vector<std::vector<double>> W(n_comp, std::vector<double>(n_comp, 1.0));
        for (int i = 0; i < n_comp; ++i) W[i][i] = 0.0;
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa(n_comp, std::vector<double>(n_comp, 0.0));
        for (int i = 0; i < n_comp; ++i) kappa[i][i] = 1.0;

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.02, 888);

        for (int s = 0; s < 30; ++s) {
            solver.step(0.015, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics diag = solver.compute_diagnostics(0.015);
        assert(diag.max_unity_deviation < 1e-12);
        if (rank == 0) {
            std::cout << "[PASS] 3D multi-component slab MPI solver executed successfully.\n";
        }
    }

    // 3. Regular Solution Thermodynamics Test
    {
        int n_comp = 3;
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        std::vector<std::vector<double>> omega = {
            {0.0, 3.5, 3.5},
            {3.5, 0.0, 3.5},
            {3.5, 3.5, 0.0}
        };
        auto fe = std::make_shared<mcch::RegularSolution>(n_comp, 1.0, omega);
        auto mob = std::make_shared<mcch::DegenerateMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa(n_comp, std::vector<double>(n_comp, 0.0));
        for (int i = 0; i < n_comp; ++i) kappa[i][i] = 1.0;

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.02, 555);

        for (int s = 0; s < 40; ++s) {
            solver.step(0.01, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics diag = solver.compute_diagnostics(0.01);
        assert(diag.max_unity_deviation < 1e-12);
        if (rank == 0) {
            std::cout << "[PASS] Regular solution model with degenerate mobility passed.\n";
        }
    }

    // 4. Multicomponent Constant Mobility Matrix Test
    {
        int n_comp = 3;
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        std::vector<std::vector<double>> W = {
            {0.0, 1.5, 1.5},
            {1.5, 0.0, 1.5},
            {1.5, 1.5, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);

        // Constant mobility matrix for multicomponent system
        std::vector<std::vector<double>> M = {
            {1.2, 0.1, 0.0},
            {0.1, 0.9, 0.1},
            {0.0, 0.1, 1.1}
        };
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, M);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.02, 333);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.01;
        for (int s = 0; s < 50; ++s) {
            solver.step(dt, mcch::TimeIntegratorType::RK2);
        }

        mcch::Diagnostics diag = solver.compute_diagnostics(dt);
        assert(diag.max_unity_deviation < 1e-12);
        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(diag.average_composition[m] - d0.average_composition[m]);
            (void)drift;
            assert(drift < 1e-12);
        }
        if (rank == 0) {
            std::cout << "[PASS] Multi-component system with constant mobility matrix passed.\n";
        }
    }

    MPI_Finalize();
    return 0;
}
