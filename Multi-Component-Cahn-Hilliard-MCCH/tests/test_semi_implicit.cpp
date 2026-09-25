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

    // =========================================================================
    // Test 1: Semi-Implicit Euler (Mass Conservation & Sum-to-One)
    // =========================================================================
    {
        int nx = 48, ny = 48, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 3;
        std::vector<std::vector<double>> W = {
            {0.0, 1.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.3, 0.3, 0.4}, 0.05, 12345);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.05;

        for (int step = 0; step < 60; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::SEMI_IMPLICIT);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);

        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            if (rank == 0) {
                std::cout << "[Semi-Implicit Euler] Component " << m 
                          << " mass drift: " << drift << "\n";
            }
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 1: Semi-Implicit Euler mass conservation and sum-to-one passed.\n";
        }
    }

    // =========================================================================
    // Test 2: Semi-Implicit BDF2 (SBDF2) Multi-step Scheme
    // =========================================================================
    {
        int nx = 48, ny = 48, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 3;
        std::vector<std::vector<double>> W = {
            {0.0, 1.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.35, 0.35, 0.30}, 0.05, 54321);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.05;

        for (int step = 0; step < 60; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::SEMI_IMPLICIT_BDF2);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);

        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            if (rank == 0) {
                std::cout << "[SBDF2] Component " << m << " mass drift: " << drift << "\n";
            }
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 2: Semi-Implicit BDF2 (SBDF2) passed.\n";
        }
    }

    // =========================================================================
    // Test 3: Large Time Step Stability (Overcoming Explicit Biharmonic CFL)
    // =========================================================================
    {
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 3;
        std::vector<std::vector<double>> W = {
            {0.0, 1.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.05, 777);

        // dt = 0.2 is 10x larger than explicit RK2 stable limit (~0.02)
        double dt_large = 0.2;

        for (int step = 0; step < 40; ++step) {
            solver.step(dt_large, mcch::TimeIntegratorType::SEMI_IMPLICIT);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt_large);
        assert(!std::isnan(d_final.total_free_energy));
        assert(d_final.max_composition <= 1.2);
        assert(d_final.min_composition >= -0.2);
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 3: Large time step stability (dt = " << dt_large 
                      << " >> explicit CFL) passed.\n";
        }
    }

    // =========================================================================
    // Test 4: 3D Multi-Component Semi-Implicit Simulation
    // =========================================================================
    {
        int nx = 24, ny = 24, nz = 24;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(3, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 3;
        std::vector<std::vector<double>> W = {
            {0.0, 1.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.02, 888);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.05;

        for (int step = 0; step < 20; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::SEMI_IMPLICIT);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);
        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            (void)drift;
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 4: 3D Multi-component Semi-Implicit simulation passed.\n";
        }
    }

    // =========================================================================
    // Test 5: 4-Component (Quaternary) Semi-Implicit Simulation
    // =========================================================================
    {
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 4;
        std::vector<std::vector<double>> W(n_comp, std::vector<double>(n_comp, 1.0));
        for (int i = 0; i < n_comp; ++i) W[i][i] = 0.0;
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

        std::vector<std::vector<double>> kappa(n_comp, std::vector<double>(n_comp, 0.0));
        for (int i = 0; i < n_comp; ++i) kappa[i][i] = 1.0;

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.25, 0.25, 0.25, 0.25}, 0.02, 999);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.05;

        for (int step = 0; step < 30; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::SEMI_IMPLICIT);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);
        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            (void)drift;
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 5: 4-Component Semi-Implicit simulation passed.\n";
        }
    }

    // =========================================================================
    // Test 6: Multicomponent Semi-Implicit with Mobility Matrix
    // =========================================================================
    {
        int nx = 32, ny = 32, nz = 1;
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, true);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        int n_comp = 3;
        std::vector<std::vector<double>> W = {
            {0.0, 1.0, 1.0},
            {1.0, 0.0, 1.0},
            {1.0, 1.0, 0.0}
        };
        auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W);
        std::vector<std::vector<double>> M = {
            {1.0, 0.1, 0.1},
            {0.1, 1.0, 0.1},
            {0.1, 0.1, 1.0}
        };
        auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, M);

        std::vector<std::vector<double>> kappa = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        mcch::Solver solver(ctx, grid, n_comp, fe, mob, kappa);
        solver.initialize_random({0.33, 0.33, 0.34}, 0.02, 555);

        mcch::Diagnostics d0 = solver.compute_diagnostics(0.0);
        double dt = 0.05;

        for (int step = 0; step < 20; ++step) {
            solver.step(dt, mcch::TimeIntegratorType::SEMI_IMPLICIT);
        }

        mcch::Diagnostics d_final = solver.compute_diagnostics(dt);
        for (int m = 0; m < n_comp; ++m) {
            double drift = std::abs(d_final.average_composition[m] - d0.average_composition[m]);
            if (rank == 0) {
                std::cout << "[Test 6] Component " << m << " drift: " << drift << "\n";
            }
            assert(drift < 1e-12);
        }
        assert(d_final.max_unity_deviation < 1e-12);

        if (rank == 0) {
            std::cout << "[PASS] Test 6: Multicomponent Semi-Implicit with Mobility Matrix passed.\n";
        }
    }

    if (rank == 0) {
        std::cout << "\n=======================================================\n";
        std::cout << "  ALL SEMI-IMPLICIT FFTW MPI TESTS PASSED (" << size << " RANKS)\n";
        std::cout << "=======================================================\n";
    }

    MPI_Finalize();
    return 0;
}
