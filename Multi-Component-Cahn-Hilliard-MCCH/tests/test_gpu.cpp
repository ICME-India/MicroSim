#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"
#include "Solver.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifndef MCCH_ENABLE_CUDA
    if (rank == 0) {
        std::cout << "[SKIP] CUDA support is not enabled in this build.\n";
    }
    MPI_Finalize();
    return 0;
#else
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        if (rank == 0) {
            std::cout << "[SKIP] No CUDA-capable GPU detected on this system.\n";
        }
        MPI_Finalize();
        return 0;
    }

    if (rank == 0) {
        std::cout << "==================================================================\n";
        std::cout << "  Testing GPU Acceleration for Explicit Euler, RK2, and RK4       \n";
        std::cout << "==================================================================\n";
    }

    int nx = 32, ny = 32, nz = 1;
    double dx = 1.0, dy = 1.0, dz = 1.0;
    int n_comp = 3;
    double dt = 0.005;

    mcch::MPIContext mpi(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1, false);
    mcch::Grid grid(2, nx, ny, nz, dx, dy, dz, mpi, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

    std::vector<std::vector<double>> W = {
        {0.0, 1.5, 1.5},
        {1.5, 0.0, 1.5},
        {1.5, 1.5, 0.0}
    };
    std::vector<double> A = {0.2, 0.2, 0.2};
    auto fe = std::make_shared<mcch::PolynomialMultiWell>(n_comp, W, A);
    auto mob = std::make_shared<mcch::ConstantMobility>(n_comp, 1.0);

    std::vector<std::vector<double>> kappa = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    std::vector<double> c_mean = {0.33, 0.33, 0.34};

    // ------------------------------------------------------------------------
    // Test 1: Explicit Euler Comparison (CPU vs GPU)
    // ------------------------------------------------------------------------
    {
        mcch::Solver solver_cpu(mpi, grid, n_comp, fe, mob, kappa);
        mcch::Solver solver_gpu(mpi, grid, n_comp, fe, mob, kappa);

        solver_cpu.initialize_random(c_mean, 0.02, 1234);
        solver_gpu.c = solver_cpu.c;
        solver_gpu.set_device(mcch::DeviceType::GPU);
        solver_gpu.gpu_solver->copy_to_device(solver_gpu.c);

        for (int step = 0; step < 10; ++step) {
            solver_cpu.step_euler(dt);
            solver_gpu.step_euler(dt);
        }

        solver_gpu.sync_from_device();

        double max_diff = 0.0;
        for (int m = 0; m < n_comp; ++m) {
            for (size_t i = 0; i < grid.local_total_size(); ++i) {
                double diff = std::abs(solver_cpu.c[m][i] - solver_gpu.c[m][i]);
                max_diff = std::max(max_diff, diff);
            }
        }
        double global_max_diff = mpi.allreduce_max(max_diff);

        if (rank == 0) {
            std::cout << "[PASS] Explicit Euler CPU vs GPU max diff = " 
                      << std::scientific << global_max_diff << "\n";
        }
        assert(global_max_diff < 1e-10 && "Explicit Euler GPU result differs from CPU");
    }

    // ------------------------------------------------------------------------
    // Test 2: RK2 Comparison (CPU vs GPU)
    // ------------------------------------------------------------------------
    {
        mcch::Solver solver_cpu(mpi, grid, n_comp, fe, mob, kappa);
        mcch::Solver solver_gpu(mpi, grid, n_comp, fe, mob, kappa);

        solver_cpu.initialize_random(c_mean, 0.02, 5678);
        solver_gpu.c = solver_cpu.c;
        solver_gpu.set_device(mcch::DeviceType::GPU);
        solver_gpu.gpu_solver->copy_to_device(solver_gpu.c);

        for (int step = 0; step < 10; ++step) {
            solver_cpu.step_rk2(dt);
            solver_gpu.step_rk2(dt);
        }

        solver_gpu.sync_from_device();

        double max_diff = 0.0;
        for (int m = 0; m < n_comp; ++m) {
            for (size_t i = 0; i < grid.local_total_size(); ++i) {
                double diff = std::abs(solver_cpu.c[m][i] - solver_gpu.c[m][i]);
                max_diff = std::max(max_diff, diff);
            }
        }
        double global_max_diff = mpi.allreduce_max(max_diff);

        if (rank == 0) {
            std::cout << "[PASS] RK2 CPU vs GPU max diff = " 
                      << std::scientific << global_max_diff << "\n";
        }
        assert(global_max_diff < 1e-10 && "RK2 GPU result differs from CPU");
    }

    // ------------------------------------------------------------------------
    // Test 3: RK4 Comparison (CPU vs GPU)
    // ------------------------------------------------------------------------
    {
        mcch::Solver solver_cpu(mpi, grid, n_comp, fe, mob, kappa);
        mcch::Solver solver_gpu(mpi, grid, n_comp, fe, mob, kappa);

        solver_cpu.initialize_random(c_mean, 0.02, 9012);
        solver_gpu.c = solver_cpu.c;
        solver_gpu.set_device(mcch::DeviceType::GPU);
        solver_gpu.gpu_solver->copy_to_device(solver_gpu.c);

        for (int step = 0; step < 10; ++step) {
            solver_cpu.step_rk4(dt);
            solver_gpu.step_rk4(dt);
        }

        solver_gpu.sync_from_device();

        double max_diff = 0.0;
        for (int m = 0; m < n_comp; ++m) {
            for (size_t i = 0; i < grid.local_total_size(); ++i) {
                double diff = std::abs(solver_cpu.c[m][i] - solver_gpu.c[m][i]);
                max_diff = std::max(max_diff, diff);
            }
        }
        double global_max_diff = mpi.allreduce_max(max_diff);

        if (rank == 0) {
            std::cout << "[PASS] RK4 CPU vs GPU max diff = " 
                      << std::scientific << global_max_diff << "\n";
        }
        assert(global_max_diff < 1e-10 && "RK4 GPU result differs from CPU");
    }

    // ------------------------------------------------------------------------
    // Test 4: Thermodynamic Consistency (dF/dt <= 0 on GPU)
    // ------------------------------------------------------------------------
    {
        mcch::Solver solver_gpu(mpi, grid, n_comp, fe, mob, kappa);
        solver_gpu.initialize_random(c_mean, 0.05, 42);
        solver_gpu.set_device(mcch::DeviceType::GPU);
        solver_gpu.gpu_solver->copy_to_device(solver_gpu.c);

        double E0 = solver_gpu.compute_total_free_energy();
        for (int step = 0; step < 50; ++step) {
            solver_gpu.step_rk4(dt);
        }
        double E_final = solver_gpu.compute_total_free_energy();

        if (rank == 0) {
            std::cout << "[PASS] GPU Free Energy Decrease: E0 = " << E0 
                      << ", E_final = " << E_final << " (dE = " << (E_final - E0) << ")\n";
        }
        assert(E_final < E0 && "Thermodynamic free energy did not decrease on GPU");
    }

    if (rank == 0) {
        std::cout << "\nALL GPU TESTS PASSED SUCCESSFULLY!\n";
    }

    MPI_Finalize();
    return 0;
#endif
}
