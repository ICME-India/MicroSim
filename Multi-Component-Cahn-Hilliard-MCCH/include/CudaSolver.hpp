#pragma once

#include "CudaCommon.hpp"
#include "Grid.hpp"
#include "MPIContext.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"

#include <vector>
#include <memory>
#include <iostream>

namespace mcch {

#ifdef MCCH_ENABLE_CUDA

class CudaSolver {
public:
    const MPIContext& mpi;
    const Grid& grid;
    int n_comp;

    // Simulation configuration POD for kernels
    CudaSimParams params;

    // Device memory pointers
    double* d_c = nullptr;         // Concentrations: n_comp * alloc_total
    double* d_mu = nullptr;        // Chemical potentials: n_comp * alloc_total
    double* d_rhs = nullptr;       // RHS time derivatives: n_comp * alloc_total
    double* d_c_stage = nullptr;   // Intermediate state for RK stages: n_comp * alloc_total
    double* d_k1 = nullptr;        // RK Stage 1 derivative: n_comp * alloc_total
    double* d_k2 = nullptr;        // RK Stage 2 derivative: n_comp * alloc_total
    double* d_k3 = nullptr;        // RK Stage 3 derivative: n_comp * alloc_total
    double* d_k4 = nullptr;        // RK Stage 4 derivative: n_comp * alloc_total

    // Device parameters
    double* d_kappa = nullptr;     // Gradient energy matrix (n_comp x n_comp)
    double* d_M = nullptr;         // Mobility matrix (n_comp x n_comp)
    double* d_W = nullptr;         // Multi-well barrier matrix (n_comp x n_comp)
    double* d_A = nullptr;         // Self-barrier vector (n_comp)
    double* d_omega = nullptr;     // Regular solution interaction matrix (n_comp x n_comp)

    // Host pinned buffers for MPI ghost exchange (zero-copy / high bandwidth)
    double* h_pinned_send_left = nullptr;
    double* h_pinned_recv_left = nullptr;
    double* h_pinned_send_right = nullptr;
    double* h_pinned_recv_right = nullptr;
    size_t ghost_bytes_per_comp = 0;

    // Execution configuration
    int threads_per_block = 256;
    int num_blocks_cells = 0;
    int num_blocks_boundary = 0;

    CudaSolver(const MPIContext& mpi_ctx, const Grid& grid_in, int num_components,
               std::shared_ptr<FreeEnergy> fe, std::shared_ptr<MobilityModel> mob,
               const std::vector<std::vector<double>>& kappa_mat);

    ~CudaSolver();

    // Disable copy
    CudaSolver(const CudaSolver&) = delete;
    CudaSolver& operator=(const CudaSolver&) = delete;

    // Memory transfers between host and device
    void copy_to_device(const std::vector<std::vector<double>>& host_c);
    void copy_to_host(std::vector<std::vector<double>>& host_c) const;
    void copy_mu_to_host(std::vector<std::vector<double>>& host_mu) const;

    // Boundary & ghost cell exchange
    void exchange_ghosts(double* d_field);

    // Compute chemical potentials on GPU
    void compute_chemical_potentials(const double* d_c_in, double* d_mu_out);

    // Compute RHS (time derivative) on GPU
    void compute_rhs(const double* d_c_in, double* d_mu_scratch, double* d_rhs_out);

    // Explicit Euler step on GPU
    void step_euler(double dt);

    // 2nd-order Runge-Kutta (Heun / Midpoint) step on GPU
    void step_rk2(double dt);

    // Classical 4th-order Runge-Kutta step on GPU
    void step_rk4(double dt);

    // Normalize partition of unity (sum c_i = 1) on GPU
    void normalize_unity();

    // Device memory allocation & deallocation
    void allocate_device_memory();
    void free_device_memory();
};

#else

// CPU fallback stub when CUDA is not compiled
class CudaSolver {
public:
    CudaSolver(const MPIContext&, const Grid&, int,
               std::shared_ptr<FreeEnergy>, std::shared_ptr<MobilityModel>,
               const std::vector<std::vector<double>>&) {
        throw std::runtime_error("MCCH was compiled without CUDA. GPU execution is unavailable.");
    }

    void copy_to_device(const std::vector<std::vector<double>>&) {}
    void copy_to_host(std::vector<std::vector<double>>&) const {}
    void copy_mu_to_host(std::vector<std::vector<double>>&) const {}
    void step_euler(double) {}
    void step_rk2(double) {}
    void step_rk4(double) {}
};

#endif

} // namespace mcch
