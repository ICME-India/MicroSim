#pragma once

#include <string>
#include <iostream>
#include <stdexcept>

#ifdef MCCH_ENABLE_CUDA
#include <cuda_runtime.h>

#define CUDA_CHECK(call)                                                      \
    do {                                                                      \
        cudaError_t err__ = (call);                                           \
        if (err__ != cudaSuccess) {                                           \
            std::cerr << "CUDA Error at " << __FILE__ << ":" << __LINE__      \
                      << " code=" << err__ << " ("                            \
                      << cudaGetErrorString(err__) << ")" << std::endl;       \
            throw std::runtime_error(std::string("CUDA Error: ") +           \
                                     cudaGetErrorString(err__));              \
        }                                                                     \
    } while (0)

#define CUDA_CHECK_LAST_ERROR()                                               \
    do {                                                                      \
        cudaError_t err__ = cudaGetLastError();                               \
        if (err__ != cudaSuccess) {                                           \
            std::cerr << "CUDA Kernel Launch Error at " << __FILE__ << ":"    \
                      << __LINE__ << " code=" << err__ << " ("                \
                      << cudaGetErrorString(err__) << ")" << std::endl;       \
            throw std::runtime_error(std::string("CUDA Kernel Error: ") +     \
                                     cudaGetErrorString(err__));              \
        }                                                                     \
    } while (0)

#endif

namespace mcch {

enum class DeviceType {
    CPU,
    GPU
};

enum CudaFreeEnergyType {
    CUDA_FE_POLYNOMIAL_MULTIWELL = 0,
    CUDA_FE_REGULAR_SOLUTION = 1,
    CUDA_FE_BINARY_DOUBLE_WELL = 2
};

enum CudaMobilityType {
    CUDA_MOB_CONSTANT = 0,
    CUDA_MOB_DEGENERATE = 1
};

enum CudaBoundaryType {
    CUDA_BC_PERIODIC = 0,
    CUDA_BC_NEUMANN = 1
};

// Simulation parameters POD structure passed by value to CUDA kernels
struct CudaSimParams {
    int dim;             // 2 or 3
    int nx, ny, nz;      // global grid dimensions
    int local_nx;        // local slab size along x
    int ghost;           // ghost layer depth (typically 1)
    int alloc_nx;        // local_nx + 2 * ghost
    size_t slice_size;   // ny * nz
    size_t alloc_total;  // alloc_nx * ny * nz
    size_t local_total;  // local_nx * ny * nz

    double dx, dy, dz;
    double idx, idy, idz;
    double idx2, idy2, idz2;

    int n_comp;
    int fe_type;         // CudaFreeEnergyType
    int mob_type;        // CudaMobilityType
    double M0;           // Mobility value
    double RT;           // For regular solution
    double eps;          // Regularization for ln(c)
    double W_binary;     // For binary double-well

    int bc_x;            // CudaBoundaryType
    int bc_y;
    int bc_z;
};

} // namespace mcch
