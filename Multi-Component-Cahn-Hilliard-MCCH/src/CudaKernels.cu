#ifdef MCCH_ENABLE_CUDA

#include "CudaSolver.hpp"
#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace mcch {

// ============================================================================
// CUDA Kernels
// ============================================================================

__global__ void kernel_compute_chemical_potentials(
    const double* __restrict__ c,
    double* __restrict__ mu,
    const double* __restrict__ kappa,
    const double* __restrict__ W,
    const double* __restrict__ A,
    const double* __restrict__ omega,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;

    // Neighbor indices along y
    int j_m, j_p;
    if (params.bc_y == CUDA_BC_PERIODIC) {
        j_m = (j == 0) ? ny - 1 : j - 1;
        j_p = (j == ny - 1) ? 0 : j + 1;
    } else {
        j_m = (j == 0) ? 1 : j - 1;
        j_p = (j == ny - 1) ? ny - 2 : j + 1;
    }

    // Neighbor indices along z
    int k_m = 0, k_p = 0;
    if (params.dim == 3) {
        if (params.bc_z == CUDA_BC_PERIODIC) {
            k_m = (k == 0) ? nz - 1 : k - 1;
            k_p = (k == nz - 1) ? 0 : k + 1;
        } else {
            k_m = (k == 0) ? 1 : k - 1;
            k_p = (k == nz - 1) ? nz - 2 : k + 1;
        }
    }

    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;
    size_t xm_idx   = ((size_t)(i_alloc - 1) * ny + j) * nz + k;
    size_t xp_idx   = ((size_t)(i_alloc + 1) * ny + j) * nz + k;
    size_t ym_idx   = ((size_t)i_alloc * ny + j_m) * nz + k;
    size_t yp_idx   = ((size_t)i_alloc * ny + j_p) * nz + k;

    size_t zm_idx = 0, zp_idx = 0;
    if (params.dim == 3) {
        zm_idx = ((size_t)i_alloc * ny + j) * nz + k_m;
        zp_idx = ((size_t)i_alloc * ny + j) * nz + k_p;
    }

    int n_comp = params.n_comp;
    double c_local[16];
    double lap_c[16];
    double df_dc[16];

    for (int m = 0; m < n_comp; ++m) {
        size_t offset = (size_t)m * params.alloc_total;
        double cm = c[offset + curr_idx];
        c_local[m] = cm;

        double lap = (c[offset + xp_idx] - 2.0 * cm + c[offset + xm_idx]) * params.idx2
                   + (c[offset + yp_idx] - 2.0 * cm + c[offset + ym_idx]) * params.idy2;
        if (params.dim == 3) {
            lap += (c[offset + zp_idx] - 2.0 * cm + c[offset + zm_idx]) * params.idz2;
        }
        lap_c[m] = lap;
    }

    // Thermodynamic bulk derivatives df0/dc_m
    if (params.fe_type == CUDA_FE_REGULAR_SOLUTION) {
        double eps = params.eps;
        double RT = params.RT;
        for (int m = 0; m < n_comp; ++m) {
            double ci = (c_local[m] > eps) ? c_local[m] : eps;
            double d_entropy = RT * (log(ci) + 1.0);
            double d_enthalpy = 0.0;
            for (int j_c = 0; j_c < n_comp; ++j_c) {
                if (m != j_c) {
                    d_enthalpy += omega[m * n_comp + j_c] * c_local[j_c];
                }
            }
            df_dc[m] = d_entropy + d_enthalpy;
        }
    } else if (params.fe_type == CUDA_FE_BINARY_DOUBLE_WELL) {
        double c0 = c_local[0];
        df_dc[0] = 2.0 * params.W_binary * c0 * (1.0 - c0) * (1.0 - 2.0 * c0);
        for (int m = 1; m < n_comp; ++m) {
            df_dc[m] = 0.0;
        }
    } else { // CUDA_FE_POLYNOMIAL_MULTIWELL
        for (int m = 0; m < n_comp; ++m) {
            double d = 0.0;
            for (int j_c = 0; j_c < n_comp; ++j_c) {
                if (m != j_c) {
                    d += 2.0 * W[m * n_comp + j_c] * c_local[m] * (c_local[j_c] * c_local[j_c]);
                }
            }
            double Am = A[m];
            if (Am != 0.0) {
                d += 2.0 * Am * c_local[m] * (1.0 - c_local[m]) * (1.0 - 2.0 * c_local[m]);
            }
            df_dc[m] = d;
        }
    }

    // Raw chemical potential: mu_m = df0/dc_m - sum_n (kappa_mn * lap_c_n)
    double mu_sum = 0.0;
    double raw_mu[16];
    for (int m = 0; m < n_comp; ++m) {
        double grad_term = 0.0;
        for (int n = 0; n < n_comp; ++n) {
            grad_term += kappa[m * n_comp + n] * lap_c[n];
        }
        double mu_val = df_dc[m] - grad_term;
        raw_mu[m] = mu_val;
        mu_sum += mu_val;
    }

    // Projected chemical potential: mu_tilde_m = mu_m - (1/N) * sum(mu_k)
    double mu_mean = mu_sum / (double)n_comp;
    for (int m = 0; m < n_comp; ++m) {
        size_t offset = (size_t)m * params.alloc_total;
        mu[offset + curr_idx] = raw_mu[m] - mu_mean;
    }
}

// In-device ghost exchange for single-GPU execution (P = 1)
__global__ void kernel_apply_x_boundary(
    double* __restrict__ field,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.slice_size) return;

    int ny = params.ny;
    int nz = params.nz;
    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    int j = (params.dim == 3) ? (int)(tid / nz) : (int)tid;

    int ghost = params.ghost;
    int local_nx = params.local_nx;

    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        if (params.bc_x == CUDA_BC_PERIODIC) {
            // Periodic:
            // Left ghost [0 .. ghost-1] receives from interior right [local_nx .. local_nx+ghost-1]
            // Right ghost [ghost+local_nx .. ghost+local_nx+ghost-1] receives from interior left [ghost .. 2*ghost-1]
            for (int g = 0; g < ghost; ++g) {
                size_t left_ghost_idx     = ((size_t)g * ny + j) * nz + k;
                size_t right_interior_idx = ((size_t)(local_nx + g) * ny + j) * nz + k;
                field[off + left_ghost_idx] = field[off + right_interior_idx];

                size_t right_ghost_idx    = ((size_t)(ghost + local_nx + g) * ny + j) * nz + k;
                size_t left_interior_idx  = ((size_t)(ghost + g) * ny + j) * nz + k;
                field[off + right_ghost_idx] = field[off + left_interior_idx];
            }
        } else {
            // Neumann reflection:
            for (int g = 0; g < ghost; ++g) {
                size_t left_ghost_idx     = ((size_t)(ghost - 1 - g) * ny + j) * nz + k;
                size_t left_interior_idx  = ((size_t)(ghost + g) * ny + j) * nz + k;
                field[off + left_ghost_idx] = field[off + left_interior_idx];

                size_t right_ghost_idx    = ((size_t)(ghost + local_nx + g) * ny + j) * nz + k;
                size_t right_interior_idx = ((size_t)(ghost + local_nx - 1 - g) * ny + j) * nz + k;
                field[off + right_ghost_idx] = field[off + right_interior_idx];
            }
        }
    }
}

// Extract interior slices to send buffers for MPI multi-rank exchange
__global__ void kernel_extract_ghost_slice(
    const double* __restrict__ field,
    double* __restrict__ send_left,
    double* __restrict__ send_right,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.slice_size) return;

    int ny = params.ny;
    int nz = params.nz;
    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    int j = (params.dim == 3) ? (int)(tid / nz) : (int)tid;

    int ghost = params.ghost;
    int local_nx = params.local_nx;

    for (int m = 0; m < params.n_comp; ++m) {
        size_t off_field = (size_t)m * params.alloc_total;
        size_t off_buf   = (size_t)m * (params.ghost * params.slice_size);

        for (int g = 0; g < ghost; ++g) {
            size_t buf_idx = ((size_t)g * ny + j) * nz + k;
            size_t left_idx  = ((size_t)(ghost + g) * ny + j) * nz + k;
            size_t right_idx = ((size_t)(local_nx + g) * ny + j) * nz + k;

            send_left[off_buf + buf_idx]  = field[off_field + left_idx];
            send_right[off_buf + buf_idx] = field[off_field + right_idx];
        }
    }
}

// Insert received ghost slices into device field
__global__ void kernel_insert_ghost_slice(
    double* __restrict__ field,
    const double* __restrict__ recv_left,
    const double* __restrict__ recv_right,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.slice_size) return;

    int ny = params.ny;
    int nz = params.nz;
    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    int j = (params.dim == 3) ? (int)(tid / nz) : (int)tid;

    int ghost = params.ghost;
    int local_nx = params.local_nx;

    for (int m = 0; m < params.n_comp; ++m) {
        size_t off_field = (size_t)m * params.alloc_total;
        size_t off_buf   = (size_t)m * (params.ghost * params.slice_size);

        for (int g = 0; g < ghost; ++g) {
            size_t buf_idx = ((size_t)g * ny + j) * nz + k;
            size_t left_ghost_idx  = ((size_t)g * ny + j) * nz + k;
            size_t right_ghost_idx = ((size_t)(ghost + local_nx + g) * ny + j) * nz + k;

            if (recv_left)  field[off_field + left_ghost_idx]  = recv_left[off_buf + buf_idx];
            if (recv_right) field[off_field + right_ghost_idx] = recv_right[off_buf + buf_idx];
        }
    }
}

// Compute flux divergence: fast path for constant isotropic mobility M_ij = M0 * delta_ij
__global__ void kernel_compute_flux_divergence_constant(
    const double* __restrict__ mu,
    double* __restrict__ rhs,
    const double* __restrict__ M,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;

    int j_m, j_p;
    if (params.bc_y == CUDA_BC_PERIODIC) {
        j_m = (j == 0) ? ny - 1 : j - 1;
        j_p = (j == ny - 1) ? 0 : j + 1;
    } else {
        j_m = (j == 0) ? 1 : j - 1;
        j_p = (j == ny - 1) ? ny - 2 : j + 1;
    }

    int k_m = 0, k_p = 0;
    if (params.dim == 3) {
        if (params.bc_z == CUDA_BC_PERIODIC) {
            k_m = (k == 0) ? nz - 1 : k - 1;
            k_p = (k == nz - 1) ? 0 : k + 1;
        } else {
            k_m = (k == 0) ? 1 : k - 1;
            k_p = (k == nz - 1) ? nz - 2 : k + 1;
        }
    }

    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;
    size_t xm_idx   = ((size_t)(i_alloc - 1) * ny + j) * nz + k;
    size_t xp_idx   = ((size_t)(i_alloc + 1) * ny + j) * nz + k;
    size_t ym_idx   = ((size_t)i_alloc * ny + j_m) * nz + k;
    size_t yp_idx   = ((size_t)i_alloc * ny + j_p) * nz + k;

    size_t zm_idx = 0, zp_idx = 0;
    if (params.dim == 3) {
        zm_idx = ((size_t)i_alloc * ny + j) * nz + k_m;
        zp_idx = ((size_t)i_alloc * ny + j) * nz + k_p;
    }

    double lap_mu[16];
    int n_comp = params.n_comp;
    for (int m = 0; m < n_comp; ++m) {
        size_t offset = (size_t)m * params.alloc_total;
        double curr_mu = mu[offset + curr_idx];
        double lap = (mu[offset + xp_idx] - 2.0 * curr_mu + mu[offset + xm_idx]) * params.idx2
                   + (mu[offset + yp_idx] - 2.0 * curr_mu + mu[offset + ym_idx]) * params.idy2;
        if (params.dim == 3) {
            lap += (mu[offset + zp_idx] - 2.0 * curr_mu + mu[offset + zm_idx]) * params.idz2;
        }
        lap_mu[m] = lap;
    }

    for (int m = 0; m < n_comp; ++m) {
        double sum = 0.0;
        for (int n = 0; n < n_comp; ++n) {
            sum += M[m * n_comp + n] * lap_mu[n];
        }
        size_t offset = (size_t)m * params.alloc_total;
        rhs[offset + curr_idx] = sum;
    }
}

// Compute flux divergence: conservative face-flux formulation for degenerate Onsager mobility
__global__ void kernel_compute_flux_divergence_degenerate(
    const double* __restrict__ c,
    const double* __restrict__ mu,
    double* __restrict__ rhs,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;
    int n_comp = params.n_comp;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;

    int j_m, j_p;
    if (params.bc_y == CUDA_BC_PERIODIC) {
        j_m = (j == 0) ? ny - 1 : j - 1;
        j_p = (j == ny - 1) ? 0 : j + 1;
    } else {
        j_m = (j == 0) ? 1 : j - 1;
        j_p = (j == ny - 1) ? ny - 2 : j + 1;
    }

    int k_m = 0, k_p = 0;
    if (params.dim == 3) {
        if (params.bc_z == CUDA_BC_PERIODIC) {
            k_m = (k == 0) ? nz - 1 : k - 1;
            k_p = (k == nz - 1) ? 0 : k + 1;
        } else {
            k_m = (k == 0) ? 1 : k - 1;
            k_p = (k == nz - 1) ? nz - 2 : k + 1;
        }
    }

    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;
    size_t xm_idx   = ((size_t)(i_alloc - 1) * ny + j) * nz + k;
    size_t xp_idx   = ((size_t)(i_alloc + 1) * ny + j) * nz + k;
    size_t ym_idx   = ((size_t)i_alloc * ny + j_m) * nz + k;
    size_t yp_idx   = ((size_t)i_alloc * ny + j_p) * nz + k;

    size_t zm_idx = 0, zp_idx = 0;
    if (params.dim == 3) {
        zm_idx = ((size_t)i_alloc * ny + j) * nz + k_m;
        zp_idx = ((size_t)i_alloc * ny + j) * nz + k_p;
    }

    double J_xp[16], J_xm[16];
    double J_yp[16], J_ym[16];
    double J_zp[16], J_zm[16];

    double M0 = params.M0;

    // 1. Face x+1/2
    {
        double cf[16], grad_mu[16];
        for (int m = 0; m < n_comp; ++m) {
            size_t off = (size_t)m * params.alloc_total;
            double c_val = 0.5 * (c[off + curr_idx] + c[off + xp_idx]);
            cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
            grad_mu[m] = (mu[off + xp_idx] - mu[off + curr_idx]) * params.idx;
        }
        for (int m = 0; m < n_comp; ++m) {
            double flux = 0.0;
            for (int n = 0; n < n_comp; ++n) {
                double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                flux -= Mij * grad_mu[n];
            }
            J_xp[m] = flux;
        }
    }

    // 2. Face x-1/2
    {
        double cf[16], grad_mu[16];
        for (int m = 0; m < n_comp; ++m) {
            size_t off = (size_t)m * params.alloc_total;
            double c_val = 0.5 * (c[off + xm_idx] + c[off + curr_idx]);
            cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
            grad_mu[m] = (mu[off + curr_idx] - mu[off + xm_idx]) * params.idx;
        }
        for (int m = 0; m < n_comp; ++m) {
            double flux = 0.0;
            for (int n = 0; n < n_comp; ++n) {
                double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                flux -= Mij * grad_mu[n];
            }
            J_xm[m] = flux;
        }
    }

    // 3. Face y+1/2
    {
        double cf[16], grad_mu[16];
        for (int m = 0; m < n_comp; ++m) {
            size_t off = (size_t)m * params.alloc_total;
            double c_val = 0.5 * (c[off + curr_idx] + c[off + yp_idx]);
            cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
            grad_mu[m] = (mu[off + yp_idx] - mu[off + curr_idx]) * params.idy;
        }
        for (int m = 0; m < n_comp; ++m) {
            double flux = 0.0;
            for (int n = 0; n < n_comp; ++n) {
                double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                flux -= Mij * grad_mu[n];
            }
            J_yp[m] = flux;
        }
    }

    // 4. Face y-1/2
    {
        double cf[16], grad_mu[16];
        for (int m = 0; m < n_comp; ++m) {
            size_t off = (size_t)m * params.alloc_total;
            double c_val = 0.5 * (c[off + ym_idx] + c[off + curr_idx]);
            cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
            grad_mu[m] = (mu[off + curr_idx] - mu[off + ym_idx]) * params.idy;
        }
        for (int m = 0; m < n_comp; ++m) {
            double flux = 0.0;
            for (int n = 0; n < n_comp; ++n) {
                double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                flux -= Mij * grad_mu[n];
            }
            J_ym[m] = flux;
        }
    }

    // 5. Z faces if 3D
    if (params.dim == 3) {
        // Face z+1/2
        {
            double cf[16], grad_mu[16];
            for (int m = 0; m < n_comp; ++m) {
                size_t off = (size_t)m * params.alloc_total;
                double c_val = 0.5 * (c[off + curr_idx] + c[off + zp_idx]);
                cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
                grad_mu[m] = (mu[off + zp_idx] - mu[off + curr_idx]) * params.idz;
            }
            for (int m = 0; m < n_comp; ++m) {
                double flux = 0.0;
                for (int n = 0; n < n_comp; ++n) {
                    double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                    flux -= Mij * grad_mu[n];
                }
                J_zp[m] = flux;
            }
        }

        // Face z-1/2
        {
            double cf[16], grad_mu[16];
            for (int m = 0; m < n_comp; ++m) {
                size_t off = (size_t)m * params.alloc_total;
                double c_val = 0.5 * (c[off + zm_idx] + c[off + curr_idx]);
                cf[m] = (c_val < 0.0) ? 0.0 : ((c_val > 1.0) ? 1.0 : c_val);
                grad_mu[m] = (mu[off + curr_idx] - mu[off + zm_idx]) * params.idz;
            }
            for (int m = 0; m < n_comp; ++m) {
                double flux = 0.0;
                for (int n = 0; n < n_comp; ++n) {
                    double Mij = (m == n) ? (M0 * cf[m] * (1.0 - cf[m])) : (-M0 * cf[m] * cf[n]);
                    flux -= Mij * grad_mu[n];
                }
                J_zm[m] = flux;
            }
        }
    }

    // Divergence: rhs = -div(J)
    for (int m = 0; m < n_comp; ++m) {
        double div_J = (J_xp[m] - J_xm[m]) * params.idx + (J_yp[m] - J_ym[m]) * params.idy;
        if (params.dim == 3) {
            div_J += (J_zp[m] - J_zm[m]) * params.idz;
        }
        size_t off = (size_t)m * params.alloc_total;
        rhs[off + curr_idx] = -div_J;
    }
}

// Explicit Euler point-wise update: c += dt * rhs
__global__ void kernel_euler_update(
    double* __restrict__ c,
    const double* __restrict__ rhs,
    double dt,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;
    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;

    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        c[off + curr_idx] += dt * rhs[off + curr_idx];
    }
}

// Stage update: c_stage = c_base + factor * k_stage
__global__ void kernel_rk_stage_update(
    double* __restrict__ c_stage,
    const double* __restrict__ c_base,
    const double* __restrict__ k_stage,
    double factor,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;
    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;

    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        c_stage[off + curr_idx] = c_base[off + curr_idx] + factor * k_stage[off + curr_idx];
    }
}

// RK2 final combination: c += 0.5 * dt * (k1 + k2)
__global__ void kernel_rk2_combine(
    double* __restrict__ c,
    const double* __restrict__ k1,
    const double* __restrict__ k2,
    double dt,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;
    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;

    double factor = 0.5 * dt;
    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        c[off + curr_idx] += factor * (k1[off + curr_idx] + k2[off + curr_idx]);
    }
}

// RK4 final combination: c += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
__global__ void kernel_rk4_combine(
    double* __restrict__ c,
    const double* __restrict__ k1,
    const double* __restrict__ k2,
    const double* __restrict__ k3,
    const double* __restrict__ k4,
    double dt,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;
    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;

    double dt6 = dt / 6.0;
    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        c[off + curr_idx] += dt6 * (k1[off + curr_idx] + 2.0 * k2[off + curr_idx] + 2.0 * k3[off + curr_idx] + k4[off + curr_idx]);
    }
}

// Normalize partition of unity locally: sum c_i = 1
__global__ void kernel_normalize_unity(
    double* __restrict__ c,
    CudaSimParams params)
{
    size_t tid = blockDim.x * (size_t)blockIdx.x + threadIdx.x;
    if (tid >= params.local_total) return;

    int nz = params.nz;
    int ny = params.ny;
    int ghost = params.ghost;

    int k = (params.dim == 3) ? (int)(tid % nz) : 0;
    size_t temp = (params.dim == 3) ? (tid / nz) : tid;
    int j = (int)(temp % ny);
    int i_loc = (int)(temp / ny);
    int i_alloc = i_loc + ghost;
    size_t curr_idx = ((size_t)i_alloc * ny + j) * nz + k;

    double sum = 0.0;
    for (int m = 0; m < params.n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        sum += c[off + curr_idx];
    }
    if (fabs(sum - 1.0) > 1e-15 && sum > 0.0) {
        double inv_sum = 1.0 / sum;
        for (int m = 0; m < params.n_comp; ++m) {
            size_t off = (size_t)m * params.alloc_total;
            c[off + curr_idx] *= inv_sum;
        }
    }
}


// ============================================================================
// CudaSolver Class Implementation
// ============================================================================

CudaSolver::CudaSolver(const MPIContext& mpi_ctx, const Grid& grid_in, int num_components,
                       std::shared_ptr<FreeEnergy> fe, std::shared_ptr<MobilityModel> mob,
                       const std::vector<std::vector<double>>& kappa_mat)
    : mpi(mpi_ctx), grid(grid_in), n_comp(num_components)
{
    // Populate parameters
    params.dim = grid.dim;
    params.nx = grid.nx;
    params.ny = grid.ny;
    params.nz = grid.nz;
    params.local_nx = grid.local_nx;
    params.ghost = grid.ghost;
    params.alloc_nx = grid.alloc_nx;
    params.slice_size = grid.slice_size();
    params.alloc_total = grid.local_total_size();
    params.local_total = grid.local_interior_size();

    params.dx = grid.dx;
    params.dy = grid.dy;
    params.dz = grid.dz;
    params.idx = 1.0 / grid.dx;
    params.idy = 1.0 / grid.dy;
    params.idz = (grid.dim == 3) ? (1.0 / grid.dz) : 0.0;
    params.idx2 = params.idx * params.idx;
    params.idy2 = params.idy * params.idy;
    params.idz2 = params.idz * params.idz;

    params.n_comp = n_comp;

    // Boundaries
    params.bc_x = (grid.bc_x == BoundaryType::NEUMANN) ? CUDA_BC_NEUMANN : CUDA_BC_PERIODIC;
    params.bc_y = (grid.bc_y == BoundaryType::NEUMANN) ? CUDA_BC_NEUMANN : CUDA_BC_PERIODIC;
    params.bc_z = (grid.bc_z == BoundaryType::NEUMANN) ? CUDA_BC_NEUMANN : CUDA_BC_PERIODIC;

    // Free Energy Type
    if (auto reg = dynamic_cast<RegularSolution*>(fe.get())) {
        params.fe_type = CUDA_FE_REGULAR_SOLUTION;
        params.RT = reg->get_RT();
        params.eps = reg->get_eps();
    } else if (auto bin = dynamic_cast<BinaryDoubleWell*>(fe.get())) {
        params.fe_type = CUDA_FE_BINARY_DOUBLE_WELL;
        params.W_binary = bin->get_barrier();
    } else {
        params.fe_type = CUDA_FE_POLYNOMIAL_MULTIWELL;
    }

    // Mobility Type
    if (mob->is_constant()) {
        params.mob_type = CUDA_MOB_CONSTANT;
        params.M0 = static_cast<ConstantMobility*>(mob.get())->value();
    } else {
        params.mob_type = CUDA_MOB_DEGENERATE;
        params.M0 = static_cast<DegenerateMobility*>(mob.get())->base_value();
    }

    // Calculate grid block configurations
    threads_per_block = 256;
    num_blocks_cells = (int)((params.local_total + threads_per_block - 1) / threads_per_block);
    num_blocks_boundary = (int)((params.slice_size + threads_per_block - 1) / threads_per_block);

    // Allocate device buffers
    allocate_device_memory();

    // Upload kappa matrix
    std::vector<double> flat_kappa(n_comp * n_comp, 0.0);
    for (int i = 0; i < n_comp; ++i) {
        for (int j = 0; j < n_comp; ++j) {
            flat_kappa[i * n_comp + j] = kappa_mat[i][j];
        }
    }
    CUDA_CHECK(cudaMemcpy(d_kappa, flat_kappa.data(), n_comp * n_comp * sizeof(double), cudaMemcpyHostToDevice));

    // Upload mobility matrix
    std::vector<double> flat_M(n_comp * n_comp, 0.0);
    std::vector<double> dummy_c(n_comp, 1.0 / n_comp);
    mob->evaluate_matrix(dummy_c, flat_M);
    CUDA_CHECK(cudaMemcpy(d_M, flat_M.data(), n_comp * n_comp * sizeof(double), cudaMemcpyHostToDevice));

    // Upload thermodynamic parameters
    if (auto poly = dynamic_cast<PolynomialMultiWell*>(fe.get())) {
        const auto& W_mat = poly->get_W();
        const auto& A_vec = poly->get_A();
        std::vector<double> flat_W(n_comp * n_comp, 0.0);
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < n_comp; ++j) {
                if (i < (int)W_mat.size() && j < (int)W_mat[i].size()) {
                    flat_W[i * n_comp + j] = W_mat[i][j];
                }
            }
        }
        std::vector<double> vec_A(n_comp, 0.0);
        for (int i = 0; i < n_comp; ++i) {
            if (i < (int)A_vec.size()) {
                vec_A[i] = A_vec[i];
            }
        }
        CUDA_CHECK(cudaMemcpy(d_W, flat_W.data(), n_comp * n_comp * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_A, vec_A.data(), n_comp * sizeof(double), cudaMemcpyHostToDevice));
    } else if (auto reg = dynamic_cast<RegularSolution*>(fe.get())) {
        const auto& omega_mat = reg->get_omega();
        std::vector<double> flat_omega(n_comp * n_comp, 0.0);
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < n_comp; ++j) {
                if (i < (int)omega_mat.size() && j < (int)omega_mat[i].size()) {
                    flat_omega[i * n_comp + j] = omega_mat[i][j];
                }
            }
        }
        CUDA_CHECK(cudaMemcpy(d_omega, flat_omega.data(), n_comp * n_comp * sizeof(double), cudaMemcpyHostToDevice));
    }
}

CudaSolver::~CudaSolver() {
    free_device_memory();
}

void CudaSolver::allocate_device_memory() {
    size_t total_field_bytes = (size_t)n_comp * params.alloc_total * sizeof(double);
    size_t mat_bytes = (size_t)n_comp * n_comp * sizeof(double);
    size_t vec_bytes = (size_t)n_comp * sizeof(double);

    CUDA_CHECK(cudaMalloc(&d_c, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_mu, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_rhs, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_c_stage, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_k1, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_k2, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_k3, total_field_bytes));
    CUDA_CHECK(cudaMalloc(&d_k4, total_field_bytes));

    CUDA_CHECK(cudaMalloc(&d_kappa, mat_bytes));
    CUDA_CHECK(cudaMalloc(&d_M, mat_bytes));
    CUDA_CHECK(cudaMalloc(&d_W, mat_bytes));
    CUDA_CHECK(cudaMalloc(&d_A, vec_bytes));
    CUDA_CHECK(cudaMalloc(&d_omega, mat_bytes));

    // Allocate pinned host buffers for MPI ghost exchange
    ghost_bytes_per_comp = (size_t)params.ghost * params.slice_size * sizeof(double);
    size_t total_ghost_bytes = (size_t)n_comp * ghost_bytes_per_comp;

    CUDA_CHECK(cudaMallocHost(&h_pinned_send_left, total_ghost_bytes));
    CUDA_CHECK(cudaMallocHost(&h_pinned_recv_left, total_ghost_bytes));
    CUDA_CHECK(cudaMallocHost(&h_pinned_send_right, total_ghost_bytes));
    CUDA_CHECK(cudaMallocHost(&h_pinned_recv_right, total_ghost_bytes));
}

void CudaSolver::free_device_memory() {
    if (d_c) { cudaFree(d_c); d_c = nullptr; }
    if (d_mu) { cudaFree(d_mu); d_mu = nullptr; }
    if (d_rhs) { cudaFree(d_rhs); d_rhs = nullptr; }
    if (d_c_stage) { cudaFree(d_c_stage); d_c_stage = nullptr; }
    if (d_k1) { cudaFree(d_k1); d_k1 = nullptr; }
    if (d_k2) { cudaFree(d_k2); d_k2 = nullptr; }
    if (d_k3) { cudaFree(d_k3); d_k3 = nullptr; }
    if (d_k4) { cudaFree(d_k4); d_k4 = nullptr; }

    if (d_kappa) { cudaFree(d_kappa); d_kappa = nullptr; }
    if (d_M) { cudaFree(d_M); d_M = nullptr; }
    if (d_W) { cudaFree(d_W); d_W = nullptr; }
    if (d_A) { cudaFree(d_A); d_A = nullptr; }
    if (d_omega) { cudaFree(d_omega); d_omega = nullptr; }
    if (h_pinned_send_left) { cudaFreeHost(h_pinned_send_left); h_pinned_send_left = nullptr; }
    if (h_pinned_recv_left) { cudaFreeHost(h_pinned_recv_left); h_pinned_recv_left = nullptr; }
    if (h_pinned_send_right) { cudaFreeHost(h_pinned_send_right); h_pinned_send_right = nullptr; }
    if (h_pinned_recv_right) { cudaFreeHost(h_pinned_recv_right); h_pinned_recv_right = nullptr; }
}

void CudaSolver::copy_to_device(const std::vector<std::vector<double>>& host_c) {
    for (int m = 0; m < n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        CUDA_CHECK(cudaMemcpy(d_c + off, host_c[m].data(), params.alloc_total * sizeof(double), cudaMemcpyHostToDevice));
    }
}

void CudaSolver::copy_to_host(std::vector<std::vector<double>>& host_c) const {
    for (int m = 0; m < n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        CUDA_CHECK(cudaMemcpy(host_c[m].data(), d_c + off, params.alloc_total * sizeof(double), cudaMemcpyDeviceToHost));
    }
}

void CudaSolver::copy_mu_to_host(std::vector<std::vector<double>>& host_mu) const {
    for (int m = 0; m < n_comp; ++m) {
        size_t off = (size_t)m * params.alloc_total;
        CUDA_CHECK(cudaMemcpy(host_mu[m].data(), d_mu + off, params.alloc_total * sizeof(double), cudaMemcpyDeviceToHost));
    }
}

void CudaSolver::exchange_ghosts(double* d_field) {
    if (mpi.size == 1) {
        // Fast in-device GPU ghost exchange
        kernel_apply_x_boundary<<<num_blocks_boundary, threads_per_block>>>(d_field, params);
        CUDA_CHECK_LAST_ERROR();
        return;
    }

    // Multi-GPU with MPI: extract boundary slices to pinned host memory
    double* d_send_left = nullptr;
    double* d_send_right = nullptr;
    double* d_recv_left = nullptr;
    double* d_recv_right = nullptr;
    size_t total_ghost_bytes = (size_t)n_comp * ghost_bytes_per_comp;

    CUDA_CHECK(cudaMalloc(&d_send_left, total_ghost_bytes));
    CUDA_CHECK(cudaMalloc(&d_send_right, total_ghost_bytes));
    CUDA_CHECK(cudaMalloc(&d_recv_left, total_ghost_bytes));
    CUDA_CHECK(cudaMalloc(&d_recv_right, total_ghost_bytes));

    kernel_extract_ghost_slice<<<num_blocks_boundary, threads_per_block>>>(d_field, d_send_left, d_send_right, params);
    CUDA_CHECK_LAST_ERROR();

    CUDA_CHECK(cudaMemcpy(h_pinned_send_left, d_send_left, total_ghost_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_pinned_send_right, d_send_right, total_ghost_bytes, cudaMemcpyDeviceToHost));

    // Post non-blocking MPI exchanges across components
    int num_requests = 0;
    MPI_Request requests[4];
    int count = (int)(params.ghost * params.slice_size);

    for (int m = 0; m < n_comp; ++m) {
        size_t off = (size_t)m * (params.ghost * params.slice_size);
        double* s_left  = h_pinned_send_left + off;
        double* r_left  = h_pinned_recv_left + off;
        double* s_right = h_pinned_send_right + off;
        double* r_right = h_pinned_recv_right + off;

        num_requests = 0;
        if (mpi.left_rank != MPI_PROC_NULL) {
            MPI_Irecv(r_left, count, MPI_DOUBLE, mpi.left_rank, 100 + m, mpi.comm, &requests[num_requests++]);
        }
        if (mpi.right_rank != MPI_PROC_NULL) {
            MPI_Irecv(r_right, count, MPI_DOUBLE, mpi.right_rank, 200 + m, mpi.comm, &requests[num_requests++]);
        }
        if (mpi.right_rank != MPI_PROC_NULL) {
            MPI_Isend(s_right, count, MPI_DOUBLE, mpi.right_rank, 100 + m, mpi.comm, &requests[num_requests++]);
        }
        if (mpi.left_rank != MPI_PROC_NULL) {
            MPI_Isend(s_left, count, MPI_DOUBLE, mpi.left_rank, 200 + m, mpi.comm, &requests[num_requests++]);
        }

        if (num_requests > 0) {
            MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
        }

        // Neumann physical boundary condition at domain ends
        if (mpi.bc_x == BoundaryType::NEUMANN) {
            if (mpi.rank == 0) {
                // Mirror left
                std::copy(s_left, s_left + count, r_left);
            }
            if (mpi.rank == mpi.size - 1) {
                // Mirror right
                std::copy(s_right, s_right + count, r_right);
            }
        }
    }

    CUDA_CHECK(cudaMemcpy(d_recv_left, h_pinned_recv_left, total_ghost_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_recv_right, h_pinned_recv_right, total_ghost_bytes, cudaMemcpyHostToDevice));

    kernel_insert_ghost_slice<<<num_blocks_boundary, threads_per_block>>>(d_field, d_recv_left, d_recv_right, params);
    CUDA_CHECK_LAST_ERROR();

    cudaFree(d_send_left);
    cudaFree(d_send_right);
    cudaFree(d_recv_left);
    cudaFree(d_recv_right);
}

void CudaSolver::compute_chemical_potentials(const double* d_c_in, double* d_mu_out) {
    kernel_compute_chemical_potentials<<<num_blocks_cells, threads_per_block>>>(
        d_c_in, d_mu_out, d_kappa, d_W, d_A, d_omega, params);
    CUDA_CHECK_LAST_ERROR();

    exchange_ghosts(d_mu_out);
}

void CudaSolver::compute_rhs(const double* d_c_in, double* d_mu_scratch, double* d_rhs_out) {
    compute_chemical_potentials(d_c_in, d_mu_scratch);

    if (params.mob_type == CUDA_MOB_CONSTANT) {
        kernel_compute_flux_divergence_constant<<<num_blocks_cells, threads_per_block>>>(
            d_mu_scratch, d_rhs_out, d_M, params);
    } else {
        kernel_compute_flux_divergence_degenerate<<<num_blocks_cells, threads_per_block>>>(
            d_c_in, d_mu_scratch, d_rhs_out, params);
    }
    CUDA_CHECK_LAST_ERROR();
}

void CudaSolver::normalize_unity() {
    kernel_normalize_unity<<<num_blocks_cells, threads_per_block>>>(d_c, params);
    CUDA_CHECK_LAST_ERROR();
}

// ----------------------------------------------------------------------------
// 1. Explicit Euler implementation on GPU
// ----------------------------------------------------------------------------
void CudaSolver::step_euler(double dt) {
    compute_rhs(d_c, d_mu, d_rhs);

    kernel_euler_update<<<num_blocks_cells, threads_per_block>>>(d_c, d_rhs, dt, params);
    CUDA_CHECK_LAST_ERROR();

    normalize_unity();
    exchange_ghosts(d_c);
}

// ----------------------------------------------------------------------------
// 2. RK2 (Heun / Midpoint) implementation on GPU
// ----------------------------------------------------------------------------
void CudaSolver::step_rk2(double dt) {
    // Stage 1: k1 = rhs(c^n)
    compute_rhs(d_c, d_mu, d_k1);

    // Intermediate state: c_stage = c^n + dt * k1
    kernel_rk_stage_update<<<num_blocks_cells, threads_per_block>>>(d_c_stage, d_c, d_k1, dt, params);
    CUDA_CHECK_LAST_ERROR();
    exchange_ghosts(d_c_stage);

    // Stage 2: k2 = rhs(c_stage)
    compute_rhs(d_c_stage, d_mu, d_k2);

    // Final combination: c^{n+1} = c^n + 0.5 * dt * (k1 + k2)
    kernel_rk2_combine<<<num_blocks_cells, threads_per_block>>>(d_c, d_k1, d_k2, dt, params);
    CUDA_CHECK_LAST_ERROR();

    normalize_unity();
    exchange_ghosts(d_c);
}

// ----------------------------------------------------------------------------
// 3. Classical 4th-order Runge-Kutta (RK4) implementation on GPU
// ----------------------------------------------------------------------------
void CudaSolver::step_rk4(double dt) {
    // Stage 1: k1 = rhs(c^n)
    compute_rhs(d_c, d_mu, d_k1);

    // Stage 2: c_stage = c^n + 0.5 * dt * k1
    kernel_rk_stage_update<<<num_blocks_cells, threads_per_block>>>(d_c_stage, d_c, d_k1, 0.5 * dt, params);
    CUDA_CHECK_LAST_ERROR();
    exchange_ghosts(d_c_stage);
    compute_rhs(d_c_stage, d_mu, d_k2);

    // Stage 3: c_stage = c^n + 0.5 * dt * k2
    kernel_rk_stage_update<<<num_blocks_cells, threads_per_block>>>(d_c_stage, d_c, d_k2, 0.5 * dt, params);
    CUDA_CHECK_LAST_ERROR();
    exchange_ghosts(d_c_stage);
    compute_rhs(d_c_stage, d_mu, d_k3);

    // Stage 4: c_stage = c^n + dt * k3
    kernel_rk_stage_update<<<num_blocks_cells, threads_per_block>>>(d_c_stage, d_c, d_k3, dt, params);
    CUDA_CHECK_LAST_ERROR();
    exchange_ghosts(d_c_stage);
    compute_rhs(d_c_stage, d_mu, d_k4);

    // Final combination: c += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    kernel_rk4_combine<<<num_blocks_cells, threads_per_block>>>(d_c, d_k1, d_k2, d_k3, d_k4, dt, params);
    CUDA_CHECK_LAST_ERROR();

    normalize_unity();
    exchange_ghosts(d_c);
}

} // namespace mcch

#endif // MCCH_ENABLE_CUDA
