#pragma once

#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"

#include <fftw3.h>
#include <fftw3-mpi.h>
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace mcch {

/**
 * @brief High-performance Parallel Semi-Implicit Fourier Spectral Solver for Multi-Component
 *        Cahn-Hilliard (MCCH) equations using FFTW with MPI 1D Slab Decomposition.
 *
 * Mathematical formulation:
 * Mass conservation: dc_i/dt = M_0 * lap(tilde_mu_i)
 * where tilde_mu_i = (df_0/dc_i - 1/N * sum_k df_0/dc_k) - sum_j tilde_kappa_ij * lap(c_j).
 *
 * In Fourier space (k^2 = |k|^2):
 * F{dc_i/dt} = -M_0 * k^2 * F{tilde_mu_{0, i}} - M_0 * k^4 * sum_j kappa_ij * c_hat_j.
 *
 * The stiff linear 4th-order biharmonic term is integrated IMPLICITLY in Fourier space,
 * unconditionally eliminating the severe explicit stability limit (dt ~ dx^4 / (M*kappa)).
 *
 * Implements:
 * 1. First-order Semi-Implicit Euler (Chen & Shen 1998)
 * 2. Second-order Semi-Implicit Backward Differentiation Formula 2 (SBDF2)
 * 3. Optional Eyre-type convex-concave stabilization parameter S >= 0
 * 4. Arbitrary number of components (N >= 2) in both 2D and 3D
 * 5. General symmetric or diagonal gradient energy tensor kappa_ij
 */
class FFTWSemiImplicitSolver {
public:
    MPIContext& mpi;
    const Grid& grid;
    int n_comp;
    std::shared_ptr<FreeEnergy> free_energy;
    std::shared_ptr<MobilityModel> mobility;
    std::vector<std::vector<double>> kappa;

    ptrdiff_t local_nx;
    ptrdiff_t local_x_start;
    ptrdiff_t alloc_local;
    size_t local_cells; // local_nx * ny * nz

    // Precomputed Fourier wavevectors
    std::vector<double> k2; // |k|^2
    std::vector<double> k4; // |k|^4

    // FFTW complex buffers for each component
    std::vector<fftw_complex*> c_hat;
    std::vector<fftw_complex*> mu_hat;

    // History buffers for 2nd-order multi-step scheme (SBDF2)
    bool has_prev_step = false;
    std::vector<fftw_complex*> c_prev;
    std::vector<fftw_complex*> mu_prev;

    // Plans for MPI DFT
    fftw_plan plan_fwd = nullptr;
    fftw_plan plan_bwd = nullptr;

    double fft_norm = 1.0;
    bool is_kappa_diagonal = true;
    bool is_M_diagonal = true;
    std::vector<double> M_matrix;

    FFTWSemiImplicitSolver(MPIContext& mpi_ctx, const Grid& grid_in, int num_components,
                           std::shared_ptr<FreeEnergy> fe, std::shared_ptr<MobilityModel> mob,
                           const std::vector<std::vector<double>>& kappa_matrix)
        : mpi(mpi_ctx), grid(grid_in), n_comp(num_components),
          free_energy(fe), mobility(mob), kappa(kappa_matrix) {

        if ((int)kappa.size() != n_comp) {
            throw std::invalid_argument("Kappa matrix row count != num_components");
        }

        // Initialize FFTW MPI
        fftw_mpi_init();

        // Query FFTW MPI slab decomposition
        if (grid.dim == 2) {
            alloc_local = fftw_mpi_local_size_2d(grid.nx, grid.ny, mpi.comm, &local_nx, &local_x_start);
            local_cells = (size_t)local_nx * grid.ny;
        } else {
            alloc_local = fftw_mpi_local_size_3d(grid.nx, grid.ny, grid.nz, mpi.comm, &local_nx, &local_x_start);
            local_cells = (size_t)local_nx * grid.ny * grid.nz;
        }

        if (local_nx != mpi.local_nx || local_x_start != mpi.local_x_start) {
            if (mpi.rank == 0) {
                std::cerr << "Warning: MPIContext slab decomposition (" << mpi.local_nx 
                          << " @ " << mpi.local_x_start << ") differs from FFTW MPI ("
                          << local_nx << " @ " << local_x_start << ").\n"
                          << "Ensure MPIContext was constructed with use_fftw_decomp=true.\n";
            }
        }

        // Normalization factor for backward DFT
        double total_points = (double)grid.nx * (double)grid.ny * (grid.dim == 3 ? (double)grid.nz : 1.0);
        fft_norm = 1.0 / total_points;

        // Check if kappa is diagonal
        is_kappa_diagonal = true;
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < n_comp; ++j) {
                if (i != j && std::abs(kappa[i][j]) > 1e-14) {
                    is_kappa_diagonal = false;
                    break;
                }
            }
        }

        // Query mobility matrix and check if diagonal
        M_matrix.resize(n_comp * n_comp, 0.0);
        std::vector<double> dummy_c(n_comp, 1.0 / n_comp);
        mobility->evaluate_matrix(dummy_c, M_matrix);

        is_M_diagonal = true;
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < n_comp; ++j) {
                if (i != j && std::abs(M_matrix[i * n_comp + j]) > 1e-14) {
                    is_M_diagonal = false;
                    break;
                }
            }
        }

        // Allocate FFTW complex buffers
        c_hat.resize(n_comp);
        mu_hat.resize(n_comp);
        c_prev.resize(n_comp);
        mu_prev.resize(n_comp);

        for (int m = 0; m < n_comp; ++m) {
            c_hat[m] = fftw_alloc_complex(alloc_local);
            mu_hat[m] = fftw_alloc_complex(alloc_local);
            c_prev[m] = fftw_alloc_complex(alloc_local);
            mu_prev[m] = fftw_alloc_complex(alloc_local);
        }

        // Create MPI plans
        if (grid.dim == 2) {
            plan_fwd = fftw_mpi_plan_dft_2d(grid.nx, grid.ny, c_hat[0], c_hat[0], mpi.comm, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bwd = fftw_mpi_plan_dft_2d(grid.nx, grid.ny, c_hat[0], c_hat[0], mpi.comm, FFTW_BACKWARD, FFTW_ESTIMATE);
        } else {
            plan_fwd = fftw_mpi_plan_dft_3d(grid.nx, grid.ny, grid.nz, c_hat[0], c_hat[0], mpi.comm, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bwd = fftw_mpi_plan_dft_3d(grid.nx, grid.ny, grid.nz, c_hat[0], c_hat[0], mpi.comm, FFTW_BACKWARD, FFTW_ESTIMATE);
        }

        if (!plan_fwd || !plan_bwd) {
            throw std::runtime_error("Failed to create FFTW MPI plans.");
        }

        // Precompute Fourier wavenumbers |k|^2 and |k|^4
        precompute_wavenumbers();
    }

    ~FFTWSemiImplicitSolver() {
        if (plan_fwd) fftw_destroy_plan(plan_fwd);
        if (plan_bwd) fftw_destroy_plan(plan_bwd);
        for (int m = 0; m < n_comp; ++m) {
            if (c_hat[m]) fftw_free(c_hat[m]);
            if (mu_hat[m]) fftw_free(mu_hat[m]);
            if (c_prev[m]) fftw_free(c_prev[m]);
            if (mu_prev[m]) fftw_free(mu_prev[m]);
        }
    }

    // Reset multi-step history (e.g. after parameter change or adaptive step restart)
    void reset_history() {
        has_prev_step = false;
    }

    bool has_history() const {
        return has_prev_step;
    }

    // Semi-implicit Euler step: 1st-order in time, spectral in space
    // c_hat^{n+1} = (c_hat^n - dt * M0 * k^2 * (mu_hat^n - S * c_hat^n)) / (1 + dt * M0 * (kappa * k^4 + S * k^2))
    void step_semi_implicit_euler(std::vector<std::vector<double>>& c_fields, double dt, double S = 0.0) {
        double M0 = get_base_mobility();
        (void)M0;
        size_t interior_offset = (size_t)grid.ghost * grid.slice_size();

        // 1. Evaluate bulk chemical potentials point-wise in physical space
        std::vector<double> c_local(n_comp);
        std::vector<double> df_dc(n_comp);

        for (size_t idx = 0; idx < local_cells; ++idx) {
            for (int m = 0; m < n_comp; ++m) {
                c_local[m] = c_fields[m][interior_offset + idx];
            }
            free_energy->chemical_derivatives(c_local, df_dc);

            // Projected chemical potential: subtract mean so sum_m tilde_mu = 0
            double sum_df = 0.0;
            for (int m = 0; m < n_comp; ++m) sum_df += df_dc[m];
            double mean_df = sum_df / n_comp;

            for (int m = 0; m < n_comp; ++m) {
                c_hat[m][idx][0] = c_local[m];
                c_hat[m][idx][1] = 0.0;
                mu_hat[m][idx][0] = df_dc[m] - mean_df;
                mu_hat[m][idx][1] = 0.0;
            }
        }

        // 2. Forward FFTs for each component
        for (int m = 0; m < n_comp; ++m) {
            fftw_mpi_execute_dft(plan_fwd, c_hat[m], c_hat[m]);
            fftw_mpi_execute_dft(plan_fwd, mu_hat[m], mu_hat[m]);
        }

        // 3. Save current state to history buffers for SBDF2
        for (int m = 0; m < n_comp; ++m) {
            for (size_t idx = 0; idx < local_cells; ++idx) {
                c_prev[m][idx][0] = c_hat[m][idx][0];
                c_prev[m][idx][1] = c_hat[m][idx][1];
                mu_prev[m][idx][0] = mu_hat[m][idx][0];
                mu_prev[m][idx][1] = mu_hat[m][idx][1];
            }
        }
        has_prev_step = true;

        // 4. Update in Fourier space
        if (is_kappa_diagonal && is_M_diagonal) {
            for (int m = 0; m < n_comp; ++m) {
                double kap = kappa[m][m];
                double M_m = M_matrix[m * n_comp + m];
                for (size_t idx = 0; idx < local_cells; ++idx) {
                    double k_sq = k2[idx];
                    double k_4 = k4[idx];
                    double denom = 1.0 + dt * M_m * (kap * k_4 + S * k_sq);
                    double dt_M_k2 = dt * M_m * k_sq;

                    double num_r = c_hat[m][idx][0] - dt_M_k2 * (mu_hat[m][idx][0] - S * c_hat[m][idx][0]);
                    double num_i = c_hat[m][idx][1] - dt_M_k2 * (mu_hat[m][idx][1] - S * c_hat[m][idx][1]);

                    c_hat[m][idx][0] = num_r / denom;
                    c_hat[m][idx][1] = num_i / denom;
                }
            }
        } else {
            // General matrix kappa and/or general matrix mobility
            std::vector<double> M_kappa(n_comp * n_comp, 0.0);
            for (int i = 0; i < n_comp; ++i) {
                for (int j = 0; j < n_comp; ++j) {
                    double sum = 0.0;
                    for (int p = 0; p < n_comp; ++p) {
                        sum += M_matrix[i * n_comp + p] * kappa[p][j];
                    }
                    M_kappa[i * n_comp + j] = sum;
                }
            }

            std::vector<double> A(n_comp * n_comp);
            std::vector<double> br(n_comp), bi(n_comp);
            std::vector<double> xr(n_comp), xi(n_comp);

            for (size_t idx = 0; idx < local_cells; ++idx) {
                double k_sq = k2[idx];
                double k_4 = k4[idx];

                for (int i = 0; i < n_comp; ++i) {
                    for (int j = 0; j < n_comp; ++j) {
                        A[i * n_comp + j] = dt * k_4 * M_kappa[i * n_comp + j]
                                          + dt * S * k_sq * M_matrix[i * n_comp + j];
                        if (i == j) {
                            A[i * n_comp + j] += 1.0;
                        }
                    }

                    double diff_mu_r = 0.0;
                    double diff_mu_i = 0.0;
                    for (int j = 0; j < n_comp; ++j) {
                        diff_mu_r += M_matrix[i * n_comp + j] * (mu_hat[j][idx][0] - S * c_hat[j][idx][0]);
                        diff_mu_i += M_matrix[i * n_comp + j] * (mu_hat[j][idx][1] - S * c_hat[j][idx][1]);
                    }

                    br[i] = c_hat[i][idx][0] - dt * k_sq * diff_mu_r;
                    bi[i] = c_hat[i][idx][1] - dt * k_sq * diff_mu_i;
                }

                solve_small_system(n_comp, A.data(), xr.data(), xi.data(), br.data(), bi.data());

                for (int m = 0; m < n_comp; ++m) {
                    c_hat[m][idx][0] = xr[m];
                    c_hat[m][idx][1] = xi[m];
                }
            }
        }

        // 5. Backward FFTs & normalize
        for (int m = 0; m < n_comp; ++m) {
            fftw_mpi_execute_dft(plan_bwd, c_hat[m], c_hat[m]);
            for (size_t idx = 0; idx < local_cells; ++idx) {
                c_fields[m][interior_offset + idx] = c_hat[m][idx][0] * fft_norm;
            }
        }
    }

    // Second-order Semi-Implicit BDF2 (SBDF2) multi-step scheme:
    void step_sbdf2(std::vector<std::vector<double>>& c_fields, double dt, double S = 0.0) {
        if (!has_prev_step) {
            // First step bootstrap via semi-implicit Euler
            step_semi_implicit_euler(c_fields, dt, S);
            return;
        }

        double M0 = get_base_mobility();
        (void)M0;
        size_t interior_offset = (size_t)grid.ghost * grid.slice_size();

        // 1. Evaluate bulk chemical potentials point-wise in physical space
        std::vector<double> c_local(n_comp);
        std::vector<double> df_dc(n_comp);

        for (size_t idx = 0; idx < local_cells; ++idx) {
            for (int m = 0; m < n_comp; ++m) {
                c_local[m] = c_fields[m][interior_offset + idx];
            }
            free_energy->chemical_derivatives(c_local, df_dc);

            double sum_df = 0.0;
            for (int m = 0; m < n_comp; ++m) sum_df += df_dc[m];
            double mean_df = sum_df / n_comp;

            for (int m = 0; m < n_comp; ++m) {
                c_hat[m][idx][0] = c_local[m];
                c_hat[m][idx][1] = 0.0;
                mu_hat[m][idx][0] = df_dc[m] - mean_df;
                mu_hat[m][idx][1] = 0.0;
            }
        }

        // 2. Forward FFTs for current state (n)
        for (int m = 0; m < n_comp; ++m) {
            fftw_mpi_execute_dft(plan_fwd, c_hat[m], c_hat[m]);
            fftw_mpi_execute_dft(plan_fwd, mu_hat[m], mu_hat[m]);
        }

        // 3. Update in Fourier space using SBDF2
        if (is_kappa_diagonal && is_M_diagonal) {
            for (int m = 0; m < n_comp; ++m) {
                double kap = kappa[m][m];
                double M_m = M_matrix[m * n_comp + m];
                for (size_t idx = 0; idx < local_cells; ++idx) {
                    double k_sq = k2[idx];
                    double k_4 = k4[idx];
                    double denom = 3.0 + 2.0 * dt * M_m * (kap * k_4 + S * k_sq);
                    double two_dt_M_k2 = 2.0 * dt * M_m * k_sq;

                    double c_n_r = c_hat[m][idx][0];
                    double c_n_i = c_hat[m][idx][1];
                    double c_prev_r = c_prev[m][idx][0];
                    double c_prev_i = c_prev[m][idx][1];

                    double mu_n_r = mu_hat[m][idx][0];
                    double mu_n_i = mu_hat[m][idx][1];
                    double mu_prev_r = mu_prev[m][idx][0];
                    double mu_prev_i = mu_prev[m][idx][1];

                    // Extrapolated explicit potentials: 2*mu^n - mu^{n-1}
                    double ext_mu_r = 2.0 * mu_n_r - mu_prev_r;
                    double ext_mu_i = 2.0 * mu_n_i - mu_prev_i;

                    // Extrapolated explicit stabilization: 2*c^n - c^{n-1}
                    double ext_c_r = 2.0 * c_n_r - c_prev_r;
                    double ext_c_i = 2.0 * c_n_i - c_prev_i;

                    double rhs_r = 4.0 * c_n_r - c_prev_r - two_dt_M_k2 * (ext_mu_r - S * ext_c_r);
                    double rhs_i = 4.0 * c_n_i - c_prev_i - two_dt_M_k2 * (ext_mu_i - S * ext_c_i);

                    // Update history for next step
                    c_prev[m][idx][0] = c_n_r;
                    c_prev[m][idx][1] = c_n_i;
                    mu_prev[m][idx][0] = mu_n_r;
                    mu_prev[m][idx][1] = mu_n_i;

                    // Store new solution
                    c_hat[m][idx][0] = rhs_r / denom;
                    c_hat[m][idx][1] = rhs_i / denom;
                }
            }
        } else {
            // General matrix kappa and/or mobility with SBDF2
            std::vector<double> M_kappa(n_comp * n_comp, 0.0);
            for (int i = 0; i < n_comp; ++i) {
                for (int j = 0; j < n_comp; ++j) {
                    double sum = 0.0;
                    for (int p = 0; p < n_comp; ++p) {
                        sum += M_matrix[i * n_comp + p] * kappa[p][j];
                    }
                    M_kappa[i * n_comp + j] = sum;
                }
            }

            std::vector<double> A(n_comp * n_comp);
            std::vector<double> br(n_comp), bi(n_comp);
            std::vector<double> xr(n_comp), xi(n_comp);

            for (size_t idx = 0; idx < local_cells; ++idx) {
                double k_sq = k2[idx];
                double k_4 = k4[idx];

                for (int i = 0; i < n_comp; ++i) {
                    for (int j = 0; j < n_comp; ++j) {
                        A[i * n_comp + j] = 2.0 * dt * k_4 * M_kappa[i * n_comp + j]
                                          + 2.0 * dt * S * k_sq * M_matrix[i * n_comp + j];
                        if (i == j) {
                            A[i * n_comp + j] += 3.0;
                        }
                    }

                    double c_n_r = c_hat[i][idx][0];
                    double c_n_i = c_hat[i][idx][1];
                    double c_prev_r = c_prev[i][idx][0];
                    double c_prev_i = c_prev[i][idx][1];

                    double ext_diff_r = 0.0;
                    double ext_diff_i = 0.0;
                    for (int j = 0; j < n_comp; ++j) {
                        double ext_mu_r = 2.0 * mu_hat[j][idx][0] - mu_prev[j][idx][0];
                        double ext_mu_i = 2.0 * mu_hat[j][idx][1] - mu_prev[j][idx][1];
                        double ext_c_r = 2.0 * c_hat[j][idx][0] - c_prev[j][idx][0];
                        double ext_c_i = 2.0 * c_hat[j][idx][1] - c_prev[j][idx][1];
                        double dr = ext_mu_r - S * ext_c_r;
                        double di = ext_mu_i - S * ext_c_i;
                        ext_diff_r += M_matrix[i * n_comp + j] * dr;
                        ext_diff_i += M_matrix[i * n_comp + j] * di;
                    }

                    br[i] = 4.0 * c_n_r - c_prev_r - 2.0 * dt * k_sq * ext_diff_r;
                    bi[i] = 4.0 * c_n_i - c_prev_i - 2.0 * dt * k_sq * ext_diff_i;

                    // Update history
                    c_prev[i][idx][0] = c_n_r;
                    c_prev[i][idx][1] = c_n_i;
                    mu_prev[i][idx][0] = mu_hat[i][idx][0];
                    mu_prev[i][idx][1] = mu_hat[i][idx][1];
                }

                solve_small_system(n_comp, A.data(), xr.data(), xi.data(), br.data(), bi.data());

                for (int m = 0; m < n_comp; ++m) {
                    c_hat[m][idx][0] = xr[m];
                    c_hat[m][idx][1] = xi[m];
                }
            }
        }

        // 4. Backward FFTs & normalize
        for (int m = 0; m < n_comp; ++m) {
            fftw_mpi_execute_dft(plan_bwd, c_hat[m], c_hat[m]);
            for (size_t idx = 0; idx < local_cells; ++idx) {
                c_fields[m][interior_offset + idx] = c_hat[m][idx][0] * fft_norm;
            }
        }
    }

private:
    double get_base_mobility() const {
        if (!M_matrix.empty()) {
            return M_matrix[0];
        }
        if (mobility->type() == MobilityType::CONSTANT) {
            return static_cast<ConstantMobility*>(mobility.get())->value();
        } else if (mobility->type() == MobilityType::DEGENERATE) {
            return static_cast<DegenerateMobility*>(mobility.get())->base_value();
        }
        return 1.0;
    }

    void precompute_wavenumbers() {
        k2.resize(local_cells);
        k4.resize(local_cells);

        double Lx = grid.nx * grid.dx;
        double Ly = grid.ny * grid.dy;
        double Lz = (grid.dim == 3) ? (grid.nz * grid.dz) : 1.0;

        if (grid.dim == 2) {
            for (ptrdiff_t i = 0; i < local_nx; ++i) {
                ptrdiff_t iglob = local_x_start + i;
                double mx = (iglob <= grid.nx / 2) ? (double)iglob : (double)(iglob - grid.nx);
                double kx = 2.0 * M_PI * mx / Lx;

                for (int j = 0; j < grid.ny; ++j) {
                    double my = (j <= grid.ny / 2) ? (double)j : (double)(j - grid.ny);
                    double ky = 2.0 * M_PI * my / Ly;

                    size_t idx = (size_t)i * grid.ny + j;
                    double k_sq = kx * kx + ky * ky;
                    k2[idx] = k_sq;
                    k4[idx] = k_sq * k_sq;
                }
            }
        } else {
            for (ptrdiff_t i = 0; i < local_nx; ++i) {
                ptrdiff_t iglob = local_x_start + i;
                double mx = (iglob <= grid.nx / 2) ? (double)iglob : (double)(iglob - grid.nx);
                double kx = 2.0 * M_PI * mx / Lx;

                for (int j = 0; j < grid.ny; ++j) {
                    double my = (j <= grid.ny / 2) ? (double)j : (double)(j - grid.ny);
                    double ky = 2.0 * M_PI * my / Ly;

                    for (int k = 0; k < grid.nz; ++k) {
                        double mz = (k <= grid.nz / 2) ? (double)k : (double)(k - grid.nz);
                        double kz = 2.0 * M_PI * mz / Lz;

                        size_t idx = ((size_t)i * grid.ny + j) * grid.nz + k;
                        double k_sq = kx * kx + ky * ky + kz * kz;
                        k2[idx] = k_sq;
                        k4[idx] = k_sq * k_sq;
                    }
                }
            }
        }
    }

    // Direct Gaussian elimination with partial pivoting for small N x N real systems with 2 RHS vectors
    static void solve_small_system(int n, const double* A_in, double* xr, double* xi,
                                  const double* br, const double* bi) {
        constexpr int MAX_N = 16;
        if (n > MAX_N) {
            throw std::runtime_error("Component count exceeds small system max (16).");
        }

        double A[MAX_N][MAX_N];
        double rhs_r[MAX_N];
        double rhs_i[MAX_N];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = A_in[i * n + j];
            }
            rhs_r[i] = br[i];
            rhs_i[i] = bi[i];
        }

        for (int p = 0; p < n; ++p) {
            int max_r = p;
            double max_val = std::abs(A[p][p]);
            for (int r = p + 1; r < n; ++r) {
                if (std::abs(A[r][p]) > max_val) {
                    max_val = std::abs(A[r][p]);
                    max_r = r;
                }
            }
            if (max_r != p) {
                for (int j = p; j < n; ++j) std::swap(A[p][j], A[max_r][j]);
                std::swap(rhs_r[p], rhs_r[max_r]);
                std::swap(rhs_i[p], rhs_i[max_r]);
            }
            double pivot = A[p][p];
            if (std::abs(pivot) < 1e-15) pivot = (pivot >= 0 ? 1e-15 : -1e-15);
            for (int r = p + 1; r < n; ++r) {
                double factor = A[r][p] / pivot;
                for (int j = p; j < n; ++j) {
                    A[r][j] -= factor * A[p][j];
                }
                rhs_r[r] -= factor * rhs_r[p];
                rhs_i[r] -= factor * rhs_i[p];
            }
        }

        for (int i = n - 1; i >= 0; --i) {
            double sum_r = rhs_r[i];
            double sum_i = rhs_i[i];
            for (int j = i + 1; j < n; ++j) {
                sum_r -= A[i][j] * xr[j];
                sum_i -= A[i][j] * xi[j];
            }
            double diag = A[i][i];
            if (std::abs(diag) < 1e-15) diag = (diag >= 0 ? 1e-15 : -1e-15);
            xr[i] = sum_r / diag;
            xi[i] = sum_i / diag;
        }
    }
};

} // namespace mcch
