#pragma once

#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"
#include "FFTWSolver.hpp"
#include "CudaCommon.hpp"
#include "CudaSolver.hpp"

#include <vector>
#include <memory>
#include <random>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace mcch {

enum class TimeIntegratorType {
    EULER,
    RK2,
    RK4,
    SEMI_IMPLICIT,
    SEMI_IMPLICIT_BDF2
};

struct Diagnostics {
    int step;
    double time;
    double total_free_energy;
    double bulk_free_energy;
    double gradient_free_energy;
    std::vector<double> average_composition; // size n_comp
    double max_unity_deviation; // max |sum(c_i) - 1|
    double min_composition;
    double max_composition;
    double dF_dt; // rate of free energy change
};

class Solver {
public:
    MPIContext& mpi;
    Grid grid;
    int n_comp;
    double T;
    std::shared_ptr<FreeEnergy> free_energy;
    std::shared_ptr<MobilityModel> mobility;

    // Gradient energy coefficient matrix kappa_ij (n_comp x n_comp)
    std::vector<std::vector<double>> kappa;

    // Concentrations: c[m] has size grid.local_total_size()
    std::vector<std::vector<double>> c;

    // Chemical potentials: mu[m] has size grid.local_total_size()
    std::vector<std::vector<double>> mu;

    // RHS (time derivatives): rhs[m] has size grid.local_total_size()
    std::vector<std::vector<double>> rhs;

    // Scratch buffers for multi-stage time integrators
    std::vector<std::vector<double>> c_stage;
    std::vector<std::vector<double>> k1;
    std::vector<std::vector<double>> k2;
    std::vector<std::vector<double>> k3;
    std::vector<std::vector<double>> k4;

    double current_time;
    int current_step;
    double prev_energy;

    std::unique_ptr<FFTWSemiImplicitSolver> fftw_solver;
    DeviceType device_type = DeviceType::CPU;
#ifdef MCCH_ENABLE_CUDA
    std::unique_ptr<CudaSolver> gpu_solver;
#endif

    void set_device(DeviceType dev) {
        device_type = dev;
        if (device_type == DeviceType::GPU) {
#ifdef MCCH_ENABLE_CUDA
            if (!gpu_solver) {
                gpu_solver = std::make_unique<CudaSolver>(mpi, grid, n_comp, free_energy, mobility, kappa);
                gpu_solver->copy_to_device(c);
            }
#else
            throw std::runtime_error("MCCH solver was built without CUDA support. Rebuild with a CUDA compiler (nvcc) to run on GPU.");
#endif
        }
    }

    void sync_from_device() {
#ifdef MCCH_ENABLE_CUDA
        if (device_type == DeviceType::GPU && gpu_solver) {
            gpu_solver->copy_to_host(c);
            gpu_solver->copy_mu_to_host(mu);
        }
#endif
    }

    Solver(MPIContext& mpi_context, const Grid& grid_in, int num_components,
           std::shared_ptr<FreeEnergy> fe, std::shared_ptr<MobilityModel> mob,
           const std::vector<std::vector<double>>& kappa_matrix)
        : mpi(mpi_context), grid(grid_in), n_comp(num_components),
          free_energy(fe), mobility(mob), kappa(kappa_matrix),
          current_time(0.0), current_step(0), prev_energy(0.0) {

        if ((int)kappa.size() != n_comp) {
            throw std::invalid_argument("Kappa matrix row count != num_components");
        }

        size_t total_local_cells = grid.local_total_size();
        c.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        mu.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        rhs.resize(n_comp, std::vector<double>(total_local_cells, 0.0));

        c_stage.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        k1.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        k2.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        k3.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
        k4.resize(n_comp, std::vector<double>(total_local_cells, 0.0));
    }

    // Initialize with uniform mean concentration plus random perturbations
    void initialize_random(const std::vector<double>& c_mean, double noise_amp, unsigned int seed = 42) {
        if ((int)c_mean.size() != n_comp) {
            throw std::invalid_argument("c_mean size != n_comp");
        }

        std::mt19937 gen(seed + mpi.rank * 1000);
        std::uniform_real_distribution<double> dist(-noise_amp, noise_amp);

        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            for (int j = 0; j < grid.ny; ++j) {
                for (int k = 0; k < grid.nz; ++k) {
                    size_t idx = grid.idx(i_loc, j, k);
                    double sum = 0.0;
                    for (int m = 0; m < n_comp - 1; ++m) {
                        c[m][idx] = std::max(1e-4, std::min(0.999, c_mean[m] + dist(gen)));
                        sum += c[m][idx];
                    }
                    // Final component enforces sum = 1
                    c[n_comp - 1][idx] = 1.0 - sum;
                    if (c[n_comp - 1][idx] < 0.0) {
                        // Normalize if boundary overshoot
                        sum += c[n_comp - 1][idx];
                        for (int m = 0; m < n_comp; ++m) {
                            c[m][idx] /= sum;
                        }
                    }
                }
            }
        }
        exchange_all_ghosts(c);
#ifdef MCCH_ENABLE_CUDA
        if (device_type == DeviceType::GPU && gpu_solver) {
            gpu_solver->copy_to_device(c);
        }
#endif
        prev_energy = compute_total_free_energy();
    }

    // Initialize droplet/inclusion geometry
    void initialize_droplet(int inside_comp, double cx, double cy, double cz, double radius, double diffuse_width) {
        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            double x = grid.coord_x(i_loc);
            for (int j = 0; j < grid.ny; ++j) {
                double y = grid.coord_y(j);
                for (int k = 0; k < grid.nz; ++k) {
                    double z = grid.coord_z(k);
                    double r = 0.0;
                    if (grid.dim == 2) {
                        r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
                    } else {
                        r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy) + (z - cz)*(z - cz));
                    }
                    // Smooth hyperbolic tangent profile
                    double phi = 0.5 * (1.0 - std::tanh((r - radius) / (2.0 * diffuse_width)));
                    size_t idx = grid.idx(i_loc, j, k);

                    double remainder = (1.0 - phi) / (n_comp - 1);
                    for (int m = 0; m < n_comp; ++m) {
                        if (m == inside_comp) {
                            c[m][idx] = phi;
                        } else {
                            c[m][idx] = remainder;
                        }
                    }
                }
            }
        }
        exchange_all_ghosts(c);
#ifdef MCCH_ENABLE_CUDA
        if (device_type == DeviceType::GPU && gpu_solver) {
            gpu_solver->copy_to_device(c);
        }
#endif
        prev_energy = compute_total_free_energy();
    }

    void exchange_all_ghosts(std::vector<std::vector<double>>& fields) {
        for (int m = 0; m < n_comp; ++m) {
            mpi.exchange_ghosts(fields[m].data());
        }
    }

    // Compute chemical potentials mu_m for all species
    void compute_chemical_potentials(const std::vector<std::vector<double>>& c_in,
                                     std::vector<std::vector<double>>& mu_out) {
        double idx2 = 1.0 / (grid.dx * grid.dx);
        double idy2 = 1.0 / (grid.dy * grid.dy);
        double idz2 = (grid.dim == 3) ? (1.0 / (grid.dz * grid.dz)) : 0.0;

        std::vector<double> c_local(n_comp);
        std::vector<double> df_dc(n_comp);
        std::vector<double> lap_c(n_comp);

        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            int i_alloc = i_loc + grid.ghost;
            for (int j = 0; j < grid.ny; ++j) {
                int j_m, j_p;
                grid.get_y_neighbors(j, j_m, j_p);

                for (int k = 0; k < grid.nz; ++k) {
                    int k_m, k_p;
                    grid.get_z_neighbors(k, k_m, k_p);

                    size_t curr_idx = grid.raw_idx(i_alloc, j, k);
                    size_t xm_idx   = grid.raw_idx(i_alloc - 1, j, k);
                    size_t xp_idx   = grid.raw_idx(i_alloc + 1, j, k);
                    size_t ym_idx   = grid.raw_idx(i_alloc, j_m, k);
                    size_t yp_idx   = grid.raw_idx(i_alloc, j_p, k);

                    for (int m = 0; m < n_comp; ++m) {
                        c_local[m] = c_in[m][curr_idx];
                        
                        double lap = (c_in[m][xp_idx] - 2.0 * c_local[m] + c_in[m][xm_idx]) * idx2
                                   + (c_in[m][yp_idx] - 2.0 * c_local[m] + c_in[m][ym_idx]) * idy2;
                        if (grid.dim == 3) {
                            size_t zm_idx = grid.raw_idx(i_alloc, j, k_m);
                            size_t zp_idx = grid.raw_idx(i_alloc, j, k_p);
                            lap += (c_in[m][zp_idx] - 2.0 * c_local[m] + c_in[m][zm_idx]) * idz2;
                        }
                        lap_c[m] = lap;
                    }

                    // Bulk thermodynamic derivatives df0/dc_m
                    free_energy->chemical_derivatives(c_local, df_dc);

                    // Compute raw chemical potential: mu_m = df0/dc_m - sum_n (kappa_mn * lap_c_n)
                    double mu_sum = 0.0;
                    for (int m = 0; m < n_comp; ++m) {
                        double grad_term = 0.0;
                        for (int n = 0; n < n_comp; ++n) {
                            grad_term += kappa[m][n] * lap_c[n];
                        }
                        mu_out[m][curr_idx] = df_dc[m] - grad_term;
                        mu_sum += mu_out[m][curr_idx];
                    }

                    // Project potentials: mu_tilde_m = mu_m - (1/N)*sum(mu_k)
                    // This subtracts the Gibbs-Duhem Lagrange multiplier ensuring exact sum(J_i) = 0
                    double mu_mean = mu_sum / n_comp;
                    for (int m = 0; m < n_comp; ++m) {
                        mu_out[m][curr_idx] -= mu_mean;
                    }
                }
            }
        }

        // Exchange ghosts for chemical potentials
        exchange_all_ghosts(mu_out);
    }

    // Compute RHS (time derivatives dc_m/dt = -div(J_m))
    void compute_rhs(const std::vector<std::vector<double>>& c_in,
                     std::vector<std::vector<double>>& mu_in,
                     std::vector<std::vector<double>>& rhs_out) {

        // Step 1: Compute chemical potentials and exchange ghosts
        compute_chemical_potentials(c_in, mu_in);

        double idx = 1.0 / grid.dx;
        double idy = 1.0 / grid.dy;
        double idz = (grid.dim == 3) ? (1.0 / grid.dz) : 0.0;

        double idx2 = idx * idx;
        double idy2 = idy * idy;
        double idz2 = idz * idz;

        bool is_const_mob = mobility->is_constant();
        if (is_const_mob) {
            std::vector<double> dummy_c(n_comp, 1.0 / n_comp);
            std::vector<double> M_mat(n_comp * n_comp);
            mobility->evaluate_matrix(dummy_c, M_mat);

            bool is_diag = true;
            bool is_scalar = true;
            double M0 = M_mat.empty() ? 1.0 : M_mat[0];
            for (int i = 0; i < n_comp; ++i) {
                for (int j = 0; j < n_comp; ++j) {
                    double v = M_mat[i * n_comp + j];
                    if (i != j && std::abs(v) > 1e-14) is_diag = false;
                    if (i == j && std::abs(v - M0) > 1e-14) is_scalar = false;
                }
            }

            if (is_diag && is_scalar) {
                // Fast path for constant isotropic scalar mobility: div(M grad mu) = M0 * lap(mu)
                for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                    int i_alloc = i_loc + grid.ghost;
                    for (int j = 0; j < grid.ny; ++j) {
                        int j_m, j_p;
                        grid.get_y_neighbors(j, j_m, j_p);

                        for (int k = 0; k < grid.nz; ++k) {
                            int k_m, k_p;
                            grid.get_z_neighbors(k, k_m, k_p);

                            size_t curr_idx = grid.raw_idx(i_alloc, j, k);
                            size_t xm_idx   = grid.raw_idx(i_alloc - 1, j, k);
                            size_t xp_idx   = grid.raw_idx(i_alloc + 1, j, k);
                            size_t ym_idx   = grid.raw_idx(i_alloc, j_m, k);
                            size_t yp_idx   = grid.raw_idx(i_alloc, j_p, k);

                            for (int m = 0; m < n_comp; ++m) {
                                double lap_mu = (mu_in[m][xp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][xm_idx]) * idx2
                                              + (mu_in[m][yp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][ym_idx]) * idy2;
                                if (grid.dim == 3) {
                                    size_t zm_idx = grid.raw_idx(i_alloc, j, k_m);
                                    size_t zp_idx = grid.raw_idx(i_alloc, j, k_p);
                                    lap_mu += (mu_in[m][zp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][zm_idx]) * idz2;
                                }
                                rhs_out[m][curr_idx] = M0 * lap_mu;
                            }
                        }
                    }
                }
            } else {
                // General constant mobility matrix: div(M grad mu) = sum_n M_mn * lap(mu_n)
                std::vector<double> lap_mu(n_comp);
                for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                    int i_alloc = i_loc + grid.ghost;
                    for (int j = 0; j < grid.ny; ++j) {
                        int j_m, j_p;
                        grid.get_y_neighbors(j, j_m, j_p);

                        for (int k = 0; k < grid.nz; ++k) {
                            int k_m, k_p;
                            grid.get_z_neighbors(k, k_m, k_p);

                            size_t curr_idx = grid.raw_idx(i_alloc, j, k);
                            size_t xm_idx   = grid.raw_idx(i_alloc - 1, j, k);
                            size_t xp_idx   = grid.raw_idx(i_alloc + 1, j, k);
                            size_t ym_idx   = grid.raw_idx(i_alloc, j_m, k);
                            size_t yp_idx   = grid.raw_idx(i_alloc, j_p, k);

                            for (int m = 0; m < n_comp; ++m) {
                                double lap = (mu_in[m][xp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][xm_idx]) * idx2
                                           + (mu_in[m][yp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][ym_idx]) * idy2;
                                if (grid.dim == 3) {
                                    size_t zm_idx = grid.raw_idx(i_alloc, j, k_m);
                                    size_t zp_idx = grid.raw_idx(i_alloc, j, k_p);
                                    lap += (mu_in[m][zp_idx] - 2.0 * mu_in[m][curr_idx] + mu_in[m][zm_idx]) * idz2;
                                }
                                lap_mu[m] = lap;
                            }

                            for (int m = 0; m < n_comp; ++m) {
                                double sum = 0.0;
                                for (int n = 0; n < n_comp; ++n) {
                                    sum += M_mat[m * n_comp + n] * lap_mu[n];
                                }
                                rhs_out[m][curr_idx] = sum;
                            }
                        }
                    }
                }
            }
        } else {
            // Conservative face-flux formulation for degenerate or matrix mobility:
            // J_{m, face} = - sum_n M_{mn}(c_face) * (mu_n,next - mu_n,curr) / dx
            // rhs_m = - (J_{m, face_p} - J_{m, face_m}) / dx
            std::vector<double> c_face(n_comp);
            std::vector<double> M_mat(n_comp * n_comp);
            std::vector<double> J_xp(n_comp), J_xm(n_comp);
            std::vector<double> J_yp(n_comp), J_ym(n_comp);
            std::vector<double> J_zp(n_comp), J_zm(n_comp);

            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                int i_alloc = i_loc + grid.ghost;
                for (int j = 0; j < grid.ny; ++j) {
                    int j_m, j_p;
                    grid.get_y_neighbors(j, j_m, j_p);

                    for (int k = 0; k < grid.nz; ++k) {
                        int k_m, k_p;
                        grid.get_z_neighbors(k, k_m, k_p);

                        size_t curr_idx = grid.raw_idx(i_alloc, j, k);
                        size_t xm_idx   = grid.raw_idx(i_alloc - 1, j, k);
                        size_t xp_idx   = grid.raw_idx(i_alloc + 1, j, k);
                        size_t ym_idx   = grid.raw_idx(i_alloc, j_m, k);
                        size_t yp_idx   = grid.raw_idx(i_alloc, j_p, k);

                        // 1. Face x+1/2
                        for (int m = 0; m < n_comp; ++m) {
                            c_face[m] = 0.5 * (c_in[m][curr_idx] + c_in[m][xp_idx]);
                        }
                        mobility->evaluate_matrix(c_face, M_mat);
                        for (int m = 0; m < n_comp; ++m) {
                            double flux = 0.0;
                            for (int n = 0; n < n_comp; ++n) {
                                flux -= M_mat[m * n_comp + n] * (mu_in[n][xp_idx] - mu_in[n][curr_idx]) * idx;
                            }
                            J_xp[m] = flux;
                        }

                        // 2. Face x-1/2
                        for (int m = 0; m < n_comp; ++m) {
                            c_face[m] = 0.5 * (c_in[m][xm_idx] + c_in[m][curr_idx]);
                        }
                        mobility->evaluate_matrix(c_face, M_mat);
                        for (int m = 0; m < n_comp; ++m) {
                            double flux = 0.0;
                            for (int n = 0; n < n_comp; ++n) {
                                flux -= M_mat[m * n_comp + n] * (mu_in[n][curr_idx] - mu_in[n][xm_idx]) * idx;
                            }
                            J_xm[m] = flux;
                        }

                        // 3. Face y+1/2
                        for (int m = 0; m < n_comp; ++m) {
                            c_face[m] = 0.5 * (c_in[m][curr_idx] + c_in[m][yp_idx]);
                        }
                        mobility->evaluate_matrix(c_face, M_mat);
                        for (int m = 0; m < n_comp; ++m) {
                            double flux = 0.0;
                            for (int n = 0; n < n_comp; ++n) {
                                flux -= M_mat[m * n_comp + n] * (mu_in[n][yp_idx] - mu_in[n][curr_idx]) * idy;
                            }
                            J_yp[m] = flux;
                        }

                        // 4. Face y-1/2
                        for (int m = 0; m < n_comp; ++m) {
                            c_face[m] = 0.5 * (c_in[m][ym_idx] + c_in[m][curr_idx]);
                        }
                        mobility->evaluate_matrix(c_face, M_mat);
                        for (int m = 0; m < n_comp; ++m) {
                            double flux = 0.0;
                            for (int n = 0; n < n_comp; ++n) {
                                flux -= M_mat[m * n_comp + n] * (mu_in[n][curr_idx] - mu_in[n][ym_idx]) * idy;
                            }
                            J_ym[m] = flux;
                        }

                        // 5. Z faces if 3D
                        if (grid.dim == 3) {
                            size_t zm_idx = grid.raw_idx(i_alloc, j, k_m);
                            size_t zp_idx = grid.raw_idx(i_alloc, j, k_p);

                            for (int m = 0; m < n_comp; ++m) {
                                c_face[m] = 0.5 * (c_in[m][curr_idx] + c_in[m][zp_idx]);
                            }
                            mobility->evaluate_matrix(c_face, M_mat);
                            for (int m = 0; m < n_comp; ++m) {
                                double flux = 0.0;
                                for (int n = 0; n < n_comp; ++n) {
                                    flux -= M_mat[m * n_comp + n] * (mu_in[n][zp_idx] - mu_in[n][curr_idx]) * idz;
                                }
                                J_zp[m] = flux;
                            }

                            for (int m = 0; m < n_comp; ++m) {
                                c_face[m] = 0.5 * (c_in[m][zm_idx] + c_in[m][curr_idx]);
                            }
                            mobility->evaluate_matrix(c_face, M_mat);
                            for (int m = 0; m < n_comp; ++m) {
                                double flux = 0.0;
                                for (int n = 0; n < n_comp; ++n) {
                                    flux -= M_mat[m * n_comp + n] * (mu_in[n][curr_idx] - mu_in[n][zm_idx]) * idz;
                                }
                                J_zm[m] = flux;
                            }
                        }

                        // Divergence: rhs = -div(J) = - [ (J_xp - J_xm)/dx + (J_yp - J_ym)/dy + (J_zp - J_zm)/dz ]
                        for (int m = 0; m < n_comp; ++m) {
                            double div_J = (J_xp[m] - J_xm[m]) * idx + (J_yp[m] - J_ym[m]) * idy;
                            if (grid.dim == 3) {
                                div_J += (J_zp[m] - J_zm[m]) * idz;
                            }
                            rhs_out[m][curr_idx] = -div_J;
                        }
                    }
                }
            }
        }
    }

    // Time integration: Forward Euler
    void step_euler(double dt) {
        if (device_type == DeviceType::GPU) {
#ifdef MCCH_ENABLE_CUDA
            if (!gpu_solver) set_device(DeviceType::GPU);
            gpu_solver->step_euler(dt);
            current_time += dt;
            current_step++;
            return;
#else
            throw std::runtime_error("MCCH solver was built without CUDA support. Rebuild with a CUDA compiler (nvcc) to run on GPU.");
#endif
        }

        compute_rhs(c, mu, rhs);

        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c[m][idx] += dt * rhs[m][idx];
                    }
                }
            }
        }
        normalize_unity();
        exchange_all_ghosts(c);
        current_time += dt;
        current_step++;
    }

    // Time integration: 2nd-order Runge-Kutta (Heun / Midpoint)
    void step_rk2(double dt) {
        if (device_type == DeviceType::GPU) {
#ifdef MCCH_ENABLE_CUDA
            if (!gpu_solver) set_device(DeviceType::GPU);
            gpu_solver->step_rk2(dt);
            current_time += dt;
            current_step++;
            return;
#else
            throw std::runtime_error("MCCH solver was built without CUDA support. Rebuild with a CUDA compiler (nvcc) to run on GPU.");
#endif
        }

        // Stage 1: k1 = rhs(c^n)
        compute_rhs(c, mu, k1);

        // Intermediate state: c_stage = c^n + dt * k1
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c_stage[m][idx] = c[m][idx] + dt * k1[m][idx];
                    }
                }
            }
        }
        exchange_all_ghosts(c_stage);

        // Stage 2: k2 = rhs(c_stage)
        compute_rhs(c_stage, mu, k2);

        // Update: c^{n+1} = c^n + 0.5 * dt * (k1 + k2)
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c[m][idx] += 0.5 * dt * (k1[m][idx] + k2[m][idx]);
                    }
                }
            }
        }
        normalize_unity();
        exchange_all_ghosts(c);
        current_time += dt;
        current_step++;
    }

    // Time integration: Classical 4th-order Runge-Kutta (RK4)
    void step_rk4(double dt) {
        if (device_type == DeviceType::GPU) {
#ifdef MCCH_ENABLE_CUDA
            if (!gpu_solver) set_device(DeviceType::GPU);
            gpu_solver->step_rk4(dt);
            current_time += dt;
            current_step++;
            return;
#else
            throw std::runtime_error("MCCH solver was built without CUDA support. Rebuild with a CUDA compiler (nvcc) to run on GPU.");
#endif
        }

        // Stage 1: k1 = rhs(c)
        compute_rhs(c, mu, k1);

        // Stage 2: c_stage = c + 0.5*dt*k1
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c_stage[m][idx] = c[m][idx] + 0.5 * dt * k1[m][idx];
                    }
                }
            }
        }
        exchange_all_ghosts(c_stage);
        compute_rhs(c_stage, mu, k2);

        // Stage 3: c_stage = c + 0.5*dt*k2
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c_stage[m][idx] = c[m][idx] + 0.5 * dt * k2[m][idx];
                    }
                }
            }
        }
        exchange_all_ghosts(c_stage);
        compute_rhs(c_stage, mu, k3);

        // Stage 4: c_stage = c + dt*k3
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c_stage[m][idx] = c[m][idx] + dt * k3[m][idx];
                    }
                }
            }
        }
        exchange_all_ghosts(c_stage);
        compute_rhs(c_stage, mu, k4);

        // Final combination: c += (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        double dt6 = dt / 6.0;
        for (int m = 0; m < n_comp; ++m) {
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        c[m][idx] += dt6 * (k1[m][idx] + 2.0*k2[m][idx] + 2.0*k3[m][idx] + k4[m][idx]);
                    }
                }
            }
        }
        normalize_unity();
        exchange_all_ghosts(c);
        current_time += dt;
        current_step++;
    }

    void init_fftw_solver() {
        if (!fftw_solver) {
            fftw_solver = std::make_unique<FFTWSemiImplicitSolver>(mpi, grid, n_comp, free_energy, mobility, kappa);
        }
    }

    // Time integration: 1st-order Semi-Implicit Euler with FFTW MPI
    void step_semi_implicit(double dt, double stabilization = 0.0) {
        init_fftw_solver();
        fftw_solver->step_semi_implicit_euler(c, dt, stabilization);
        normalize_unity();
        exchange_all_ghosts(c);
        compute_chemical_potentials(c, mu);
        current_time += dt;
        current_step++;
    }

    // Time integration: 2nd-order Semi-Implicit BDF2 (SBDF2) with FFTW MPI
    void step_sbdf2(double dt, double stabilization = 0.0) {
        init_fftw_solver();
        fftw_solver->step_sbdf2(c, dt, stabilization);
        normalize_unity();
        exchange_all_ghosts(c);
        compute_chemical_potentials(c, mu);
        current_time += dt;
        current_step++;
    }

    // Step with chosen integrator
    void step(double dt, TimeIntegratorType method = TimeIntegratorType::RK2, double stabilization = 0.0) {
        if (method == TimeIntegratorType::EULER) {
            step_euler(dt);
        } else if (method == TimeIntegratorType::RK2) {
            step_rk2(dt);
        } else if (method == TimeIntegratorType::RK4) {
            step_rk4(dt);
        } else if (method == TimeIntegratorType::SEMI_IMPLICIT) {
            step_semi_implicit(dt, stabilization);
        } else if (method == TimeIntegratorType::SEMI_IMPLICIT_BDF2) {
            step_sbdf2(dt, stabilization);
        }
    }

    // Enforce sum_i c_i = 1 to prevent machine-precision round-off accumulation
    void normalize_unity() {
        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            for (int j = 0; j < grid.ny; ++j) {
                for (int k = 0; k < grid.nz; ++k) {
                    size_t idx = grid.idx(i_loc, j, k);
                    double sum = 0.0;
                    for (int m = 0; m < n_comp; ++m) {
                        sum += c[m][idx];
                    }
                    if (std::abs(sum - 1.0) > 1e-15) {
                        for (int m = 0; m < n_comp; ++m) {
                            c[m][idx] /= sum;
                        }
                    }
                }
            }
        }
    }

    // Calculate total Free Energy F = int (f0 + 0.5*sum(kappa_mn * grad(c_m) . grad(c_n))) dV
    double compute_total_free_energy() const {
        const_cast<Solver*>(this)->sync_from_device();
        double dV = grid.cell_volume();
        double idx = 0.5 / grid.dx;
        double idy = 0.5 / grid.dy;
        double idz = (grid.dim == 3) ? (0.5 / grid.dz) : 0.0;

        double local_energy = 0.0;
        std::vector<double> c_local(n_comp);
        std::vector<double> grad_x(n_comp), grad_y(n_comp), grad_z(n_comp);

        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            int i_alloc = i_loc + grid.ghost;
            for (int j = 0; j < grid.ny; ++j) {
                int j_m, j_p;
                grid.get_y_neighbors(j, j_m, j_p);

                for (int k = 0; k < grid.nz; ++k) {
                    int k_m, k_p;
                    grid.get_z_neighbors(k, k_m, k_p);

                    size_t curr_idx = grid.raw_idx(i_alloc, j, k);
                    size_t xm_idx   = grid.raw_idx(i_alloc - 1, j, k);
                    size_t xp_idx   = grid.raw_idx(i_alloc + 1, j, k);
                    size_t ym_idx   = grid.raw_idx(i_alloc, j_m, k);
                    size_t yp_idx   = grid.raw_idx(i_alloc, j_p, k);

                    for (int m = 0; m < n_comp; ++m) {
                        c_local[m] = c[m][curr_idx];
                        grad_x[m] = (c[m][xp_idx] - c[m][xm_idx]) * idx;
                        grad_y[m] = (c[m][yp_idx] - c[m][ym_idx]) * idy;
                        if (grid.dim == 3) {
                            size_t zm_idx = grid.raw_idx(i_alloc, j, k_m);
                            size_t zp_idx = grid.raw_idx(i_alloc, j, k_p);
                            grad_z[m] = (c[m][zp_idx] - c[m][zm_idx]) * idz;
                        } else {
                            grad_z[m] = 0.0;
                        }
                    }

                    // Bulk energy density
                    double f_bulk = free_energy->density(c_local);

                    // Gradient energy density
                    double f_grad = 0.0;
                    for (int m = 0; m < n_comp; ++m) {
                        for (int n = 0; n < n_comp; ++n) {
                            double dot = grad_x[m] * grad_x[n] + grad_y[m] * grad_y[n] + grad_z[m] * grad_z[n];
                            f_grad += 0.5 * kappa[m][n] * dot;
                        }
                    }

                    local_energy += (f_bulk + f_grad) * dV;
                }
            }
        }

        return mpi.allreduce_sum(local_energy);
    }

    // Diagnostics calculation
    Diagnostics compute_diagnostics(double dt) {
        sync_from_device();
        Diagnostics diag;
        diag.step = current_step;
        diag.time = current_time;

        double dV = grid.cell_volume();
        double total_V = grid.volume();

        // 1. Average compositions: (1/V) * int c_m dV
        diag.average_composition.resize(n_comp, 0.0);
        for (int m = 0; m < n_comp; ++m) {
            double local_mass = 0.0;
            for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
                for (int j = 0; j < grid.ny; ++j) {
                    for (int k = 0; k < grid.nz; ++k) {
                        size_t idx = grid.idx(i_loc, j, k);
                        local_mass += c[m][idx] * dV;
                    }
                }
            }
            diag.average_composition[m] = mpi.allreduce_sum(local_mass) / total_V;
        }

        // 2. Max unity deviation and min/max composition bounds
        double local_max_dev = 0.0;
        double local_min_c = 1e9;
        double local_max_c = -1e9;

        for (int i_loc = 0; i_loc < grid.local_nx; ++i_loc) {
            for (int j = 0; j < grid.ny; ++j) {
                for (int k = 0; k < grid.nz; ++k) {
                    size_t idx = grid.idx(i_loc, j, k);
                    double sum = 0.0;
                    for (int m = 0; m < n_comp; ++m) {
                        double val = c[m][idx];
                        sum += val;
                        local_min_c = std::min(local_min_c, val);
                        local_max_c = std::max(local_max_c, val);
                    }
                    local_max_dev = std::max(local_max_dev, std::abs(sum - 1.0));
                }
            }
        }

        diag.max_unity_deviation = mpi.allreduce_max(local_max_dev);
        diag.min_composition = mpi.allreduce_min(local_min_c);
        diag.max_composition = mpi.allreduce_max(local_max_c);

        // 3. Free energy
        diag.total_free_energy = compute_total_free_energy();
        if (current_step > 0 && dt > 0.0) {
            diag.dF_dt = (diag.total_free_energy - prev_energy) / dt;
        } else {
            diag.dF_dt = 0.0;
        }
        prev_energy = diag.total_free_energy;

        return diag;
    }

    // Estimate conservative stable time step:
    // - For explicit methods: 4th-order biharmonic CFL (dt ~ h^4 / (M * kappa))
    // - For semi-implicit methods: relaxed 2nd-order diffusion CFL (dt ~ h^2 / (dim * M))
    double estimate_stable_dt(TimeIntegratorType method = TimeIntegratorType::RK2, double safety_factor = 0.4) const {
        double min_h = grid.dx;
        min_h = std::min(min_h, grid.dy);
        if (grid.dim == 3) {
            min_h = std::min(min_h, grid.dz);
        }

        double M_val = 1.0;
        if (mobility->type() == MobilityType::CONSTANT) {
            M_val = static_cast<ConstantMobility*>(mobility.get())->value();
        } else if (mobility->type() == MobilityType::DEGENERATE) {
            M_val = static_cast<DegenerateMobility*>(mobility.get())->base_value();
        } else {
            std::vector<double> dummy_c(n_comp, 1.0 / n_comp);
            std::vector<double> M_mat(n_comp * n_comp);
            mobility->evaluate_matrix(dummy_c, M_mat);
            double max_m = 1e-6;
            for (double v : M_mat) max_m = std::max(max_m, std::abs(v));
            M_val = max_m;
        }
        if (M_val <= 0.0) M_val = 1.0;

        if (method == TimeIntegratorType::SEMI_IMPLICIT || method == TimeIntegratorType::SEMI_IMPLICIT_BDF2) {
            return safety_factor * (min_h * min_h) / (2.0 * grid.dim * M_val);
        }

        double max_kappa = 1e-6;
        for (int m = 0; m < n_comp; ++m) {
            for (int n = 0; n < n_comp; ++n) {
                max_kappa = std::max(max_kappa, std::abs(kappa[m][n]));
            }
        }

        double h4 = min_h * min_h * min_h * min_h;
        double lambda_max = (4.0 * grid.dim) * (4.0 * grid.dim) / h4;
        return safety_factor * 2.0 / (lambda_max * M_val * max_kappa);
    }
};

} // namespace mcch
