/**
 * @file Grand_Potential.cpp
 * @brief SYCL Grand Potential Phase Field - Multi-Pass Architecture (AMReX-style)
 * 
 * KEY CHANGE: Uses separate kernel passes like AMReX:
 *   Pass 1: Compute gradient energy term (dAdphi, divdAdgradphi) -> term1
 *   Pass 2: Compute double-well potential (dwdphi) -> term2
 *   Pass 3: Compute thermodynamic driving force (dpsi) -> term3
 *   Pass 4: Combine terms and update phi
 *   Pass 5: Diffusion evolution
 * 
 * Each pass has a .wait() barrier ensuring all cells complete before next pass.
 * This matches AMReX's multi-kernel approach that is proven to work with AlZn.
 */

#include "physics/PF_model.hpp"
#include "physics/Anisotropy.hpp"
#include "utilities/helper_functions.hpp"
#include "thermo/Function_F.hpp"
#include "thermo/Function_F4.hpp"

#include <cmath>
#include <iostream>

namespace {
    constexpr double L_PI = 3.14159265358979323846;
    constexpr double INV_PI_SQ = 1.0 / (L_PI * L_PI);
    constexpr double EPS_SMALL = 1.0e-12;
    constexpr double COMPOSITION_MIN = 0.0;
    constexpr double COMPOSITION_MAX = 1.0;
    constexpr double INTERFACE_THRESHOLD = 1.0e-10;
}

// ============================================================================
// Interpolation functions
// ============================================================================
SYCL_EXTERNAL inline double hphi(const double* phi, int a, int nP) {
    double sum1 = 0.0;
    double sum = 3.0 * phi[a] * phi[a] - 2.0 * phi[a] * phi[a] * phi[a];
    for (int b = 0; b < nP; b++) {
        for (int c = 0; c < nP; c++) {
            if (b != a && c != a && b < c) {
                sum1 += phi[b] * phi[c];
            }
        }
    }
    sum1 *= 2.0 * phi[a];
    return sum + sum1;
}

SYCL_EXTERNAL inline double dhphi(const double* phi, int b, int a, int nP) {
    double sum = 0.0;
    if (a == b) {
        sum = 6.0 * phi[a] * (1.0 - phi[a]);
        for (int e = 0; e < nP; e++) {
            for (int f = 0; f < nP; f++) {
                if (e != b && f != b && e < f) {
                    sum += 2.0 * phi[e] * phi[f];
                }
            }
        }
    } else {
        for (int e = 0; e < nP; e++) {
            if (e != b && e != a) {
                sum += 2.0 * phi[e];
            }
        }
        sum *= phi[b];
    }
    return sum;
}

SYCL_EXTERNAL inline double FunctionTau(const double* phi, double tau_default, int nP) {
    double sum = 0.0, sum1 = 0.0;
    for (int a = 0; a < nP; a++) {
        for (int b = 0; b < nP; b++) {
            if (a < b) {
                sum += tau_default * phi[a] * phi[b];
                sum1 += phi[a] * phi[b];
            }
        }
    }
    return (sum1 > 0.0) ? sum / sum1 : tau_default;
}

// ============================================================================
// Interface detection with threshold
// ============================================================================
SYCL_EXTERNAL inline bool is_at_interface(double divphi) {
    return sycl::fabs(divphi) > INTERFACE_THRESHOLD;
}

// ============================================================================
// Simplex projection
// ============================================================================
SYCL_EXTERNAL inline void projection_on_simplex(
    const double* phi, double* deltaphi, const double* divphi, int nP) 
{
    double Deltaphi = 0.0;
    int count_phases = 0;
    double sum_phib = 0.0;
    
    // Count active phases
    for (int a = 0; a < nP; a++) {
        double phi_new = phi[a] + deltaphi[a];
        if (is_at_interface(divphi[a]) && (phi_new > 0.0) && (phi_new < 1.0)) {
            count_phases++;
            sum_phib += phi[a];
        }
    }
    
    // Fix negative values
    for (int a = 0; a < nP; a++) {
        double phi_new = phi[a] + deltaphi[a];
        if (phi_new < 0.0) {
            double Deltaphi_alpha = sycl::fabs(deltaphi[a] + phi[a]);
            deltaphi[a] += Deltaphi_alpha;
            Deltaphi += Deltaphi_alpha;
        }
    }
    
    // Redistribute using phi-weighted average
    if (Deltaphi > 0.0) {
        for (int b = 0; b < nP; b++) {
            double phi_new = phi[b] + deltaphi[b];
            if (is_at_interface(divphi[b]) && (phi_new > 0.0) && (phi_new < 1.0)) {
                if (sycl::fabs(sum_phib) > EPS_SMALL) {
                    deltaphi[b] -= Deltaphi * phi[b] / sum_phib;
                } else if (count_phases > 0) {
                    deltaphi[b] -= Deltaphi / count_phases;
                }
            }
        }
    }
    
    // Fix values > 1
    for (int a = 0; a < nP; a++) {
        if (is_at_interface(divphi[a])) {
            if ((phi[a] + deltaphi[a]) > 1.0) {
                deltaphi[a] = 1.0 - phi[a];
                for (int b = 0; b < nP; b++) {
                    if (is_at_interface(divphi[b]) && (b != a)) {
                        deltaphi[b] = -phi[b];
                    }
                }
                break;
            }
        }
    }
}

// ============================================================================
// Constructor
// ============================================================================
microsim::Solver::Solver(microsim::Info *global_fd_, mpi::mpi_comm *comm_, sycl::queue &q_)
    : global_fd(global_fd_), comm(comm_), q(q_) {
    nPhase = global_fd->parameters->NUM_PHASES;
    nComp = global_fd->parameters->NUM_COMPONENTS;
    current_step = 0;
    
    // Allocate intermediate term buffers (like AMReX MultiFabs)
    workers *wk = global_fd->workers_mpi;
    size_t total_cells = wk->ghost_d * wk->ghost_h * wk->ghost_w;
    
    term1_buffer = sycl::malloc_device<double>(total_cells * nPhase, q);
    term2_buffer = sycl::malloc_device<double>(total_cells * nPhase, q);
    term3_buffer = sycl::malloc_device<double>(total_cells * nPhase, q);
    divphi_buffer = sycl::malloc_device<double>(total_cells * nPhase, q);
    phase_comp_buffer = sycl::malloc_device<double>(total_cells * nPhase * (nComp - 1), q);
    
    std::cout << "\n=== Grand Potential SYCL (Multi-Pass AMReX-style) ===" << std::endl;
    std::cout << "Phases: " << nPhase << ", Components: " << nComp << std::endl;
    std::cout << "Nsmooth: " << global_fd->parameters->Nsmooth << std::endl;
    std::cout << "Allocated intermediate buffers for " << total_cells << " cells" << std::endl;
}

microsim::Solver::~Solver() {
    if (term1_buffer) sycl::free(term1_buffer, q);
    if (term2_buffer) sycl::free(term2_buffer, q);
    if (term3_buffer) sycl::free(term3_buffer, q);
    if (divphi_buffer) sycl::free(divphi_buffer, q);
    if (phase_comp_buffer) sycl::free(phase_comp_buffer, q);
}

// ============================================================================
// PASS 0: Compute divergence and pre-compute phase compositions
// ============================================================================
template <int Dim>
void kernel_pass0_divphi_and_phase_comp(
    fields *gf, workers *wk,
    double *divphi_buf, double *phase_comp_buf,
    const microsim::ThermoDynamics &thermo,
    double dx, double dy, double T, int nP, int nV,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            int front  = idx(z, y, x+1, wk->ghost_h, wk->ghost_w);
            int back   = idx(z, y, x-1, wk->ghost_h, wk->ghost_w);
            int right  = idx(z, y+1, x, wk->ghost_h, wk->ghost_w);
            int left   = idx(z, y-1, x, wk->ghost_h, wk->ghost_w);
            
            // Compute divergence (Laplacian) of phi - same as AMReX
            for (int a = 0; a < nP; a++) {
                double div = (gf[front].phia[a] - 2.0*gf[center].phia[a] + gf[back].phia[a]) / (dx*dx)
                           + (gf[right].phia[a] - 2.0*gf[center].phia[a] + gf[left].phia[a]) / (dy*dy);
                divphi_buf[center * nP + a] = div;
            }
            
            // Check bulk vs interface
            bool is_interface = true;
            int bulk_phase = -1;
            for (int a = 0; a < nP; a++) {
                if (gf[center].phia[a] == 1.0) {
                    bulk_phase = a;
                    is_interface = false;
                    break;
                }
            }
            
            // Pre-compute phase compositions (like AMReX c_mu call)
            for (int a = 0; a < nP; a++) {
                double c_temp[MAX_COMP - 1];
                
                if (is_interface || a != bulk_phase) {
                    thermo.Function_C_mu(gf[center].mu, c_temp, T, a);
                } else {
                    // Bulk phase uses actual composition
                    for (int k = 0; k < nV; k++) {
                        c_temp[k] = gf[center].C[k];
                    }
                }
                
                // Store with clamping
                for (int k = 0; k < nV; k++) {
                    double c_clamped = sycl::fmax(COMPOSITION_MIN, sycl::fmin(COMPOSITION_MAX, c_temp[k]));
                    phase_comp_buf[center * nP * nV + a * nV + k] = c_clamped;
                }
            }
        });
    }).wait();  // BARRIER - all cells must complete before next pass
}

// ============================================================================
// PASS 1: Compute gradient energy term (dAdphi + divdAdgradphi) -> term1
// ============================================================================
template <int Dim>
void kernel_pass1_gradient_energy(
    fields *gf, workers *wk,
    double *term1_buf, double *gamma,
    double dx, double dy, int nP,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            int front  = idx(z, y, x+1, wk->ghost_h, wk->ghost_w);
            int back   = idx(z, y, x-1, wk->ghost_h, wk->ghost_w);
            int right  = idx(z, y+1, x, wk->ghost_h, wk->ghost_w);
            int left   = idx(z, y-1, x, wk->ghost_h, wk->ghost_w);
            
            // Compute staggered gradients (AMReX style)
            double gradphi_X[MAX_PHASES], gradphi_Y[MAX_PHASES];
            double phistagg_X[MAX_PHASES], phistagg_Y[MAX_PHASES];
            double gradphi_X_back[MAX_PHASES], gradphi_Y_left[MAX_PHASES];
            double phistagg_X_back[MAX_PHASES], phistagg_Y_left[MAX_PHASES];
            
            for (int a = 0; a < nP; a++) {
                // Forward gradients at i+1/2, j+1/2
                gradphi_X[a] = (gf[front].phia[a] - gf[center].phia[a]) / dx;
                gradphi_Y[a] = (gf[right].phia[a] - gf[center].phia[a]) / dy;
                phistagg_X[a] = 0.5 * (gf[front].phia[a] + gf[center].phia[a]);
                phistagg_Y[a] = 0.5 * (gf[right].phia[a] + gf[center].phia[a]);
                
                // Backward gradients at i-1/2, j-1/2
                gradphi_X_back[a] = (gf[center].phia[a] - gf[back].phia[a]) / dx;
                gradphi_Y_left[a] = (gf[center].phia[a] - gf[left].phia[a]) / dy;
                phistagg_X_back[a] = 0.5 * (gf[center].phia[a] + gf[back].phia[a]);
                phistagg_Y_left[a] = 0.5 * (gf[center].phia[a] + gf[left].phia[a]);
            }
            
            // Central gradients (AMReX uses central difference)
            double gradphi_c_X[MAX_PHASES], gradphi_c_Y[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                gradphi_c_X[a] = (gf[front].phia[a] - gf[back].phia[a]) / (2.0 * dx);
                gradphi_c_Y[a] = (gf[right].phia[a] - gf[left].phia[a]) / (2.0 * dy);
            }
            
            // Compute dAdphi and divdAdgradphi for each phase
            for (int a = 0; a < nP; a++) {
                double dAdphi = 0.0;
                double divdAdgradphi = 0.0;
                
                for (int b = 0; b < nP; b++) {
                    if (b != a) {
                        double gab = gamma[a * nP + b];
                        
                        // qab at center (AMReX Function_Q)
                        double qab_c_x = gf[center].phia[a] * gradphi_c_X[b] - gf[center].phia[b] * gradphi_c_X[a];
                        double qab_c_y = gf[center].phia[a] * gradphi_c_Y[b] - gf[center].phia[b] * gradphi_c_Y[a];
                        
                        // dAdphi = 2 * gamma * (qab · grad_phi_b)
                        dAdphi += 2.0 * gab * (qab_c_x * gradphi_c_X[b] + qab_c_y * gradphi_c_Y[b]);
                        
                        // qab at staggered positions for divergence
                        double qab_iph_x = phistagg_X[a] * gradphi_X[b] - phistagg_X[b] * gradphi_X[a];
                        double qab_imh_x = phistagg_X_back[a] * gradphi_X_back[b] - phistagg_X_back[b] * gradphi_X_back[a];
                        double qab_jph_y = phistagg_Y[a] * gradphi_Y[b] - phistagg_Y[b] * gradphi_Y[a];
                        double qab_jmh_y = phistagg_Y_left[a] * gradphi_Y_left[b] - phistagg_Y_left[b] * gradphi_Y_left[a];
                        
                        // Flux = qab * (-phi_b_stag)
                        double flux_iph = qab_iph_x * (-phistagg_X[b]);
                        double flux_imh = qab_imh_x * (-phistagg_X_back[b]);
                        double flux_jph = qab_jph_y * (-phistagg_Y[b]);
                        double flux_jmh = qab_jmh_y * (-phistagg_Y_left[b]);
                        
                        // Divergence of flux
                        divdAdgradphi += 2.0 * gab * ((flux_iph - flux_imh) / dx + (flux_jph - flux_jmh) / dy);
                    }
                }
                
                // Store term1 = -dAdphi + divdAdgradphi (AMReX convention)
                term1_buf[center * nP + a] = -dAdphi + divdAdgradphi;
            }
        });
    }).wait();  // BARRIER
}

// ============================================================================
// PASS 2: Compute double-well potential (dwdphi) -> term2
// ============================================================================
template <int Dim>
void kernel_pass2_double_well(
    fields *gf, workers *wk,
    double *term2_buf, double *divphi_buf,
    double *gamma, double *gamma_abc, int nP,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            
            for (int a = 0; a < nP; a++) {
                double sum = 0.0;
                
                // Binary interaction term (AMReX style with sign check)
                for (int b = 0; b < nP; b++) {
                    if (b != a) {
                        double divphi_b = divphi_buf[center * nP + b];
                        if (sycl::fabs(divphi_b) > INTERFACE_THRESHOLD) {
                            double phi_a = gf[center].phia[a];
                            double phi_b = gf[center].phia[b];
                            
                            // Sign handling (AMReX function_W_01)
                            if (phi_a * phi_b >= 0.0) {
                                sum += gamma[a * nP + b] * phi_b;
                            } else {
                                sum -= gamma[a * nP + b] * phi_b;
                            }
                        }
                    }
                }
                sum *= 16.0 * INV_PI_SQ;
                
                // Triple junction term (for nP > 2)
                if (nP > 2) {
                    for (int b = 0; b < nP; b++) {
                        for (int cc = b + 1; cc < nP; cc++) {
                            if (b != a && cc != a) {
                                double divphi_b = divphi_buf[center * nP + b];
                                double divphi_c = divphi_buf[center * nP + cc];
                                
                                if (sycl::fabs(divphi_b) > INTERFACE_THRESHOLD && 
                                    sycl::fabs(divphi_c) > INTERFACE_THRESHOLD) {
                                    sum += gamma_abc[a * nP * nP + b * nP + cc] * 
                                           gf[center].phia[b] * gf[center].phia[cc];
                                }
                            }
                        }
                    }
                }
                
                term2_buf[center * nP + a] = sum;
            }
        });
    }).wait();  // BARRIER
}

// ============================================================================
// PASS 3: Compute thermodynamic driving force (dpsi) -> term3
// ============================================================================
template <int Dim>
void kernel_pass3_dpsi(
    fields *gf, workers *wk,
    double *term3_buf, double *phase_comp_buf,
    const microsim::ThermoDynamics &thermo,
    double T, int nP, int nV,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            
            double psi[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                // Get pre-computed phase composition
                double c_phase[MAX_COMP - 1];
                for (int k = 0; k < nV; k++) {
                    c_phase[k] = phase_comp_buf[center * nP * nV + a * nV + k];
                }
                
                // psi_a = f_a(c_a) - mu · c_a
                double fe = thermo.Free_energy(c_phase, T, a);
                double mu_dot_c = 0.0;
                for (int k = 0; k < nV; k++) {
                    mu_dot_c += gf[center].mu[k] * c_phase[k];
                }
                psi[a] = fe - mu_dot_c;
            }
            
            // Compute dpsi/dphi_a = sum_b psi_b * dhphi(b,a)
            for (int a = 0; a < nP; a++) {
                double sum = 0.0;
                for (int b = 0; b < nP; b++) {
                    sum += psi[b] * dhphi(gf[center].phia, b, a, nP);
                }
                term3_buf[center * nP + a] = sum;
            }
        });
    }).wait();  // BARRIER
}

// ============================================================================
// PASS 4: Combine terms and update phi
// ============================================================================
template <int Dim>
void kernel_pass4_update_phi(
    fields *gf, fields *gfn, workers *wk,
    double *term1_buf, double *term2_buf, double *term3_buf, double *divphi_buf,
    double eps, double Vm, double tau0, double dt, int nP,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            
            // Read stored divphi
            double divphi[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                divphi[a] = divphi_buf[center * nP + a];
            }
            
            // Compute lambda_phi
            double lambda_phi[MAX_PHASES];
            double sum_lambda = 0.0;
            double active_phases = 0.0;
            
            for (int a = 0; a < nP; a++) {
                if (sycl::fabs(divphi[a]) > INTERFACE_THRESHOLD) {
                    // lambda = eps * term1 - term2/eps - term3/Vm
                    lambda_phi[a] = eps * term1_buf[center * nP + a]
                                  - term2_buf[center * nP + a] / eps
                                  - term3_buf[center * nP + a] / Vm;
                    sum_lambda += lambda_phi[a];
                    active_phases += 1.0;
                } else {
                    lambda_phi[a] = 0.0;
                }
            }
            
            if (active_phases > 0.0) {
                sum_lambda /= active_phases;
            }
            
            // Compute deltaphi
            double tau_val = FunctionTau(gf[center].phia, tau0, nP);
            double deltaphi[MAX_PHASES];
            
            for (int a = 0; a < nP; a++) {
                if (sycl::fabs(divphi[a]) > INTERFACE_THRESHOLD) {
                    deltaphi[a] = dt * (lambda_phi[a] - sum_lambda) / (tau_val * eps);
                } else {
                    deltaphi[a] = 0.0;
                }
            }
            
            // Simplex projection
            projection_on_simplex(gf[center].phia, deltaphi, divphi, nP);
            
            // Update phi
            for (int a = 0; a < nP; a++) {
                double phi_new = gf[center].phia[a] + deltaphi[a];
                phi_new = sycl::fmax(0.0, sycl::fmin(1.0, phi_new));
                gfn[center].phia[a] = phi_new;
                gfn[center].deltaPhi[a] = deltaphi[a] / dt;
            }
        });
    }).wait();  // BARRIER
}

// ============================================================================
// PASS 5: Diffusion evolution
// ============================================================================
template <int Dim>
void kernel_pass5_diffusion(
    fields *gf, fields *gfn, workers *wk,
    double *phase_comp_buf,
    const microsim::ThermoDynamics &thermo,
    Phasefield_matraix *pm,
    double dx, double dy, double dt, double T, double eps, int nP, int nV,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            int front  = idx(z, y, x+1, wk->ghost_h, wk->ghost_w);
            int back   = idx(z, y, x-1, wk->ghost_h, wk->ghost_w);
            int right  = idx(z, y+1, x, wk->ghost_h, wk->ghost_w);
            int left   = idx(z, y-1, x, wk->ghost_h, wk->ghost_w);
            
            // Check bulk vs interface
            bool is_interface = true;
            int bulk_phase = -1;
            for (int a = 0; a < nP; a++) {
                if (gf[center].phia[a] == 1.0) {
                    bulk_phase = a;
                    is_interface = false;
                    break;
                }
            }
            

            // Compute h values for susceptibility
            double h_center[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                h_center[a] = hphi(gf[center].phia, a, nP);
            }
            
            // phi at staggered positions
            double phi_center[MAX_PHASES], phi_front[MAX_PHASES], phi_back[MAX_PHASES];
            double phi_right[MAX_PHASES], phi_left[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                phi_center[a] = gf[center].phia[a];
                phi_front[a]  = gf[front].phia[a];
                phi_back[a]   = gf[back].phia[a];
                phi_right[a]  = gf[right].phia[a];
                phi_left[a]   = gf[left].phia[a];
            }
            
            // Compute mu gradients
            double gradmu_X[MAX_COMP-1], gradmu_Y[MAX_COMP-1];
            double gradmu_X_back[MAX_COMP-1], gradmu_Y_left[MAX_COMP-1];
            
            for (int k = 0; k < nV; k++) {
                gradmu_X[k] = (gf[front].mu[k] - gf[center].mu[k]) / dx;
                gradmu_Y[k] = (gf[right].mu[k] - gf[center].mu[k]) / dy;
                gradmu_X_back[k] = (gf[center].mu[k] - gf[back].mu[k]) / dx;
                gradmu_Y_left[k] = (gf[center].mu[k] - gf[left].mu[k]) / dy;
            }
            
            // Compute fluxes with phase-weighted mobility
            double flux_X[MAX_COMP-1] = {0}, flux_X_back[MAX_COMP-1] = {0};
            double flux_Y[MAX_COMP-1] = {0}, flux_Y_left[MAX_COMP-1] = {0};
            
            for (int k = 0; k < nV; k++) {
                for (int l = 0; l < nV; l++) {
                    double D_front = 0.0, D_back = 0.0, D_right = 0.0, D_left = 0.0;
                    
                    for (int a = 0; a < nP; a++) {
                        double D_a = pm->Diffusivity[a * nV * nV + k * nV + l];
                        double dcdmu_a = thermo.cmu[thermo.a_ind(a, k, l)];
                        double M_a = D_a * dcdmu_a;
                        
                        // Use phi directly for mobility interpolation (like AMReX/MPI)
                        D_front += 0.5 * (phi_center[a] + phi_front[a]) * M_a;
                        D_back  += 0.5 * (phi_center[a] + phi_back[a]) * M_a;
                        D_right += 0.5 * (phi_center[a] + phi_right[a]) * M_a;
                        D_left  += 0.5 * (phi_center[a] + phi_left[a]) * M_a;
                    }
                    
                    flux_X[k]      += D_front * gradmu_X[l];
                    flux_X_back[k] += D_back * gradmu_X_back[l];
                    flux_Y[k]      += D_right * gradmu_Y[l];
                    flux_Y_left[k] += D_left * gradmu_Y_left[l];
                }
            }
            
            // Divergence of flux * dt
            double divflux[MAX_COMP-1];
            for (int k = 0; k < nV; k++) {
                divflux[k] = ((flux_X[k] - flux_X_back[k]) / dx + 
                              (flux_Y[k] - flux_Y_left[k]) / dy) * dt;
            }
            
            // ================================================================
            // ANTI-TRAPPING CURRENT (Karma-Rappel correction)
            // ================================================================
            // This corrects for solute trapping at the diffuse interface
            // j_at = (1 - D_s/D_l) * |n_s · n_l| * (c_l - c_s) * (dphi/dt) * (grad_phi/|grad_phi|)
            // ================================================================
            
            double divjat[MAX_COMP-1] = {0};
            
            if (is_interface) {
                // Need additional neighbor indices for 4-point averaging
                int front_right = idx(z, y+1, x+1, wk->ghost_h, wk->ghost_w);
                int front_left  = idx(z, y-1, x+1, wk->ghost_h, wk->ghost_w);
                int back_right  = idx(z, y+1, x-1, wk->ghost_h, wk->ghost_w);
                int back_left   = idx(z, y-1, x-1, wk->ghost_h, wk->ghost_w);
                
                // Compute gradients of phi at staggered positions for all phases
                // At i+1/2 (iph): gradphi_x = (phi[i+1] - phi[i])/dx
                //                 gradphi_y = avg of central differences
                double norm_x[MAX_PHASES][5] = {{0}};  // [phase][position: 0=cent, 1=iph, 2=imh, 3=jph, 4=jmh]
                double norm_y[MAX_PHASES][5] = {{0}};
                
                constexpr int CENT = 0, IPH = 1, IMH = 2, JPH = 3, JMH = 4;
                
                for (int a = 0; a < nP; a++) {
                    // At i+1/2 (between center and front)
                    double gx_iph = (gf[front].phia[a] - gf[center].phia[a]) / dx;
                    double gy_iph = (gf[front_right].phia[a] - gf[front_left].phia[a] + 
                                     gf[right].phia[a] - gf[left].phia[a]) / (4.0 * dy);
                    double mod_iph = sycl::sqrt(gx_iph*gx_iph + gy_iph*gy_iph);
                    if (mod_iph > 1e-15) {
                        norm_x[a][IPH] = gx_iph / mod_iph;
                        norm_y[a][IPH] = gy_iph / mod_iph;
                    }
                    
                    // At i-1/2 (between back and center)
                    double gx_imh = (gf[center].phia[a] - gf[back].phia[a]) / dx;
                    double gy_imh = (gf[right].phia[a] - gf[left].phia[a] + 
                                     gf[back_right].phia[a] - gf[back_left].phia[a]) / (4.0 * dy);
                    double mod_imh = sycl::sqrt(gx_imh*gx_imh + gy_imh*gy_imh);
                    if (mod_imh > 1e-15) {
                        norm_x[a][IMH] = gx_imh / mod_imh;
                        norm_y[a][IMH] = gy_imh / mod_imh;
                    }
                    
                    // At j+1/2 (between center and right)
                    double gx_jph = (gf[front_right].phia[a] - gf[back_right].phia[a] + 
                                     gf[front].phia[a] - gf[back].phia[a]) / (4.0 * dx);
                    double gy_jph = (gf[right].phia[a] - gf[center].phia[a]) / dy;
                    double mod_jph = sycl::sqrt(gx_jph*gx_jph + gy_jph*gy_jph);
                    if (mod_jph > 1e-15) {
                        norm_x[a][JPH] = gx_jph / mod_jph;
                        norm_y[a][JPH] = gy_jph / mod_jph;
                    }
                    
                    // At j-1/2 (between left and center)
                    double gx_jmh = (gf[front].phia[a] - gf[back].phia[a] + 
                                     gf[front_left].phia[a] - gf[back_left].phia[a]) / (4.0 * dx);
                    double gy_jmh = (gf[center].phia[a] - gf[left].phia[a]) / dy;
                    double mod_jmh = sycl::sqrt(gx_jmh*gx_jmh + gy_jmh*gy_jmh);
                    if (mod_jmh > 1e-15) {
                        norm_x[a][JMH] = gx_jmh / mod_jmh;
                        norm_y[a][JMH] = gy_jmh / mod_jmh;
                    }
                }
                
                // Anti-trapping flux at each staggered position
                double jat_iph[MAX_COMP-1] = {0}, jat_imh[MAX_COMP-1] = {0};
                double jat_jph[MAX_COMP-1] = {0}, jat_jmh[MAX_COMP-1] = {0};
                
                // Liquid phase index (last phase)
                const int liq = nP - 1;
                
                // Get phase compositions at neighbors for liquid
                double c_liq_center[MAX_COMP-1], c_liq_front[MAX_COMP-1], c_liq_back[MAX_COMP-1];
                double c_liq_right[MAX_COMP-1], c_liq_left[MAX_COMP-1];
                for (int k = 0; k < nV; k++) {
                    c_liq_center[k] = phase_comp_buf[center * nP * nV + liq * nV + k];
                    c_liq_front[k]  = phase_comp_buf[front * nP * nV + liq * nV + k];
                    c_liq_back[k]   = phase_comp_buf[back * nP * nV + liq * nV + k];
                    c_liq_right[k]  = phase_comp_buf[right * nP * nV + liq * nV + k];
                    c_liq_left[k]   = phase_comp_buf[left * nP * nV + liq * nV + k];
                }
                
                // Loop over solid phases (all except liquid)
                for (int a = 0; a < nP - 1; a++) {
                    // Get phase compositions for solid phase a
                    double c_sol_center[MAX_COMP-1], c_sol_front[MAX_COMP-1], c_sol_back[MAX_COMP-1];
                    double c_sol_right[MAX_COMP-1], c_sol_left[MAX_COMP-1];
                    for (int k = 0; k < nV; k++) {
                        c_sol_center[k] = phase_comp_buf[center * nP * nV + a * nV + k];
                        c_sol_front[k]  = phase_comp_buf[front * nP * nV + a * nV + k];
                        c_sol_back[k]   = phase_comp_buf[back * nP * nV + a * nV + k];
                        c_sol_right[k]  = phase_comp_buf[right * nP * nV + a * nV + k];
                        c_sol_left[k]   = phase_comp_buf[left * nP * nV + a * nV + k];
                    }
                    
                    // delta_phi = phi_new - phi_old (from Pass 4)
                    double dphi_center = gfn[center].deltaPhi[a] * dt;
                    double dphi_front  = gfn[front].deltaPhi[a] * dt;
                    double dphi_back   = gfn[back].deltaPhi[a] * dt;
                    double dphi_right  = gfn[right].deltaPhi[a] * dt;
                    double dphi_left   = gfn[left].deltaPhi[a] * dt;
                    
                    // Scalar product: n_solid · n_liquid at each staggered position
                    double scalprod_iph = norm_x[a][IPH] * norm_x[liq][IPH] + norm_y[a][IPH] * norm_y[liq][IPH];
                    double scalprod_imh = norm_x[a][IMH] * norm_x[liq][IMH] + norm_y[a][IMH] * norm_y[liq][IMH];
                    double scalprod_jph = norm_x[a][JPH] * norm_x[liq][JPH] + norm_y[a][JPH] * norm_y[liq][JPH];
                    double scalprod_jmh = norm_x[a][JMH] * norm_x[liq][JMH] + norm_y[a][JMH] * norm_y[liq][JMH];
                    
                    for (int k = 0; k < nV; k++) {
                        // Diffusivity ratio: (1 - D_solid/D_liquid)
                        double D_sol = pm->Diffusivity[a * nV * nV + k * nV + k];
                        double D_liq = pm->Diffusivity[liq * nV * nV + k * nV + k];
                        double diff_ratio = (sycl::fabs(D_liq) > EPS_SMALL) ? (1.0 - D_sol / D_liq) : 0.0;
                        
                        // Anti-trapping coefficient for double-well: 1/(2*sqrt(2))
                        double at_coef = 0.5 / sycl::sqrt(2.0);
                        
                        // jat at i+1/2 (x-direction flux)
                        double dc_front = c_liq_front[k] - c_sol_front[k];
                        double dc_center = c_liq_center[k] - c_sol_center[k];
                        jat_iph[k] += diff_ratio * at_coef * sycl::fabs(scalprod_iph) * 
                                      (dc_front * dphi_front + dc_center * dphi_center) * norm_x[a][IPH];
                        
                        // jat at i-1/2
                        double dc_back = c_liq_back[k] - c_sol_back[k];
                        jat_imh[k] += diff_ratio * at_coef * sycl::fabs(scalprod_imh) * 
                                      (dc_center * dphi_center + dc_back * dphi_back) * norm_x[a][IMH];
                        
                        // jat at j+1/2 (y-direction flux)
                        double dc_right = c_liq_right[k] - c_sol_right[k];
                        jat_jph[k] += diff_ratio * at_coef * sycl::fabs(scalprod_jph) * 
                                      (dc_right * dphi_right + dc_center * dphi_center) * norm_y[a][JPH];
                        
                        // jat at j-1/2
                        double dc_left = c_liq_left[k] - c_sol_left[k];
                        jat_jmh[k] += diff_ratio * at_coef * sycl::fabs(scalprod_jmh) * 
                                      (dc_center * dphi_center + dc_left * dphi_left) * norm_y[a][JMH];
                    }
                }
                
                // Divergence of anti-trapping current
                // For double-well potential: coefficient is -epsilon (not -pi*epsilon/4)
                double at_div_coef = -eps;
                for (int k = 0; k < nV; k++) {
                    divjat[k] = at_div_coef * ((jat_iph[k] - jat_imh[k]) / dx + 
                                               (jat_jph[k] - jat_jmh[k]) / dy);
                }
            }
            
            if (!is_interface) {
                // Bulk: update composition, then compute mu
                // No anti-trapping in bulk (divjat = 0)
                for (int k = 0; k < nV; k++) {
                    double c_new = gf[center].C[k] + divflux[k] - divjat[k];
                    gfn[center].C[k] = sycl::fmax(COMPOSITION_MIN, sycl::fmin(COMPOSITION_MAX, c_new));
                }
                thermo.Function_Mu(gfn[center].C, T, bulk_phase, gfn[center].mu);
            } else {
                // Interface: update mu directly
                
                // Source term from phase change
                double source[MAX_COMP-1] = {0};
                for (int a = 0; a < nP; a++) {
                    double dhdt = 0.0;
                    for (int b = 0; b < nP; b++) {
                        dhdt += dhphi(gf[center].phia, a, b, nP) * gfn[center].deltaPhi[b] * dt;
                    }
                    for (int k = 0; k < nV; k++) {
                        double c_a_k = phase_comp_buf[center * nP * nV + a * nV + k];
                        source[k] += c_a_k * dhdt;
                    }
                }
                
                // Susceptibility matrix chi
                double chi[MAX_COMP-1][MAX_COMP-1] = {{0}};
                for (int a = 0; a < nP; a++) {
                    for (int k = 0; k < nV; k++) {
                        for (int l = 0; l < nV; l++) {
                            chi[k][l] += h_center[a] * thermo.cmu[thermo.a_ind(a, k, l)];
                        }
                    }
                }
                
                // Invert chi and update mu
                // delta_mu = chi^{-1} * (divflux - source - divjat)
                // delta_c = divflux - divjat
                if (nV == 1) {
                    double chi_inv = (sycl::fabs(chi[0][0]) > EPS_SMALL) ? 1.0 / chi[0][0] : 0.0;
                    gfn[center].mu[0] = gf[center].mu[0] + (divflux[0] - source[0] - divjat[0]) * chi_inv;
                    gfn[center].C[0] = gf[center].C[0] + divflux[0] - divjat[0];
                } else if (nV == 2) {
                    double det = chi[0][0] * chi[1][1] - chi[0][1] * chi[1][0];
                    if (sycl::fabs(det) > EPS_SMALL) {
                        double inv_det = 1.0 / det;
                        double chi_inv[2][2] = {
                            {chi[1][1] * inv_det, -chi[0][1] * inv_det},
                            {-chi[1][0] * inv_det, chi[0][0] * inv_det}
                        };
                        for (int k = 0; k < 2; k++) {
                            double dmu = 0.0;
                            for (int l = 0; l < 2; l++) {
                                dmu += chi_inv[k][l] * (divflux[l] - source[l] - divjat[l]);
                            }
                            gfn[center].mu[k] = gf[center].mu[k] + dmu;
                            gfn[center].C[k] = gf[center].C[k] + divflux[k] - divjat[k];
                        }
                    } else {
                        // Singular - diagonal approximation
                        for (int k = 0; k < 2; k++) {
                            if (sycl::fabs(chi[k][k]) > EPS_SMALL) {
                                gfn[center].mu[k] = gf[center].mu[k] + (divflux[k] - source[k] - divjat[k]) / chi[k][k];
                            }
                            gfn[center].C[k] = gf[center].C[k] + divflux[k] - divjat[k];
                        }
                    }
                } else {
                    // General case - diagonal
                    for (int k = 0; k < nV; k++) {
                        double chi_inv = (sycl::fabs(chi[k][k]) > EPS_SMALL) ? 1.0 / chi[k][k] : 0.0;
                        gfn[center].mu[k] = gf[center].mu[k] + (divflux[k] - source[k] - divjat[k]) * chi_inv;
                        gfn[center].C[k] = gf[center].C[k] + divflux[k] - divjat[k];
                    }
                }
                
                // Clamp compositions
                for (int k = 0; k < nV; k++) {
                    gfn[center].C[k] = sycl::fmax(COMPOSITION_MIN, sycl::fmin(COMPOSITION_MAX, gfn[center].C[k]));
                }
            }
        });
    }).wait();  // BARRIER
}

// ============================================================================
// SMOOTHING KERNEL - Simplified physics (no thermodynamics)
// ============================================================================
template <int Dim>
void kernel_smooth(
    fields *gf, fields *gfn, workers *wk,
    double *gamma, double dx, double dy, double dt, double eps, double tau0, int nP,
    sycl::queue &q)
{
    auto range = Coords<Dim>::to_range(wk->ghost_d, wk->ghost_h, wk->ghost_w);
    
    q.submit([&](sycl::handler &h) {
        h.parallel_for(range, [=](sycl::id<Dim> id) {
            auto c = Coords<Dim>::from_id(id);
            if (c.x < 1 || c.x > wk->local_w || c.y < 1 || c.y > wk->local_h) return;
            if constexpr (Dim == 3) { if (c.z < 1 || c.z > wk->local_d) return; }

            const int z = c.z, y = c.y, x = c.x;
            int center = idx(z, y, x, wk->ghost_h, wk->ghost_w);
            int front  = idx(z, y, x+1, wk->ghost_h, wk->ghost_w);
            int back   = idx(z, y, x-1, wk->ghost_h, wk->ghost_w);
            int right  = idx(z, y+1, x, wk->ghost_h, wk->ghost_w);
            int left   = idx(z, y-1, x, wk->ghost_h, wk->ghost_w);
            
            // Compute divphi (Laplacian)
            double divphi[MAX_PHASES];
            for (int a = 0; a < nP; a++) {
                divphi[a] = (gf[front].phia[a] - 2.0*gf[center].phia[a] + gf[back].phia[a]) / (dx*dx)
                          + (gf[right].phia[a] - 2.0*gf[center].phia[a] + gf[left].phia[a]) / (dy*dy);
            }
            
            // Simplified gradient energy (isotropic)
            double gradphi_c_X[MAX_PHASES], gradphi_c_Y[MAX_PHASES];
            double gradphi_X[MAX_PHASES], gradphi_Y[MAX_PHASES];
            double gradphi_X_back[MAX_PHASES], gradphi_Y_left[MAX_PHASES];
            double phistagg_X[MAX_PHASES], phistagg_Y[MAX_PHASES];
            double phistagg_X_back[MAX_PHASES], phistagg_Y_left[MAX_PHASES];
            
            for (int a = 0; a < nP; a++) {
                gradphi_c_X[a] = (gf[front].phia[a] - gf[back].phia[a]) / (2.0 * dx);
                gradphi_c_Y[a] = (gf[right].phia[a] - gf[left].phia[a]) / (2.0 * dy);
                gradphi_X[a] = (gf[front].phia[a] - gf[center].phia[a]) / dx;
                gradphi_Y[a] = (gf[right].phia[a] - gf[center].phia[a]) / dy;
                gradphi_X_back[a] = (gf[center].phia[a] - gf[back].phia[a]) / dx;
                gradphi_Y_left[a] = (gf[center].phia[a] - gf[left].phia[a]) / dy;
                phistagg_X[a] = 0.5 * (gf[front].phia[a] + gf[center].phia[a]);
                phistagg_Y[a] = 0.5 * (gf[right].phia[a] + gf[center].phia[a]);
                phistagg_X_back[a] = 0.5 * (gf[center].phia[a] + gf[back].phia[a]);
                phistagg_Y_left[a] = 0.5 * (gf[center].phia[a] + gf[left].phia[a]);
            }
            
            // Compute terms
            double term1[MAX_PHASES], term2[MAX_PHASES];
            
            for (int a = 0; a < nP; a++) {
                double dAdphi = 0.0, divdAdgradphi = 0.0, dwdphi = 0.0;
                
                for (int b = 0; b < nP; b++) {
                    if (b != a) {
                        double gab = gamma[a * nP + b];
                        
                        // dAdphi
                        double qab_x = gf[center].phia[a] * gradphi_c_X[b] - gf[center].phia[b] * gradphi_c_X[a];
                        double qab_y = gf[center].phia[a] * gradphi_c_Y[b] - gf[center].phia[b] * gradphi_c_Y[a];
                        dAdphi += 2.0 * gab * (qab_x * gradphi_c_X[b] + qab_y * gradphi_c_Y[b]);
                        
                        // divdAdgradphi
                        double qab_iph = phistagg_X[a] * gradphi_X[b] - phistagg_X[b] * gradphi_X[a];
                        double qab_imh = phistagg_X_back[a] * gradphi_X_back[b] - phistagg_X_back[b] * gradphi_X_back[a];
                        double qab_jph = phistagg_Y[a] * gradphi_Y[b] - phistagg_Y[b] * gradphi_Y[a];
                        double qab_jmh = phistagg_Y_left[a] * gradphi_Y_left[b] - phistagg_Y_left[b] * gradphi_Y_left[a];
                        
                        double flux_iph = qab_iph * (-phistagg_X[b]);
                        double flux_imh = qab_imh * (-phistagg_X_back[b]);
                        double flux_jph = qab_jph * (-phistagg_Y[b]);
                        double flux_jmh = qab_jmh * (-phistagg_Y_left[b]);
                        
                        divdAdgradphi += 2.0 * gab * ((flux_iph - flux_imh) / dx + (flux_jph - flux_jmh) / dy);
                        
                        // dwdphi (simplified - no sign check for smoothing)
                        if (sycl::fabs(divphi[b]) > INTERFACE_THRESHOLD) {
                            dwdphi += gab * gf[center].phia[b];
                        }
                    }
                }
                dwdphi *= 16.0 * INV_PI_SQ;
                
                term1[a] = -dAdphi + divdAdgradphi;
                term2[a] = dwdphi;
            }
            
            // Compute lambda and update
            double lambda_phi[MAX_PHASES];
            double sum_lambda = 0.0, active_phases = 0.0;
            
            for (int a = 0; a < nP; a++) {
                if (sycl::fabs(divphi[a]) > INTERFACE_THRESHOLD) {
                    lambda_phi[a] = eps * term1[a] - term2[a] / eps;  // NO dpsi!
                    sum_lambda += lambda_phi[a];
                    active_phases += 1.0;
                } else {
                    lambda_phi[a] = 0.0;
                }
            }
            if (active_phases > 0.0) sum_lambda /= active_phases;
            
            double tau_val = FunctionTau(gf[center].phia, tau0, nP);
            double deltaphi[MAX_PHASES];
            
            for (int a = 0; a < nP; a++) {
                if (sycl::fabs(divphi[a]) > INTERFACE_THRESHOLD) {
                    deltaphi[a] = dt * (lambda_phi[a] - sum_lambda) / (tau_val * eps);
                } else {
                    deltaphi[a] = 0.0;
                }
            }
            
            projection_on_simplex(gf[center].phia, deltaphi, divphi, nP);
            
            for (int a = 0; a < nP; a++) {
                double phi_new = gf[center].phia[a] + deltaphi[a];
                gfn[center].phia[a] = sycl::fmax(0.0, sycl::fmin(1.0, phi_new));
            }
        });
    }).wait();
}

// ============================================================================
// Main Evolution - Multi-Pass Architecture
// ============================================================================
template <int Dim>
void microsim::Solver::evolve_system() {
    OFFSETS *off = global_fd->offsets;
    fields *gf = global_fd->gridfields;
    fields *gfn = global_fd->gridfields_new;
    workers *wk = global_fd->workers_mpi;
    Input_parameters *inp = global_fd->parameters;
    Domain *dm = global_fd->domainInfo;
    Phasefield_matraix *pm = global_fd->PFMAT;

    const int nP = this->nPhase;
    const int nC = this->nComp;
    const int nV = nC - 1;
    const double dx = dm->DELTA_X, dy = dm->DELTA_Y;
    const double dt = dm->DELTA_T;
    const double eps = inp->epsilon;
    const double Vm = inp->Volm;
    const double tau0 = inp->tau;
    const double T = inp->T;

    microsim::ThermoDynamics thermo(global_fd, true);

    // ========================================================================
    // SMOOTHING LOOP - First timestep only
    // ========================================================================
    if (current_step == 0) {
        int nsmooth = global_fd->parameters->Nsmooth;
        if (nsmooth > 1) {
            std::cout << "Running " << nsmooth << " smoothing iterations..." << std::endl;
            
            for (int t = 1; t < nsmooth; t++) {
                comm->exchange_boundary<Dim>({{off->OFF_PHI, nPhase}});
                
                kernel_smooth<Dim>(gf, gfn, wk,
                    pm->gamma, dx, dy, dt, eps, tau0, nP, q);
                
                std::swap(global_fd->gridfields, global_fd->gridfields_new);
                gf = global_fd->gridfields;
                gfn = global_fd->gridfields_new;
            }
            std::cout << "Smoothing complete." << std::endl;
        }
    }
    current_step++;

    // ========================================================================
    // MAIN EVOLUTION - Multi-Pass as in AMREX
    // ========================================================================
    
    // Exchange boundaries first (like AMReX FillBoundary)
    comm->exchange_boundary<Dim>({{off->OFF_PHI, nPhase},
                                  {off->OFF_MU, nComp - 1},
                                  {off->OFF_COMP, nComp - 1}});

    // PASS 0: Compute divphi and pre-compute phase compositions
    kernel_pass0_divphi_and_phase_comp<Dim>(
        gf, wk, divphi_buffer, phase_comp_buffer,
        thermo, dx, dy, T, nP, nV, q);

    // PASS 1: Compute gradient energy term -> term1
    kernel_pass1_gradient_energy<Dim>(
        gf, wk, term1_buffer, pm->gamma, dx, dy, nP, q);

    // PASS 2: Compute double-well potential -> term2
    kernel_pass2_double_well<Dim>(
        gf, wk, term2_buffer, divphi_buffer,
        pm->gamma, pm->gamma_abc, nP, q);

    // PASS 3: Compute thermodynamic driving force -> term3
    kernel_pass3_dpsi<Dim>(
        gf, wk, term3_buffer, phase_comp_buffer,
        thermo, T, nP, nV, q);

    // PASS 4: Combine terms and update phi
    kernel_pass4_update_phi<Dim>(
        gf, gfn, wk,
        term1_buffer, term2_buffer, term3_buffer, divphi_buffer,
        eps, Vm, tau0, dt, nP, q);

    // Exchange phi after update
    comm->exchange_boundary<Dim>({{off->OFF_PHI, nPhase}});

    // PASS 5: Diffusion evolution
    kernel_pass5_diffusion<Dim>(
        gf, gfn, wk, phase_comp_buffer,
        thermo, pm, dx, dy, dt, T, eps, nP, nV, q);

    // Exchange mu and comp
    comm->exchange_boundary<Dim>({{off->OFF_MU, nComp - 1},
                                  {off->OFF_COMP, nComp - 1}});

    // Swap buffers
    std::swap(global_fd->grid_data, global_fd->grid_data_new);
    std::swap(global_fd->gridfields, global_fd->gridfields_new);
}

// Legacy wrappers (kept for compatibility)
template <int Dim>
void microsim::Solver::solve_Phasefield() {
    // Now handled in multi-pass evolve_system
}

template <int Dim>
void microsim::Solver::solve_diffusion() {
    // Now handled in multi-pass evolve_system
}

template <int Dim>
void microsim::Solver::solve_flux() {}

// Template instantiations
template void microsim::Solver::evolve_system<2>();
template void microsim::Solver::evolve_system<3>();
template void microsim::Solver::solve_Phasefield<2>();
template void microsim::Solver::solve_Phasefield<3>();
template void microsim::Solver::solve_diffusion<2>();
template void microsim::Solver::solve_diffusion<3>();
template void microsim::Solver::solve_flux<2>();
template void microsim::Solver::solve_flux<3>();
