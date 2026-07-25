/**
 * @file Function_F_GP.hpp
 * @brief Thermodynamic Functions for Grand Potential Phase-Field Model
 * 
 * Implementation based on:
 * Choudhury, A. "Quantitative Phase-Field Model for Phase Transformations 
 * in Multi-Component Alloys" (PhD Thesis, KIT 2012)
 * 
 * Key functions:
 * - Free energy f_α(c_α, T)
 * - Chemical potential μ = ∂f/∂c
 * - Phase composition c_α(μ, T) from μ = ∂f_α/∂c_α
 * - Grand potential ψ_α = f_α - μ·c_α
 * - Susceptibility χ_α = ∂c_α/∂μ
 * - Interpolation functions h_α(φ), ∂h_α/∂φ_β
 */

#pragma once
#include <sycl/sycl.hpp>
#include "utilities/helper_functions.hpp"
#include "utilities/Input_parameters.hpp"
#include "utilities/microsim.hpp"

namespace microsim {

struct ThermoDynamics {
    // Data members (Copies of pointers and values)
    int nPhase;
    int nComp;
    double T_eq;
    double T;
    bool isothermal;

    double *A;      // Hessian matrix A[a][i][j] for phase a
    double *B;      // Linear term B[a][i]
    double *C;      // Constant term C[a]
    double *Beq;    // B at equilibrium temperature
    double *dBbdT;  // Temperature derivative of B
    double *dcbdT;  // Temperature derivative of equilibrium composition
    double *ceq;    // Equilibrium compositions
    double *DELTA_C;
    double *DELTA_T;
    double *cmu;    // ∂c/∂μ matrix (inverse of Hessian)
    double *muc;    // Hessian ∂μ/∂c = 2A (diagonal) or A (off-diagonal)
    double *slopes;
    double *c_guess;

    // Constructor
    ThermoDynamics(microsim::Info *global_fd, bool isothermal_) {
        nPhase = global_fd->parameters->NUM_PHASES;
        nComp = global_fd->parameters->NUM_COMPONENTS;
        T_eq = global_fd->parameters->T_eq;
        T = global_fd->parameters->T;
        isothermal = isothermal_;

        Phasefield_matraix *pfm = global_fd->PFMAT;
        A = pfm->A;
        B = pfm->B;
        C = pfm->C;
        Beq = pfm->Beq;
        dBbdT = pfm->dBbdT;
        dcbdT = pfm->dcbdT;
        ceq = pfm->ceq;
        DELTA_C = pfm->DELTA_C;
        DELTA_T = pfm->DELTA_T;
        cmu = pfm->cmu;
        muc = pfm->muc;
        slopes = pfm->slopes;
        c_guess = pfm->c_guess;
    }

    // ========================================================================
    // Index helper functions
    // ========================================================================
    
    /// Index for A matrix: A[phase][comp1][comp2]
    SYCL_EXTERNAL __always_inline int a_ind(int p, int c1, int c2) const {
        return p * (nComp - 1) * (nComp - 1) + c1 * (nComp - 1) + c2;
    }
    
    /// Index for equilibrium composition: ceq[phase1][phase2][component]
    SYCL_EXTERNAL __always_inline int ceq_index(int p1, int p2, int a) const {
        return p1 * nPhase * (nComp - 1) + p2 * (nComp - 1) + a;
    }
    
    /// Index for B and DELTA_C: [phase][component]
    SYCL_EXTERNAL __always_inline int delC_ind(int p, int c) const {
        return p * (nComp - 1) + c;
    }
    
    /// Index for DELTA_T: [phase1][phase2]
    SYCL_EXTERNAL __always_inline int delT_ind(int p1, int p2) const {
        return p1 * nPhase + p2;
    }

    // ========================================================================
    // Core thermodynamic functions
    // ========================================================================
    
    /**
     * @brief Free energy density f_α(c, T) for phase α
     * f_α = Σ_{i≤j} A_α^{ij} c_i c_j + Σ_i B_α^i c_i + C_α
     * (Parabolic approximation, Eq. 3.3 in thesis)
     */
    SYCL_EXTERNAL __always_inline double Free_energy(const double *c, double T, long a) const;
    
    /**
     * @brief Temperature-dependent B coefficient
     */
    SYCL_EXTERNAL __always_inline double functionB(double T, long i, long a) const;
    
    /**
     * @brief Temperature-dependent C coefficient
     */
    SYCL_EXTERNAL __always_inline double functionC(double T, long a) const;
    
    /**
     * @brief Chemical potential μ = ∂f_α/∂c from composition
     * μ_i = 2 A_α^{ii} c_i + Σ_{j≠i} A_α^{ij} c_j + B_α^i
     */
    SYCL_EXTERNAL __always_inline void Function_Mu(double *c, double T, long a, double *Mu) const;
    
    /**
     * @brief Phase composition c_α from chemical potential μ
     * c_α = (∂²f_α/∂c²)^{-1} · (μ - B_α)
     * This is the key GP relation for computing phase concentrations
     */
    SYCL_EXTERNAL __always_inline void Function_C_mu(double *mu, double *c, double T, long a) const;
    
    /**
     * @brief Susceptibility ∂c_α/∂μ (stored in cmu matrix)
     * For parabolic free energy: ∂c/∂μ = (2A)^{-1}
     */
    SYCL_EXTERNAL __always_inline void Function_dcdmu(double *mu, double *phase_comp, double T, long a, double *dcdmu) const;
    
    /**
     * @brief Grand potential derivative ∂Ψ/∂φ_α
     * Ψ = Σ_β h_β(φ) ψ_β where ψ_β = f_β - μ·c_β
     */
    SYCL_EXTERNAL __always_inline double Function_dpsi(double *mu, double *phase_comp, double T, double *phi, long a) const;
    
    /**
     * @brief Interpolated relaxation time τ(φ)
     * τ = Σ_{α<β} τ_αβ φ_α φ_β / Σ_{α<β} φ_α φ_β
     */
    SYCL_EXTERNAL __always_inline double FunctionTau(double Tau, const double *phi) const;

    /**
     * @brief Grand potential density ψ_α = f_α - μ·c_α
     * (Eq. 3.20 in thesis)
     */
    SYCL_EXTERNAL __always_inline double GrandPotential(const double *c, const double *mu, double T, long a) const;

    // ========================================================================
    // Interpolation functions for multi-phase
    // ========================================================================
    
    /**
     * @brief Derivative of interpolation function ∂h_β/∂φ_α
     * For multi-phase: h_α = φ_α² Σ_{β≠α} φ_β² + φ_α³ (normalized)
     * (Eq. 3.10-3.11 in thesis)
     */
    static SYCL_EXTERNAL __always_inline double dhphi(double *phi, long b, long a, int nPhase);
    
    /**
     * @brief Interpolation function h_α(φ)
     * Satisfies: h_α(1,0,...) = 1, h_α(0,...) = 0, Σh_α = 1
     */
    static SYCL_EXTERNAL __always_inline double hphi(double *phi, long a, int nPhase);

    /**
     * @brief Standard smooth interpolation h(φ) = φ³(10 - 15φ + 6φ²)
     * Used for two-phase systems
     */
    static SYCL_EXTERNAL __always_inline double h_smooth(double phi);
    
    /**
     * @brief Derivative dh/dφ = 30φ²(1-φ)²
     */
    static SYCL_EXTERNAL __always_inline double dh_smooth(double phi);

    // Initialize from host
    void Initialize_host();
};

} // namespace microsim

// ============================================================================
// Implementation
// ============================================================================

inline double microsim::ThermoDynamics::Free_energy(const double *c, double T, long a) const {
    double sum = 0.0;
    
    // Quadratic terms: Σ_{i≤j} A_α^{ij} c_i c_j
    for (int i = 0; i < nComp - 1; i++) {
        for (int j = 0; j < nComp - 1; j++) {
            if (i <= j) {
                sum += A[a_ind(a, i, j)] * c[i] * c[j];
            }
        }
        
        // Update B if non-isothermal
        if (!isothermal) {
            B[delC_ind(a, i)] = Beq[delC_ind(a, i)] + dBbdT[delC_ind(a, i)] * (T - T_eq);
        }
        
        // Linear term: B_α^i c_i
        sum += B[delC_ind(a, i)] * c[i];
    }
    
    // Update C if non-isothermal
    if (!isothermal) {
        C[a] = functionC(T, a);
    }
    
    // Constant term
    sum += C[a];
    
    return sum;
}

inline double microsim::ThermoDynamics::GrandPotential(const double *c, const double *mu, double T, long a) const {
    // ψ_α = f_α(c_α) - μ·c_α
    double f_alpha = Free_energy(c, T, a);
    double mu_dot_c = 0.0;
    
    for (int k = 0; k < nComp - 1; ++k) {
        mu_dot_c += mu[k] * c[k];
    }
    
    return f_alpha - mu_dot_c;
}

inline double microsim::ThermoDynamics::functionB(double T, long i, long a) const {
    double c_liq[MAX_COMP - 1];
    double c_sol[MAX_COMP - 1];
    double sum_c = 0.0;
    double B_ai = 0.0;

    if (a != (nPhase - 1)) {
        for (int k = 0; k < nComp - 1; ++k) {
            c_liq[k] = ceq[ceq_index(i, nPhase - 1, k)] - 
                       (DELTA_C[delC_ind(i, k)] * (T_eq - T) / DELTA_T[delT_ind(a, nPhase - 1)]);
            c_sol[k] = ceq[ceq_index(a, a, k)] - 
                       (DELTA_C[delC_ind(a, k)] * (T_eq - T) / DELTA_T[delT_ind(a, a)]);

            if (k != i) {
                sum_c += A[a_ind(nPhase - 1, k, i)] * c_liq[k] - A[a_ind(a, k, i)] * c_sol[k];
            }
        }
        B_ai = 2.0 * (A[a_ind(nPhase - 1, i, i)] * c_liq[i] - A[a_ind(a, i, i)] * c_sol[i]) + sum_c;
    }

    return B_ai;
}

inline double microsim::ThermoDynamics::functionC(double T, long a) const {
    double c_liq[MAX_COMP - 1];
    double c_sol[MAX_COMP - 1];
    double sum_c = 0.0;

    if (a != (nPhase - 1)) {
        for (int k = 0; k < nComp - 1; k++) {
            c_liq[k] = ceq[ceq_index(a, nPhase - 1, k)] - 
                       (DELTA_C[delC_ind(a, k)]) * (T_eq - T) / DELTA_T[delT_ind(a, nPhase - 1)];
            c_sol[k] = ceq[ceq_index(a, a, k)] - 
                       (DELTA_C[delC_ind(a, k)]) * (T_eq - T) / DELTA_T[delT_ind(a, a)];
        }

        for (int i = 0; i < nComp - 1; i++) {
            for (int j = 0; j < nComp - 1; j++) {
                if (i <= j) {
                    sum_c += A[a_ind(a, i, j)] * c_sol[i] * c_sol[j] - 
                             A[a_ind(nPhase - 1, i, j)] * c_liq[i] * c_liq[j];
                }
            }
        }
    }

    return sum_c;
}

inline void microsim::ThermoDynamics::Function_Mu(double *c, double T, long a, double *Mu) const {
    // μ_k = ∂f_α/∂c_k = 2 A_α^{kk} c_k + Σ_{j≠k} A_α^{kj} c_j + B_α^k
    for (int k = 0; k < nComp - 1; k++) {
        double sum = 2.0 * A[a_ind(a, k, k)] * c[k] + 
                     (Beq[delC_ind(a, k)] + dBbdT[delC_ind(a, k)] * (T - T_eq));
        
        for (int j = 0; j < nComp - 1; j++) {
            if (k != j) {
                sum += A[a_ind(a, k, j)] * c[j];
            }
        }
        Mu[k] = sum;
    }
}

inline void microsim::ThermoDynamics::Function_C_mu(double *mu, double *c, double T, long a) const {
    // c_α = (∂²f_α/∂c²)^{-1} · (μ - B_α)
    // cmu stores the inverse Hessian
    for (int l = 0; l < nComp - 1; l++) {
        double sum = 0.0;
        for (int k = 0; k < nComp - 1; k++) {
            sum += cmu[a_ind(a, l, k)] * (mu[k] - (Beq[delC_ind(a, k)] + dBbdT[delC_ind(a, k)] * (T - T_eq)));
        }
        c[l] = sum;
    }
}

inline void microsim::ThermoDynamics::Function_dcdmu(double *mu, double *phase_comp, double T, long a, double *dcdmu) const {
    // ∂c_α/∂μ = (∂²f_α/∂c²)^{-1} = cmu
    for (int i = 0; i < nComp - 1; i++) {
        for (int j = 0; j < nComp - 1; j++) {
            dcdmu[i * (nComp - 1) + j] = cmu[a_ind(a, i, j)];
        }
    }
}

inline double microsim::ThermoDynamics::Function_dpsi(double *mu, double *phase_comp, double T, double *phi, long a) const {
    // ∂Ψ/∂φ_α = Σ_β (∂h_β/∂φ_α) ψ_β
    // where ψ_β = f_β - μ·c_β
    double sum = 0.0;
    double c[MAX_COMP - 1];

    for (int b = 0; b < nPhase; b++) {
        // Get phase concentration from storage
        double psi = 0.0;
        for (int k = 0; k < nComp - 1; k++) {
            c[k] = phase_comp[b * (nComp - 1) + k];
            psi -= mu[k] * c[k];
        }
        psi += Free_energy(c, T, b);
        sum += psi * dhphi(phi, b, a, nPhase);
    }
    
    return sum;
}

inline double microsim::ThermoDynamics::dhphi(double *phi, long b, long a, int nPhase) {
    // Multi-phase interpolation derivative (Eq. 3.11 in thesis)
    // ∂h_β/∂φ_α for the multi-obstacle potential formulation
    double sum = 0.0;

    if (a == b) {
        // Diagonal term: ∂h_α/∂φ_α = 6φ_α(1-φ_α) + Σ_{e<f, e,f≠α} 2φ_eφ_f
        sum = 6.0 * phi[a] * (1.0 - phi[a]);
        for (int e = 0; e < nPhase; e++) {
            for (int f = 0; f < nPhase; f++) {
                if (e != b && f != b && e < f) {
                    sum += 2.0 * phi[e] * phi[f];
                }
            }
        }
    } else {
        // Off-diagonal term: ∂h_β/∂φ_α = 2φ_β Σ_{e≠α,β} φ_e
        for (int e = 0; e < nPhase; e++) {
            if (e != b && e != a) {
                sum += 2.0 * phi[e];
            }
        }
        sum *= phi[b];
    }
    
    return sum;
}

inline double microsim::ThermoDynamics::hphi(double *phi, long a, int nPhase) {
    // Multi-phase interpolation function (Eq. 3.10 in thesis)
    // h_α = 3φ_α² - 2φ_α³ + 2φ_α Σ_{β<γ, β,γ≠α} φ_βφ_γ
    double sum1 = 0.0;
    double sum = 3.0 * phi[a] * phi[a] - 2.0 * phi[a] * phi[a] * phi[a];
    
    for (int b = 0; b < nPhase; b++) {
        for (int c = 0; c < nPhase; c++) {
            if (b != a && c != a && b < c) {
                sum1 += phi[b] * phi[c];
            }
        }
    }
    sum1 *= 2.0 * phi[a];
    
    return sum + sum1;
}

inline double microsim::ThermoDynamics::h_smooth(double phi) {
    // Standard smooth interpolation h(φ) = φ³(10 - 15φ + 6φ²)
    // Satisfies h(0) = 0, h(1) = 1, h'(0) = h'(1) = 0
    return phi * phi * phi * (10.0 - 15.0 * phi + 6.0 * phi * phi);
}

inline double microsim::ThermoDynamics::dh_smooth(double phi) {
    // dh/dφ = 30φ²(1-φ)²
    return 30.0 * phi * phi * (1.0 - phi) * (1.0 - phi);
}

inline double microsim::ThermoDynamics::FunctionTau(double Tau, const double *phi) const {
    // Interpolated relaxation time (Eq. 3.32 in thesis)
    // τ(φ) = Σ_{α<β} τ_αβ φ_α φ_β / Σ_{α<β} φ_α φ_β
    double sum = 0.0;
    double sum1 = 0.0;
    
    for (int a = 0; a < nPhase; ++a) {
        for (int b = 0; b < nPhase; ++b) {
            if (a < b) {
                sum += Tau * phi[a] * phi[b];
                sum1 += phi[a] * phi[b];
            }
        }
    }

    if (sum1 > 1.0e-12) {
        return sum / sum1;
    } else {
        return Tau;
    }
}

__attribute__((weak)) void microsim::ThermoDynamics::Initialize_host() {
    // Build the Hessian matrix muc = ∂μ/∂c and its inverse cmu = ∂c/∂μ
    for (int a = 0; a < nPhase; ++a) {
        for (int i = 0; i < nComp - 1; ++i) {
            for (int j = 0; j < nComp - 1; ++j) {
                if (i == j) {
                    muc[a_ind(a, i, j)] = 2.0 * A[a_ind(a, i, j)];
                } else {
                    muc[a_ind(a, i, j)] = A[a_ind(a, i, j)];
                }
            }
        }
        // Invert muc to get cmu
        matrix_inversion_GSL(muc, cmu, a, nComp - 1);
    }

    // Compute temperature-dependent parameters
    for (int a = 0; a < nPhase; ++a) {
        DELTA_T[delT_ind(a, nPhase - 1)] = 0.0;
        DELTA_T[delT_ind(a, a)] = 0.0;

        for (int k = 0; k < nComp - 1; ++k) {
            DELTA_T[delT_ind(a, a)] += slopes[ceq_index(a, a, k)] * 
                                       (ceq[ceq_index(a, nPhase - 1, k)] - ceq[ceq_index(a, a, k)]);
            DELTA_T[delT_ind(a, nPhase - 1)] += slopes[ceq_index(a, nPhase - 1, k)] * 
                                                (ceq[ceq_index(a, nPhase - 1, k)] - ceq[ceq_index(a, a, k)]);
            DELTA_C[delC_ind(a, k)] = ceq[ceq_index(a, nPhase - 1, k)] - ceq[ceq_index(a, a, k)];
        }

        for (int k = 0; k < nComp - 1; ++k) {
            dcbdT[ceq_index(a, a, k)] = DELTA_C[delC_ind(a, k)] / DELTA_T[delT_ind(a, a)];
            dcbdT[ceq_index(a, nPhase - 1, k)] = DELTA_C[delC_ind(a, k)] / DELTA_T[delT_ind(a, nPhase - 1)];

            dBbdT[delC_ind(a, k)] = 2.0 * (A[a_ind(nPhase - 1, k, k)] * dcbdT[ceq_index(a, nPhase - 1, k)] - 
                                          A[a_ind(a, k, k)] * dcbdT[ceq_index(a, a, k)]);
            for (int i = 0; i < nComp - 1; i++) {
                if (k != i) {
                    dBbdT[delC_ind(a, k)] += (A[a_ind(nPhase - 1, k, i)] * dcbdT[ceq_index(a, nPhase - 1, i)] - 
                                             A[a_ind(a, k, i)] * dcbdT[ceq_index(a, a, i)]);
                }
            }
        }
    }

    // Zero out liquid phase dBbdT
    for (int k = 0; k < nComp - 1; ++k) {
        dBbdT[delC_ind(nPhase - 1, k)] = 0.0;
    }

    // Compute B at equilibrium
    for (int a = 0; a < nPhase; a++) {
        for (int i = 0; i < nComp - 1; i++) {
            Beq[delC_ind(a, i)] = functionB(T_eq, i, a);
        }
    }

    // Compute B and C at operating temperature
    for (int a = 0; a < nPhase; ++a) {
        for (int i = 0; i < nComp - 1; ++i) {
            B[delC_ind(a, i)] = Beq[delC_ind(a, i)] + dBbdT[delC_ind(a, i)] * (T - T_eq);
        }
        C[a] = functionC(T, a);
    }
}