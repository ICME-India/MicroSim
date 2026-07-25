#pragma once

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>
#include <iostream>
#include <vector>

#include <sycl/sycl.hpp>
#include "utilities/helper_functions.hpp"
#include "utilities/Input_parameters.hpp"
#include "utilities/microsim.hpp"


namespace microsim {
    
    struct ThermoDynamicsF4 {
        // Member variables
        int nPhase;
        int nComp;
        double T_eq;
        double T;
        bool isothermal;
        int *thermo_phase;  // Mapping from phase to thermodynamic phase
        
        // Thermodynamic matrices (flat arrays)
        double *A;           // [nPhase * (nComp-1) * (nComp-1)]
        double *B;           // [nPhase * (nComp-1)]
        double *C;           // [nPhase]
        double *Beq;         // [nPhase * (nComp-1)]
        double *dBbdT;       // [nPhase * (nComp-1)]
        double *dcbdT;       // [nPhase * nPhase * (nComp-1)]
        double *ceq;         // [nPhase * nPhase * (nComp-1)]
        double *DELTA_C;     // [nPhase * (nComp-1)]
        double *DELTA_T;     // [nPhase * nPhase]
        double *cmu;         // [nPhase * (nComp-1) * (nComp-1)]
        double *muc;         // [nPhase * (nComp-1) * (nComp-1)]
        double *slopes;      // [nPhase * nPhase * (nComp-1)]
        double *c_guess;     // [nPhase * nPhase * (nComp-1)]
        
        // Spline coefficients (host and device)
        SplineCoefficients *ES_coeffs_flat;   // Host
        SplineCoefficients *Thf_coeffs_flat;  // Host
        SplineCoefficients *d_ES_coeffs_flat; // Device
        SplineCoefficients *d_Thf_coeffs_flat;// Device
        
        // Device spline coefficient arrays (contiguous memory)
        double *d_ES_c0, *d_ES_c1, *d_ES_c2, *d_ES_c3, *d_ES_x;
        double *d_Thf_c0, *d_Thf_c1, *d_Thf_c2, *d_Thf_c3, *d_Thf_x;
        
        // Temporary arrays for non-isothermal calculations
        double *temp_muc;    // Temporary for matrix operations
        double *temp_cmu;    // Temporary for matrix inversion
        
        // SYCL queue pointer (optional, for device operations)
        sycl::queue *queue_ptr;
        
        // Constructor
        ThermoDynamicsF4(microsim::Info *global_fd, bool isothermal_) {
            nPhase = global_fd->parameters->NUM_PHASES;
            nComp = global_fd->parameters->NUM_COMPONENTS;
            T_eq = global_fd->parameters->T_eq;
            T = global_fd->parameters->T;
            isothermal = isothermal_;

            Phasefield_matraix *pfm = global_fd->PFMAT;
            
            // Assign pointers to pre-allocated arrays
            A = pfm->A;
            B = pfm->B;
            C = pfm->C;
            muc = pfm->muc;
            cmu = pfm->cmu;
            slopes = pfm->slopes;
            c_guess = pfm->c_guess;
            ceq = pfm->ceq;
            Beq = pfm->Beq;
            dBbdT = pfm->dBbdT;
            dcbdT = pfm->dcbdT;
            DELTA_C = pfm->DELTA_C;
            DELTA_T = pfm->DELTA_T;
            
            // Initialize spline pointers to nullptr
            ES_coeffs_flat = nullptr;
            Thf_coeffs_flat = nullptr;
            d_ES_coeffs_flat = nullptr;
            d_Thf_coeffs_flat = nullptr;
            
            d_ES_c0 = d_ES_c1 = d_ES_c2 = d_ES_c3 = d_ES_x = nullptr;
            d_Thf_c0 = d_Thf_c1 = d_Thf_c2 = d_Thf_c3 = d_Thf_x = nullptr;
            
            temp_muc = nullptr;
            temp_cmu = nullptr;
            thermo_phase = nullptr;
            queue_ptr = nullptr;
        }
        
        // ====================================================================
        // INDEXING HELPER FUNCTIONS
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline int a_ind(int p, int c1, int c2) const {
            return p * (nComp - 1) * (nComp - 1) + c1 * (nComp - 1) + c2;
        }
        
        SYCL_EXTERNAL __always_inline int b_ind(int p, int c) const {
            return p * (nComp - 1) + c;
        }
        
        SYCL_EXTERNAL __always_inline int c_ind(int p) const {
            return p;
        }
        
        SYCL_EXTERNAL __always_inline int ceq_index(int p1, int p2, int a) const {
            return p1 * nPhase * (nComp - 1) + p2 * (nComp - 1) + a;
        }
        
        SYCL_EXTERNAL __always_inline int delC_ind(int p, int c) const {
            return p * (nComp - 1) + c;
        }
        
        SYCL_EXTERNAL __always_inline int delT_ind(int p1, int p2) const {
            return p1 * nPhase + p2;
        }
        
        SYCL_EXTERNAL __always_inline int get_ES_flat_index(int phase, int comp, int idx) const {
            return phase * (nComp - 1) * 2 + comp * 2 + idx;
        }
        
        SYCL_EXTERNAL __always_inline int get_Thf_flat_index(int phase, int comp1, int comp2) const {
            return phase * (nComp - 1) * (nComp - 1) + comp1 * (nComp - 1) + comp2;
        }
        
        SYCL_EXTERNAL __always_inline int muc_ind(int p, int i, int j) const {
            return p * (nComp - 1) * (nComp - 1) + i * (nComp - 1) + j;
        }
        
        // ====================================================================
        // SPLINE EVALUATION (STATIC - DEVICE COMPATIBLE)
        // ====================================================================
        
        static SYCL_EXTERNAL __always_inline double evaluate_spline_sycl(
            const SplineCoefficients &coeffs,
            double x
        ) {
            // Check for invalid coefficients
            if (coeffs.n_points <= 0 || !coeffs.c0 || !coeffs.x) {
                return 0.0; // Return safe default instead of NaN
            }

            // Binary search to find the interval
            int low = 0, high = coeffs.n_points;

            if (x <= coeffs.x[0]) {
                low = 0;
            } else if (x >= coeffs.x[coeffs.n_points - 1]) {  // Fixed: n_points-1 is the last index
                low = coeffs.n_points - 1;
            } else {
                while (high - low > 1) {
                    int mid = (low + high) / 2;
                    if (coeffs.x[mid] > x) {
                        high = mid;
                    } else {
                        low = mid;
                    }
                }
            }

            // Evaluate cubic polynomial with bounds checking
            double dx = x - coeffs.x[low];
            double result = coeffs.c0[low] + dx * (coeffs.c1[low] +
                   dx * (coeffs.c2[low] + dx * coeffs.c3[low]));

            // Check for NaN or Inf in result
            #ifdef __SYCL_DEVICE_ONLY__
            if (sycl::isnan(result) || sycl::isinf(result)) {
                return 0.0;
            }
            #else
            if (std::isnan(result) || std::isinf(result)) {
                return 0.0;
            }
            #endif

            return result;
        }
        
        // ====================================================================
        // FUNCTION A: Populate A matrix from splines
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline void function_A(double T_val) const {
            for(int a = 0; a < nPhase; ++a) {
                int thermo_ph = thermo_phase[a];
                
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int thf_idx = get_Thf_flat_index(thermo_ph, i, j);
                        double val = evaluate_spline_sycl(d_Thf_coeffs_flat[thf_idx], T_val);
                        
                        int a_idx = a_ind(a, i, j);
                        if(i == j) {
                            A[a_idx] = 0.5 * val;
                        } else {
                            A[a_idx] = val;
                        }
                    }
                }
            }
        }
        
        // ====================================================================
        // FUNCTION B: Calculate B vector
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline double functionB(double T_val, long i, long a) const {
            if(a == nPhase - 1) {
                return 0.0; // Liquid phase
            }
            
            int thermo_ph = thermo_phase[a];
            double sum_c = 0.0;
            
            // Get equilibrium compositions
            double c_liq[20]; // Max components
            double c_sol[20];
            
            for(int k = 0; k < nComp - 1; ++k) {
                int es_liq_idx = get_ES_flat_index(thermo_ph, k, 1);
                int es_sol_idx = get_ES_flat_index(thermo_ph, k, 0);
                
                c_liq[k] = evaluate_spline_sycl(d_ES_coeffs_flat[es_liq_idx], T_val);
                c_sol[k] = evaluate_spline_sycl(d_ES_coeffs_flat[es_sol_idx], T_val);
            }
            
            // Calculate sum_c
            for(int k = 0; k < nComp - 1; ++k) {
                if(k != i) {
                    int a_liq_idx = a_ind(nPhase - 1, k, i);
                    int a_sol_idx = a_ind(a, k, i);
                    
                    sum_c += A[a_liq_idx] * c_liq[k] - A[a_sol_idx] * c_sol[k];
                }
            }
            
            // Calculate B_ai
            int a_liq_ii = a_ind(nPhase - 1, i, i);
            int a_sol_ii = a_ind(a, i, i);
            
            return 2.0 * (A[a_liq_ii] * c_liq[i] - A[a_sol_ii] * c_sol[i]) + sum_c;
        }
        
        // ====================================================================
        // FUNCTION C: Calculate C scalar
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline double functionC(double T_val, long a) const {
            if(a == nPhase - 1) {
                return 0.0; // Liquid phase
            }
            
            int thermo_ph = thermo_phase[a];
            double sum_c = 0.0;
            
            // Get equilibrium compositions
            double c_liq[20];
            double c_sol[20];
            
            for(int k = 0; k < nComp - 1; ++k) {
                int es_liq_idx = get_ES_flat_index(thermo_ph, k, 1);
                int es_sol_idx = get_ES_flat_index(thermo_ph, k, 0);
                
                c_liq[k] = evaluate_spline_sycl(d_ES_coeffs_flat[es_liq_idx], T_val);
                c_sol[k] = evaluate_spline_sycl(d_ES_coeffs_flat[es_sol_idx], T_val);
            }
            
            // Calculate sum
            for(int i = 0; i < nComp - 1; ++i) {
                for(int j = 0; j < nComp - 1; ++j) {
                    if(i <= j) {
                        int a_sol_idx = a_ind(a, i, j);
                        int a_liq_idx = a_ind(nPhase - 1, i, j);
                        
                        sum_c += A[a_sol_idx] * c_sol[i] * c_sol[j] - 
                                 A[a_liq_idx] * c_liq[i] * c_liq[j];
                    }
                }
            }
            
            return sum_c;
        }
        
        // ====================================================================
        // FREE ENERGY CALCULATION
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline double Free_energy(
            const double *c, 
            double T_val, 
            long a
        ) const {
            double sum = 0.0;
            
            // Update A, B, C if not isothermal (this modifies the arrays)
            // Note: For true const, you'd need to pass these as parameters
            if(!isothermal) {
                // Update A matrix
                int thermo_ph = thermo_phase[a];
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int thf_idx = get_Thf_flat_index(thermo_ph, i, j);
                        double val = evaluate_spline_sycl(d_Thf_coeffs_flat[thf_idx], T_val);
                        
                        int a_idx = a_ind(a, i, j);
                        if(i == j) {
                            A[a_idx] = 0.5 * val;
                        } else {
                            A[a_idx] = val;
                        }
                    }
                }
            }
            
            // Calculate quadratic term
            for(int i = 0; i < nComp - 1; ++i) {
                for(int j = 0; j < nComp - 1; ++j) {
                    if(i <= j) {
                        int a_idx = a_ind(a, i, j);
                        sum += A[a_idx] * c[i] * c[j];
                    }
                }
                
                // Update B if not isothermal
                if(!isothermal) {
                    int b_idx = b_ind(a, i);
                    B[b_idx] = functionB(T_val, i, a);
                }
                
                // Add linear term
                int b_idx = b_ind(a, i);
                sum += B[b_idx] * c[i];
            }
            
            // Update C if not isothermal
            if(!isothermal) {
                int c_idx = c_ind(a);
                C[c_idx] = functionC(T_val, a);
            }
            
            // Add constant term
            sum += C[c_ind(a)];
            
            return sum;
        }
        
        // ====================================================================
        // CHEMICAL POTENTIAL CALCULATION
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline void Function_Mu(
            double *c, 
            double T_val, 
            long a, 
            double *Mu
        ) const {
            // Update A if not isothermal
            if(!isothermal) {
                int thermo_ph = thermo_phase[a];
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int thf_idx = get_Thf_flat_index(thermo_ph, i, j);
                        double val = evaluate_spline_sycl(d_Thf_coeffs_flat[thf_idx], T_val);
                        
                        int a_idx = a_ind(a, i, j);
                        if(i == j) {
                            A[a_idx] = 0.5 * val;
                        } else {
                            A[a_idx] = val;
                        }
                    }
                }
            }
            
            // Update B if not isothermal
            if(!isothermal) {
                for(int k = 0; k < nComp - 1; ++k) {
                    int b_idx = b_ind(a, k);
                    B[b_idx] = functionB(T_val, k, a);
                }
            }
            
            // Calculate chemical potential: Mu_k = 2*A_kk*c_k + sum(A_kj*c_j for j!=k) + B_k
            for(int k = 0; k < nComp - 1; ++k) {
                int a_kk_idx = a_ind(a, k, k);
                int b_k_idx = b_ind(a, k);
                
                double sum = 2.0 * A[a_kk_idx] * c[k] + B[b_k_idx];
                
                for(int j = 0; j < nComp - 1; ++j) {
                    if(j != k) {
                        int a_kj_idx = a_ind(a, k, j);
                        sum += A[a_kj_idx] * c[j];
                    }
                }
                
                Mu[k] = sum;
            }
        }
        
        // ====================================================================
        // COMPOSITION FROM CHEMICAL POTENTIAL
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline void Function_C_mu(
            double *mu, 
            double *c, 
            double T_val, 
            long a
        ) const {
            // Update A if not isothermal
            if(!isothermal) {
                int thermo_ph = thermo_phase[a];
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int thf_idx = get_Thf_flat_index(thermo_ph, i, j);
                        double val = evaluate_spline_sycl(d_Thf_coeffs_flat[thf_idx], T_val);
                        
                        int a_idx = a_ind(a, i, j);
                        if(i == j) {
                            A[a_idx] = 0.5 * val;
                        } else {
                            A[a_idx] = val;
                        }
                    }
                }
                
                // Build muc matrix from A
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int a_idx = a_ind(a, i, j);
                        int muc_idx = muc_ind(a, i, j);
                        
                        if(i == j) {
                            muc[muc_idx] = 2.0 * A[a_idx];
                        } else {
                            muc[muc_idx] = A[a_idx];
                        }
                    }
                }
                
                // Invert muc to get cmu
                // Note: You need a device-compatible matrix inversion
                MatrixInversion_device(&muc[a * (nComp-1) * (nComp-1)], 
                                      &cmu[a * (nComp-1) * (nComp-1)], 
                                      nComp - 1);
            }
            
            // Update B if not isothermal
            if(!isothermal) {
                for(int k = 0; k < nComp - 1; ++k) {
                    int b_idx = b_ind(a, k);
                    B[b_idx] = functionB(T_val, k, a);
                }
            }
            
            // Calculate composition: c = cmu * (mu - B)
            for(int l = 0; l < nComp - 1; ++l) {
                double sum = 0.0;
                for(int k = 0; k < nComp - 1; ++k) {
                    int cmu_idx = muc_ind(a, l, k); // Using same indexing for cmu
                    int b_k_idx = b_ind(a, k);
                    sum += cmu[cmu_idx] * (mu[k] - B[b_k_idx]);
                }
                c[l] = sum;
            }
        }
        
        // ====================================================================
        // DERIVATIVE OF COMPOSITION W.R.T. CHEMICAL POTENTIAL
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline void Function_dcdmu(
            double *mu, 
            double *phase_comp, 
            double T_val, 
            long a, 
            double *dcdmu
        ) const {
            // Update A if not isothermal
            if(!isothermal) {
                int thermo_ph = thermo_phase[a];
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int thf_idx = get_Thf_flat_index(thermo_ph, i, j);
                        double val = evaluate_spline_sycl(d_Thf_coeffs_flat[thf_idx], T_val);
                        
                        int a_idx = a_ind(a, i, j);
                        if(i == j) {
                            A[a_idx] = 0.5 * val;
                        } else {
                            A[a_idx] = val;
                        }
                    }
                }
                
                // Build muc matrix from A
                for(int i = 0; i < nComp - 1; ++i) {
                    for(int j = 0; j < nComp - 1; ++j) {
                        int a_idx = a_ind(a, i, j);
                        int muc_idx = muc_ind(a, i, j);
                        
                        if(i == j) {
                            muc[muc_idx] = 2.0 * A[a_idx];
                        } else {
                            muc[muc_idx] = A[a_idx];
                        }
                    }
                }
                
                // Invert muc to get cmu
                MatrixInversion_device(&muc[a * (nComp-1) * (nComp-1)], 
                                      &cmu[a * (nComp-1) * (nComp-1)], 
                                      nComp - 1);
            }
            
            // Copy cmu to dcdmu (dcdmu is the derivative matrix)
            for(int i = 0; i < nComp - 1; ++i) {
                for(int j = 0; j < nComp - 1; ++j) {
                    int idx = i * (nComp - 1) + j;
                    int cmu_idx = muc_ind(a, i, j);
                    dcdmu[idx] = cmu[cmu_idx];
                }
            }
        }
        
        // ====================================================================
        // DRIVING FORCE
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline double Function_dpsi(
            double *mu, 
            double *phase_comp, 
            double T_val, 
            double *phi, 
            long a
        ) const {
            double sum = 0.0;
            double c[20]; // Temporary array
            
            for(int i = 0; i < nComp - 1; ++i) {
                c[i] = 0.0;
                double psi = 0.0;
                
                for(int j = 0; j < nPhase; ++j) {
                    int comp_idx = ceq_index(i, j, 0); // Adjust indexing as needed
                    c[i] += phase_comp[comp_idx] * phi[j];
                }
                
                for(int j = 0; j < nPhase; ++j) {
                    psi -= mu[j] * c[j];
                }
                
                // Add free energy
                psi += Free_energy(c, T_val, i);
                
                // Multiply by derivative of interpolation function
                sum += psi * dhphi(phi, i, a, nPhase);
            }
            
            return sum;
        }
        
        // ====================================================================
        // TAU FUNCTION (for phase field evolution)
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline double FunctionTau(
            double Tau, 
            const double *phi
        ) const {
            double sum = 0.0; double sum1 = 0.0;
            for(int a=0; a<nPhase; ++a){
                for(int b=0; b<nPhase; ++b){
                    if(a<b){
                        sum += Tau * phi[a] * phi[b];
                        sum1 += phi[a] * phi[b];
                    }
                }
            }
        
            if(sum1){
                return sum/sum1;
            }else{
                return Tau;
            }
        }
        
        // ====================================================================
        // STATIC HELPER FUNCTIONS
        // ====================================================================
        
        static SYCL_EXTERNAL __always_inline double dhphi(
            double *phi, 
            long b, 
            long a, 
            int nPhase
        ) {
            double sum = 0.0;
            if (a==b) {
                sum = 6.0*phi[a]*(1.0-phi[a]);
                for (int e=0; e < nPhase; e++) {
                    for (int f=0; f < nPhase; f++) {
                        if (e!=b && f!=b && e<f) {
                            sum += 2.0*phi[e]*phi[f];
                        }
                    }
                }
            } else {
                for (int e=0; e < nPhase; e++) {
                    if (e!=b && e!=a) {
                        sum += 2.0*phi[e];
                    }
                }
                sum *= phi[b];
            }
            return sum;
        }
        
        static SYCL_EXTERNAL __always_inline double hphi(
            double *phi, 
            long a, 
            int nPhase
        ) {
            double sum1 = 0.0;
            double sum  = 3.0*phi[a]*phi[a] - 2.0*phi[a]*phi[a]*phi[a];
            for (int b=0; b < nPhase; b++) {
                for (int c=0; c < nPhase; c++) {
                    if (b!=a && c!=a && b < c) {
                        sum1 += phi[b]*phi[c];
                    }
                }
            }
            sum1 *= 2.0*phi[a];
            return sum + sum1;
        }
        
        // ====================================================================
        // DEVICE-COMPATIBLE MATRIX INVERSION
        // ====================================================================
        
        SYCL_EXTERNAL __always_inline void MatrixInversion_device(
            double *input,
            double *output,
            int n
        ) const {
            // Simple Gauss-Jordan elimination for small matrices
            // Initialize output as identity matrix
            for(int i = 0; i < n * n; ++i) {
                output[i] = (i % (n + 1) == 0) ? 1.0 : 0.0; // Identity matrix
            }

            double temp[400]; // Temporary matrix (max 20x20)
            for(int i = 0; i < n * n; ++i) {
                temp[i] = input[i];
            }

            // Check for NaN/Inf in input
            bool input_invalid = false;
            for(int i = 0; i < n * n; ++i) {
                #ifdef __SYCL_DEVICE_ONLY__
                if (sycl::isnan(temp[i]) || sycl::isinf(temp[i])) {
                    input_invalid = true;
                    break;
                }
                #else
                if (std::isnan(temp[i]) || std::isinf(temp[i])) {
                    input_invalid = true;
                    break;
                }
                #endif
            }

            if (input_invalid) {
                // Return identity matrix for invalid input
                return;
            }

            // Gauss-Jordan elimination with partial pivoting
            for(int i = 0; i < n; ++i) {
                // Find best pivot in column i
                int pivot_row = i;
                double max_pivot = 0.0;
                #ifdef __SYCL_DEVICE_ONLY__
                max_pivot = sycl::fabs(temp[i * n + i]);
                #else
                max_pivot = std::fabs(temp[i * n + i]);
                #endif

                for(int k = i + 1; k < n; ++k) {
                    #ifdef __SYCL_DEVICE_ONLY__
                    double abs_val = sycl::fabs(temp[k * n + i]);
                    #else
                    double abs_val = std::fabs(temp[k * n + i]);
                    #endif
                    if (abs_val > max_pivot) {
                        max_pivot = abs_val;
                        pivot_row = k;
                    }
                }

                // Check for singular matrix
                if(max_pivot < 1e-10) {
                    // Matrix is singular - return identity
                    for(int ii = 0; ii < n * n; ++ii) {
                        output[ii] = (ii % (n + 1) == 0) ? 1.0 : 0.0;
                    }
                    return;
                }

                // Swap rows if needed
                if (pivot_row != i) {
                    for(int j = 0; j < n; ++j) {
                        double tmp = temp[i * n + j];
                        temp[i * n + j] = temp[pivot_row * n + j];
                        temp[pivot_row * n + j] = tmp;

                        tmp = output[i * n + j];
                        output[i * n + j] = output[pivot_row * n + j];
                        output[pivot_row * n + j] = tmp;
                    }
                }

                double pivot = temp[i * n + i];

                // Scale row i
                for(int j = 0; j < n; ++j) {
                    temp[i * n + j] /= pivot;
                    output[i * n + j] /= pivot;
                }

                // Eliminate column i in other rows
                for(int k = 0; k < n; ++k) {
                    if(k != i) {
                        double factor = temp[k * n + i];
                        for(int j = 0; j < n; ++j) {
                            temp[k * n + j] -= factor * temp[i * n + j];
                            output[k * n + j] -= factor * output[i * n + j];
                        }
                    }
                }
            }

            // Final check for NaN/Inf in output
            for(int i = 0; i < n * n; ++i) {
                #ifdef __SYCL_DEVICE_ONLY__
                if (sycl::isnan(output[i]) || sycl::isinf(output[i])) {
                    // Replace with identity matrix
                    for(int ii = 0; ii < n * n; ++ii) {
                        output[ii] = (ii % (n + 1) == 0) ? 1.0 : 0.0;
                    }
                    return;
                }
                #else
                if (std::isnan(output[i]) || std::isinf(output[i])) {
                    // Replace with identity matrix
                    for(int ii = 0; ii < n * n; ++ii) {
                        output[ii] = (ii % (n + 1) == 0) ? 1.0 : 0.0;
                    }
                    return;
                }
                #endif
            }
        }
        
        // ====================================================================
        // HOST-SIDE INITIALIZATION FUNCTIONS
        // ====================================================================
        
        void Initialize_host(phases_components *ref);
        void cleanup_spline_coefficients();
        void extract_spline_coefficients(gsl_spline *spline, SplineCoefficients &coeffs, int n_points);
        void copy_to_device(sycl::queue &q);
    };
    
} 