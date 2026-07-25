#pragma once

#include <sycl/sycl.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

/// @attention This is a Abstraction layer to handle multi-dimensional
///
/// @if DIMENSIONS 3
/// id[0] -> z axis, id[1] -> y axis, id[2] -> x axis
/// @elseif DIMENSIONS 2
/// id[0] -> y axis, id[1] -> x axis
/// @endif
/// @note Coords struct provides conversion between sycl::id and (z,y,x)
///
/// @brief idx converts 3D coordinates to linear index

template <int Dim> struct Coords;
template <> struct Coords<2> {
    int z, y, x;
    static Coords from_id(sycl::id<2> id) { return {0, (int)id[0], (int)id[1]}; }
    static sycl::range<2> to_range(int d, int h, int w) { return { (size_t)h, (size_t)w }; }
    static sycl::range<2> local_range(int x){return {(size_t(std::sqrt(x))), (size_t(std::sqrt(x)))};}
};
template <> struct Coords<3> {
    int z, y, x;
    static Coords from_id(sycl::id<3> id) { return {(int)id[0], (int)id[1], (int)id[2]}; }
    static sycl::range<3> to_range(int d, int h, int w) { return { (size_t)d, (size_t)h, (size_t)w }; }
    static sycl::range<3> local_range(int x) {
        // Find factors that are as close as possible to each other
        size_t z_dim = (size_t)std::cbrt(x);
        size_t y_dim = z_dim;
        size_t x_dim = z_dim;
        size_t product = z_dim * y_dim * x_dim;
        if (product < (size_t)x) {
            x_dim = (size_t)x / (z_dim * y_dim);
        }
        
        return {z_dim, y_dim, x_dim};
    }
};

inline int idx(int z, int y, int x, int h_stride, int w_stride) {
    return z * (h_stride * w_stride) + y * w_stride + x;
}


/* // Simplex projection to keep the phi in bounds [0,1] and sum to 1
SYCL_EXTERNAL __always_inline
void projection_on_simplex(const double *phi, const double div[] ,double *deltaphi, int nPhase) {
    double Deltaphi = 0.0;
    double count_phases = 0;
    double Deltaphi_alpha = 0.0;
    double sum_phib = 0.0;
  
    for (int a=0; a < nPhase; a++) {
      if ((sycl::fabs(div[a]) > 0.0) && (phi[a]+deltaphi[a]) > 0.0  && (phi[a]+deltaphi[a]) < 1.0) {
        count_phases++;
        sum_phib += phi[a];
      }
    }

    for (int a=0; a < nPhase; a++) {
        if ((phi[a]+deltaphi[a]) < 0.0) {
            Deltaphi_alpha = sycl::fabs(phi[a]+deltaphi[a]);
            deltaphi[a]   += Deltaphi_alpha;
            Deltaphi      += Deltaphi_alpha;
        }
    }
    
    if (Deltaphi > 0.0) {
        for (int b=0; b < nPhase; b++) {
            if (((sycl::fabs(div[b]) > 0.0)) && ((phi[b]+deltaphi[b]) > 0.0)  && ((phi[b]+deltaphi[b]) < 1.0)) {
                if (sycl::fabs(sum_phib) > 0.0) {
                    deltaphi[b] -= Deltaphi*(phi[b])/(sum_phib);
                } else {
                    deltaphi[b] -= Deltaphi/(count_phases);
                }
            }
        }
    }
    
    for (int a=0; a < nPhase; a++) {
        if (sycl::fabs(div[a]) > 0.0) {
            if ((phi[a]+deltaphi[a]) > 1.0) {
                deltaphi[a]  = (1.0-phi[a]);
                //Correct all the other phases due to this correction,
                //If you are bulk, all other phases must go to zero.
                for (int b=0; b < nPhase; b++) {
                    if ((sycl::fabs(div[b]) > 0.0) && b!=a) {
                        deltaphi[b] = -phi[b];
                    }
                }
                Deltaphi = 0.0;
                break;
            }
        }
    }
}
 */

__attribute__((weak)) void matrix_inversion_GSL(double *Mat, double *InvMat, int a, int N){
    int s;
    // Allocate memory
    gsl_matrix *matrix      = gsl_matrix_alloc(N, N);
    gsl_matrix *Inv_matrix  = gsl_matrix_alloc(N, N);
    gsl_permutation *perm = gsl_permutation_alloc(N);

    gsl_matrix_set_all(matrix, 0.0);
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            gsl_matrix_set(matrix, i, j, Mat[a*N*N + i*N + j]);
        }
    }

    int status =gsl_linalg_LU_decomp(matrix, perm, &s);

    if(status != GSL_SUCCESS){
        std::__throw_runtime_error("GSL matrix inversion failed");
    }else{
        gsl_linalg_LU_invert(matrix, perm, Inv_matrix);
    }

    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            InvMat[a*N*N + i*N + j] = gsl_matrix_get(Inv_matrix, i, j);
        }
    }

    // Free the allocated memory
    gsl_permutation_free(perm);
    gsl_matrix_free(matrix);
    gsl_matrix_free(Inv_matrix);
}