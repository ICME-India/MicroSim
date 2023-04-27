#ifndef FUNCTIONA_02_CUH_
#define FUNCTIONA_02_CUH_

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "structures.h"
#include "anisotropy_01.cuh"

#define ix 0
#define iy 1
#define iz 2

#define centre 0
#define left   1
#define right  2
#define top    3
#define bottom 4
#define front  5
#define back   6

/*
 * __device__ double calcAnisotropy_02
 *
 * Calculate anisotropy
 *
 * Arguments:
 *              1. double phi[MAX_NUM_PHASES][3][3][3] - Values of phi at 26 nearest-neighbours and second-nearest neighbours,
 *                                                       arranged phase-by-phase (flattened). [0:N-1][i-1:i+1][j-1:j+1][k-1][k+1]
 *              2. double *dab - Anisotropy strength for every pair of phases.
 *              3. double *eps_ab - kappa
 *              4. double *Rotation_matrix - [N][N][3][3]
 *              5. double *Inv_rotation_matrix - [N][N][3][3]
 *              6. long phase - Master phase
 *              7. long NUMPHASES - Number of phases
 *              8. double DELTA_X
 *              9. double DELTA_Y
 *              10. double DELTA_Z
 *
 */
extern __device__
double calcAnisotropy_02(double phi[MAX_NUM_PHASES][3][3][3],
                         double *dab, double *eps_ab,
                         double *Rotation_matrix, double *Inv_rotation_matrix,
                         long phase, long NUMPHASES, long DIMENSION,
                         double DELTA_X, double DELTA_Y, double DELTA_Z);

#endif
