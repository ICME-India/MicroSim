#ifndef UPDATEPHI_CUH_
#define UPDATEPHI_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "Thermo.cuh"
#include "utilityKernels.cuh"
#include "matrix.cuh"
#include "functionA_01.cuh"
#include "functionA_02.cuh"
#include "functionTau.cuh"
#include "boundary.cuh"

/*
 * Kernel that solves d(\phi_{i})/dt = -L/N \sum_{j=1, i \neq j}^{N} (df/dphi_{i} - df/dphi_[j})
 */
__global__
void __updatePhi__(double **phi, double **dfdphi, double **phiNew,
                   double *relaxCoeff, double *kappaPhi,
                   double *dab, int FUNCTION_ANISOTROPY,
                   long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ,
                   long yStep, long zStep, long padding,
                   double DELTA_X, double DELTA_Y, double DELTA_Z,
                   double DELTA_t);

/*
 * Host-side wrapper function for __updatePhi__
 */
#ifdef __cplusplus
extern "C"
#endif
void updatePhi(double **phi, double **dfdphi, double **phiNew, double **phaseComp,
               domainInfo* simDomain, controls* simControls,
               simParameters* simParams, subdomainInfo* subdomain,
               dim3 gridSize, dim3 blockSize);

#endif
