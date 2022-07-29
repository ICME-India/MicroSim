#ifndef UPDATEPHI_CUH_
#define UPDATEPHI_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"

#ifndef MAX_NUM_PHASES
#define MAX_NUM_PHASES 5
#endif

#ifndef MAX_NUM_COMP
#define MAX_NUM_COMP 5
#endif

#ifndef MAX_NUM_PHASE_COMP
#define MAX_NUM_PHASE_COMP 16
#endif

/*
 * Kernel that solves d(\phi_{i})/dt = -L/N \sum_{j=1, i \neq j}^{N} (df/dphi_{i} - df/dphi_[j})
 */
__global__
void __updatePhi__(double **phi, double **dfdphi, double **phiNew,
                   double *relaxCoeff,
                   int NUMPHASES, int NUMCOMPONENTS,
                   int sizeX, int sizeY, int sizeZ,
                   double DELTA_t);

__global__
void __updatePhiBinary__(double **phi, double **dfdphi, double **phiNew,
                         double *relaxCoeff, double kappaPhi,
                         int NUMPHASES, int NUMCOMPONENTS,
                         int sizeX, int sizeY, int sizeZ,
                         double DELTA_X, double DELTA_Y, double DELTA_Z,
                         double DELTA_t);

/*
 * Host-side wrapper function for __updatePhi__
 */
#ifdef __cplusplus
extern "C"
#endif
void updatePhi(double **phi, double **dfdphi, double **phiNew,
               domainInfo* simDomain, controls* simControls,
               simParameters* simParams, subdomainInfo* subdomain,
               dim3 gridSize, dim3 blockSize);

#endif
