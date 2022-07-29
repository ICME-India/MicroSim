#ifndef UPDATECOMPOSITION_CUH_
#define UPDATECOMPOSITION_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "Thermo.cuh"

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
 * Kernel that solves dc_{j}/dt = div(M grad(mu)) using a fourth order FD stencil and forward Euler
 */
__global__
void __updateComposition__(double **phi,
                           double **comp, double **compNew,
                           double **phaseComp,
                           double *F0_A, double *F0_B,
                           double *diffusivity,
                           int NUMPHASES, int NUMCOMPONENTS,
                           int sizeX, int sizeY, int sizeZ,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t);

__global__
void __updateCompositionBinary__(double **phi,
                                 double **comp, double **compNew,
                                 double **phaseComp,
                                 double *F0_A, double *F0_B,
                                 double *diffusivity,
                                 int NUMPHASES, int NUMCOMPONENTS,
                                 int sizeX, int sizeY, int sizeZ,
                                 double DELTA_X, double DELTA_Y, double DELTA_Z,
                                 double DELTA_t);

/*
 * Host-side wrapper function for __updateComposition__
 */
#ifdef __cplusplus
extern "C"
#endif
void updateComposition(double **phi, double **comp, double **phiNew, double **compNew,
                       double **phaseComp, double **mu,
                       domainInfo* simDomain, controls* simControls,
                       simParameters* simParams, subdomainInfo* subdomain,
                       dim3 gridSize, dim3 blockSize);

#endif
