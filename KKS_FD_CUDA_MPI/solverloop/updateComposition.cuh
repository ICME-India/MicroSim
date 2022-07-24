#ifndef UPDATECOMPOSITION_CUH_
#define UPDATECOMPOSITION_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"

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
void updateComposition(double **phi, double **comp, double **compNew,
                       double **phaseComp,
                       domainInfo* simDomain, controls* simControls,
                       simParameters* simParams, subdomainInfo* subdomain,
                       dim3 gridSize, dim3 blockSize);

#endif
