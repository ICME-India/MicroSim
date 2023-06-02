#ifndef UPDATECOMPOSITION_CUH_
#define UPDATECOMPOSITION_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "Thermo.cuh"
#include "utilityKernels.cuh"
#include "matrix.cuh"
#include "functionH.cuh"
#include "functionF.cuh"
#include "boundary.cuh"

/*
 * Kernel that solves dc_{j}/dt = div(M grad(mu)) using a fourth order FD stencil and forward Euler
 */

// For functions F_01, F_03, F_04
__global__
void __updateComposition__(double **phi, double **phiNew,
                           double **comp, double **compNew,
                           double **phaseComp,
                           double *F0_A, double *F0_B,
                           double *mobility, double *diffusivity, double *kappaPhi, double *theta_ij,
                           long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                           long sizeX, long sizeY, long sizeZ,
                           long xStep, long yStep, long padding, int antiTrapping,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t);

// For function F_02
__global__
void __updateComposition_02__(double **phi, double **phiNew,
                              double **comp, double **compNew, double **mu,
                              double **phaseComp, long *thermo_phase,
                              double *diffusivity, double *kappaPhi, double *theta_ij,
                              double temperature, double molarVolume,
                              long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                              long sizeX, long sizeY, long sizeZ,
                              long xStep, long yStep, long padding,
                              double DELTA_X, double DELTA_Y, double DELTA_Z,
                              double DELTA_t);

__global__
void __updateMu_02__(double **phi, double **comp,
                     double **phiNew, double **compNew,
                     double **phaseComp, double **mu,
                     long *thermo_phase, double temperature, double molarVolume,
                     long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                     long sizeX, long sizeY, long sizeZ,
                     long xStep, long yStep, long padding,
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
                       domainInfo simDomain, controls simControls,
                       simParameters simParams, subdomainInfo subdomain,
                       dim3 gridSize, dim3 blockSize);

#endif
