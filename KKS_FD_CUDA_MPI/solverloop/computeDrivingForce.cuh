#ifndef COMPUTEDRIVINGFORCE_CUH_
#define COMPUTEDRIVINGFORCE_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "utilityKernels.cuh"
#include "Thermo.cuh"
#include "matrix.cuh"
#include "functionW_01.cuh"
#include "functionW_02.cuh"
#include "functionH.cuh"

/*
 * Explicit calculation of the right-hand side of the Allen-Cahn equation.
 * Evaluation of the mobility function in the Cahn-Hilliard equation.
 */
__global__
void __computeDrivingForce__(double **phi, double **comp,
                             double **dfdphi,
                             double **phaseComp,
                             double *F0_A, double *F0_B, double *F0_C,
                             double molarVolume,
                             double *theta_i, double *theta_ij, double *theta_ijk,
                             long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                             long sizeX, long sizeY, long sizeZ,
                             long yStep, long zStep, long padding);

__global__
void __computeDrivingForce_02__(double **phi, double **comp,
                                double **dfdphi, double **phaseComp,
                                double **mu,
                                double molarVolume,
                                double *theta_i, double *theta_ij, double *theta_ijk,
                                double temperature, long *thermo_phase,
                                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                                long sizeX, long sizeY, long sizeZ,
                                long yStep, long zStep, long padding);

/*
 * Wrapper function for __computeDrivingForce__
 */
#ifdef __cplusplus
extern "C"
#endif
void computeDrivingForce(double **phi, double **comp,
                         double **dfdphi,
                         double **phaseComp ,double **mu,
                         domainInfo* simDomain, controls* simControls,
                         simParameters* simParams, subdomainInfo* subdomain,
                         dim3 gridSize, dim3 blockSize);
#endif
