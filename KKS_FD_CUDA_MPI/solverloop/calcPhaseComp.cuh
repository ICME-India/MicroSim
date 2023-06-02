#ifndef CALCPHASECOMP_CUH_
#define CALCPHASECOMP_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "utilityKernels.cuh"
#include "Thermo.cuh"
#include "matrix.cuh"
#include "functionF.cuh"
#include "functionH.cuh"


/*
 * Calculate phase compositions for Function_F != 2
 */
__global__
void __calcPhaseComp__(double **phi, double **comp,
                       double **phaseComp,
                       double *F0_A, double *F0_B, double *F0_C,
                       long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                       long sizeX, long sizeY, long sizeZ, long padding);

/*
 *  Initialise diffusion potentials for Function_F == 2
 */
__global__
void __initMu__(double **phi, double **comp, double **phaseComp, double **mu,
                long *thermo_phase, double temperature,
                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                long xStep, long yStep, long padding);

/*
 * Calculate phase compositions for Function_F == 2
 */
__global__
void __calcPhaseComp_02__(double **phi, double **comp,
                          double **phaseComp, double **mu, double *cguess,
                          double temperature, long *thermo_phase,
                          long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                          long sizeX, long sizeY, long sizeZ, long padding);

/*
 * Wrapper function for __calcPhaseComp__
 */
#ifdef __cplusplus
extern "C"
#endif
void calcPhaseComp(double **phi, double **comp,
                   double **phaseComp, double **mu,
                   domainInfo* simDomain, controls* simControls,
                   simParameters* simParams, subdomainInfo* subdomain,
                   dim3 gridSize, dim3 blockSize);

#endif
