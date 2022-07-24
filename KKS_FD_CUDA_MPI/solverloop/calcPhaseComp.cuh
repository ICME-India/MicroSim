#ifndef CALCPHASECOMP_CUH_
#define CALCPHASECOMP_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "utilityKernels.cuh"
#include "Thermo.cuh"

#ifndef NUM_PHASE_COMP
#define NUM_PHASE_COMP 12
#endif

/*
 * Calculate phase compositions
 */
__global__
void __calcPhaseComp_01_03__(double **phi, double **comp,
                             double **phaseComp,
                             double *F0_A, double *F0_B, double *F0_C,
                             int NUMPHASES, int NUMCOMPONENTS,
                             int sizeX, int sizeY, int sizeZ);

__global__
void __calcPhaseCompBinary_01_03__(double **phi, double **comp,
                                   double **phaseComp,
                                   double *F0_A, double *F0_B, double *F0_C,
                                   int NUMPHASES, int NUMCOMPONENTS,
                                   int sizeX, int sizeY, int sizeZ);


__global__
void initPhaseComp_02(double **phi, double **comp,
                      double **phaseComp, double *cguess,
                      double temperature, int *thermo_phase,
                      int sizeX, int sizeY, int sizeZ);

__global__
void __calcPhaseCompBinary_02__(double **phi, double **comp,
                                double **phaseComp,
                                double temperature, int *thermo_phase,
                                int sizeX, int sizeY, int sizeZ);

/*
 * Wrapper function for __calcPhaseComp__
 */
#ifdef __cplusplus
extern "C"
#endif
void calcPhaseComp(double **phi, double **comp,
                   double **phaseComp,
                   domainInfo* simDomain, controls* simControls,
                   simParameters* simParams, subdomainInfo* subdomain,
                   dim3 gridSize, dim3 blockSize);

#endif
