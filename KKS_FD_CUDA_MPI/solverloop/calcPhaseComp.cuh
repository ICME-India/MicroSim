#ifndef CALCPHASECOMP_CUH_
#define CALCPHASECOMP_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "utilityKernels.cuh"
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
 * Calculate phase compositions
 */
__global__
void __calcPhaseComp___(double **phi, double **comp,
                             double **phaseComp,
                             double *F0_A, double *F0_B, double *F0_C,
                             int NUMPHASES, int NUMCOMPONENTS,
                             int sizeX, int sizeY, int sizeZ);

__global__
void __calcPhaseCompBinary___(double **phi, double **comp,
                                   double **phaseComp,
                                   double *F0_A, double *F0_B, double *F0_C,
                                   int NUMPHASES, int NUMCOMPONENTS,
                                   int sizeX, int sizeY, int sizeZ);


__global__
void __initMu__(double **comp, double **mu,
                double molarVolume, int *thermo_phase, double temperature,
                int NUMPHASES, int NUMCOMPONENTS,
                int sizeX, int sizeY, int sizeZ);

__global__
void __calcPhaseComp_02__(double **phi, double **comp,
                          double **phaseComp, double **mu, double *cguess,
                          double molarVolume,
                          double temperature, int *thermo_phase,
                          int NUMPHASES, int NUMCOMPONENTS,
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
                   double **phaseComp, double **mu,
                   domainInfo* simDomain, controls* simControls,
                   simParameters* simParams, subdomainInfo* subdomain,
                   dim3 gridSize, dim3 blockSize);

#endif
