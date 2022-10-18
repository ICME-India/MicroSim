#ifndef UTILITYKERNELS_CUH_
#define UTILITYKERNELS_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include "helper_string.h"
#include "helper_cuda.h"

#include "Thermo.cuh"
#include "structures.h"

extern __device__
double calcPhaseEnergy(double **phaseComp, long phase,
                       double *F0_A, double *F0_B, double *F0_C,
                       long idx,
                       long NUMPHASES, long NUMCOMPONENTS);

extern __device__
double calcDiffusionPotential(double **phaseComp,
                              long phase, long component,
                              double *F0_A, double *F0_B,
                              long idx,
                              long NUMPHASES, long NUMCOMPONENTS);

__global__
void computeChange(double *A, double *B, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ);

void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize);

#endif
