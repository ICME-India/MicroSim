#ifndef UTILITYKERNELS_CUH_
#define UTILITYKERNELS_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include "helper_string.h"
#include "helper_cuda.h"

#include "Thermo.cuh"
#include "structures.h"

__global__
void __computeChange__(double *A, double *B, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ, long padding);

__global__
void __resetArray__(double **arr, long numArr,
                    long sizeX, long sizeY, long sizeZ);

void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize);

void resetArray(double **arr, long numArr,
                long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                dim3 gridSize, dim3 blockSize);

#endif
