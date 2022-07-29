#ifndef UTILITYKERNELS_CUH_
#define UTILITYKERNELS_CUH_

#include <cuda.h>
#include <cuda_runtime.h>

#include "Thermo.cuh"
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

// __device__ __host__
// double evalFunc(void f(double, double*, double*), double x, double temperature);


__device__ __host__
double spline_eval(double x, double *controlPoints,
                   double *a, double *b, double *c, double *d,
                   int numControlPoints);

__global__
void computeChange(double *A, double *B,
                   int sizeX, int sizeY, int sizeZ);

void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize);

#endif
