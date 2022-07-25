#ifndef UTILITYKERNELS_CUH_
#define UTILITYKERNELS_CUH_

#include <cuda.h>
#include <cuda_runtime.h>

#include "Thermo.cuh"
#include "structures.h"

/*
 * __device__ __host__ double evalFunc
 *
 * Evaluates thermo writer functions of the specified form at the given composition and temperature
 * Can be called from both host-side code and device-side code
 *
 * Arguments:
 *              1. void f(double, double*, double*) -> free-energy function and derivatives as written by SymPy codegen
 *              2. double x -> specifies c[1], as Al-Zn test simulations are being run using X_{ZN}
 *              3. double temperature -> temperature
 *
 * Return:
 *          the value of the function, as a double
 *
 */
__device__ __host__
double evalFunc(void f(double, double*, double*), double x, double temperature);


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
