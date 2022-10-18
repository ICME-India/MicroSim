#ifndef ANISOTROPY_01_CUH_
#define ANISOTROPY_01_CUH_

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

extern __device__
void anisotropy_01_dAdq(double *qab, double *dadq, long a, long b, double *dab, long NUMPHASES);

extern __device__
double anisotropy_01_function_ac(double *qab, long a, long b, double *dab, long NUMPHASES);

#endif
