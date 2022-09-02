#ifndef FUNCTIONH_CUH_
#define FUNCTIONH_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

extern __device__
double calcInterp5th(double **phi, long a, long idx, long NUMPHASES);

extern __device__
double calcInterp5thDiff(double **phi, long a, long b, long idx, long NUMPHASES);

extern __device__
double calcInterp3rd(double **phi, long a, long idx, long NUMPHASES);

extern __device__
double calcInterp3rdDiff(double **phi, long a, long b, long idx, long NUMPHASES);

#endif
