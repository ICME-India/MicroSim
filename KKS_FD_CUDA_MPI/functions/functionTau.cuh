#ifndef FUNCTIONTAU_CUH_
#define FUNCTIONTAU_CUH_

#include "structures.h"
#include "matrix.cuh"
#include "utilityFunctions.h"

extern __device__ __host__
double FunctionTau(double **phi, double *relaxCoeff, long idx, long NUMPHASES);

void calculateTau(domainInfo *simDomain, controls *simControls, simParameters *simParams);

#endif
