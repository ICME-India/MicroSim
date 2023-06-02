#ifndef SMOOTH_CUH_
#define SMOOTH_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "Thermo.cuh"
// #include "functionA_02.cuh"
#include "functionF.cuh"
#include "functionTau.cuh"
#include "utilityKernels.cuh"
#include "functionA_01.cuh"

__global__
void __smooth__(double **phi, double **phiNew,
                double *relaxCoeff, double *kappaPhi,
                double *dab, double *Rotation_matrix, double *Inv_rotation_matrix, int FUNCTION_ANISOTROPY,
                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                long xStep, long yStep, long padding,
                double DELTA_X, double DELTA_Y, double DELTA_Z,
                double DELTA_t);

void smooth(double **phi, double **phiNew,
            domainInfo* simDomain, controls* simControls,
            simParameters* simParams, subdomainInfo* subdomain,
            dim3 gridSize, dim3 blockSize);

#endif
