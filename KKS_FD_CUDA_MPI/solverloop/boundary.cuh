#ifndef BOUNDARY_CUH_
#define BOUNDARY_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"

#define xm 1
#define xp 2
#define ym 3
#define yp 4
#define zm 5
#define zp 6

/*
 *  1. double **field - phi/comp/mu
 *  2. long face - Encode boundaries from 0 to 5, X+/X-/Y+/Y-/Z+/Z-
 *  3. long numFields - numPhases/numComponents-1/(numComponents-1)*numPhases [phi/comp,mu/phaseComp]
 *  4. sizeX, sizeY, sizeZ - Dimensions of the subdomain, including the boundary buffer layers
 */

__global__
void applyNeumann(double **field,
                  long face, long numFields, long DIMENSION,
                  long sizeX, long sizeY, long sizeZ, long padding);

__global__
void applyPeriodic(double **field,
                   long face, long numFields, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ, long padding);

#ifdef __cplusplus
extern "C"
#endif
void applyBoundaryCondition(double **field, long fieldCode, long numFields,
                            domainInfo* simDomain, controls* simControls,
                            simParameters* simParams, subdomainInfo* subdomain,
                            dim3 gridSize, dim3 blockSize);

#endif
