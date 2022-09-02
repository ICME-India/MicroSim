#ifndef INITIALIZE_VARIABLES_H_
#define INITIALIZE_VARIABLES_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "structures.h"
#include "utilityFunctions.h"

/*
 * Allocate memory on the CPU and GPU for all the phase-field variables
 * and other required variables.
 */
void decomposeDomain(domainInfo simDomain, controls *simControls, subdomainInfo *subdomain,
                     int rank, int size);

/*
 * Initialize variables and move to GPU
 */
void moveParamsToGPU(domainInfo *simDomain, controls *simControls, simParameters *simParams);

/*
 * Calculate grid size and block size for kernel calls
 */
void calcKernelParams(dim3 *gridSize, dim3 *blockSize, domainInfo simDomain, controls simControls, subdomainInfo *subdomain);

#endif
