#ifndef FUNCTIONF_H_
#define FUNCTIONF_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <sys/stat.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "utilityFunctions.h"
#include "structures.h"
#include "utilityKernels.cuh"
#include "matrix.cuh"
#include "Thermo.cuh"


void function_F_04_init_propertymatrices(domainInfo *simDomain, controls *simControls, simParameters *simParams);

void calcFreeEnergyCoeffs(domainInfo *simDomain, controls *simControls, simParameters *simParams);

#endif
