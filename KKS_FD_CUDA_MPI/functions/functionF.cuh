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

// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>

#include "utilityFunctions.h"
#include "structures.h"
#include "utilityKernels.cuh"
#include "matrix.cuh"
#include "Thermo.cuh"

/*
 * __device__ double calcPhaseEnergy
 *
 * Calculate f_{p}(c^{p}) = \sum_{i = 1}^{K-1}\sum_{j = 1}^{K-1} A_{ij}^{p}c_{i}^{p}c_{j}^{p}
 *                          + \sum_{i = 1}^{K-1} c_{i}^{p}
 *                          + C^{p}
 *
 * Arguments:
 *              1. double **phaseComp -> all the phase compositions, ordered phase-by-phase
 *              2. long phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. long idx -> position of cell in 1D
 *              7. long NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
extern __device__
double calcPhaseEnergy(double **phaseComp, long phase,
                       double *F0_A, double *F0_B, double *F0_C,
                       long idx,
                       long NUMPHASES, long NUMCOMPONENTS);

/*
 * __device__ double calcDiffusionPotential
 *
 * Calculate \frac{df_{p}}{dc^{p}_{c}}
 *
 * Arguments:
 *              1. double **phaseComp -> all the phase compositions, ordered phase-by-phase
 *              2. long phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. long idx -> position of cell in 1D
 *              7. long NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
extern __device__
double calcDiffusionPotential(double **phaseComp,
                              long phase, long component,
                              double *F0_A, double *F0_B,
                              long idx,
                              long NUMPHASES, long NUMCOMPONENTS);

void function_F_04_init_propertymatrices(domainInfo *simDomain, controls *simControls, simParameters *simParams);

void calcFreeEnergyCoeffs(domainInfo *simDomain, controls *simControls, simParameters *simParams);

#endif
