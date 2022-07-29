#ifndef COMPUTEDRIVINGFORCE_CUH_
#define COMPUTEDRIVINGFORCE_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "structures.h"
#include "utilityKernels.cuh"
#include "Thermo.cuh"

#ifndef MAX_NUM_PHASES
#define MAX_NUM_PHASES 5
#endif

#ifndef MAX_NUM_COMP
#define MAX_NUM_COMP 5
#endif

#ifndef MAX_NUM_PHASE_COMP
#define MAX_NUM_PHASE_COMP 16
#endif

/*
 * Explicit calculation of the right-hand side of the Allen-Cahn equation.
 * Evaluation of the mobility function in the Cahn-Hilliard equation.
 */
__global__ void __computeDrivingForce__(double **phi, double **comp,
                                        double **dfdphi,
                                        double **phaseComp,
                                        double *F0_A, double *F0_B, double *F0_C,
                                        double kappaPhi, double molarVolume,
                                        double *theta_i, double *theta_ij, double *theta_ijk,
                                        int NUMPHASES, int NUMCOMPONENTS,
                                        int sizeX, int sizeY, int sizeZ,
                                        double DELTA_X, double DELTA_Y, double DELTA_Z);

__global__
void __computeDrivingForce_02__(double **phi, double **comp,
                                 double **dfdphi,
                                 double **phaseComp,
                                 double molarVolume, double *kappaPhi,
                                 double *theta_i, double *theta_ij, double *theta_ijk,
                                 double temperature, int *thermo_phase,
                                 int NUMPHASES, int NUMCOMPONENTS,
                                 int sizeX, int sizeY, int sizeZ,
                                 double DELTA_X, double DELTA_Y, double DELTA_Z);

__global__
void __computeDrivingForceBinary_02__(double **phi, double **comp,
                                      double **dfdphi, double **mu,
                                      double **phaseComp,
                                      double molarVolume,
                                      double *theta_i, double *theta_ij, double *theta_ijk,
                                      double temperature, int *thermo_phase,
                                      int NUMPHASES, int NUMCOMPONENTS,
                                      int sizeX, int sizeY, int sizeZ);

__global__ void __computeDrivingForceBinary__(double **phi, double **comp,
                                                    double **dfdphi,
                                                    double **phaseComp,
                                                    double *F0_A, double *F0_B, double *F0_C,
                                                    double *diffusivity,
                                                    double *theta_i, double *theta_ij, double *theta_ijk,
                                                    int NUMPHASES, int NUMCOMPONENTS,
                                                    int sizeX, int sizeY, int sizeZ);

/*
 * Wrapper function for __computeDrivingForce__
 */
#ifdef __cplusplus
extern "C"
#endif
void computeDrivingForce(double **phi, double **comp,
                         double **dfdphi, double **mu,
                         double **phaseComp,
                         domainInfo* simDomain, controls* simControls,
                         simParameters* simParams, subdomainInfo* subdomain,
                         dim3 gridSize, dim3 blockSize);
#endif
