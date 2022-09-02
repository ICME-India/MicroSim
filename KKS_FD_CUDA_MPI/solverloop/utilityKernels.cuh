#ifndef UTILITYKERNELS_CUH_
#define UTILITYKERNELS_CUH_

#include <cuda.h>
#include <cuda_runtime.h>

#include "Thermo.cuh"
#include "structures.h"

#ifndef MAX_NUM_PHASES
#define MAX_NUM_PHASES 5
#endif

#ifndef MAX_NUM_COMP
#define MAX_NUM_COMP 5
#endif

#ifndef MAX_NUM_PHASE_COMP
#define MAX_NUM_PHASE_COMP 16
#endif

extern __device__
int LUPDecomposeC1(double A[][MAX_NUM_COMP], long N, double Tol, int *P);

extern __device__
int LUPDecomposeC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], long N, double Tol, int *P);

extern __device__
void LUPSolveC1(double A[][MAX_NUM_COMP], int *P, double *b, long N, double *x);

extern __device__
void LUPSolveC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, double *b, long N, double *x);

extern __device__
void LUPInvertC1(double A[][MAX_NUM_COMP], int *P, long N, double IA[][MAX_NUM_COMP]);

extern __device__
void LUPInvertC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, long N, double IA[][MAX_NUM_COMP]);

extern __device__
int LUPDecomposePC1(double A[][MAX_NUM_PHASE_COMP], long N, double Tol, int *P);

extern __device__
int LUPDecomposePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], long N, double Tol, int *P);

extern __device__
void LUPSolvePC1(double A[][MAX_NUM_PHASE_COMP], int *P, double *b, long N, double *x);

extern __device__
void LUPSolvePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, double *b, long N, double *x);

extern __device__
void LUPInvertPC1(double A[][MAX_NUM_PHASE_COMP], int *P, long N, double IA[][MAX_NUM_PHASE_COMP]);

extern __device__
void LUPInvertPC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, long N, double IA[][MAX_NUM_PHASE_COMP]);

extern __device__
double calcPhaseEnergy(double **phaseComp, long phase,
                       double *F0_A, double *F0_B, double *F0_C,
                       long idx,
                       long NUMPHASES, long NUMCOMPONENTS);

extern __device__
double calcDiffusionPotential(double **phaseComp,
                              long phase, long component,
                              double *F0_A, double *F0_B,
                              long idx,
                              long NUMPHASES, long NUMCOMPONENTS);

extern __device__
double FunctionTau(double **phi, double *relaxCoeff, long idx, long NUMPHASES);

__global__
void computeChange(double *A, double *B,
                   long sizeX, long sizeY, long sizeZ);

void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize);

#endif
