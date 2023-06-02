#ifndef UTILITYFUNCTIONS_H_
#define UTILITYFUNCTIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <string.h>

#ifndef THERMO
#define THERMO 1
#endif

#if THERMO == 1
#include "Thermo.cuh"
#endif
#include "structures.h"

#include "matrix.cuh"

void testThermoFuncs(domainInfo simDomain, simParameters simParams);

void get_Rotation_Matrix(double **R, double theta, int axis);

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_matrix3M(double ***Mat, char *tmpstr, long NUMPHASES);
void populate_thetaij_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_thetai_matrix(double *Mat, char *tmpstr, long NUMPHASES);
void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_A_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_string_array(char **string, char *tmpstr, long size);
void populate_rotation_matrix(double ****Mat, double ****Mat_Inv, char *tmpstr);
void populate_symmetric_tensor(symmetric_tensor *Mat, char *tmpstr, long NUMPHASES);
void populate_cubic_stiffness(Stiffness_cubic *Mat, char *tmpstr);

/*
 * Generate random number using user-specified seed.
 */
double ran(long *idum);

/*
 *  Allocate a 1d double array
 */
double* MallocV(long m);

/*
 * Allocate a 2d double array
 */
double** MallocM(long a, long b);

/*
 * Allocate a 3D double array
 */
double*** Malloc3M(long a, long b, long c);

/*
 * Allocate a 4D double array
 */
double**** Malloc4M(long a, long b, long c, long d);

/*
 * Free a 2D double array
 */
void FreeM(double **Mat, long a);

/*
 * Free a 3D double array
 */
void Free3M(double ***Mat, long a, long b);

/*
 * Free a 4D double array
 */
void Free4M(double ****Mat, long a, long b, long c);

/*
 * Allocate memory on GPU and create pointers to make access to each phase/component simpler
 * arr is a continguous block containing information for every phase/component in a linear order
 * arr2 is an array of pointers pointing to the start of each phase/comp block in device memory
 * N is the number of phases/components, and stride is the size of each block (numCompCells)
 */
void allocOnDev(double **arr, double ***arr2, long N, long stride);

/*
 * De-allocate memory on GPU and memory taken by related helper pointers
 */
void freeOnDev(double **arr, double ***arr2);

/*
 * Aggregation of all free() calls.
 */
void freeVars(domainInfo *simDomain, controls *simControls, simParameters *simParams);

#endif
