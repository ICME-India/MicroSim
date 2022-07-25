#ifndef UTILITYFUNCTIONS_H_
#define UTILITYFUNCTIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <string.h>

#include "structures.h"

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_matrix3M(double ***Mat, char *tmpstr, long NUMPHASES);
void populate_thetaij_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_thetai_matrix(double *Mat, char *tmpstr, long NUMPHASES);
void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_A_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_string_array(char **string, char *tmpstr, long size);

/*
 * Generate random number using user-specified seed.
 */
double ran(long *idum);

/*
 * Allocate a 2d double array
 */
double** malloc2M(int a, int b);

/*
 * Allocate a 3D double array
 */
double*** malloc3M(int a, int b, int c);

/*
 * Allocate a 4D double array
 */
double**** malloc4M(int a, int b, int c, int d);

/*
 * Free a 2D double array
 */
void free2M(double **Mat, int a);

/*
 * Free a 3D double array
 */
void free3M(double ***Mat, int a, int b);

/*
 * Free a 4D double array
 */
void free4M(double ****Mat, int a, int b, int c);

/*
 * Allocate memory on GPU and create pointers to make access to each phase/component simpler
 * arr is a continguous block containing information for every phase/component in a linear order
 * arr2 is an array of pointers pointing to the start of each phase/comp block in device memory
 * N is the number of phases/components, and stride is the size of each block (numCompCells)
 */
void allocOnDev(double **arr, double ***arr2, int N, int stride);

/*
 * De-allocate memory on GPU and memory taken by related helper pointers
 */
void freeOnDev(double **arr, double ***arr2);

/*
 * Aggregation of all free() calls.
 */
void freeVars(domainInfo *simDomain, simParameters *simParams);

#endif
