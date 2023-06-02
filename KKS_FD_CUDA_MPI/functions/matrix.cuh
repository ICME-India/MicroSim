#ifndef MATRIX_CUH_
#define MATRIX_CUH_

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "structures.h"
#include "utilityFunctions.h"

/*
 * Matrix operations, taken from the MicroSim GP solver's matrix.h
 */
void matinvnew(double **coeffmatrix,double **inv, long size);
void multiply(double **inv,double *y,double *prod, long size);
void multiply2d(double **m1,double **m2,double **prod, long size);
void vectorsum(double *y1,double *y2, double *sum, long size);
void substituteb(double **fac,double *y,double *vec, long size);
void substitutef(double **fac,double **y1,int index,double *vec, long size);
void pivot(double **coeffmatrix,double **factor,int k,int *tag, long size);
void colswap(double **m1,double **m2,int *tag, long size);
void rowswap(double **m1,double **m2,int *tag, long size);

/*
 *  The following are host-side functions
 */
int LUPDecompose(double **A, int N, double Tol, int *P);
void LUPInvert(double **A, int *P, int N, double **IA);
void matrixMultiply(double **A, double **B, double **C, int N);

/*
 *  The following are device-side functions
 */
extern __device__
void multiply(double *A, double *B, double *C, long ip1, long ip2, long NUMPHASES, int N);

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

#endif
