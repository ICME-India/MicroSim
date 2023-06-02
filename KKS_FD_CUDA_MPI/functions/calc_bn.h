#ifndef CALC_BN_HPP
#define CALC_BN_HPP

#ifndef ENABLE_CUFFTMP
#define ENABLE_CUFFTMP 0
#endif

#if ENABLE_CUFFTMP == 1

#include "structures.h"
#include <math.h>
#include <stdio.h>

/***********************************************************************
 *       Eigenstrain
 ***********************************************************************/
void calc_strn(double eigen_strn[3][3], double eps11, double eps22,
               double eps33, double eps12, double eps13, double eps23);

/*************************************************************
 *       Function to translate C[6][6] to lambda[3][3][3][3] *
 *************************************************************/
void eval_lambda(double C[6][6], double lambda[3][3][3][3]);

/***********************************************************************
 * Function to calculate inverse of omega, which is normalized G_inverse
 * Khachaturyan calls G as Fourier transform of Green function of anisotropic
 * elasticity
 ***********************************************************************/
void calc_inv_omega(double kx, double ky, double kz,
                    double lambda[3][3][3][3], double inv_omega[3][3]);

/*******************************************************************
 * This function just inverts a three by three matrix
 *******************************************************************/
void calc_inverse(double a[3][3], double ainv[3][3]);

/************************************************************
 * This is linear elastic Hooke's law. Here, we are
 * interested in calculating the stress-like quantity calculated
 * from the eigen strain.
 ************************************************************/
void calc_stress(double eigen_strn[3][3], double lambda[3][3][3][3],
                 double stress[3][3]);

/****************************************************************
 * Given the components of the Fourier vector, the following function
 * evaluates B, which is the Fourier transform of the strain induced
 * energy
 ****************************************************************/
double elast(double kx, double ky, double kz, double omega[3][3], double
eigen_strn1[3][3], double eigen_strn2[3][3], double stress1[3][3],
double stress2[3][3], double lambda[3][3][3][3]);

void calculate_Bn(double *B, double *kx, double *ky, double *kz, domainInfo simDomain, simParameters simParams, subdomainInfo subdomain, long p, long q);

#endif
#endif
