#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "solverloop/defines.h"
#include "solverloop/defines1.h"

void F3_c_mu(double *mu, double *cp, double T, int phas, double *c_gues, struct propmatf3 *propf3); 

void F3_c_mu(double *mu, double *cp, double T, int phas, double *c_gues, struct propmatf3 *propf3) { 
  long k,l;
  double sum=0.0;
  for(l=0; l < NUMCOMPONENTS-1; l++) {
    sum = 0.0;
    for (k=0;k < NUMCOMPONENTS-1;k++) {
      sum += cmu[a][l][k]*(mu[k]-(Beq[a][k] + dBbdT[a][k]*(T-Teq)));
    }
    c[l] = sum;
  }
}

double function_C(double T, int ip);

double function_C(double T, int a) {

  double c_liq[nsol-1], c_sol[nsol-1];
  
  int i, j, k;
  double sum_s=0.0, sum_l=0.0, sum_c=0.0;
  
  if (a != (npha-1)) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      c_liq[k] = ceq[a][NUMPHASES-1][k] - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][NUMPHASES-1]);
      c_sol[k] = ceq[a][a][k]           - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][a]);
    }
    for (i=0; i < NUMCOMPONENTS-1; i++) {
      for (j=0; j < NUMCOMPONENTS-1; j++) {
        if (i <= j) {
          sum_c += (A[a][i][j]*c_sol[i]*c_sol[j] - A[NUMPHASES-1][i][j]*c_liq[i]*c_liq[j]);
        }
      }
    }
  }
  return sum_c;
}

