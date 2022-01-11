#ifndef FUNCTION_F_01_H_
#define FUNCTION_F_01_H_

#include "global_vars.h"
double function_F_01_free_energy(double *c, double T, long a) {
  long i,j;
  double sum=0.0;
  for (i=0;i<NUMCOMPONENTS-1;i++) {
    for (j=0;j<NUMCOMPONENTS-1;j++) {
      if (i<=j) {
        sum += A[a][i][j]*c[i]*c[j];
      }
    }
    sum += B[a][i]*c[i];
  }
  sum += C[a];
  return sum;
}
double function_F_01_Mu(double *c, double T, long a, long i) {
  long j;
  double sum=0.0;
  sum  = 2.0*A[a][i][i]*c[i] + (Beq[a][i] + dBbdT[a][i]*(T-Teq));
  for (j=0; j<NUMCOMPONENTS-1;j++) {
    if (i!=j) {
      sum += A[a][i][j]*c[j];
    }
  }
  return sum;
}
double function_F_01_c_mu(double *mu, double T, long a, long i) {
  long k;
  double sum=0.0;
  for (k=0;k < NUMCOMPONENTS-1;k++) {
    sum += cmu[a][i][k]*(mu[k]-(Beq[a][k] + dBbdT[a][k]*(T-Teq)));
  }
  return sum;
}
double function_F_01_dc_dmu(double *mu, double T, long a, long i, long j) {
  return cmu[a][i][j];
}
double function_F_01_dpsi(double *mu, double T, double *phi, long a) {
  double psi=0.0;
  double c[NUMCOMPONENTS-1];
  long k=0, b;
  double sum=0.0;
  
  for (b=0; b < NUMPHASES; b++) {
    psi = 0.0;
    for (k=0;k <NUMCOMPONENTS-1; k++) {
      c[k] = c_mu(mu, T, b, k);
      psi -= mu[k]*c[k];
    }
    psi += free_energy(c,T,b);
    sum += psi*dhphi(phi, b, a);
  }
  return sum;
}
double function_F_01_function_A(double T1, long i, long j, long a) {
  return A[a][i][j];
}

double function_F_01_function_B(double T, long i, long a) {
  double c_liq[NUMCOMPONENTS-1];
  double c_sol[NUMCOMPONENTS-1];
  
  long k;
  double sum_s=0.0, sum_l=0.0, sum_c=0.0;
  double B_ai=0.0;
  
  if (a != (NUMPHASES-1)) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      c_liq[k] = ceq[a][NUMPHASES-1][k] - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][NUMPHASES-1]);
      c_sol[k] = ceq[a][a][k]           - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][a]);
      if (k!=i) {
        sum_c += A[NUMPHASES-1][k][i]*c_liq[k] - A[a][k][i]*c_sol[k];
      }
    }
    B_ai = (2.0*(A[NUMPHASES-1][i][i]*c_liq[i] - A[a][i][i]*c_sol[i]) + sum_c);
  }
  return B_ai;
}

double function_F_01_function_C(double T, long a) {
  double c_liq[NUMCOMPONENTS-1];
  double c_sol[NUMCOMPONENTS-1];
  
  long i, j, k;
  double sum_s=0.0, sum_l=0.0, sum_c=0.0;
  
  if (a != (NUMPHASES-1)) {
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
void function_F_01_compute_chemicalpotential(struct fields* gridinfo) {
  long x, y;
  long gidy;
  double c[NUMCOMPONENTS]; 
  long b, k, l;
  long a;
  long bulk_phase;
  double **dcdmu, **inv_dcdmu, *chempot, *sum;
    
  dcdmu     = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
  inv_dcdmu = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  chempot   = MallocV(NUMCOMPONENTS-1);
  sum       = MallocV(NUMCOMPONENTS-1);
  
  for(x=0;x<=(MESH_X-1);x++) {
    for (y=0; y<=(MESH_Y-1); y++) {
      gidy = x*MESH_Y + y;
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        sum[k] = 0.0;
        c[k]   = gridinfo[gidy].compi[k];
        for (l=0; l < NUMCOMPONENTS-1; l++) {
          dcdmu[k][l] = 0.0;
        }
      }
      interface = 1;
      for (a=0; a < NUMPHASES; a++) {
        if (gridinfo[gidy].phia[a] == 1.0) {
          bulk_phase=a;
          interface = 0;
          break;
        }
      }
      if (interface) {
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          sum[k] = 0.0;
          for (l=0; l < NUMCOMPONENTS-1; l++) {
            for (a=0; a < NUMPHASES; a++) {
              sum[k] += cmu[a][k][l]*B[a][l]*hphi(gridinfo[gidy].phia, a);
            }
          }
          sum[k] += c[k];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          for (l=0; l < NUMCOMPONENTS-1; l++) {
            for (a=0; a < NUMPHASES; a++) {
              dcdmu[k][l] += dc_dmu(c, T, a, k, l)*hphi(gridinfo[gidy].phia, a);
            }
          }
        }
        matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
        multiply(inv_dcdmu,sum,gridinfo[gidy].compi,NUMCOMPONENTS-1);
      } else {
        for (k=0; k < NUMCOMPONENTS-1; k++) {
         gridinfo[gidy].compi[k] = Mu(c, T, bulk_phase, k);
        }
      }
    }
  }
  FreeM(dcdmu,NUMCOMPONENTS-1);
  FreeM(inv_dcdmu,NUMCOMPONENTS-1);
  free(chempot);
  free(sum);
}

void init_propertymatrices(double T) {
  //Initialize property matrices
  long a, b, i, j, k;
    
//   for (a=0; a < NUMPHASES; a++) {
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       B[a][i] = function_B(T, i, a);
//     }
//     C[a] = function_C(T,a);
//   }  
#ifndef PHASE_DIAGRAM_PROP
#define PHASE_DIAGRAM_PROP
  for (a=0;a<NUMPHASES;a++) {
    for (i=0;i<NUMCOMPONENTS-1;i++) {
      for (j=0;j<NUMCOMPONENTS-1;j++) {
	if (i==j) {
	  muc[a][i][j]=2.0*A[a][i][j];
	} else {
	  muc[a][i][j]=A[a][i][j];
	}
      }
    }
    matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
  }

  for (a=0; a < NUMPHASES-1; a++) {
    DELTA_T[a][NUMPHASES-1] = 0.0;
    DELTA_T[a][a]           = 0.0;
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      DELTA_T[a][a]           += slopes[a][a][k]*(ceq[a][NUMPHASES-1][k]              - ceq[a][a][k]);
      DELTA_T[a][NUMPHASES-1] += slopes[a][NUMPHASES-1][k]*(ceq[a][NUMPHASES-1][k]    - ceq[a][a][k]);
      DELTA_C[a][k]            = ceq[a][NUMPHASES-1][k]   - ceq[a][a][k];
    }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
     dcbdT[a][a][k]           = DELTA_C[a][k]/DELTA_T[a][a];
     dcbdT[a][NUMPHASES-1][k] = DELTA_C[a][k]/DELTA_T[a][NUMPHASES-1];
    }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      dBbdT[a][k]      = 2.0*(A[NUMPHASES-1][k][k]*dcbdT[a][NUMPHASES-1][k] - A[a][k][k]*dcbdT[a][a][k]);
      for (i=0; i < NUMCOMPONENTS-1; i++) {
        if (k!=i) {
          dBbdT[a][k] += (A[NUMPHASES-1][k][i]*dcbdT[a][NUMPHASES-1][i] - A[a][k][i]*dcbdT[a][a][i]);
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    dBbdT[NUMPHASES-1][k] = 0.0;
  }
  for (a=0; a < NUMPHASES; a++) {
    for (i=0; i < NUMCOMPONENTS-1; i++) {
      Beq[a][i] = function_B(Teq, i, a); 
    }
  }
#endif
  for (a=0; a < NUMPHASES; a++) {
    for (i=0; i < NUMCOMPONENTS-1; i++) {
      B[a][i] = function_B(T, i, a); 
    }
    C[a] = function_C(T,a);
  }
}

#endif