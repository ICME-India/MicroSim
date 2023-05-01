#ifndef FUNCTION_F_05_H_
#define FUNCTION_F_05_H_

#include "global_vars.h"

// double function_F_03_ComputeSlopes(double T, long a);

// double function_F_05_free_energy(double *c, double T, long a) {
//   long i,j;
//   double sum=0.0;
//   for (i=0;i<NUMCOMPONENTS-1;i++) {
//     for (j=0;j<NUMCOMPONENTS-1;j++) {
//       if (i<=j) {
//         sum += A[a][i][j]*c[i]*c[j];
//       }
//     }
//     if (!ISOTHERMAL) {
//       B[a][i] = (Beq[a][i] + dBbdT[a][i]*(T-Teq));
//     }
//     sum += B[a][i]*c[i];
//   }
//   if (!ISOTHERMAL) {
//     C[a] = function_C(T,a);
//   }
//   sum += C[a];
//   return sum;
// }
// void function_F_03_Mu(double *c, double T, long a, double *Mu) {
//   long j,k;
//   double sum=0.0;
//   for(k=0; k < NUMCOMPONENTS-1; k++) {
//     sum  = 2.0*A[a][k][k]*c[k] + (Beq[a][k] + dBbdT[a][k]*(T-Teq));
//     for (j=0; j<NUMCOMPONENTS-1;j++) {
//       if (k!=j) {
//         sum += A[a][k][j]*c[j];
//       }
//     }
//     Mu[k] = sum;
//   }
// }
// void function_F_03_c_mu(double *mu, double *c, double T, long a, double *c_guess) {
//   long k,l;
//   double sum=0.0;
//   for(l=0; l < NUMCOMPONENTS-1; l++) {
//     sum = 0.0;
//     for (k=0;k < NUMCOMPONENTS-1;k++) {
//       sum += cmu[a][l][k]*(mu[k]-(Beq[a][k] + dBbdT[a][k]*(T-Teq)));
//     }
//     c[l] = sum;
//   }
// }
// void function_F_03_dc_dmu(double *mu, double *phase_comp, double T, long a, double **dcdmu) {
//   long i, j;
//   for (i=0; i<NUMCOMPONENTS-1; i++) {
//     for (j=0; j<NUMCOMPONENTS-1; j++) {
//       dcdmu[i][j] = cmu[a][i][j];
//     }
//   }
// }
double function_F_05_dpsi(double *mu, double **phase_comp, double T, double *phi, long a) {
  long k=0, b;
  double sum=0.0;
  sum = -(Lf*(T-Teq)/Teq)*dhphi(phi, NUMPHASES-1, a);
  return sum;
}
// void function_F_03_function_A(double T, double ***c) {
//   double dmudc[(NUMCOMPONENTS-1)*(NUMCOMPONENTS-1)];
//   long i,j,k, index;
//   double **dmudc_matrix;
//   double sum=0.0;
//   
//   double y[NUMCOMPONENTS];
//   dmudc_matrix = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
//   
//   for (a=0; a < NUMPHASES; a++) {
//     sum = 0.0;
//     for (k=0; k<NUMCOMPONENTS-1; k++) {
//       y[k] = c[a][a][k];
//       sum += y[k];
//     }
//     y[NUMCOMPONENTS-1] = 1.0 - sum;
//     
//     (*dmudc_tdb[thermo_phase[a]])(T, y, dmudc);
//     
//     for(i=0; i<NUMCOMPONENTS-1; i++) {
//       for(j=0; j<NUMCOMPONENTS-1; j++) {
//         index = i*(NUMCOMPONENTS-1) + j;
//         dmudc_matrix[i][j] = dmudc[index];
//       }
//     }
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       for (j=0; j < NUMCOMPONENTS-1; j++) {
//         if (i==j) {
//           A[a][i][j] = 0.5*dmudc_matrix[i][j];
//         } else {
//           A[a][i][j] = dmudc_matrix[i][j];
//         }
//       }
//     }
//   }
//   FreeM(dmudc_matrix, NUMCOMPONENTS-1);
// }
// 
// double function_F_03_ComputeSlopes(double T, long a) {
// //   printf("I am coming to compute slopes\n");
//   double dmudc[(NUMCOMPONENTS-1)*(NUMCOMPONENTS-1)];
//   long k,l,index;
//   double **dmudc_matrix;
//   double sum=0.0;
//   double sum_cs[NUMCOMPONENTS-1];
//   double sum_cl[NUMCOMPONENTS-1];
//   double mu_s[NUMCOMPONENTS-1];
//   double dmu_s[NUMCOMPONENTS-1];
//   double mu_l[NUMCOMPONENTS-1];
//   double dmu_l[NUMCOMPONENTS-1];
//   double f_s, f_l;
//   double df_s, df_l;
//   double DT = 1;
//   double dS;
//   double dmuS, dmuL;
//   
//   double y[NUMCOMPONENTS];
//   
//   for (k=0; k<NUMCOMPONENTS-1; k++) {
//     y[k] = ceq[a][NUMPHASES-1][k];
//     sum += y[k];
//     sum_cl[k] = 0.0;
//   }
//   y[NUMCOMPONENTS-1] = 1.0 - sum;
//   (*free_energy_tdb[thermo_phase[NUMPHASES-1]])(T - DT, y, &f_l);
//   (*free_energy_tdb[thermo_phase[NUMPHASES-1]])(T + DT, y, &df_l);
//   
//   (*Mu_tdb[thermo_phase[NUMPHASES-1]])(T-DT, y, mu_l);
//   (*Mu_tdb[thermo_phase[NUMPHASES-1]])(T+DT, y, dmu_l);
//   
//   
//   sum = 0.0;
//   for (k=0; k<NUMCOMPONENTS-1; k++) {
//     y[k] = ceq[a][a][k];
//     sum += y[k];
//     sum_cs[k] = 0.0;
//   }
//   y[NUMCOMPONENTS-1] = 1.0 - sum;
//   
//   (*free_energy_tdb[thermo_phase[a]])(T - DT, y, &f_s);
//   (*free_energy_tdb[thermo_phase[a]])(T + DT, y, &df_s);
//   
//   (*Mu_tdb[thermo_phase[a]])(T-DT, y, mu_s);
//   (*Mu_tdb[thermo_phase[a]])(T+DT, y, dmu_s);
//   
//   dS   = 0.5*(df_s - f_s)/DT - 0.5*(df_l - f_l)/DT;
//   
//   dmuS = 0.0;
//   dmuL = 0.0;
//   
//   for (k=0; k<NUMCOMPONENTS-1; k++) {
//     dmuS += 0.5*(dmu_s[k] - mu_s[k])*(ceq[a][a][k] - ceq[a][NUMPHASES-1][k])/DT;
//     dmuL += 0.5*(dmu_l[k] - mu_l[k])*(ceq[a][a][k] - ceq[a][NUMPHASES-1][k])/DT;
//   }
//   
//   for (k=0; k<NUMCOMPONENTS-1; k++) {
//     for (l=0; l<NUMCOMPONENTS-1; l++) {
//       if (k==l) {
//         sum_cs[k] += 2.0*(ceq[a][a][l] - ceq[a][NUMPHASES-1][l])*A[a][l][k];
//         sum_cl[k] += 2.0*(ceq[a][a][l] - ceq[a][NUMPHASES-1][l])*A[NUMPHASES-1][l][k];
//       } else {
//         sum_cs[k] += (ceq[a][a][l] - ceq[a][NUMPHASES-1][l])*A[a][l][k];
//         sum_cl[k] += (ceq[a][a][l] - ceq[a][NUMPHASES-1][l])*A[NUMPHASES-1][l][k];
//       }
//     }
//     slopes[a][a][k]           = sum_cs[k]/(dS - dmuS);
//     slopes[a][NUMPHASES-1][k] = sum_cl[k]/(dS - dmuL);
//     slopes[NUMPHASES-1][a][k] = slopes[a][NUMPHASES-1][k];
//   }
// }
// 
// double function_F_03_function_B(double T, long i, long a) {
//   double c_liq[NUMCOMPONENTS-1];
//   double c_sol[NUMCOMPONENTS-1];
//   
//   long k;
//   double sum_s=0.0, sum_l=0.0, sum_c=0.0;
//   double B_ai=0.0;
//   
//   if (a != (NUMPHASES-1)) {
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       c_liq[k] = ceq[a][NUMPHASES-1][k] - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][NUMPHASES-1]);
//       c_sol[k] = ceq[a][a][k]           - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][a]);
//       if (k!=i) {
//         sum_c += A[NUMPHASES-1][k][i]*c_liq[k] - A[a][k][i]*c_sol[k];
//       }
//     }
//     B_ai = (2.0*(A[NUMPHASES-1][i][i]*c_liq[i] - A[a][i][i]*c_sol[i]) + sum_c);
//   }
//   return B_ai;
// }
// 
// double function_F_03_function_C(double T, long a) {
//   double c_liq[NUMCOMPONENTS-1];
//   double c_sol[NUMCOMPONENTS-1];
//   
//   long i, j, k;
//   double sum_s=0.0, sum_l=0.0, sum_c=0.0;
//   
//   if (a != (NUMPHASES-1)) {
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       c_liq[k] = ceq[a][NUMPHASES-1][k] - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][NUMPHASES-1]);
//       c_sol[k] = ceq[a][a][k]           - (DELTA_C[a][k])*(Teq-T)/(DELTA_T[a][a]);
//     }
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       for (j=0; j < NUMCOMPONENTS-1; j++) {
//         if (i <= j) {
//           sum_c += (A[a][i][j]*c_sol[i]*c_sol[j] - A[NUMPHASES-1][i][j]*c_liq[i]*c_liq[j]);
//         }
//       }
//     }
//   }
//   return sum_c;
// }
// void function_F_03_init_propertymatrices(double T) {
//   //Initialize property matrices
//   long a, b, i, j, k;
//   function_A((T+Teq)*0.5, ceq);
//   for (a=0;a<NUMPHASES;a++) {
//     if (a < (NUMPHASES-1)) {
//       function_F_03_ComputeSlopes((T+Teq)*0.5, a);
//     }
//   }
//   function_A(T, c_guess);
//   for (a=0;a<NUMPHASES;a++) {
//     for (i=0;i<NUMCOMPONENTS-1;i++) {
//       for (j=0;j<NUMCOMPONENTS-1;j++) {
// 	if (i==j) {
// 	  muc[a][i][j]=2.0*A[a][i][j];
// 	} else {
// 	  muc[a][i][j]=A[a][i][j];
// 	}
//       }
//     }
//     matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
//   }
//   for (a=0; a < NUMPHASES-1; a++) {
//     DELTA_T[a][NUMPHASES-1] = 0.0;
//     DELTA_T[a][a]           = 0.0;
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       DELTA_T[a][a]           += slopes[a][a][k]*(ceq[a][NUMPHASES-1][k]              - ceq[a][a][k]);
//       DELTA_T[a][NUMPHASES-1] += slopes[a][NUMPHASES-1][k]*(ceq[a][NUMPHASES-1][k]    - ceq[a][a][k]);
//       DELTA_C[a][k]            = ceq[a][NUMPHASES-1][k]   - ceq[a][a][k];
//     }
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//      dcbdT[a][a][k]           = DELTA_C[a][k]/DELTA_T[a][a];
//      dcbdT[a][NUMPHASES-1][k] = DELTA_C[a][k]/DELTA_T[a][NUMPHASES-1];
//     }
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       dBbdT[a][k]      = 2.0*(A[NUMPHASES-1][k][k]*dcbdT[a][NUMPHASES-1][k] - A[a][k][k]*dcbdT[a][a][k]);
//       for (i=0; i < NUMCOMPONENTS-1; i++) {
//         if (k!=i) {
//           dBbdT[a][k] += (A[NUMPHASES-1][k][i]*dcbdT[a][NUMPHASES-1][i] - A[a][k][i]*dcbdT[a][a][i]);
//         }
//       }
//     }
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     dBbdT[NUMPHASES-1][k] = 0.0;
//   }
//   for (a=0; a < NUMPHASES; a++) {
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       Beq[a][i] = function_B(Teq, i, a); 
//     }
//   }
//   for (a=0; a < NUMPHASES; a++) {
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       B[a][i] = function_B(T, i, a); 
//     }
//     C[a] = function_C(T,a);
//   }
// }

// void function_F_01_compute_chemicalpotential(struct fields* gridinfo) {
//   long x, y;
//   long gidy;
//   double c[NUMCOMPONENTS]; 
//   long b, k, l;
//   long a;
//   long bulk_phase;
//   double **dcdmu, **inv_dcdmu, *chempot, *sum;
//     
//   dcdmu     = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
//   inv_dcdmu = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
//   chempot   = MallocV(NUMCOMPONENTS-1);
//   sum       = MallocV(NUMCOMPONENTS-1);
//   
//   for(x=0;x<=(MESH_X-1);x++) {
//     for (y=0; y<=(MESH_Y-1); y++) {
//       gidy = x*MESH_Y + y;
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
//         sum[k] = 0.0;
//         c[k]   = gridinfo[gidy].compi[k];
//         for (l=0; l < NUMCOMPONENTS-1; l++) {
//           dcdmu[k][l] = 0.0;
//         }
//       }
//       interface = 1;
//       for (a=0; a < NUMPHASES; a++) {
//         if (gridinfo[gidy].phia[a] == 1.0) {
//           bulk_phase=a;
//           interface = 0;
//           break;
//         }
//       }
//       if (interface) {
//         for (k=0; k < NUMCOMPONENTS-1; k++) {
//           sum[k] = 0.0;
//           for (l=0; l < NUMCOMPONENTS-1; l++) {
//             for (a=0; a < NUMPHASES; a++) {
//               sum[k] += cmu[a][k][l]*B[a][l]*hphi(gridinfo[gidy].phia, a);
//             }
//           }
//           sum[k] += c[k];
//         }
//         for (k=0; k < NUMCOMPONENTS-1; k++) {
//           for (l=0; l < NUMCOMPONENTS-1; l++) {
//             for (a=0; a < NUMPHASES; a++) {
//               dc_dmu(double *mu, double *phase_comp, double T, long a, double **dcdmu)
//               dcdmu[k][l] += dc_dmu(c, T, a, k, l)*hphi(gridinfo[gidy].phia, a);
//             }
//           }
//         }
//         matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
//         multiply(inv_dcdmu,sum,gridinfo[gidy].compi,NUMCOMPONENTS-1);
//       } else {
//         for (k=0; k < NUMCOMPONENTS-1; k++) {
//          gridinfo[gidy].compi[k] = Mu(c, T, bulk_phase, k);
//         }
//       }
//     }
//   }
//   FreeM(dcdmu,NUMCOMPONENTS-1);
//   FreeM(inv_dcdmu,NUMCOMPONENTS-1);
//   free(chempot);
//   free(sum);
// }

// void init_propertymatrices(double T) {
//   //Initialize property matrices
//   long a, b, i, j, k;
//     
// //   for (a=0; a < NUMPHASES; a++) {
// //     for (i=0; i < NUMCOMPONENTS-1; i++) {
// //       B[a][i] = function_B(T, i, a);
// //     }
// //     C[a] = function_C(T,a);
// //   }  
// #ifndef PHASE_DIAGRAM_PROP
// #define PHASE_DIAGRAM_PROP
//   for (a=0;a<NUMPHASES;a++) {
//     for (i=0;i<NUMCOMPONENTS-1;i++) {
//       for (j=0;j<NUMCOMPONENTS-1;j++) {
// 	if (i==j) {
// 	  muc[a][i][j]=2.0*A[a][i][j];
// 	} else {
// 	  muc[a][i][j]=A[a][i][j];
// 	}
//       }
//     }
//     matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
//   }
// 
//   for (a=0; a < NUMPHASES-1; a++) {
//     DELTA_T[a][NUMPHASES-1] = 0.0;
//     DELTA_T[a][a]           = 0.0;
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       DELTA_T[a][a]           += slopes[a][a][k]*(ceq[a][NUMPHASES-1][k]              - ceq[a][a][k]);
//       DELTA_T[a][NUMPHASES-1] += slopes[a][NUMPHASES-1][k]*(ceq[a][NUMPHASES-1][k]    - ceq[a][a][k]);
//       DELTA_C[a][k]            = ceq[a][NUMPHASES-1][k]   - ceq[a][a][k];
//     }
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//      dcbdT[a][a][k]           = DELTA_C[a][k]/DELTA_T[a][a];
//      dcbdT[a][NUMPHASES-1][k] = DELTA_C[a][k]/DELTA_T[a][NUMPHASES-1];
//     }
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       dBbdT[a][k]      = 2.0*(A[NUMPHASES-1][k][k]*dcbdT[a][NUMPHASES-1][k] - A[a][k][k]*dcbdT[a][a][k]);
//       for (i=0; i < NUMCOMPONENTS-1; i++) {
//         if (k!=i) {
//           dBbdT[a][k] += (A[NUMPHASES-1][k][i]*dcbdT[a][NUMPHASES-1][i] - A[a][k][i]*dcbdT[a][a][i]);
//         }
//       }
//     }
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     dBbdT[NUMPHASES-1][k] = 0.0;
//   }
//   for (a=0; a < NUMPHASES; a++) {
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       Beq[a][i] = function_B(Teq, i, a); 
//     }
//   }
// #endif
//   for (a=0; a < NUMPHASES; a++) {
//     for (i=0; i < NUMCOMPONENTS-1; i++) {
//       B[a][i] = function_B(T, i, a); 
//     }
//     C[a] = function_C(T,a);
//   }
// }

#endif
