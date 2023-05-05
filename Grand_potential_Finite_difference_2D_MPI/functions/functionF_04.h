#ifndef FUNCTION_F_04_H_
#define FUNCTION_F_04_H_

#include "global_vars.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// double function_F_03_ComputeSlopes(double T, long a);

double function_F_04_free_energy(double *c, double T, long a) {
  long i,j;
  double sum=0.0;
  
  if(!ISOTHERMAL) {
    for(i=0; i<NUMCOMPONENTS-1; i++) {
      for(j=0; j<NUMCOMPONENTS-1; j++) {
        if (i==j) {
          A[a][i][j] = 0.5*gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        } else {
          A[a][i][j] = gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        }
      }
    }
  }
  for (i=0;i<NUMCOMPONENTS-1;i++) {
    for (j=0;j<NUMCOMPONENTS-1;j++) {
      if (i<=j) {
        sum += A[a][i][j]*c[i]*c[j];
      }
    }
    if (!ISOTHERMAL) {
      B[a][i] = function_B(T,i,a);
    }
    sum += B[a][i]*c[i];
  }
  if (!ISOTHERMAL) {
    C[a] = function_C(T,a);
  }
  sum += C[a];
  return sum;
}
void function_F_04_Mu(double *c, double T, long a, double *Mu) {
  long j,k;
  double sum=0.0;
  
  if(!ISOTHERMAL) {
   for(i=0; i<NUMCOMPONENTS-1; i++) {
      for(j=0; j<NUMCOMPONENTS-1; j++) {
        if (i==j) {
          A[a][i][j] = 0.5*gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        } else {
          A[a][i][j] = gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        }
      }
    }
  }
  
  if(!ISOTHERMAL) {
    for (k=0;k < NUMCOMPONENTS-1;k++) {
      B[a][k] = function_B(T,k,a);
    }
  }
  for(k=0; k < NUMCOMPONENTS-1; k++) {
    sum  = 2.0*A[a][k][k]*c[k] + B[a][k];
    for (j=0; j<NUMCOMPONENTS-1;j++) {
      if (k!=j) {
        sum += A[a][k][j]*c[j];
      }
    }
    Mu[k] = sum;
  }
}
void function_F_04_c_mu(double *mu, double *c, double T, long a, double *c_guess) {
  long k,l;
  double sum=0.0;
  if(!ISOTHERMAL) {
    for(i=0; i<NUMCOMPONENTS-1; i++) {
      for(j=0; j<NUMCOMPONENTS-1; j++) {
        if (i==j) {
          A[a][i][j] = 0.5*gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        } else {
          A[a][i][j] = gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        }
      }
    }
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
  if(!ISOTHERMAL) {
    for (k=0;k < NUMCOMPONENTS-1;k++) {
       B[a][k] = function_B(T,k,a);
    }
  }
  for(l=0; l < NUMCOMPONENTS-1; l++) {
    sum = 0.0;
    for (k=0;k < NUMCOMPONENTS-1;k++) {
      sum += cmu[a][l][k]*(mu[k]-B[a][k]);
    }
    c[l] = sum;
  }
}
void function_F_04_dc_dmu(double *mu, double *phase_comp, double T, long a, double **dcdmu) {
  long i, j;
  if(!ISOTHERMAL) {
    for(i=0; i<NUMCOMPONENTS-1; i++) {
      for(j=0; j<NUMCOMPONENTS-1; j++) {
        if (i==j) {
          A[a][i][j] = 0.5*gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        } else {
          A[a][i][j] = gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        }
      }
    }
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
  for (i=0; i<NUMCOMPONENTS-1; i++) {
    for (j=0; j<NUMCOMPONENTS-1; j++) {
      dcdmu[i][j] = cmu[a][i][j];
    }
  }
}
double function_F_04_dpsi(double *mu, double **phase_comp, double T, double *phi, long a) {
  double psi=0.0;
  double c[NUMCOMPONENTS-1];
  long k=0, b;
  double sum=0.0;
  
  for (b=0; b < NUMPHASES; b++) {
    psi = 0.0;
    for (k=0;k <NUMCOMPONENTS-1; k++) {
      c[k] = phase_comp[b][k];
      psi -= mu[k]*c[k];
    }
    psi += free_energy(c,T,b);
    sum += psi*dhphi(phi, b, a);
  }
  return sum;
}
void function_F_04_function_A(double T, double ***c) {
  long i,j,k, index;
  double **dmudc_matrix;
  double sum=0.0;
  
  double y[NUMCOMPONENTS];  
  for (a=0; a < NUMPHASES; a++) {
    for(i=0; i<NUMCOMPONENTS-1; i++) {
      for(j=0; j<NUMCOMPONENTS-1; j++) {
        if (i==j) {
          A[a][i][j] = 0.5*gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        } else {
          A[a][i][j] = gsl_spline_eval (spline_ThF[thermo_phase[a]][i][j], T, acc_ThF[thermo_phase[a]][i][j]);
        }
      }
    }
  }
}

// double function_F_04_ComputeSlopes(double T, long a) {
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

double function_F_04_function_B(double T, long i, long a) {
  double c_liq[NUMCOMPONENTS-1];
  double c_sol[NUMCOMPONENTS-1];
  
  long k;
  double sum_s=0.0, sum_l=0.0, sum_c=0.0;
  double B_ai=0.0;
  if (a != (NUMPHASES-1)) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      c_liq[k] = gsl_spline_eval (spline_ES[thermo_phase[a]][k][1], T, acc_ES[thermo_phase[a]][k][1]);
      c_sol[k] = gsl_spline_eval (spline_ES[thermo_phase[a]][k][0], T, acc_ES[thermo_phase[a]][k][0]);
      if (k!=i) {
        sum_c += A[NUMPHASES-1][k][i]*c_liq[k] - A[a][k][i]*c_sol[k];
      }
    }
    B_ai = (2.0*(A[NUMPHASES-1][i][i]*c_liq[i] - A[a][i][i]*c_sol[i]) + sum_c);
  }
  return B_ai;
}

double function_F_04_function_C(double T, long a) {
  double c_liq[NUMCOMPONENTS-1];
  double c_sol[NUMCOMPONENTS-1];
  
  long i, j, k;
  double sum_s=0.0, sum_l=0.0, sum_c=0.0;
  
  if (a != (NUMPHASES-1)) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      c_liq[k] = gsl_spline_eval (spline_ES[thermo_phase[a]][k][1], T, acc_ES[thermo_phase[a]][k][1]);
      c_sol[k] = gsl_spline_eval (spline_ES[thermo_phase[a]][k][0], T, acc_ES[thermo_phase[a]][k][0]);
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
void function_F_04_init_propertymatrices(double T) {
  //Initialize property matrices
  FILE *fp;
  long a;
  long k,j,i;
  char filename[1000];
  long numlines,lines;
  double composition_solid[NUMCOMPONENTS-1];
  double composition_liquid[NUMCOMPONENTS-1];
  char *file_contents;
  
  struct stat sb;
  
  
//   for (a=0; a < NUM_THERMO_PHASES-1; a++) {
  sprintf(filename,"tdbs_encrypted/Composition_%s.csv",Phases_tdb[0]);
  fp = fopen(filename, "r");
  if (stat(filename, &sb) == -1) {
    perror("stat");
    exit(EXIT_FAILURE);
  }
  
  file_contents = malloc(sb.st_size);
  
  numlines = 0;
  while (fscanf(fp, "%[^\n] ", file_contents) != EOF) {
    numlines++;
  }
  fclose(fp);
//   }
  comp_ES = Malloc4M(NUM_THERMO_PHASES, NUMCOMPONENTS-1, 2,              numlines);
  ThF     = Malloc4M(NUM_THERMO_PHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1, numlines);
  T_ES    = MallocM(NUM_THERMO_PHASES, numlines);
  T_ThF   = MallocM(NUM_THERMO_PHASES, numlines);
  
  for (a=0; a < NUM_THERMO_PHASES-1; a++) {
    sprintf(filename,"tdbs_encrypted/Composition_%s.csv",Phases_tdb[a]);
    fp = fopen(filename, "r");
    fscanf(fp, "%*[^\n]\n");
    for (lines=0; lines < numlines-1; lines++) {
      fscanf(fp, "%le,",&T_ES[a][lines]);
//       printf("%le\n", T_ES[a][lines]);
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        fscanf(fp, "%le,",&comp_ES[a][k][0][lines]);
      }
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        fscanf(fp, "%le,",&comp_ES[a][k][1][lines]);
//         comp_ES[a][k][lines] = composition_liquid[k] - composition_solid[k];
      }
    }
    fclose(fp);
  }
  for (a=0; a < NUM_THERMO_PHASES; a++) {
    printf("%s\n",Phases_tdb[a]);
    sprintf(filename,"tdbs_encrypted/HSN_%s.csv",Phases_tdb[a]);
    fp = fopen(filename, "r");
    fscanf(fp, "%*[^\n]\n");
    for (lines=0; lines < numlines-1; lines++) {
      fscanf(fp, "%le,",&T_ThF[a][lines]);
//       printf("%le\n", T_ThF[a][lines]);
      for (j=0; j < NUMCOMPONENTS-1; j++) {
        fscanf(fp, "%le,",&ThF[a][j][j][lines]);
      }
      for (i=0; i < NUMCOMPONENTS-1; i++) {
        for (k=i+1; k < NUMCOMPONENTS-1; k++) {
           fscanf(fp, "%le,",&ThF[a][i][k][lines]);
           ThF[a][k][i][lines] = ThF[a][i][k][lines];
        }
      }
    }
    fclose(fp);
  }
  
  
  acc_ES     = (gsl_interp_accel****)malloc((NUM_THERMO_PHASES-1)*sizeof(gsl_interp_accel***));
  spline_ES  = (gsl_spline ****)malloc((NUM_THERMO_PHASES-1)*sizeof(gsl_spline***));
  acc_ThF    = (gsl_interp_accel****)malloc(NUM_THERMO_PHASES*sizeof(gsl_interp_accel***));
  spline_ThF = (gsl_spline****)malloc(NUM_THERMO_PHASES*sizeof(gsl_spline***));
  
  for (a=0; a<NUM_THERMO_PHASES; a++) {
    acc_ThF[a]    = (gsl_interp_accel***)malloc((NUMCOMPONENTS-1)*sizeof(gsl_interp_accel**));
    spline_ThF[a] = (gsl_spline***)malloc((NUMCOMPONENTS-1)*sizeof(gsl_spline**));
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      acc_ThF[a][k]    = (gsl_interp_accel**)malloc((NUMCOMPONENTS-1)*sizeof(gsl_interp_accel*));
      spline_ThF[a][k] = (gsl_spline**)malloc((NUMCOMPONENTS-1)*sizeof(gsl_spline*));
      for (j=0; j<NUMCOMPONENTS-1; j++) {
        acc_ThF[a][k][j]    = gsl_interp_accel_alloc ();
        spline_ThF[a][k][j] = gsl_spline_alloc (gsl_interp_cspline, numlines-1);
        gsl_spline_init (spline_ThF[a][k][j], T_ThF[a], ThF[a][k][j], numlines-1);
      }
    }
  }
  for (a=0; a<NUM_THERMO_PHASES-1; a++) {
    acc_ES[a]     = (gsl_interp_accel***)malloc((NUMCOMPONENTS-1)*sizeof(gsl_interp_accel**));
    spline_ES[a]  = (gsl_spline***)malloc((NUMCOMPONENTS-1)*sizeof(gsl_spline**));
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      acc_ES[a][k]     = (gsl_interp_accel**)malloc((NUMCOMPONENTS-1)*sizeof(gsl_interp_accel*));
      spline_ES[a][k]  = (gsl_spline**)malloc((NUMCOMPONENTS-1)*sizeof(gsl_spline*));
      for (i=0; i<2; i++) {
        acc_ES[a][k][i]    = gsl_interp_accel_alloc ();
        spline_ES[a][k][i] = gsl_spline_alloc (gsl_interp_cspline, numlines-1);
        gsl_spline_init (spline_ES[a][k][i],  T_ES[a], comp_ES[a][k][i], numlines-1);
      }
    }
  }
  
  Free4M(comp_ES, NUM_THERMO_PHASES, NUMCOMPONENTS-1, 2);
  Free4M(ThF, NUM_THERMO_PHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
  FreeM(T_ES, NUM_THERMO_PHASES);
  FreeM(T_ThF, NUM_THERMO_PHASES);
  
  function_A(T, c_guess);
  
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
  for(a=0; a<NUMPHASES; a++) {
    for (k=0;k < NUMCOMPONENTS-1;k++) {
        B[a][k] = function_B(T,k,a);
    }
    C[a] = function_C(T,a);
  }
}



#endif
