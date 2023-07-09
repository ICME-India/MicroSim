#ifndef FUNCTION_F_02_H_
#define FUNCTION_F_02_H_

// #include "global_vars.h"
double L2norm(double *deltac);
// void tdb_dc_dmu(double *mu, double *phase_comp, double T, long a, double **dcdmu);
// void tdb_Mu(double *c, double T, long a, double *Mu);
// double tdb_dpsi(double *mu, double **phase_comp, double T, double *phi, long a);

struct rparams {
  double *mu;
  double T;
  long phase;
};

double function_F_02_free_energy(double *c, double T, long a) {
  double free_energy;
  double y[NUMCOMPONENTS];
  double sum=0.0;
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    y[k] = c[k];
    sum += y[k];
  }
  y[NUMCOMPONENTS-1] = 1.0 - sum;
  (*free_energy_tdb[thermo_phase[a]])(T, y, &free_energy);
  return free_energy;
}

double f_Mu(const gsl_vector *c, void *p, gsl_vector *f) {
  double *mu = ((struct rparams *) p)->mu;
  double T   = ((struct rparams *) p)->T;
  long a     = ((struct rparams *) p)->phase;
  
  double c0[NUMCOMPONENTS];
  double y[NUMCOMPONENTS-1];
  double mu0[NUMCOMPONENTS-1];
  double sum=0.0;
   
  long i;

  for (i=0; i<NUMCOMPONENTS-1; i++) {
    c0[i] = gsl_vector_get(c, i);
    sum += c0[i];
  }
  c0[NUMCOMPONENTS-1] = 1.0-sum;
  (*Mu_tdb[thermo_phase[a]])(T, c0, mu0);
  
  for (i=0; i<NUMCOMPONENTS-1; i++) {
    y[i]  = (mu0[i] - mu[i])/(R*T);
    gsl_vector_set (f, i, y[i]);
  }
  
  return GSL_SUCCESS;
}
int df_Mu (const gsl_vector *c, void *p, gsl_matrix *J) {

  double *mu = ((struct rparams *) p)->mu;
  double T   = ((struct rparams *) p)->T;
  long a     = ((struct rparams *) p)->phase;
  
  double c0[NUMCOMPONENTS];
  double y[NUMCOMPONENTS-1];
  double mu0[NUMCOMPONENTS-1];
  double Dmudc[(NUMCOMPONENTS-1)*(NUMCOMPONENTS-1)];
  long index;
  double sum=0.0;
  
  for (i=0; i<NUMCOMPONENTS-1; i++) {
    c0[i] = gsl_vector_get(c, i);
    sum  += c0[i];
  }
  c0[NUMCOMPONENTS-1] = 1.0-sum;
  (*dmudc_tdb[thermo_phase[a]])(T, c0, Dmudc);
  
  for (i=0; i<NUMCOMPONENTS-1; i++) {
    for (j=0; j<NUMCOMPONENTS-1; j++) {
      index = i*(NUMCOMPONENTS-1) + j;
      gsl_matrix_set (J, i, j, Dmudc[index]/(R*T));
    }
  }

  return GSL_SUCCESS;
}
int fdf_Mu(const gsl_vector *c, void *params, gsl_vector *f, gsl_matrix *J) {
  f_Mu(c, params, f);
  df_Mu(c, params, J);
  return GSL_SUCCESS;
}
void function_F_02_c_mu(double *mu, double *c, double T, long a, double *c_guess) {
  const gsl_multiroot_fdfsolver_type *solvertype;
  gsl_multiroot_fdfsolver *solver;

  int status;
  size_t i, iter = 0;
  long k;
  
  struct rparams params;
  
  params.mu    = mu;
  params.T     = T;
  params.phase = a;
  
  gsl_multiroot_function_fdf f;

  f.f      = (void*)&f_Mu;
  f.df     = (void*)&df_Mu;
  f.fdf    = (void*)&fdf_Mu;
  f.n      = NUMCOMPONENTS-1;
  f.params = (void*)&params;
  
  gsl_vector *x = gsl_vector_alloc(NUMCOMPONENTS-1);

  for (k=0; k<NUMCOMPONENTS-1; k++) {
    gsl_vector_set (x, k, c_guess[k]);
  }
  
  //solvertype = gsl_multiroot_fdfsolver_gnewton;
  solvertype = gsl_multiroot_fdfsolver_hybridj;
  solver     = gsl_multiroot_fdfsolver_alloc (solvertype, NUMCOMPONENTS-1);
  
  gsl_multiroot_fdfsolver_set (solver, &f, x);
   
  do {
      
    iter++;

    status = gsl_multiroot_fdfsolver_iterate (solver);

//     print_state (iter, s);

    if (status)
      break;

    status = gsl_multiroot_test_residual (solver->f, 1e-8);
  } while (status == GSL_CONTINUE && iter < 1000);

//   printf ("status = %s, solution=%le\n", gsl_strerror (status), gsl_vector_get(s->x, 0));

  for (k=0; k<NUMCOMPONENTS-1; k++) {
    c[k] = gsl_vector_get(solver->x,k);
  }
  
  gsl_multiroot_fdfsolver_free (solver);
  gsl_vector_free (x);
}
void function_F_02_getMu(double *mu, double *c, double *phi, double T, struct gradlayer *grad) {
  long a, l;
  double comp[NUMCOMPONENTS-1];
  double y[NUMCOMPONENTS];
  double sum_c;
  int iter=0;
  do {
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      comp[k] = 0.0;
      for(l=0; l<NUMCOMPONENTS-1; l++) {
        dcdmu[k][l] = 0.0;
      }
    }
    for (a=0; a<NUMPHASES; a++) {
      c_mu(mu, grad->phase_comp[a], T, a, ceq[a][a]);
      dc_dmu(mu, grad->phase_comp[a], T, a, grad->dcdmu_phase[a]);
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        comp[k] += grad->phase_comp[a][k]*hphi(phi, a);
      }
    }
    for(k=0; k<NUMCOMPONENTS-1; k++) {
      deltac[k] = c[k] - comp[k];
      for(l=0; l<NUMCOMPONENTS-1; l++) {
        for (a=0; a<NUMPHASES; a++) {
          dcdmu[k][l] += grad->dcdmu_phase[a][k][l]*hphi(phi,a); 
        }
      }
    }
    matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
//     inv_dcdmu[0][0] = 1.0/dcdmu[0][0];
//     deltamu[0]      = deltac[0]*inv_dcdmu[0][0];
//     mu[0]          += deltamu[0];
    multiply(inv_dcdmu,deltac,deltamu,NUMCOMPONENTS-1);
    vectorsum(deltamu,mu,mu,NUMCOMPONENTS-1);
    iter++;
  } while(L2norm(deltac) > 1e-12);
}
double L2norm(double *deltac){
  long k;
  double norm=0.0;
  for(k=0; k<NUMCOMPONENTS-1; k++) {
    norm += deltac[k]*deltac[k];
  }
  return norm;
}
void function_F_02_dc_dmu(double *mu, double *phase_comp, double T, long a, double **dcdmu) {
  double dmudc[(NUMCOMPONENTS-1)*(NUMCOMPONENTS-1)];
//   double c[NUMCOMPONENTS-1];
  long i,j,index;
  double **dmudc_matrix;
  double sum=0.0;
  
  double y[NUMCOMPONENTS];
    
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    y[k] = phase_comp[k];
    sum += y[k];
  }
  y[NUMCOMPONENTS-1] = 1.0 - sum;
  dmudc_matrix = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
  
  (*dmudc_tdb[thermo_phase[a]])(T, y, dmudc);
   
  for(i=0; i<NUMCOMPONENTS-1; i++) {
    for(j=0; j<NUMCOMPONENTS-1; j++) {
      index = i*(NUMCOMPONENTS-1) + j;
      dmudc_matrix[i][j] = dmudc[index];
    }
  }
  
//   dcdmu[0][0] = 1.0/dmudc_matrix[0][0];
  matinvnew(dmudc_matrix,dcdmu,NUMCOMPONENTS-1);
  FreeM(dmudc_matrix, NUMCOMPONENTS-1);
}
void function_F_02_Mu(double *c, double T, long a, double *Mu) {
  double y[NUMCOMPONENTS];
  double sum=0.0;
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    y[k] = c[k];
    sum += y[k];
  }
  y[NUMCOMPONENTS-1] = 1.0 - sum;
  (*Mu_tdb[thermo_phase[a]])(T, y, Mu);
}
double function_F_02_dpsi(double *mu, double **phase_comp, double T, double *phi, long a) {
  double psi=0.0;
  long k=0, b;
  double sum=0.0;
  
  double y[NUMCOMPONENTS];
  double sum_c;
  
  for (b=0; b < NUMPHASES; b++) {
    psi   = 0.0;
    sum_c = 0.0;
//     function_F_02_getPhaseComposition(mu, c, T, b);
    for (k=0;k <NUMCOMPONENTS-1; k++) {
//       c[k] = c_mu(mu, T, b, k);
      psi -= mu[k]*phase_comp[b][k];
      y[k] = phase_comp[b][k];
      sum_c += y[k];
    }
    y[NUMCOMPONENTS-1] = 1.0 -sum_c;
    psi += free_energy(y,T,b);
    sum += psi*dhphi(phi, b, a);
  }
  return sum;
}

#endif
