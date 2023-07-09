double function_F_03_free_energy_CL(double *c, double T, int a, double *F3A, double *F3B, double *F3C);
void function_F_03_Mu_CL(double *c, double T, int a, double *Mu, double *F3A, double *F3Beq, double *F3BbdT, double Teq);
void function_F_03_c_mu_CL(double *mu, double *c, double T, int a, double *c_guess, double *F3cmu, double *F3Beq, double *F3dBbdT, double Teq);
void function_F_03_dc_dmu_CL(double *mu, double *phase_comp, double T, int a, double *cmu, double *dcdmu);


double function_F_03_free_energy_CL(double *c, double T, int a, double *F3A, double *F3B, double *F3C) {
  int i,j;
  double sum=0.0;
  for (i=0;i<nsol;i++) {
    for (j=0;j<nsol;j++) {
      if (i<=j) {
        sum += F3A[(a*nsol+i)*nsol+j]*c[i]*c[j];
      }
    }
    sum += F3B[a*nsol+i]*c[i];
  }
  sum += F3C[a];
  return sum;
}

void function_F_03_Mu_CL(double *c, double T, int a, double *Mu, double *F3A, double *F3Beq, double *F3dBbdT, double Teq) {
  int j,k;
  double sum=0.0;
  for(k=0; k < nsol; k++) {
    sum  = 2.0*F3A[(a*nsol+k)*nsol+k]*c[k] + (F3Beq[a*nsol+k] + F3dBbdT[a*nsol+k]*(T-Teq));
    for (j=0; j<nsol;j++) {
      if (k!=j) {
        sum += F3A[(a*nsol+k)*nsol+j]*c[j];
      }
    }
    Mu[k] = sum;
    //printf("%le, %le, %le, %le, %le, %le, %le\n", Mu[k], A[a][k][j], c[j], Beq[a][k], dBbdT[a][k], T, Teq);
  }
}

void function_F_03_c_mu_CL(double *mu, double *c, double T, int a, double *c_guess, double *F3cmu, double *F3Beq, double *F3dBbdT, double Teq) {
  int k,l;
  double sum=0.0;
  for(l=0; l < nsol; l++) {
    sum = 0.0;
    for (k=0;k < nsol;k++) {
      sum += F3cmu[(a*nsol+l)*nsol+k]*(mu[k]-(F3Beq[a*nsol+k] + F3dBbdT[a*nsol+k]*(T-Teq)));
    }
    c[l] = sum;
  }
}

void function_F_03_dc_dmu_CL(double *mu, double *phase_comp, double T, int a, double *cmu, double *dcdmu) {
  int i, j;
  for (i=0; i<nsol; i++) {
    for (j=0; j<nsol; j++) {
      dcdmu[i*nsol+j] = cmu[(a*nsol+i)*nsol+j];
    }
  }
}

// double function_F_03_dpsi(double *mu, double **phase_comp, double T, double *phi, int a) {
//   double psi=0.0;
//   double c[nsol];
//   int k=0, b;
//   double sum=0.0;
  
//   for (b=0; b < NUMPHASES; b++) {
//     psi = 0.0;
//     for (k=0;k <nsol; k++) {
//       c[k] = phase_comp[b][k];
//       psi -= mu[k]*c[k];
//     }
//     psi += free_energy(c,T,b);
//     sum += psi*dhphi(phi, b, a);
//   }
//   return sum;
// }

