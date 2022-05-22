double FunctionTau(double *phi) {
  double sum=0.0, sum1=0.0;
  long a, b;
  for (a=0; a<NUMPHASES; a++) {
    for (b=0; b<NUMPHASES; b++) {
      if (a<b) {
       sum  += tau_ab[a][b]*phi[a]*phi[b];
       sum1 += phi[a]*phi[b];
      }
    }
  }
  if (sum1) {
    return sum/sum1;
  } else {
    return tau;
  }
}
void Calculate_Tau() {
  //Initialize property matrices
  long a, b, i, j, k;
  double mu[NUMCOMPONENTS-1];
  Mu(ceq[NUMPHASES-1][NUMPHASES-1], T, NUMPHASES-1, mu);
  dc_dmu(mu, ceq[NUMPHASES-1][NUMPHASES-1], T, NUMPHASES-1, dcdmu);
  double min_tau=0.0;
//   min_tau = tau_ab[0][NUMPHASES-1];
  for (a=0; a<NUMPHASES-1; a++) {
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      deltac[k] = ceq[a][NUMPHASES-1][k] - ceq[a][a][k];
    }
    multiply2d(Diffusivity[NUMPHASES-1],dcdmu,Ddcdmu,NUMCOMPONENTS-1);
    matinvnew(Ddcdmu,inv_dcdmu,NUMCOMPONENTS-1);
    multiply(inv_dcdmu,deltac,deltamu,NUMCOMPONENTS-1);
    double sum=0.0;
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      sum += deltamu[k]*deltac[k];
    }
    tau_ab[a][NUMPHASES-1] = sum*epsilon*(0.2222)/V;
    tau_ab[NUMPHASES-1][a] = tau_ab[a][NUMPHASES-1];
    
    if (a==0) {
      min_tau = tau_ab[a][NUMPHASES-1];
    }
    if (tau_ab[a][NUMPHASES-1] < min_tau) {
      min_tau = tau_ab[a][NUMPHASES-1];
    }
  }
  for (a=0; a<NUMPHASES-1; a++) {
    for (b=0; b<NUMPHASES-1; b++) {
      tau_ab[a][b] = min_tau;
    }
  }
  tau = min_tau;
  printf("tau_ab[0][1]=%le, tau=%le, dcdmu[0][0]=%le\n", tau_ab[0][1], tau, dcdmu[0][0]);
}
