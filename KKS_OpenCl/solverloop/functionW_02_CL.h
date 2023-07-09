
double F_W_02_dwdphi(double *phi, double *divphi, double *WH,  double *Gamma_abc, int a);

double F_W_02_dwdphi(double *phi, double *divphi, double *WH,  double *Gamma_abc, int a) {
  int b,c;
  double sum=0.0;
  double phibphic;
  
  //The well potential
  for (b=0; b < npha; b++) {
    if (b!=a) {
        sum += 2.0*WH[a*npha+b]*(phi[b]*phi[b]*phi[a]);
    }
  }
  
  //The third order potential
  for (b=0; b < npha; b++) {
    for (c=0; c < npha; c++) {
      if (b!=a && c!=a && b < c) {
          phibphic = phi[b]*phi[c];
            sum += Gamma_abc[(a*npha+b)*npha+c]*phibphic;
      }
    }
  }
  return sum;
}
