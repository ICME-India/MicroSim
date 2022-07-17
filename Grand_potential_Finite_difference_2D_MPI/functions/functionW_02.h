#ifndef FUNCTION_W_02_H_
#define FUNCTION_W_02_H_

double function_W_02_dwdphi(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a) {
  long b,c;
  double sum=0.0;
  double phibphic;
  //The well potential
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      if (fabs(divphi[b]) > 0.0) {
        sum += 2.0*Gamma[a][b]*(phi[b]*phi[b]*phi[a]);
      }
    }
  }
  sum *= 9.0;
  //The third order potential
  for (b=0; b < NUMPHASES; b++) {
    for (c=0; c < NUMPHASES; c++) {
      if (b!=a && c!=a && b < c) {
        if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
          phibphic = phi[b]*phi[c];
          sum += Gamma_abc[a][b][c]*phibphic;
//           if (phi[a]*phibphic >= 0.00000) {
//             sum += Gamma_abc[a][b][c]*phibphic;
//           } else {
//             sum += -Gamma_abc[a][b][c]*phibphic;
//           }
        }
      }
    }
  }
  return sum;
}

double function_W_02_dwdphi_smooth(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a) {
  long b,c;
  double sum=0.0;
  double phibphic;
  //The well potential
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      if (fabs(divphi[b]) > 0.0) {
        sum += 2.0*Gamma[a][b]*(phi[b]*phi[b]*phi[a]);
      }
    }
  }
  sum *= 9.0;
  //The third order potential
  for (b=0; b < NUMPHASES; b++) {
    for (c=0; c < NUMPHASES; c++) {
      if (b!=a && c!=a && b < c) {
        if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
          sum += Gamma_abc[a][b][c]*phi[b]*phi[c];
        }
      }
    }
  }
  return sum;
}

#endif