double dhphi(double *phi, long b, long a) {
  long e, f;
  double sum=0.0;
  
  if (a==b) {
    sum = 6.0*phi[a]*(1.0-phi[a]);
    for (e=0; e < NUMPHASES; e++) {
      for (f=0; f < NUMPHASES; f++) {
        if (e!=b && f!=b && e<f) {
          sum += 2.0*phi[e]*phi[f];
        }
      }
    }
  } else {
    for (e=0; e < NUMPHASES; e++) {
      if (e!=b && e!=a) {
        sum += 2.0*phi[e];
      }
    }
    sum *= phi[b];
  }
  return sum;
}
double hphi(double *phi, long a) {
  long b,c;
  double sum1 = 0.0;
  double sum  = 3.0*phi[a]*phi[a] - 2.0*phi[a]*phi[a]*phi[a];
  for (b=0; b < NUMPHASES; b++) {
    for (c=0; c < NUMPHASES; c++) {
      if (b!=a && c!=a && b < c) {
        sum1 += phi[b]*phi[c];
      }
    }
  }
  sum1 *= 2.0*phi[a];
  return sum + sum1;
}