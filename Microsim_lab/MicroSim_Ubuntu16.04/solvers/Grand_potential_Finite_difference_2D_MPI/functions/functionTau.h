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