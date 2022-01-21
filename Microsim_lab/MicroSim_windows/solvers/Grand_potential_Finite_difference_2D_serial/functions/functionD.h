double D(struct fields *gridinfo, double T, long index, long i, long j) {
  if (!ISOTHERMAL) {
//     #ifndef CONSTANT_SUSCEPTIBILITY
    init_propertymatrices(T);
//     #endif
  }
  long b, k;
  double sum=0.0;
  
  for (b=0; b < NUMPHASES; b++) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sum += Diffusivity[b][i][k]*dc_dmu(gridinfo[index].compi, T, b, k, j)*gridinfo[index].phia[b];
    }
  }
  return sum;
}