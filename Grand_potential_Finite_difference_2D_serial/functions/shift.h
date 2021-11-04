void apply_shiftY(struct variables* gridinfo1, long INTERFACE_POS_GLOBAL, int *offset, int *rows) {
  //Shift by one cell in the negative y-direction
  long x, y;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  for(x=0;x < rows_x; x++) {
    for (y=0; y <=(rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
      gidy = x*rows_y + y;
      gridinfo1[gidy] = gridinfo1[gidy+(INTERFACE_POS_GLOBAL-shiftj)];
    }
    if (workers_mpi.lasty==1) {
      for (y=(rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(rows_y-1); y++) {
        gidy = x*rows_y + y;
        for (b=0; b < NUMPHASES-1; b++) {
          gridinfo1[gidy].phia[b] = 0.0;
        }
        gridinfo1[gidy].phia[NUMPHASES-1] = 1.0;
        for (k=0; k < NUMCOMPONENTS-1; k++) {
//          c[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
          c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          chemical_potential       = Mu(c, Tfill, NUMPHASES-1, k);
          gridinfo1[gidy].compi[k] = chemical_potential;
        }
      }
    }
  }
  init_propertymatrices(T);
}
long check_SHIFT(long x) {
  long gidy,center;
  long INTERFACE_POS=0;
  for (gidy=1; gidy <=(rows_y-1); gidy++) {
//     center =  gidy   + (x)*numy[levels];
    center = x*rows_y + gidy; 
//     printf("center=%ld\n",center);
    if ((gridinfo1[center-1].phia[NUMPHASES-1]-(1.0-gridinfo1[center-1].phia[NUMPHASES-1]) < 0.0) 
      && (gridinfo1[center].phia[NUMPHASES-1]-(1.0-gridinfo1[center].phia[NUMPHASES-1]) > 0.0) ) {
      if (gidy > INTERFACE_POS) {
        INTERFACE_POS = gidy;
      }
    }
  }
  if (INTERFACE_POS > 0) {
    return INTERFACE_POS - offset_y + offset[Y];
  } else {
    return 0;
  }
}