#ifndef APPLY_SHIFTY_H
#define APPLY_SHIFTY_H

void apply_shiftY(struct fields* gridinfo_w, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  for(x=0; x < workers_mpi.rows_x; x++) {
    for(z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y <= (workers_mpi.rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
        gidy = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
        for (b=0; b < NUMPHASES; b++) {
          gridinfo_w[gidy].phia[b] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].phia[b];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo_w[gidy].compi[k] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].compi[k];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo_w[gidy].composition[k] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].composition[k];
        }
        gridinfo_w[gidy].temperature = gridinfo_w[gidy + (INTERFACE_POS_GLOBAL-shiftj)].temperature;
      }
      if (workers_mpi.lasty==1) {
        for (y=(workers_mpi.rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(workers_mpi.rows_y-1); y++) {
          gidy = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
          for (b=0; b < NUMPHASES-1; b++) {
            gridinfo_w[gidy].phia[b] = 0.0;
          }
          gridinfo_w[gidy].phia[NUMPHASES-1] = 1.0;
          for (k=0; k < NUMCOMPONENTS-1; k++) {
  //          c[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
            c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
          }
          Mu(c, Tfill, NUMPHASES-1, gridinfo_w[gidy].compi); 
//           for (k=0; k < NUMCOMPONENTS-1; k++) {
//             chemical_potential         = Mu(c, Tfill, NUMPHASES-1, k);
//             gridinfo_w[gidy].compi[k]  = chemical_potential;
//           }
        }
      }
    }
  }
//   init_propertymatrices(T);
}
long check_SHIFT(long x) {
  long center;
  long y, z;
  long INTERFACE_POS_MAX = 0;
  for (z=0; z < workers_mpi.rows_z; z++) {
    for (y=1; y <=(workers_mpi.rows_y-1); y++) {
  //     center =  gidy   + (x)*numy[levels];
      center = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y; 
  //     printf("center=%ld\n",center);
      if ((gridinfo_w[center-1].phia[NUMPHASES-1]-(1.0-gridinfo_w[center-1].phia[NUMPHASES-1]) < 0.0) 
        && (gridinfo_w[center].phia[NUMPHASES-1]-(1.0-gridinfo_w[center].phia[NUMPHASES-1]) > 0.0) ) {
        if (y > INTERFACE_POS_MAX) {
          INTERFACE_POS_MAX = y;
        }
      }
    }
  }
  if (INTERFACE_POS_MAX > 0) {
    return INTERFACE_POS_MAX - workers_mpi.offset_y + workers_mpi.offset[Y];
  } else {
    return 0;
  }
}
#endif