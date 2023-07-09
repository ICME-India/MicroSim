#ifndef APPLY_SHIFTY_H
#define APPLY_SHIFTY_H

void apply_shiftY(struct fields* gridinfo, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  for(x=0; x < rows_x; x++) {
    for(z=0; z < rows_z; z++) {
      for (y=0; y <= (rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
        gidy = x*layer_size + z*rows_y + y;
        for (b=0; b < NUMPHASES; b++) {
          gridinfo[gidy].phia[b] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].phia[b];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo[gidy].compi[k] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].compi[k];
        }
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          gridinfo[gidy].composition[k] = gridinfo[gidy+(INTERFACE_POS_GLOBAL-shiftj)].composition[k];
        }
        gridinfo[gidy].temperature = gridinfo[gidy + (INTERFACE_POS_GLOBAL-shiftj)].temperature;
      }
//       if (workers_mpi.lasty==1) {
      for (y=(rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(rows_y-1); y++) {
        gidy = x*layer_size + z*rows_y + y;
        for (b=0; b < NUMPHASES-1; b++) {
          gridinfo[gidy].phia[b] = 0.0;
        }
        gridinfo[gidy].phia[NUMPHASES-1] = 1.0;
        for (k=0; k < NUMCOMPONENTS-1; k++) {
//          c[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
          //c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
          c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
        }
        Mu(c, Teq, NUMPHASES-1, gridinfo[gidy].compi); 
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          //chemical_potential         = Mu(c, Tfill, NUMPHASES-1, k);
          gridinfo[gidy].composition[k]    = c[k]; //chemical_potential;
        }
      }
//       }
    }
  }
//   init_propertymatrices(T);
}
long check_SHIFT(long x) {
  long center;
  long y, z;
  long INTERFACE_POS_MAX = 0;
  for (z=0; z < rows_z; z++) {
    for (y=1; y <=(rows_y-1); y++) {
  //     center =  gidy   + (x)*numy[levels];
      center = x*layer_size + z*rows_y + y; 
  //     printf("center=%ld\n",center);
      if ((gridinfo[center-1].phia[NUMPHASES-1]-(1.0-gridinfo[center-1].phia[NUMPHASES-1]) < 0.0) 
        && (gridinfo[center].phia[NUMPHASES-1]-(1.0-gridinfo[center].phia[NUMPHASES-1]) > 0.0) ) {
        if (y > INTERFACE_POS_MAX) {
          INTERFACE_POS_MAX = y;
        }
      }
    }
  }
  if (INTERFACE_POS_MAX > 0) {
    return INTERFACE_POS_MAX;
  } else {
    return 0;
  }
}

void apply_shiftY_cscl(struct csle *cscl, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  //for(x=0; x < rows_x; x++) {
  //  for(z=0; z < rows_z; z++) {
  //    for (y=0; y <= (rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
  //      gidy = x*layer_size + z*rows_y + y;
  //      for (b=0; b < NUMPHASES; b++) {
  //        for (k=0; k < NUMCOMPONENTS-1; k++) {
  //          cscl[gidy].comie[b][k] = cscl[gidy+(INTERFACE_POS_GLOBAL-shiftj)].comie[b][k];
  //        }
  //      }
  //    }
  //    
  //    for (y=(rows_y-(INTERFACE_POS_GLOBAL-shiftj)); y<=(rows_y-1); y++) {
  //      gidy = x*layer_size + z*rows_y + y;
  //      for (b=0; b < NUMPHASES; b++) {
  //        for (k=0; k < NUMCOMPONENTS-1; k++) {
  //          cscl[gidy].comie[b][k] = ceq[b][b][k];
  //        }
  //      }
  //    }
  //  }
  //}
}

#endif
