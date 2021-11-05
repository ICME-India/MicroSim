#include<complex.h>
#include <gsl/gsl_rng.h>
#include "fftw3.h"

void serialinfo_xy() {
  
  long a;
  long k;
  //long index_count;
  long index;
  start[X]   = 0;
  start[Y]   = 0;
  start[Z]   = 0;
  
  rows_x     = MESH_X + 0;
  rows_y     = MESH_Y + 0;
  rows_z     = MESH_Z + 0;
  end[X]     = rows_x - 1;
  end[Y]     = rows_y - 1;
  end[Z]     = rows_z - 1;
  
  layer_size = rows_y*rows_z;
  
  if (DIMENSION == 2) {
    rows_z     = 1;
    start[Z]   = 0; 
    end[Z]     = 0;
    layer_size = rows_y;
  }
  
  index_count = layer_size*rows_x;  
  
  gridinfo = (struct fields* )malloc((index_count)*sizeof(*gridinfo));
  
  for (index=0; index < index_count; index++) {
    allocate_memory_fields(&gridinfo[index]);
  }

  global_max_min.phi_max        = (double*)malloc(NUMPHASES*sizeof(double));
  global_max_min.phi_min        = (double*)malloc(NUMPHASES*sizeof(double));
  global_max_min.com_max         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  global_max_min.com_min         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  global_max_min.rel_change_phi = (double*)malloc((NUMPHASES)*sizeof(double));
  global_max_min.rel_change_com  = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  
  for (a=0; a<NUMPHASES; a++) {
    global_max_min.phi_max[a] = 1.0;
    global_max_min.phi_min[a] = 0.0;
  }
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    global_max_min.com_max[k] = 1.0;
    global_max_min.com_min[k] = 1.0;
  }

  
}