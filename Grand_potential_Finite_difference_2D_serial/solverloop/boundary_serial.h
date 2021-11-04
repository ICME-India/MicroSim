#ifndef BOUNDARY_01_H_
#define BOUNDARY_01_H_

void copyXZ(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo1, char *field_type);
void copyXY(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo1, char *field_type);
void copyYZ(struct bc_scalars *boundary, struct fields* gridinfo1, char *field_type);

void copyXZ(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo1, char *field_type) {
  long gidy_from, gidy_to, y, a, k;
  long copy_from, copy_to;
  long j;
  long x, z;
  if (strcmp(field_type, "PHI") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if((boundary[0].type ==1) || (boundary[0].type == 3)) {
        copy_from = boundary[0].proxy[j];
        copy_to   = boundary[0].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < rows_z; z++) {
            gidy_from          = x*layer_size + z*rows_y + copy_from;
            gidy_to            = x*layer_size + z*rows_y + copy_to;
            for (a=0; a < NUMPHASES; a++) {
              gridinfo1[gidy_to].phia[a] = gridinfo1[gidy_from].phia[a];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "MU") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if((boundary[1].type == 1) || (boundary[1].type == 3)) {
        copy_from = boundary[1].proxy[j];
        copy_to   = boundary[1].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < rows_z; z++) {
            gidy_from          = x*layer_size + z*rows_y + copy_from;
            gidy_to            = x*layer_size + z*rows_y + copy_to;
            for (k=0; k< NUMCOMPONENTS-1; k++) {
              gridinfo1[gidy_to].compi[k] = gridinfo1[gidy_from].compi[k];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "T") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if((boundary[2].type == 1) || (boundary[2].type == 3)) {
        copy_from = boundary[2].proxy[j];
        copy_to   = boundary[2].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < rows_z; z++) {
            gidy_from                      = x*layer_size + z*rows_y + copy_from;
            gidy_to                        = x*layer_size + z*rows_y + copy_to;
            gridinfo1[gidy_to].temperature = gridinfo1[gidy_from].temperature;
          }
        }
      }
    }
  }
}
void copyYZ(struct bc_scalars *boundary, struct fields* gridinfo, char *field_type) {
  long gidy_from, gidy_to, y, a, k;
  long copy_from, copy_to;
  long j;
  long x, z;
  if (strcmp(field_type, "PHI") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if((boundary[0].type ==1) || (boundary[0].type == 3)) {
        copy_from = boundary[0].proxy[j];
        copy_to   = boundary[0].points[j];
        for (y=0; y < rows_y; y++) {
          for (z=0; z < rows_z; z++) {
            gidy_from          = copy_from*layer_size + z*rows_y  + y;
            gidy_to            = copy_to*layer_size   + z*rows_y  + y;
            for (a=0; a < NUMPHASES; a++) {
              gridinfo[gidy_to].phia[a] = gridinfo[gidy_from].phia[a];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "MU") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if((boundary[1].type ==1) || (boundary[1].type == 3)) {
        copy_from = boundary[1].proxy[j];
        copy_to   = boundary[1].points[j];
        for (y=0; y < rows_y; y++) {
          for (z=0; z < rows_z; z++) {
            gidy_from          = copy_from*layer_size  +  z*rows_y + y;
            gidy_to            = copy_to*layer_size    +  z*rows_y + y;
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              gridinfo[gidy_to].compi[k] = gridinfo[gidy_from].compi[k];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "T") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if((boundary[2].type == 1) || (boundary[2].type == 3)) {
        copy_from = boundary[2].proxy[j];
        copy_to   = boundary[2].points[j];
        for (y=0; y < rows_y; y++) {
          for (z=0; z < rows_z; z++) {
            gidy_from          = copy_from*layer_size + z*rows_y + y;
            gidy_to            = copy_to*layer_size   + z*rows_y + y;
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              gridinfo[gidy_to].temperature = gridinfo[gidy_from].temperature;
            }
          }
        }
      }
    }
  }
}
void copyXY(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo1, char *field_type) {
  long gidy_from, gidy_to, y, a, k;
  long copy_from, copy_to;
  long j;
  long x, z;
  if (strcmp(field_type, "PHI") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if((boundary[0].type ==1) || (boundary[0].type == 3)) {
          copy_from = boundary[0].proxy[j];
          copy_to   = boundary[0].points[j];
          for (x=x_start; x < x_end; x++) {
            for (y=0; y < rows_y; y++) {
              gidy_from          = x*layer_size + copy_from*rows_y + y;
              gidy_to            = x*layer_size + copy_to*rows_y   + y;
              for (a=0; a < NUMPHASES; a++) {
                gridinfo1[gidy_to].phia[a]     = gridinfo1[gidy_from].phia[a];
                gridinfo1[gidy_to].deltaphi[a] = gridinfo1[gidy_from].deltaphi[a];
              }
            }
          }
        }
      }
   }
   if (strcmp(field_type, "MU") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
        if((boundary[1].type ==1) || (boundary[1].type == 3)) {
          copy_from = boundary[1].proxy[j];
          copy_to   = boundary[1].points[j];
          for (x=x_start; x < x_end; x++) {
            for (y=0; y < rows_y; y++) {
              gidy_from          = x*layer_size + copy_from*rows_y + y;
              gidy_to            = x*layer_size + copy_to*rows_y   + y;
              for(k=0; k<NUMCOMPONENTS-1; k++) {
                gridinfo1[gidy_to].compi[k] = gridinfo1[gidy_from].compi[k];
              }
            }
          }
        }
      }
    }
   if (strcmp(field_type, "T") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if((boundary[2].type ==1) || (boundary[2].type == 3)) {
          copy_from = boundary[2].proxy[j];
          copy_to   = boundary[2].points[j];
          for (x=x_start; x < x_end; x++) {
            for (y=0; y < rows_y; y++) {
              gidy_from                      = x*layer_size + copy_from*rows_y + y;
              gidy_to                        = x*layer_size + copy_to*rows_y   + y;
              gridinfo1[gidy_to].temperature = gridinfo1[gidy_from].temperature;
            }
          }
        }
      }
   }
}
void apply_boundary_conditions(long taskid) {
  int i, j, field_num;
  for (i=0; i<6; i++) {
    if ((i==0) || (i==1)) {
      copyYZ(boundary[i], gridinfo, "PHI");
      copyYZ(boundary[i], gridinfo, "MU");
      if (!ISOTHERMAL) {
        copyYZ(boundary[i], gridinfo, "T");
      }
    }
    if ((i==2) || (i==3)) {
      copyXZ(boundary[i], 0, rows_x-1, gridinfo, "PHI");
      copyXZ(boundary[i], 0, rows_x-1, gridinfo, "MU");
      if (!ISOTHERMAL) {
        copyXZ(boundary[i], 0, rows_x-1, gridinfo, "T");
      }
    }
    if (DIMENSION != 2) {
      if ((i==4)||(i==5)) {
        copyXY(boundary[i], 0, rows_x-1, gridinfo, "PHI");
        copyXY(boundary[i], 0, rows_x-1, gridinfo, "MU");
        if (!ISOTHERMAL) {
          copyXY(boundary[i], 0, rows_x-1, gridinfo, "T");
        }
      }
    }
  }
}




// void apply_boundary_conditions(long taskid){
// #ifdef PERIODIC
//   if (!((workers_mpi.firstx ==1) && (workers_mpi.lastx ==1))) {
//     mpiboundary_left_right(taskid);
//   }
// #endif
// #ifndef PERIODIC
//   if (workers_mpi.firstx ==1) {
//     copyYZ(0,5,gridinfo1);
//     copyYZ(1,4,gridinfo1);
//     copyYZ(2,3,gridinfo1);
//   }
//   if (workers_mpi.lastx ==1) {
//     copyYZ(end[X]+1,end[X],gridinfo1);
//     copyYZ(end[X]+2,end[X]-1,gridinfo1);
//     copyYZ(end[X]+3,end[X]-2,gridinfo1);
//   }
// #endif
// #ifdef PERIODIC_Y
//   if (!((workers_mpi.firsty ==1) && (workers_mpi.lasty ==1))) {
//      mpiboundary_top_bottom(taskid);
//   }
// //   copyXZ(2,MESH_Y-4,start,end,gridinfo1);
// //   copyXZ(1,MESH_Y-5,start,end,gridinfo1);
// //   copyXZ(0,MESH_Y-6,start,end,gridinfo1);
// //   
// //   copyXZ(MESH_Y-3,3,start,end,gridinfo1);
// //   copyXZ(MESH_Y-2,4,start,end,gridinfo1);
// //   copyXZ(MESH_Y-1,5,start,end,gridinfo1);
// //   
// #endif
// #ifdef ISOLATE_Y
//  if (workers_mpi.firsty ==1) {
//     copyXZ(2, 3, start[X], end[X], gridinfo1);
//     copyXZ(1, 4, start[X], end[X], gridinfo1);
//     copyXZ(0, 5, start[X], end[X], gridinfo1);
//   }
//   if(workers_mpi.lasty ==1) {
//     copyXZ(rows_y-3, rows_y-4, start[X], end[X], gridinfo1);
//     copyXZ(rows_y-2, rows_y-5, start[X], end[X], gridinfo1);
//     copyXZ(rows_y-1, rows_y-6, start[X], end[X], gridinfo1);
//   }
// #endif
// }
#endif
