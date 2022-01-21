#ifndef BOUNDARY_H_
#define BOUNDARY_H_


void copyXZ(long copy_to, long copy_from, long x_start, long x_end, struct variables* gridinfo1) {
  long x, gidy_from, gidy_to;
  for (x=x_start; x<=x_end; x++) {
    gidy_from = x*rows_y + copy_from;
    gidy_to   = x*rows_y + copy_to;
    gridinfo1[gidy_to]= gridinfo1[gidy_from];
  }
}
void copyYZ(long copy_to, long copy_from, struct variables* gridinfo1) {
  long gidy_from, gidy_to, y;
  for (y=0; y < rows_y; y++) {
    gidy_from          = copy_from*rows_y + y;
    gidy_to            = copy_to*rows_y   + y;
    gridinfo1[gidy_to] = gridinfo1[gidy_from];
  }
}
void apply_boundary_conditions(long taskid){
#ifdef PERIODIC
  if (!((workers_mpi.firstx ==1) && (workers_mpi.lastx ==1))) {
    mpiboundary_left_right(taskid);
  }
#endif
#ifndef PERIODIC
  if (workers_mpi.firstx ==1) {
    copyYZ(0,5,gridinfo1);
    copyYZ(1,4,gridinfo1);
    copyYZ(2,3,gridinfo1);
  }
  if (workers_mpi.lastx ==1) {
    copyYZ(end[X]+1,end[X],gridinfo1);
    copyYZ(end[X]+2,end[X]-1,gridinfo1);
    copyYZ(end[X]+3,end[X]-2,gridinfo1);
  }
#endif
#ifdef PERIODIC_Y
  if (!((workers_mpi.firsty ==1) && (workers_mpi.lasty ==1))) {
     mpiboundary_top_bottom(taskid);
  }
//   copyXZ(2,MESH_Y-4,start,end,gridinfo1);
//   copyXZ(1,MESH_Y-5,start,end,gridinfo1);
//   copyXZ(0,MESH_Y-6,start,end,gridinfo1);
//   
//   copyXZ(MESH_Y-3,3,start,end,gridinfo1);
//   copyXZ(MESH_Y-2,4,start,end,gridinfo1);
//   copyXZ(MESH_Y-1,5,start,end,gridinfo1);
//   
#endif
#ifdef ISOLATE_Y
 if (workers_mpi.firsty ==1) {
    copyXZ(2, 3, start[X], end[X], gridinfo1);
    copyXZ(1, 4, start[X], end[X], gridinfo1);
    copyXZ(0, 5, start[X], end[X], gridinfo1);
  }
  if(workers_mpi.lasty ==1) {
    copyXZ(rows_y-3, rows_y-4, start[X], end[X], gridinfo1);
    copyXZ(rows_y-2, rows_y-5, start[X], end[X], gridinfo1);
    copyXZ(rows_y-1, rows_y-6, start[X], end[X], gridinfo1);
  }
#endif
}
#endif