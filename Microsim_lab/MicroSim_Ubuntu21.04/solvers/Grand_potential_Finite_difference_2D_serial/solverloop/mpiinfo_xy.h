void Mpiinfo(long taskid) {
  long averow[DIMENSION];
  int rank;
  long i, j;
//   long averow_y;
  workers_mpi.lastx=0;
  workers_mpi.firstx=0;
  workers_mpi.lasty=0;
  workers_mpi.firsty=0;
  
  if (taskid == MASTER) {
    /* Distribute work to workers.  Must first figure out how many rows to */
    /* send and what to do with extra rows.  */
    
    averow[X] = (MESH_X)/numworkers_x;
    averow[Y] = (MESH_Y)/numworkers_y;
    
    extra[X]  = (MESH_X)%numworkers_x;
    extra[Y]  = (MESH_Y)%numworkers_y;
    
//     printf("numworkers_x=%d, numworkers_y=%d, extra[X]=%d, extra[Y]=%d\n",numworkers_x,  numworkers_y, extra[X], extra[Y]);
//     exit(0);
    
    offset[X] = 0;
    offset[Y] = 0;
    for (rank=1; rank<=(numworkers_x*numworkers_y); rank++) {
      if (extra[X] > 0) {
        rows[X] = ((rank-1)/(numworkers_y) < extra[X]) ? averow[X]+1 : averow[X];
      } else {
        rows[X] = averow[X];
      }
      
      if (extra[Y] > 0) {
       rows[Y]  = ((rank-1)%(numworkers_y) < extra[Y]) ? averow[Y]+1 : averow[Y];
      } else {
        rows[Y] = averow[Y];
      }
      /* Tell each worker who its neighbors are, since they must exchange */
      /* data with each other. */  
      if ((rank-1)/(numworkers_y) == 0) {
#ifdef PERIODIC
	      left_node = rank + (numworkers_x-1)*numworkers_y;
#else
	      left_node = NONE;
#endif
              workers_mpi.firstx  = 1;
      } else {
	      left_node = rank - numworkers_y;
      }
      if (((rank-1)/numworkers_y) == (numworkers_x-1)) {
#ifdef PERIODIC
	      right_node =  rank - (numworkers_x-1)*numworkers_y;
#else
	      right_node = NONE;
#endif
              workers_mpi.lastx  = 1;
      } else {
	      right_node = rank + numworkers_y;
      }
      if ((rank-1)%(numworkers_y) == 0) {
#ifdef PERIODIC_Y
	      bottom_node = rank + (numworkers_y-1);
#else
	      bottom_node = NONE;
#endif
              workers_mpi.firsty  = 1;
      } else {
	      bottom_node = rank - 1;
      }
      if (((rank-1)%numworkers_y) == (numworkers_y-1)) {
#ifdef PERIODIC_Y
	      top_node   =  rank - (numworkers_y-1);
#else
	      top_node   = NONE;
#endif
              workers_mpi.lasty  = 1;
      } else {
	      top_node   = rank + 1;
      }
      workers_mpi.rank_x = (rank-1)/numworkers_y;
      workers_mpi.rank_y = (rank-1)%numworkers_y;
      
      /*  Now send startup information to each worker  */
      dest = rank;
      MPI_Send(offset,       2, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(rows,         2, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&left_node,   1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&right_node,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&top_node,    1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&bottom_node, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      
      MPI_Send(&workers_mpi.firstx, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&workers_mpi.firsty, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&workers_mpi.lastx,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&workers_mpi.lasty,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&workers_mpi.rank_x, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&workers_mpi.rank_y, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      
      
      for (i=0; i < rows[X]; i++) {
        MPI_Send(gridinfo+offset[X]*MESH_Y+ offset[Y]+ i*(MESH_Y), rows[Y], MPI_gridinfo, dest, BEGIN, MPI_COMM_WORLD);
      }
//       printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
//       printf("left= %d right= %d, rows=%d\n",left_node,right_node,rows);
      if (((offset[Y] + rows[Y])%MESH_Y) == 0) {
        offset[X] = (offset[X] + rows[X]);
      }
      offset[Y] = ((offset[Y] + rows[Y])%MESH_Y);
      
      workers_mpi.firstx=0;
      workers_mpi.firsty=0;
      workers_mpi.lasty=0;
      workers_mpi.lastx=0;
    }
  } else {
    source = MASTER;
    msgtype = BEGIN;
    
    MPI_Recv(offset, 2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, offset_x=%d, offset_y=%d\n",taskid, offset[X], offset[Y]);
    MPI_Recv(rows, 2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, rows[X]=%d, rows[Y]=%d\n",taskid,rows[X], rows[Y]);
    MPI_Recv(&left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, left=%d\n",taskid, left_node);
    MPI_Recv(&right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, right=%d\n",taskid, right_node);
    
    MPI_Recv(&top_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, top=%d\n",taskid, top_node);
    MPI_Recv(&bottom_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, bottom=%d\n", taskid,bottom_node);
    
    MPI_Recv(&workers_mpi.firstx,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.firsty,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.lastx,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.lasty,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.rank_x,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.rank_y,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    start[X] = 3;
    if ((offset[X]==0) || (offset[X]+rows[X] == MESH_X) || (offset[Y]==0) || (offset[Y]+rows[Y] == MESH_Y)) {
      if ((offset[X] == 0) || (offset[X]+rows[X] == MESH_X)) {
        rows_x     = rows[X] + 3;
        end[X]     = rows_x  - 4;
        if ((offset[X]+rows[X]) == MESH_X) {
          offset_x = 3; 
        } else {
          offset_x = 0;
        }
        if ((offset[X] == 0) && (offset[X]+rows[X] == MESH_X)) { //Just one worker in the x-direction
          offset_x   = 0;
          rows_x     = rows[X];
          end[X]     = rows_x  - 4;
        }
      } else {
        rows_x   = rows[X] + 6;
        end[X]   = rows_x  - 4;
        offset_x = 3; 
      }
      if ((offset[Y] == 0) || (offset[Y]+rows[Y] == MESH_Y)) {
        rows_y     = rows[Y] + 3;
        end[Y]     = rows_y  - 4;
        if (offset[Y]+rows[Y] == MESH_Y) {
          offset_y = 3; 
        } else {
          offset_y = 0;
        }
        if ((offset[Y] == 0) && (offset[Y]+rows[Y] == MESH_Y)) { //Just one worker in y-direction
         offset_y = 0;
         rows_y   = rows[Y];
         end[Y]   = rows_y - 4;
        }
      } else {
        rows_y     = rows[Y] + 6;
        end[Y]     = rows_y  - 4;
        offset_y   = 3;
      }
    } else {
      rows_x   = rows[X] + 6;
      rows_y   = rows[Y] + 6;
      end[X]   = rows_x  - 4;
      end[Y]   = rows_y  - 4;
      offset_x = 3;
      offset_y = 3;
    }
    gridinfo1 = (struct variables *)malloc((rows_x*rows_y)*(sizeof(*gridinfo1)));
      
    for (i=0; i < rows[X]; i++) {
      MPI_Recv(gridinfo1 + (i+offset_x)*rows_y + offset_y, rows[Y], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
    printf("taskid=%ld, rows_x=%d, rows_y=%d, firstx=%d, firsty=%d, lastx=%d, lasty=%d,\n",
           taskid, rows_x, rows_y,workers_mpi.firstx, workers_mpi.firsty, workers_mpi.lastx, workers_mpi.lasty);
  }
}
void sendtomaster() {
  long i;
  MPI_Send(offset, 2, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
  MPI_Send(rows,   2, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
    
  for (i=0; i < rows[X]; i++) {
    MPI_Send(gridinfo1 + (i+offset_x)*rows_y + offset_y, rows[Y], MPI_gridinfo, MASTER, DONE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  long i;
  for (rank=1; rank <=numworkers; rank++) {
    source = rank;
    msgtype = DONE;
    MPI_Recv(offset, 2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(rows,   2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    for (i=0; i < rows[X]; i++) {
      MPI_Recv(gridinfo + offset[X]*MESH_Y + offset[Y] + i*MESH_Y,  rows[Y], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
  }
}
void mpiexchange_left_right(long taskid) {
 if (workers_mpi.rank_x%2) {
    if (workers_mpi.firstx ==0) {
      MPI_Send(&gridinfo1[3*rows_y], 3*rows_y, MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    
      source = left_node;
      msgtype = LTAG;
      MPI_Recv(&gridinfo1[0*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
    if (workers_mpi.lastx == 0) {
      MPI_Send(&gridinfo1[(end[X]-2)*rows_y], 3*rows_y, MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
      
      source = right_node;
      msgtype = RTAG;
      
      MPI_Recv(&gridinfo1[(end[X]+1)*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
  } else {
    if (workers_mpi.lastx == 0) {
      source = right_node;
      msgtype = RTAG;
      MPI_Recv(&gridinfo1[(end[X]+1)*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      MPI_Send(&gridinfo1[(end[X]-2)*rows_y], 3*rows_y, MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
    }
    if(workers_mpi.firstx ==0) {
      source = left_node;
      msgtype = LTAG;
      MPI_Recv(&gridinfo1[0*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      MPI_Send(&gridinfo1[3*rows_y], 3*rows_y, MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    }
  }
}
void mpiexchange_top_bottom(long taskid) {
 long i;
 if (workers_mpi.rank_y%2) {
    if (workers_mpi.firsty ==0) {
      for (i=0; i < rows_x; i++) {
        MPI_Send(&gridinfo1[i*rows_y + 3],          3, MPI_gridinfo, bottom_node, TTAG, MPI_COMM_WORLD);
        source = bottom_node;
        msgtype = BTAG;
      
        MPI_Recv(&gridinfo1[i*rows_y + 0],          3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      }
    }
    if (workers_mpi.lasty == 0) {
      for (i=0; i < rows_x; i++) {
        MPI_Send(&gridinfo1[i*rows_y + (end[Y]-2)], 3, MPI_gridinfo, top_node, BTAG, MPI_COMM_WORLD);
        
        source = top_node;
        msgtype = TTAG;
        
        MPI_Recv(&gridinfo1[i*rows_y + (end[Y]+1)], 3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      }
    }
  } else {
    if (workers_mpi.lasty == 0) {
      for (i=0; i < rows_x; i++) {
        source = top_node;
        msgtype = TTAG;
        MPI_Recv(&gridinfo1[i*rows_y + (end[Y]+1)], 3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
        
        MPI_Send(&gridinfo1[i*rows_y + (end[Y]-2)], 3, MPI_gridinfo, top_node, BTAG, MPI_COMM_WORLD);
      }
    }
    if(workers_mpi.firsty ==0) {
      for (i=0; i < rows_x; i++) {
        source = bottom_node;
        msgtype = BTAG;
        MPI_Recv(&gridinfo1[i*rows_y + 0],          3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
        
        MPI_Send(&gridinfo1[i*rows_y + 3],          3, MPI_gridinfo, bottom_node, TTAG, MPI_COMM_WORLD);
      }
    }
  }
}
void mpiboundary_left_right(long taskid) {
  if (workers_mpi.firstx ==1) {
    MPI_Send(&gridinfo1[3*rows_y], 3*rows_y, MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    source = left_node;
    msgtype = LTAG;
    MPI_Recv(&gridinfo1[0*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
  }
  if(workers_mpi.lastx ==1) {
    source = right_node;
    msgtype = RTAG;
    MPI_Recv(&gridinfo1[(end[X]+1)*rows_y], 3*rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Send(&gridinfo1[(end[X]-2)*rows_y], 3*rows_y, MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
  }
}
void mpiboundary_top_bottom(long taskid) {
  long i;
  if (workers_mpi.firsty ==1) {
    for (i=0; i < rows_x; i++) {
      MPI_Send(&gridinfo1[i*rows_y + 3],          3, MPI_gridinfo, bottom_node, TTAG, MPI_COMM_WORLD);
      source = bottom_node;
      msgtype = BTAG;
      MPI_Recv(&gridinfo1[i*rows_y + 0],          3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
  }
  if(workers_mpi.lasty ==1) {
    for (i=0; i < rows_x; i++) {
      source = top_node;
      msgtype = TTAG;
      MPI_Recv(&gridinfo1[i*rows_y + (end[Y]+1)], 3, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Send(&gridinfo1[i*rows_y + (end[Y]-2)], 3, MPI_gridinfo, top_node, BTAG, MPI_COMM_WORLD);
    }
  }
}