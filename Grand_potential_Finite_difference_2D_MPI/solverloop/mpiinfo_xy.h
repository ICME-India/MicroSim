#ifndef MPIINFO_XY_H_
#define MPIINFO_XY_H_


void fill_buffer_x(long start_j, long x_start, long x_end);
void fill_buffer_y(long start_j, long y_start, long y_end);
void fill_gridinfo_x(long start_j, long x_start, long x_end);
void fill_gridinfo_y(long start_j, long y_start, long y_end);

void Mpiinfo(long taskid) {
//   long averow[DIMENSION];
  long rank;
  long i, j;
  long x, y;
//   long averow_y;
  long lastx;
  long firstx;
  long lasty;
  long firsty;
  long rank_x;
  long rank_y;
  long index;
  long index_count;
  long index_count_1;
  long index_w;
  
  
  lastx=0;
  firstx=0;
  lasty=0;
  firsty=0;
  
  if (taskid == MASTER) {
    /* Distribute work to workers.  Must first figure out how many rows to */
    /* send and what to do with extra rows.  */
    
    averow[X] = (rows_x)/numworkers_x;
    averow[Y] = (rows_y)/numworkers_y;
    
    extra[X]  = (rows_x)%numworkers_x;
    extra[Y]  = (rows_y)%numworkers_y;
    
    offset[X] = 0;
    offset[Y] = 0;
    
    for (rank=0; rank < (numworkers_x*numworkers_y); rank++) {
      if (extra[X] > 0) {
        rows[X] = ((rank)/(numworkers_y) < extra[X]) ? averow[X]+1 : averow[X];
      } else {
        rows[X] = averow[X];
      }
      
      if (extra[Y] > 0) {
        rows[Y]  = ((rank)%(numworkers_y) < extra[Y]) ? averow[Y]+1 : averow[Y];
      } else {
        rows[Y] = averow[Y];
      }
      /* Tell each worker who its neighbors are, since they must exchange */
      /* data with each other. */  
      if ((rank)/(numworkers_y) == 0) {
        //By default the neighboring will be set to periodic
// #ifdef PERIODIC
	      left_node = rank + (numworkers_x-1)*numworkers_y;
// #else
// 	      left_node = NONE;
// #endif
              firstx  = 1;
      } else {
	      left_node = rank - numworkers_y;
      }
      if (((rank)/numworkers_y) == (numworkers_x-1)) {
        //By default the neighboring will be set to periodic
// #ifdef PERIODIC
	      right_node =  rank - (numworkers_x-1)*numworkers_y;
// #else
// 	      right_node = NONE;
// #endif
              lastx  = 1;
      } else {
	      right_node = rank + numworkers_y;
      }
      if ((rank)%(numworkers_y) == 0) {
// #ifdef PERIODIC_Y
	      bottom_node = rank + (numworkers_y-1);
// #else
// 	      bottom_node = NONE;
// #endif
             firsty  = 1;
      } else {
	      bottom_node = rank - 1;
      }
      if (((rank)%numworkers_y) == (numworkers_y-1)) {
// #ifdef PERIODIC_Y
	      top_node   =  rank - (numworkers_y-1);
// #else
// 	      top_node   = NONE;
// #endif
             lasty  = 1;
      } else {
	      top_node   = rank + 1;
      }
      rank_x = (rank)/numworkers_y;
      rank_y = (rank)%numworkers_y;
      
      
      /*  Now send startup information to each worker  */
      if (rank != MASTER) {
        if (rank==1) {
          printf("sending, taskid=%ld, rows[X]=%ld, rows[Y]=%ld\n",rank, rows[X], rows[Y]);
        }
        dest = rank;
        MPI_Send(offset,       2, MPI_LONG, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(rows,         2, MPI_LONG, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&left_node,   1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&right_node,  1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&top_node,    1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&bottom_node, 1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        
        MPI_Send(&firstx, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&firsty, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&lastx,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&lasty,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rank_x, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rank_y, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        

        index_count = (rows[Y])*(rows[X]);
        
        buffer = (double *)malloc((index_count)*(SIZE_STRUCT_FIELDS)*sizeof(double));
        
        j=0;
        for (x=0; x < rows[X]; x++) {
          for (y=0; y < rows[Y]; y++) {
            index = (x + offset[X])*layer_size + (y+offset[Y]);
            for (a=0; a<NUMPHASES; a++) {
              buffer[j] = gridinfo[index].phia[a];
              j++;
            }
            for (k=0; k<NUMCOMPONENTS-1; k++) {
              buffer[j] = gridinfo[index].compi[k];
              j++;
            }
            for (k=0; k<NUMCOMPONENTS-1; k++) {
              buffer[j] = gridinfo[index].composition[k];
              j++;
            }
            for (a=0; a<NUMPHASES; a++) {
              buffer[j] = gridinfo[index].deltaphi[a];
              j++;
            }
            buffer[j] = gridinfo[index].temperature;
            j++;
          }
        }
        
        MPI_Datatype MPI_gridinfo_vector;
    
        MPI_Type_vector(index_count, SIZE_STRUCT_FIELDS, SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector);
        MPI_Type_commit(&MPI_gridinfo_vector);

        MPI_Send(buffer, 1, MPI_gridinfo_vector, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Type_free(&MPI_gridinfo_vector);
        free(buffer);
      } else {
        start[X] = 3;
        end[X]   = 3;
        workers_mpi.firstx      = firstx;
        workers_mpi.lastx       = lastx;
        workers_mpi.firsty      = firsty;
        workers_mpi.lasty       = lasty;
        workers_mpi.rank_x      = rank_x;
        workers_mpi.rank_y      = rank_y;
        workers_mpi.offset[X]   = offset[X];
        workers_mpi.offset[Y]   = offset[Y];
        workers_mpi.rows[X]     = rows[X];
        workers_mpi.rows[Y]     = rows[Y];
        workers_mpi.left_node   = left_node;
        workers_mpi.right_node  = right_node;
        workers_mpi.top_node    = top_node;
        workers_mpi.bottom_node = bottom_node;
        workers_mpi.rows_x      = rows[X] + 3;
        workers_mpi.rows_y      = rows[Y] + 3; 
        workers_mpi.rows_z      = rows[Z] + 3;
        workers_mpi.start[X]    = start[X];
        workers_mpi.end[X]      = workers_mpi.rows_x - 4;
        workers_mpi.start[Y]    = start[Y];
        workers_mpi.end[Y]      = workers_mpi.rows_y - 4;
        workers_mpi.layer_size  = workers_mpi.rows_y*workers_mpi.rows_z; 
        workers_mpi.offset_x    = 0;
        workers_mpi.offset_y    = 0;
        workers_mpi.offset_z    = 0;
        
        if ((workers_mpi.firstx==1) || (workers_mpi.lastx==1) || (workers_mpi.firsty==1) || (workers_mpi.lasty==1)) {
          if((workers_mpi.firstx==1) || (workers_mpi.lastx==1)) {
            workers_mpi.rows_x = workers_mpi.rows[X] + 3;
            workers_mpi.end[X] = workers_mpi.rows_x  - 4;
            if (workers_mpi.lastx==1) {
              workers_mpi.offset_x = 3; 
            } else {
              workers_mpi.offset_x = 0;
            }
            if (workers_mpi.firstx && workers_mpi.lastx) { //Just one worker in the x-direction
              workers_mpi.offset_x   = 0;
              workers_mpi.rows_x     = workers_mpi.rows[X];
              workers_mpi.end[X]     = workers_mpi.rows_x - 4;
            }
          } else {
            workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
            workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
            workers_mpi.offset_x = 3;
          }
          if ((workers_mpi.firsty==1) || (workers_mpi.lasty==1)) {
            workers_mpi.rows_y     = workers_mpi.rows[Y] + 3;
            workers_mpi.end[Y]     = workers_mpi.rows_y  - 4;
            
            if (workers_mpi.lasty==1) {
              workers_mpi.offset_y = 3; 
            } else {
              workers_mpi.offset_y = 0;
            }
            if (workers_mpi.firsty && workers_mpi.lasty) { //Just one worker in y-direction
              workers_mpi.offset_y = 0;
              workers_mpi.rows_y   = workers_mpi.rows[Y];
              workers_mpi.end[Y]   = workers_mpi.rows_y - 4;
            }
          } else {
            workers_mpi.rows_y    = workers_mpi.rows[Y] + 6;
            workers_mpi.end[Y]    = workers_mpi.rows_y  - 4;
            workers_mpi.offset_y   = 3;
          }
        } else {
          workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
          workers_mpi.rows_y   = workers_mpi.rows[Y] + 6;
          workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
          workers_mpi.end[Y]   = workers_mpi.rows_y  - 4;
          workers_mpi.offset_x = 3;
          workers_mpi.offset_y = 3;
        }
        
        workers_mpi.layer_size  = workers_mpi.rows_y*workers_mpi.rows_z;
        
        if (DIMENSION == 2) {
          workers_mpi.rows_z     = 1;
          workers_mpi.rows[Z]    = 1;
          workers_mpi.start[Z]   = 0; 
          workers_mpi.end[Z]     = 0;
          workers_mpi.layer_size = workers_mpi.rows_y;
        }
        
        if(workers_mpi.firstx || workers_mpi.firsty || workers_mpi.lastx || workers_mpi.lasty) {
          boundary_worker = 1;
          assign_boundary_points_mpi();
        }
        
        index_count = workers_mpi.layer_size*workers_mpi.rows_x;
        
        gridinfo_w = (struct fields* )malloc((index_count)*sizeof(struct fields));
        
        for (index=0; index < index_count; index++) {
          allocate_memory_fields(&gridinfo_w[index]);
        }
        
        int layer;
        
        for (layer=0; layer < 4; layer++) {
          gradient1[layer] = (struct gradlayer *)malloc((workers_mpi.layer_size)*(sizeof(*gradient1[layer])));
          for (index=0; index < workers_mpi.layer_size; index++) {
            allocate_memory_gradlayer(&gradient1[layer][index]);
          }
        }
        
        gradient = gradient1+1;
        
        if (RESTART == 0) {
          for (x=0; x < rows[X]; x++) {
            for (y=0; y < rows[Y]; y++) {
              index_w             = (x+workers_mpi.offset_x)*workers_mpi.layer_size + (y+workers_mpi.offset_y);
              index               = (x+offset[X])*layer_size + (y+offset[Y]);
              for (a=0; a<NUMPHASES; a++) {
                gridinfo_w[index_w].phia[a]     = gridinfo[index].phia[a];
              }
              for (k=0; k<(NUMCOMPONENTS-1);k++) {
                gridinfo_w[index_w].compi[k]    = gridinfo[index].compi[k];
              }
              for (k=0; k<(NUMCOMPONENTS-1);k++) {
                gridinfo_w[index_w].composition[k]    = gridinfo[index].composition[k];
              }
              for (a=0; a<(NUMPHASES);a++) {
                gridinfo_w[index_w].deltaphi[a] = gridinfo[index].deltaphi[a];
              }
              gridinfo_w[index_w].temperature   = gridinfo[index].temperature;
            }
          }
        }
        printf("taskid=%ld, rows_x=%ld, rows_y=%ld, firstx=%d, firsty=%d, lastx=%d, lasty=%d,offset[X]=%ld, offset[Y]=%ld, offset_x=%ld, offset_y=%ld\n",
           taskid, workers_mpi.rows_x, workers_mpi.rows_y, workers_mpi.firstx, workers_mpi.firsty, workers_mpi.lastx, workers_mpi.lasty, workers_mpi.offset[X], workers_mpi.offset[Y], offset_x, offset_y);
        
      
        buffer_boundary_x = (double *)malloc((workers_mpi.layer_size*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
        buffer_boundary_y = (double *)malloc((workers_mpi.rows_x*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    
        MPI_Type_vector(workers_mpi.rows_x, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_b);
        MPI_Type_commit(&MPI_gridinfo_vector_b);
        
        workers_max_min.phi_max           = (double*)malloc(NUMPHASES*sizeof(double));
        workers_max_min.phi_min           = (double*)malloc(NUMPHASES*sizeof(double));
        workers_max_min.mu_max            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        workers_max_min.mu_min            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        workers_max_min.rel_change_phi    = (double*)malloc((NUMPHASES)*sizeof(double));
        workers_max_min.rel_change_mu     = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        
        workers_max_min.INTERFACE_POS_MAX = 0;
        workers_max_min.INTERFACE_POS_MIN = 0;
        
        for (a=0; a<NUMPHASES; a++) {
          workers_max_min.phi_max[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]].phia[a];
          workers_max_min.phi_min[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]].phia[a];
        }
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          workers_max_min.mu_max[k] = 1.0;
          workers_max_min.mu_min[k] = 1.0;
        }
      }
      if (((offset[Y] + rows[Y])%rows_y) == 0) {
        offset[X] = (offset[X] + rows[X]);
      }
      offset[Y] = ((offset[Y] + rows[Y])%rows_y);
      
      firstx=0;
      firsty=0;
      lasty=0;
      lastx=0;
    }
  } else {
    source = MASTER;
    msgtype = BEGIN;
    
    MPI_Recv(workers_mpi.offset, 2, MPI_LONG, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, offset_x=%ld, offset_y=%ld\n",taskid, workers_mpi.offset[X], workers_mpi.offset[Y]);
    MPI_Recv(workers_mpi.rows, 2, MPI_LONG, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, rows[X]=%ld, rows[Y]=%ld\n",taskid, workers_mpi.rows[X], workers_mpi.rows[Y]);
    MPI_Recv(&workers_mpi.left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, left=%d\n",taskid, workers_mpi.left_node);
    MPI_Recv(&workers_mpi.right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, right=%d\n",taskid, workers_mpi.right_node);
    
    MPI_Recv(&workers_mpi.top_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, top=%d\n",taskid, workers_mpi.top_node);
    MPI_Recv(&workers_mpi.bottom_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, bottom=%d\n", taskid, workers_mpi.bottom_node);
    
    MPI_Recv(&workers_mpi.firstx,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.firsty,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.lastx,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.lasty,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.rank_x,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.rank_y,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    workers_mpi.start[X] = 3;
    workers_mpi.start[Y] = 3;
    workers_mpi.start[Z] = 3;
    
    if ((workers_mpi.firstx==1) || (workers_mpi.lastx==1) || (workers_mpi.firsty==1) || (workers_mpi.lasty==1)) {
      if((workers_mpi.firstx==1) || (workers_mpi.lastx==1)) {
        workers_mpi.rows_x = workers_mpi.rows[X] + 3;
        workers_mpi.end[X] = workers_mpi.rows_x  - 4;
        if (workers_mpi.lastx==1) {
          workers_mpi.offset_x = 3; 
        } else {
          workers_mpi.offset_x = 0;
        }
        if (workers_mpi.firstx && workers_mpi.lastx) { //Just one worker in the x-direction
          workers_mpi.offset_x   = 0;
          workers_mpi.rows_x     = workers_mpi.rows[X];
          workers_mpi.end[X]     = workers_mpi.rows_x - 4;
        }
      } else {
        workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
        workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
        workers_mpi.offset_x = 3;
      }
      if ((workers_mpi.firsty==1) || (workers_mpi.lasty==1)) {
        workers_mpi.rows_y     = workers_mpi.rows[Y] + 3;
        workers_mpi.end[Y]     = workers_mpi.rows_y  - 4;
        if (workers_mpi.lasty==1) {
          workers_mpi.offset_y = 3; 
        } else {
          workers_mpi.offset_y = 0;
        }
        if (workers_mpi.firsty && workers_mpi.lasty) { //Just one worker in y-direction
         workers_mpi.offset_y = 0;
         workers_mpi.rows_y   = workers_mpi.rows[Y];
         workers_mpi.end[Y]   = workers_mpi.rows_y - 4;
        }
      } else {
        workers_mpi.rows_y    = workers_mpi.rows[Y] + 6;
        workers_mpi.end[Y]    = workers_mpi.rows_y  - 4;
        workers_mpi.offset_y   = 3;
      }
    } else {
      workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
      workers_mpi.rows_y   = workers_mpi.rows[Y] + 6;
      workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
      workers_mpi.end[Y]   = workers_mpi.rows_y  - 4;
      workers_mpi.offset_x = 3;
      workers_mpi.offset_y = 3;
    }
    
    workers_mpi.layer_size = workers_mpi.rows_y*workers_mpi.rows_z;
  
    if (DIMENSION == 2) {
      workers_mpi.rows_z     = 1;
      workers_mpi.rows[Z]    = 1;
      workers_mpi.start[Z]   = 0; 
      workers_mpi.end[Z]     = 0;
      workers_mpi.layer_size = workers_mpi.rows_y;
    }
    
    if(workers_mpi.firstx || workers_mpi.firsty || workers_mpi.lastx || workers_mpi.lasty) {
      boundary_worker =1;
      assign_boundary_points_mpi();
    }
    
    index_count = workers_mpi.layer_size*workers_mpi.rows_x;  
    
    gridinfo_w = (struct fields* )malloc((index_count)*sizeof(struct fields));
    
    for (index=0; index < index_count; index++) {
      allocate_memory_fields(&gridinfo_w[index]);
    }
    
    buffer = (double *)malloc((workers_mpi.rows[X]*workers_mpi.rows[Y])*(SIZE_STRUCT_FIELDS)*sizeof(double));

    MPI_Datatype MPI_gridinfo_vector;
    
    MPI_Type_vector(workers_mpi.rows[X]*workers_mpi.rows[Y], SIZE_STRUCT_FIELDS, SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector);
    MPI_Type_commit(&MPI_gridinfo_vector);
    

    MPI_Recv(buffer, 1, MPI_gridinfo_vector, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Type_free(&MPI_gridinfo_vector);
    
    printf("Received data, status=%d,taskid=%ld, test=%d, workers_mpi.offset_x=%ld, workers_mpi.offset_y=%ld, workers_mpi.rows[X]=%ld,workers_mpi.rows_y=%ld,index_count=%ld,workers_mpi.layer_size=%ld\n",status.MPI_ERROR,taskid, MPI_SUCCESS, workers_mpi.offset_x, workers_mpi.offset_y, workers_mpi.rows[X],workers_mpi.rows_y,index_count,workers_mpi.layer_size);
    
    
    if (RESTART == 0) {
      for (index=0; index < index_count; index++) {
        for (a=0; a<NUMPHASES; a++) {
          gridinfo_w[index].phia[a] = 0.0;
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index].compi[k] = 0.0;
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index].composition[k] = 0.0;
        }
        for (a=0; a<NUMPHASES; a++) {
          gridinfo_w[index].deltaphi[a] = 0.0;
        }
        gridinfo_w[index].temperature = 0.0;
      }
      j=0;
      for (i=0; i < workers_mpi.rows[X]*workers_mpi.rows[Y]; i++) {
        index = (i/workers_mpi.rows[Y] + workers_mpi.offset_x)*workers_mpi.layer_size + (i%workers_mpi.rows[Y] + workers_mpi.offset_y);
        for (a=0; a<NUMPHASES; a++) {
          gridinfo_w[index].phia[a] = buffer[j];
          j++;
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index].compi[k] = buffer[j];
          j++;
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index].composition[k] = buffer[j];
          j++;
        }
        for (a=0; a<NUMPHASES; a++) {
          gridinfo_w[index].deltaphi[a] = buffer[j];
          j++;
        }
        gridinfo_w[index].temperature = buffer[j];
        j++;
      }
    }
    free(buffer);
    
    buffer_boundary_x = (double *)malloc((workers_mpi.layer_size*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    buffer_boundary_y = (double *)malloc((workers_mpi.rows_x*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    
    MPI_Type_vector(workers_mpi.rows_x, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_b);
    MPI_Type_commit(&MPI_gridinfo_vector_b);
    
    
    int layer;
    
    for (layer=0; layer < 4; layer++) {
      gradient1[layer] = (struct gradlayer *)malloc((workers_mpi.layer_size)*(sizeof(*gradient1[layer])));
      for (index=0; index < workers_mpi.layer_size; index++) {
        allocate_memory_gradlayer(&gradient1[layer][index]);
      }
    }
    
    gradient = gradient1+1;
    
    workers_max_min.phi_max        = (double*)malloc(NUMPHASES*sizeof(double));
    workers_max_min.phi_min        = (double*)malloc(NUMPHASES*sizeof(double));
    workers_max_min.mu_max         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    workers_max_min.mu_min         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    workers_max_min.rel_change_phi = (double*)malloc((NUMPHASES)*sizeof(double));
    workers_max_min.rel_change_mu  = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    
    for (a=0; a<NUMPHASES; a++) {
      workers_max_min.phi_max[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]].phia[a];
      workers_max_min.phi_min[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]].phia[a];
    }
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      workers_max_min.mu_max[k] = 1.0;
      workers_max_min.mu_min[k] = 1.0;
    }
    
    printf("taskid=%ld, rows_x=%ld, rows_y=%ld, firstx=%d, firsty=%d, lastx=%d, lasty=%d,offset[X]=%ld, offset[Y]=%ld, offset_x=%ld, offset_y=%ld\n",
           taskid, workers_mpi.rows_x, workers_mpi.rows_y, workers_mpi.firstx, workers_mpi.firsty, workers_mpi.lastx, workers_mpi.lasty, workers_mpi.offset[X], workers_mpi.offset[Y], workers_mpi.offset_x, workers_mpi.offset_y);
  }
}
// void sendtomaster() {
//   long i;
//   MPI_Send(workers_mpi.offset, 2, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
//   MPI_Send(workers_mpi.rows,   2, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
//     
//   for (i=0; i < workers_mpi.rows[X]; i++) {
//     MPI_Send(gridinfo_w + (i+workers_mpi.offset_x)*workers_mpi.rows_y + workers_mpi.offset_y, workers_mpi.rows[Y], MPI_gridinfo, MASTER, DONE, MPI_COMM_WORLD);
//   }
// }
// void receivefrmworker() {
//   int rank;
//   long i;
//   for (rank=1; rank <numworkers; rank++) {
//     source = rank;
//     msgtype = DONE;
//     MPI_Recv(offset, 2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//     MPI_Recv(rows,   2, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//     
//     for (i=0; i < rows[X]; i++) {
//       MPI_Recv(gridinfo + offset[X]*rows_y + offset[Y] + i*rows_y,  rows[Y], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
//     }
//   }
// }
void mpiexchange_left_right(long taskid) {
 if (workers_mpi.rank_x%2) {
    if (workers_mpi.firstx ==0) {
      fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);

      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(0, 0, 2);
    }
    if (workers_mpi.lastx == 0) {
      fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
  
      MPI_Send(buffer_boundary_x + 6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.right_node, LTAG,   MPI_COMM_WORLD);
      
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      
      MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    }
  } else {
    if (workers_mpi.lastx == 0) {
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo,                 source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
      
      fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
      
      MPI_Send(buffer_boundary_x + 6*workers_mpi.rows_y*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.right_node,    LTAG, MPI_COMM_WORLD);
    }
    if(workers_mpi.firstx ==0) {
      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x + 0*workers_mpi.rows_y*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(0, 0, 2);
      
      fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_x + 3*workers_mpi.rows_y*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.left_node, RTAG,    MPI_COMM_WORLD);
    }
  }
}
void mpiexchange_top_bottom(long taskid) {
 long i;
 if (workers_mpi.rank_y%2) {
    if (workers_mpi.firsty ==0) {
     
      fill_buffer_y(3*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_y + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
      source = workers_mpi.bottom_node;
      msgtype = BTAG;
    
      MPI_Recv(buffer_boundary_y + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y(0, 0, 2);
      
    }
    if (workers_mpi.lasty == 0) {
      
      fill_buffer_y(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
      
      MPI_Send(buffer_boundary_y + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
      
      source  = workers_mpi.top_node;
      msgtype = TTAG;
      
      MPI_Recv(buffer_boundary_y + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);
    }
  } else {
    if (workers_mpi.lasty == 0) {
      
      source = workers_mpi.top_node;
      msgtype = TTAG;

      MPI_Recv(buffer_boundary_y + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);
      
      fill_buffer_y(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
      
      MPI_Send(buffer_boundary_y + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
      
    }
    if(workers_mpi.firsty ==0) {
      source = workers_mpi.bottom_node;
      msgtype = BTAG;
      
      MPI_Recv(buffer_boundary_y + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y(0, 0, 2);
      
      fill_buffer_y(3*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_y + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
    }
  }
}
void mpiboundary_left_right(long taskid) {
  if (workers_mpi.firstx ==1) {
    
    fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
    
    MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);
    source  = workers_mpi.left_node;
    msgtype = LTAG;
    MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x(0, 0, 2);
  }
  if(workers_mpi.lastx ==1) {
    source  = workers_mpi.right_node;
    msgtype = RTAG;
    MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    
    fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
    
    MPI_Send(buffer_boundary_x + 6*workers_mpi.rows_y*SIZE_STRUCT_FIELDS, 3*workers_mpi.rows_y, MPI_gridinfo, workers_mpi.right_node, LTAG, MPI_COMM_WORLD);
  }
}
void mpiboundary_top_bottom(long taskid) {
  long i;
  if (workers_mpi.firsty ==1) {
    fill_buffer_y(3*SIZE_STRUCT_FIELDS, 3, 5);
        
    MPI_Send(buffer_boundary_y + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
    source = workers_mpi.bottom_node;
    msgtype = BTAG;
  
    MPI_Recv(buffer_boundary_y + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_y(0, 0, 2);
    
  }
  if(workers_mpi.lasty ==1) {
    source = workers_mpi.top_node;
    msgtype = TTAG;
    MPI_Recv(buffer_boundary_y + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, source, msgtype, MPI_COMM_WORLD, &status);

    fill_gridinfo_y(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);  
          
    fill_buffer_y(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
    
    MPI_Send(buffer_boundary_y + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_b, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
  }
}
void fill_buffer_x(long start_j, long x_start, long x_end) {
  long j, x, y;
  long index;
  j=start_j;
  for (x=x_start; x <= x_end; x++) {
    for (y=0; y < workers_mpi.rows_y; y++) {
      index = x*workers_mpi.layer_size + y;
      for (a=0; a<NUMPHASES; a++) {
        buffer_boundary_x[j] = gridinfo_w[index].phia[a];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_x[j] = gridinfo_w[index].compi[k];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_x[j] = gridinfo_w[index].composition[k];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        buffer_boundary_x[j] = gridinfo_w[index].deltaphi[a];
        j++;
      }
      buffer_boundary_x[j] = gridinfo_w[index].temperature;
      j++;
    }
  }
}

void fill_buffer_y(long start_j, long y_start, long y_end) {
  long j, x, y;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
    for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + y;
      for (a=0; a<NUMPHASES; a++) {
        buffer_boundary_y[j] = gridinfo_w[index].phia[a];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_y[j] = gridinfo_w[index].compi[k];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_y[j] = gridinfo_w[index].composition[k];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        buffer_boundary_y[j] = gridinfo_w[index].deltaphi[a];
        j++;
      }
      buffer_boundary_y[j] = gridinfo_w[index].temperature;
      j++;
    }
    j += 9*SIZE_STRUCT_FIELDS;
  }
}

void fill_gridinfo_x(long start_j, long x_start, long x_end) {
  long j, x, y;
  long index;
  j = start_j;
  for (x=x_start; x <= x_end; x++) {
    for (y=0; y < workers_mpi.rows_y; y++) {
      index = x*workers_mpi.layer_size + y;
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].phia[a] = buffer_boundary_x[j];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_w[index].compi[k] = buffer_boundary_x[j];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_w[index].composition[k] = buffer_boundary_x[j];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].deltaphi[a] = buffer_boundary_x[j];
        j++;
      }
      gridinfo_w[index].temperature = buffer_boundary_x[j];
      j++;
    }
  }
}

void fill_gridinfo_y(long start_j, long y_start, long y_end) {
  long j, x, y;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
    for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + y;
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].phia[a] = buffer_boundary_y[j];
        j++;
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        gridinfo_w[index].compi[k] = buffer_boundary_y[j];
        j++;
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        gridinfo_w[index].composition[k] = buffer_boundary_y[j];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].deltaphi[a] = buffer_boundary_y[j];
        j++;
      }
      gridinfo_w[index].temperature = buffer_boundary_y[j];
      j++;
    }
    j += 9*SIZE_STRUCT_FIELDS;
  }
}

#endif

