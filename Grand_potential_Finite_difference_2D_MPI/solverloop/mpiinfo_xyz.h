#ifndef MPIINFO_XYZ_H_
#define MPIINFO_XYZ_H_


void fill_buffer_x(long start_j, long x_start, long x_end);
void fill_buffer_y(long start_j, long y_start, long y_end);
void fill_buffer_z(long start_j, long z_start, long z_end);

void fill_gridinfo_x(long start_j, long x_start, long x_end);
void fill_gridinfo_y(long start_j, long y_start, long y_end);
void fill_gridinfo_z(long start_j, long z_start, long z_end);

void fill_buffer_x_stress(long start_j, long x_start, long x_end);
void fill_buffer_y_stress(long start_j, long y_start, long y_end);
void fill_buffer_z_stress(long start_j, long z_start, long z_end);

void fill_gridinfo_x_stress(long start_j, long x_start, long x_end);
void fill_gridinfo_y_stress(long start_j, long y_start, long y_end);
void fill_gridinfo_z_stress(long start_j, long z_start, long z_end);

void Mpiinfo(long taskid) {
//   long averow[DIMENSION];
  long rank;
  long i, j;
  long x, y, z;
//   long averow_y;
  long lastx;
  long firstx;
  long lasty;
  long firsty;
  long lastz;
  long firstz;
  long rank_x;
  long rank_y;
  long rank_z;
  long index;
  long index_count;
  long index_count_1;
  long index_w;
  
  
  lastx=0;
  firstx=0;
  lasty=0;
  firsty=0;
  lastz=0;
  firstz=0;
  
  if (taskid == MASTER) {
    /* Distribute work to workers.  Must first figure out how many rows to */
    /* send and what to do with extra rows.  */
    
    averow[X] = (rows_x)/numworkers_x;
    averow[Y] = (rows_y)/numworkers_y;
    averow[Z] = (rows_z)/numworkers_z;
     
    extra[X]  = (rows_x)%numworkers_x;
    extra[Y]  = (rows_y)%numworkers_y;
    extra[Z]  = (rows_z)%numworkers_z;
    
    offset[X] = 0;
    offset[Y] = 0;
    offset[Z] = 0;
    
    for (rank=0; rank < (numworkers_x*numworkers_y*numworkers_z); rank++) {
      if (extra[X] > 0) {
        rows[X] = (((rank)/(numworkers_y*numworkers_z)) < extra[X]) ? averow[X]+1 : averow[X];
      } else {
        rows[X] = averow[X];
      }
      
      if (extra[Z] > 0) {
       rows[Z]  = (((rank)%(numworkers_y*numworkers_z))/(numworkers_y) < extra[Z]) ? averow[Z]+1 : averow[Z];
      } else {
       rows[Z] = averow[Z];
      }
      
      if (extra[Y] > 0) {
       rows[Y]  = (((rank)%(numworkers_y*numworkers_z))%(numworkers_y) < extra[Y]) ? averow[Y]+1 : averow[Y];
      } else {
       rows[Y] = averow[Y];
      }
      
      /* Tell each worker who its neighbors are, since they must exchange */
      /* data with each other. */  
      if ((rank)/(numworkers_y*numworkers_z) == 0) {
	      left_node = rank + (numworkers_x-1)*(numworkers_y*numworkers_z);
        
        firstx  = 1;
      } else {
	      left_node = rank - numworkers_y*numworkers_z;
      }
      if (((rank)/(numworkers_y*numworkers_z)) == (numworkers_x-1)) {
	      right_node =  rank - (numworkers_x-1)*(numworkers_y*numworkers_z);
        
        lastx  = 1;
      } else {
	      right_node = rank + numworkers_y*numworkers_z;
      }
      //finished initializing neighbors in the x-direction............................................................
      
      //Begun initializing neighbors in the y-direction...............................................................
      if (((rank)%(numworkers_y*numworkers_z))%numworkers_y == 0) {
	      bottom_node = rank + (numworkers_y-1);
        
        firsty  = 1;
      } else {
	      bottom_node = rank - 1;
      }
      if (((rank)%(numworkers_y*numworkers_z))%numworkers_y == (numworkers_y-1)) {
	      top_node   =  rank - (numworkers_y-1);
        
        lasty  = 1;
      } else {
	      top_node   = rank + 1;
      }
      //Finished initializing neighbors in the y-direction..............................................................
      
      //Begun initializing neighbors in the z-direction
      if (((rank)%(numworkers_y*numworkers_z))/numworkers_y == 0) {
        back_node = rank + (numworkers_z-1)*(numworkers_y);

        firstz  = 1;
      } else {
        back_node = rank - numworkers_y;
      }
      if (((rank)%(numworkers_y*numworkers_z))/numworkers_y == (numworkers_z-1)) {

        front_node   =  rank - (numworkers_z-1)*(numworkers_y);

        lastz  = 1;
      } else {
        front_node   = rank  + numworkers_y;
      }
      //Finished initializing neighbors in the z-direction
      
      
       rank_x = (rank)/(numworkers_y*numworkers_z);
       rank_z = ((rank)%(numworkers_y*numworkers_z))/numworkers_y;
       rank_y = ((rank)%(numworkers_y*numworkers_z))%numworkers_y;
      
      
      
      /*  Now send startup information to each worker  */
      if (rank != MASTER) {
        if (rank==1) {
          printf("sending, taskid=%ld, rows[X]=%ld, rows[Y]=%ld\n",rank, rows[X], rows[Y]);
        }
        dest = rank;
        MPI_Send(offset,       3, MPI_LONG, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(rows,         3, MPI_LONG, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&left_node,   1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&right_node,  1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&top_node,    1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&bottom_node, 1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&front_node,  1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&back_node,   1, MPI_INT,  dest, BEGIN, MPI_COMM_WORLD);
        
        MPI_Send(&firstx, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&firsty, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&firstz, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&lastx,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&lasty,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&lastz,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rank_x, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rank_y, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rank_z, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      
        index_count = (rows[Y])*(rows[X])*(rows[Z]);
        
        buffer = (double *)malloc((index_count)*(SIZE_STRUCT_FIELDS)*sizeof(double));
        
        j=0;
        for (x=0; x < rows[X]; x++) {
         for (z=0; z < rows[Z]; z++) {
          for (y=0; y < rows[Y]; y++) {
            index = (x + offset[X])*layer_size + (z+offset[Z])*rows_y + (y+offset[Y]);
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
        start[Y] = 3;
        end[Y]   = 3;
        start[Z] = 3;
        end[Z]   = 3;
        workers_mpi.firstx      = firstx;
        workers_mpi.lastx       = lastx;
        workers_mpi.firsty      = firsty;
        workers_mpi.lasty       = lasty;
        workers_mpi.firstz      = firstz;
        workers_mpi.lastz       = lastz;
        workers_mpi.rank_x      = rank_x;
        workers_mpi.rank_y      = rank_y;
        workers_mpi.rank_z      = rank_z;
        
        workers_mpi.offset[X]   = offset[X];
        workers_mpi.offset[Y]   = offset[Y];
        workers_mpi.offset[Z]   = offset[Z];
        
        workers_mpi.rows[X]     = rows[X];
        workers_mpi.rows[Y]     = rows[Y];
        workers_mpi.rows[Z]     = rows[Z];
        
        workers_mpi.left_node   = left_node;
        workers_mpi.right_node  = right_node;
        workers_mpi.top_node    = top_node;
        workers_mpi.bottom_node = bottom_node;
        workers_mpi.front_node  = front_node;
        workers_mpi.back_node   = back_node;
        
        workers_mpi.rows_x      = rows[X] + 3;
        workers_mpi.rows_y      = rows[Y] + 3; 
        workers_mpi.rows_z      = rows[Z] + 3;
        workers_mpi.start[X]    = start[X];
        workers_mpi.end[X]      = workers_mpi.rows_x - 4;
        workers_mpi.start[Y]    = start[Y];
        workers_mpi.end[Y]      = workers_mpi.rows_y - 4;
        workers_mpi.start[Z]    = start[Z];
        workers_mpi.end[Z]      = workers_mpi.rows_z - 4;
        
        workers_mpi.layer_size  = workers_mpi.rows_y*workers_mpi.rows_z; 
        workers_mpi.offset_x    = 0;
        workers_mpi.offset_y    = 0;
        workers_mpi.offset_z    = 0;
        
        if ((workers_mpi.firstx==1) || (workers_mpi.lastx==1) || (workers_mpi.firsty==1) || (workers_mpi.lasty==1)||(workers_mpi.firstz==1) || (workers_mpi.lastz==1)) {
          
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
            workers_mpi.offset_y  = 3;
          }
          if ((workers_mpi.firstz==1) || (workers_mpi.lastz==1)) {
            workers_mpi.rows_z     = workers_mpi.rows[Z] + 3;
            workers_mpi.end[Z]     = workers_mpi.rows_z  - 4;
            
            if (workers_mpi.lastz==1) {
              workers_mpi.offset_z = 3; 
            } else {
              workers_mpi.offset_z = 0;
            }
            if (workers_mpi.firstz && workers_mpi.lastz) { //Just one worker in z-direction
              workers_mpi.offset_z = 0;
              workers_mpi.rows_z   = workers_mpi.rows[Z];
              workers_mpi.end[Z]   = workers_mpi.rows_z - 4;
            }
          } else {
            workers_mpi.rows_z    = workers_mpi.rows[Z] + 6;
            workers_mpi.end[Z]    = workers_mpi.rows_z  - 4;
            workers_mpi.offset_z   = 3;
          }
        } else {
          workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
          workers_mpi.rows_y   = workers_mpi.rows[Y] + 6;
          workers_mpi.rows_z   = workers_mpi.rows[Z] + 6;
          workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
          workers_mpi.end[Y]   = workers_mpi.rows_y  - 4;
          workers_mpi.end[Z]   = workers_mpi.rows_z  - 4;
          workers_mpi.offset_x = 3;
          workers_mpi.offset_y = 3;
          workers_mpi.offset_z = 3;
        }
        
        workers_mpi.layer_size  = workers_mpi.rows_y*workers_mpi.rows_z;
        
        if (DIMENSION == 2) {
          workers_mpi.rows_z     = 1;
          workers_mpi.rows[Z]    = 1;
          workers_mpi.start[Z]   = 0; 
          workers_mpi.end[Z]     = 0;
          workers_mpi.offset[Z]  = 0;
          workers_mpi.offset_z   = 0;
          workers_mpi.firstz     = 0;
          workers_mpi.lastz      = 0;
          workers_mpi.layer_size = workers_mpi.rows_y;
        }
        
        if(workers_mpi.firstx || workers_mpi.firsty || workers_mpi.lastx || workers_mpi.lasty || workers_mpi.firstz || workers_mpi.lastz) {
          boundary_worker = 1;
          assign_boundary_points_mpi();
        }
        
        index_count = workers_mpi.layer_size*workers_mpi.rows_x;
        
        gridinfo_w = (struct fields* )malloc((index_count)*sizeof(struct fields));
        if (ELASTICITY) {
          iter_gridinfo_w = (struct iter_variables*)malloc((index_count)*(sizeof(*iter_gridinfo_w)));
        }
        
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
           for (z=0; z < rows[Z]; z++) {
            for (y=0; y < rows[Y]; y++) {
              index_w             = (x+workers_mpi.offset_x)*workers_mpi.layer_size + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y+workers_mpi.offset_y);
              index               = (x+offset[X])*layer_size + (z + offset[Z])*rows_y + (y+offset[Y]);
              for (a=0; a<NUMPHASES; a++) {
                gridinfo_w[index_w].phia[a]     = gridinfo[index].phia[a];
              }
              for (k=0; k<(NUMCOMPONENTS-1);k++) {
                gridinfo_w[index_w].compi[k]    = gridinfo[index].compi[k];
              }
              for (k=0; k<(NUMCOMPONENTS-1);k++) {
                gridinfo_w[index_w].composition[k] = gridinfo[index].composition[k];
              }
              for (a=0; a<(NUMPHASES);a++) {
                gridinfo_w[index_w].deltaphi[a] = gridinfo[index].deltaphi[a];
              }
              gridinfo_w[index_w].temperature   = gridinfo[index].temperature;
            }
          }
         }
        }
        printf("taskid=%ld, rows_x=%ld, rows_y=%ld, firstx=%d, firsty=%d, lastx=%d, lasty=%d,offset[X]=%ld, offset[Y]=%ld, offset_x=%ld, offset_y=%ld\n",
           taskid, workers_mpi.rows_x, workers_mpi.rows_y, workers_mpi.firstx, workers_mpi.firsty, workers_mpi.lastx, workers_mpi.lasty, workers_mpi.offset[X], workers_mpi.offset[Y], offset_x, offset_y);
        
      
        buffer_boundary_x = (double *)malloc((workers_mpi.rows_y*workers_mpi.rows_z*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
        buffer_boundary_y = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_z*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
        buffer_boundary_z = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_y*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    
        MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_z, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_b);
        MPI_Type_commit(&MPI_gridinfo_vector_b);
        
        MPI_Type_vector(workers_mpi.rows_y*workers_mpi.rows_x, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_c);
        MPI_Type_commit(&MPI_gridinfo_vector_c);
        
        if (ELASTICITY) {
          buffer_boundary_x_stress = (double *)malloc((workers_mpi.rows_z*workers_mpi.rows_y*12)*9*sizeof(double));
          buffer_boundary_y_stress = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_z*12)*9*sizeof(double));
          buffer_boundary_z_stress = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_y*12)*9*sizeof(double));
          
          MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_z, 3*9, 12*9, MPI_DOUBLE, &MPI_gridinfo_vector_b_stress);
          MPI_Type_commit(&MPI_gridinfo_vector_b_stress);
          
          MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_y, 3*9, 12*9, MPI_DOUBLE, &MPI_gridinfo_vector_c_stress);
          MPI_Type_commit(&MPI_gridinfo_vector_c_stress);
        }
        
        workers_max_min.phi_max           = (double*)malloc(NUMPHASES*sizeof(double));
        workers_max_min.phi_min           = (double*)malloc(NUMPHASES*sizeof(double));
        workers_max_min.mu_max            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        workers_max_min.mu_min            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        workers_max_min.rel_change_phi    = (double*)malloc((NUMPHASES)*sizeof(double));
        workers_max_min.rel_change_mu     = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
        
        workers_max_min.INTERFACE_POS_MAX = 0;
        workers_max_min.INTERFACE_POS_MIN = 0;
        
        for (a=0; a<NUMPHASES; a++) {
          workers_max_min.phi_max[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].phia[a];
          workers_max_min.phi_min[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].phia[a];
        }
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          workers_max_min.mu_max[k] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].compi[k];
          workers_max_min.mu_min[k] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].compi[k];
        }
        for (a=0; a<NUMPHASES; a++) {
          global_max_min.phi_max[a] = workers_max_min.phi_max[a];
          global_max_min.phi_min[a] = workers_max_min.phi_min[a];
        }
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          global_max_min.mu_max[k] = workers_max_min.mu_max[k];
          global_max_min.mu_min[k] = workers_max_min.mu_min[k];
        }
      }
      if ((((offset[Z] + rows[Z])%rows_z) == 0) && (((offset[Y] + rows[Y])%rows_y) == 0)) {
        offset[X] = (offset[X] + rows[X]);
      }
      
      if (((offset[Y] + rows[Y])%rows_y) == 0) {
        offset[Z] = (offset[Z] + rows[Z])%rows_z;
      }
      
      offset[Y] = ((offset[Y] + rows[Y])%rows_y);
      
      firstx=0;
      firsty=0;
      firstz=0;
      lasty=0;
      lastx=0;
      lastz=0;
    }
  } else {
    source = MASTER;
    msgtype = BEGIN;
    
    MPI_Recv(workers_mpi.offset, 3, MPI_LONG, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, offset_x=%ld, offset_y=%ld\n",taskid, workers_mpi.offset[X], workers_mpi.offset[Y]);
    MPI_Recv(workers_mpi.rows, 3, MPI_LONG, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, rows[X]=%ld, rows[Y]=%ld\n",taskid, workers_mpi.rows[X], workers_mpi.rows[Y]);
    MPI_Recv(&workers_mpi.left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, left=%d\n",taskid, workers_mpi.left_node);
    MPI_Recv(&workers_mpi.right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, right=%d\n",taskid, workers_mpi.right_node);
    
    MPI_Recv(&workers_mpi.top_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, top=%d\n",taskid, workers_mpi.top_node);
    MPI_Recv(&workers_mpi.bottom_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, bottom=%d\n", taskid, workers_mpi.bottom_node);
    
    MPI_Recv(&workers_mpi.front_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, front=%d\n",taskid, workers_mpi.front_node);
    MPI_Recv(&workers_mpi.back_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("taskid=%ld, recieved, back=%d\n", taskid, workers_mpi.back_node);
    
    MPI_Recv(&workers_mpi.firstx,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.firsty,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.firstz,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.lastx,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.lasty,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.lastz,    1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Recv(&workers_mpi.rank_x,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.rank_y,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&workers_mpi.rank_z,   1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    
    
    workers_mpi.start[X] = 3;
    workers_mpi.start[Y] = 3;
    workers_mpi.start[Z] = 3;
    
    if ((workers_mpi.firstx==1) || (workers_mpi.lastx==1) || (workers_mpi.firsty==1) || (workers_mpi.lasty==1)||(workers_mpi.firstz==1) 
      || (workers_mpi.lastz==1)) {
          
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
      if ((workers_mpi.firstz==1) || (workers_mpi.lastz==1)) {
        workers_mpi.rows_z     = workers_mpi.rows[Z] + 3;
        workers_mpi.end[Z]     = workers_mpi.rows_z  - 4;
        
        if (workers_mpi.lastz==1) {
          workers_mpi.offset_z = 3; 
        } else {
          workers_mpi.offset_z = 0;
        }
        if (workers_mpi.firstz && workers_mpi.lastz) { //Just one worker in y-direction
          workers_mpi.offset_z = 0;
          workers_mpi.rows_z   = workers_mpi.rows[Z];
          workers_mpi.end[Z]   = workers_mpi.rows_z - 4;
        }
      } else {
        workers_mpi.rows_z    = workers_mpi.rows[Z] + 6;
        workers_mpi.end[Z]    = workers_mpi.rows_z  - 4;
        workers_mpi.offset_z  = 3;
      }
    } else {
      workers_mpi.rows_x   = workers_mpi.rows[X] + 6;
      workers_mpi.rows_y   = workers_mpi.rows[Y] + 6;
      workers_mpi.rows_z   = workers_mpi.rows[Z] + 6;
      workers_mpi.end[X]   = workers_mpi.rows_x  - 4;
      workers_mpi.end[Y]   = workers_mpi.rows_y  - 4;
      workers_mpi.end[Z]   = workers_mpi.rows_z  - 4;
      workers_mpi.offset_x = 3;
      workers_mpi.offset_y = 3;
      workers_mpi.offset_z = 3;
    }
    
    workers_mpi.layer_size  = workers_mpi.rows_y*workers_mpi.rows_z;
    
    if (DIMENSION == 2) {
      workers_mpi.rows_z     = 1;
      workers_mpi.rows[Z]    = 1;
      workers_mpi.start[Z]   = 0; 
      workers_mpi.end[Z]     = 0;
      workers_mpi.offset[Z]  = 0;
      workers_mpi.offset_z   = 0;
      workers_mpi.firstz     = 0;
      workers_mpi.lastz      = 0;
      workers_mpi.layer_size = workers_mpi.rows_y;
    }
    
    if(workers_mpi.firstx || workers_mpi.firsty || workers_mpi.lastx || workers_mpi.lasty || workers_mpi.firstz || workers_mpi.lastz) {
      boundary_worker =1;
      assign_boundary_points_mpi();
    }
    
    index_count = workers_mpi.layer_size*workers_mpi.rows_x;  
    
    gridinfo_w = (struct fields* )malloc((index_count)*sizeof(struct fields));
    if (ELASTICITY) {
      iter_gridinfo_w = (struct iter_variables*)malloc((index_count)*(sizeof(*iter_gridinfo_w)));
    }
    
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
    
    buffer = (double *)malloc((workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z])*(SIZE_STRUCT_FIELDS)*sizeof(double));

    MPI_Datatype MPI_gridinfo_vector;
    
    MPI_Type_vector(workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z], SIZE_STRUCT_FIELDS, SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector);
    MPI_Type_commit(&MPI_gridinfo_vector);
    

    MPI_Recv(buffer, 1, MPI_gridinfo_vector, source, msgtype, MPI_COMM_WORLD, &status);
    
    MPI_Type_free(&MPI_gridinfo_vector);
    
//     printf("Received data, status=%d,taskid=%ld, test=%d, workers_mpi.offset_x=%ld, workers_mpi.offset_y=%ld, workers_mpi.rows[X]=%ld,workers_mpi.rows_y=%ld,index_count=%ld,workers_mpi.layer_size=%ld\n",status.MPI_ERROR,taskid, MPI_SUCCESS, workers_mpi.offset_x, workers_mpi.offset_y, workers_mpi.rows[X],workers_mpi.rows_y,index_count,workers_mpi.layer_size);
//     
    
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
      for (i=0; i < workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z]; i++) {
        if (DIMENSION == 3) {
          index = (i/(workers_mpi.rows[Y]*workers_mpi.rows[Z]) + workers_mpi.offset_x)*workers_mpi.layer_size + 
          ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z])/workers_mpi.rows[Y]) + workers_mpi.offset_z)*workers_mpi.rows_y + 
          ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z]))%workers_mpi.rows[Y] + workers_mpi.offset_y);
        } else {
          index = (i/workers_mpi.rows[Y] + workers_mpi.offset_x)*workers_mpi.layer_size + (i%workers_mpi.rows[Y] + workers_mpi.offset_y);
        }
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
    
    buffer_boundary_x = (double *)malloc((workers_mpi.rows_z*workers_mpi.rows_y*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    buffer_boundary_y = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_z*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    buffer_boundary_z = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_y*12)*(SIZE_STRUCT_FIELDS)*sizeof(double));
    
    MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_z, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_b);
    MPI_Type_commit(&MPI_gridinfo_vector_b);
    
    MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_y, 3*SIZE_STRUCT_FIELDS, 12*SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_c);
    MPI_Type_commit(&MPI_gridinfo_vector_c);
    
    if (ELASTICITY) {
      buffer_boundary_x_stress = (double *)malloc((workers_mpi.rows_z*workers_mpi.rows_y*12)*9*sizeof(double));
      buffer_boundary_y_stress = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_z*12)*9*sizeof(double));
      buffer_boundary_z_stress = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_y*12)*9*sizeof(double));
      
      MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_z, 3*9, 12*9, MPI_DOUBLE, &MPI_gridinfo_vector_b_stress);
      MPI_Type_commit(&MPI_gridinfo_vector_b_stress);
      
      MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_y, 3*9, 12*9, MPI_DOUBLE, &MPI_gridinfo_vector_c_stress);
      MPI_Type_commit(&MPI_gridinfo_vector_c_stress);
    }
    
    
    workers_max_min.phi_max        = (double*)malloc(NUMPHASES*sizeof(double));
    workers_max_min.phi_min        = (double*)malloc(NUMPHASES*sizeof(double));
    workers_max_min.mu_max         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    workers_max_min.mu_min         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    workers_max_min.rel_change_phi = (double*)malloc((NUMPHASES)*sizeof(double));
    workers_max_min.rel_change_mu  = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    
    for (a=0; a<NUMPHASES; a++) {
      workers_max_min.phi_max[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].phia[a];
      workers_max_min.phi_min[a] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].phia[a];
    }
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      workers_max_min.mu_max[k] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].compi[k];
      workers_max_min.mu_min[k] = gridinfo_w[3*workers_mpi.rows_y + workers_mpi.rows[Y]*workers_mpi.rows[Z]].compi[k];
    }
    
    printf("Received data, taskid=%ld, rows_x=%ld, rows_y=%ld, firstx=%d, firsty=%d, firstz=%d, lastx=%d, lasty=%d, lastz=%d, offset[X]=%ld, offset[Y]=%ld, offset[Z]=%ld, offset_x=%ld, offset_y=%ld, offset_z=%ld, workers_mpi.layer_size=%ld\n",
           taskid, workers_mpi.rows_x, workers_mpi.rows_y, workers_mpi.firstx, workers_mpi.firsty, workers_mpi.firstz, workers_mpi.lastx, workers_mpi.lasty, workers_mpi.lastz, workers_mpi.offset[X], workers_mpi.offset[Y], workers_mpi.offset[Z], workers_mpi.offset_x, workers_mpi.offset_y, workers_mpi.offset_z, workers_mpi.layer_size);
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
      
      MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);

      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(0, 0, 2);
    }
    if (workers_mpi.lastx == 0) {
      fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
  
      MPI_Send(buffer_boundary_x + 6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.right_node, LTAG,   MPI_COMM_WORLD);
      
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      
      MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    }
  } else {
    if (workers_mpi.lastx == 0) {
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
      
      fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
      
      MPI_Send(buffer_boundary_x + 6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.right_node,    LTAG, MPI_COMM_WORLD);
    }
    if(workers_mpi.firstx ==0) {
      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x(0, 0, 2);
      
      fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.left_node, RTAG,    MPI_COMM_WORLD);
    }
  }
}
void mpiexchange_left_right_stress(long taskid) {
 if (workers_mpi.rank_x%2) {
    if (workers_mpi.firstx ==0) {
      fill_buffer_x_stress(3*workers_mpi.layer_size*9, 3, 5);
      
      MPI_Send(buffer_boundary_x_stress + 3*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);

      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x_stress + 0*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x_stress(0, 0, 2);
    }
    if (workers_mpi.lastx == 0) {
      fill_buffer_x_stress(6*workers_mpi.layer_size*9, workers_mpi.end[X]-2, workers_mpi.end[X]);
  
      MPI_Send(buffer_boundary_x_stress + 6*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.right_node, LTAG,   MPI_COMM_WORLD);
      
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      
      MPI_Recv(buffer_boundary_x_stress + 9*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x_stress(9*workers_mpi.layer_size*9, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    }
  } else {
    if (workers_mpi.lastx == 0) {
      source  = workers_mpi.right_node;
      msgtype = RTAG;
      MPI_Recv(buffer_boundary_x_stress + 9*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x_stress(9*workers_mpi.layer_size*9, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
      
      fill_buffer_x_stress(6*workers_mpi.layer_size*9, workers_mpi.end[X]-2, workers_mpi.end[X]);
      
      MPI_Send(buffer_boundary_x_stress + 6*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.right_node,    LTAG, MPI_COMM_WORLD);
    }
    if(workers_mpi.firstx ==0) {
      source  = workers_mpi.left_node;
      msgtype = LTAG;
      MPI_Recv(buffer_boundary_x_stress + 0*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_x_stress(0, 0, 2);
      
      fill_buffer_x_stress(3*workers_mpi.layer_size*9, 3, 5);
      
      MPI_Send(buffer_boundary_x_stress + 3*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.left_node, RTAG,    MPI_COMM_WORLD);
    }
  }
}
// void mpiexchange_left_right(long taskid) {
//  if (workers_mpi.rank_x%2) {
//     if (workers_mpi.firstx ==0) {
//       fill_buffer_x(3*SIZE_STRUCT_FIELDS, 3, 5);
//       
//       MPI_Send(buffer_boundary_x + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_a, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);
// 
//       source  = workers_mpi.left_node;
//       msgtype = LTAG;
//       MPI_Recv(buffer_boundary_x + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_a, source, msgtype, MPI_COMM_WORLD, &status);
//       
//       fill_gridinfo_x(0, 0, 2);
//     }
//     if (workers_mpi.lastx == 0) {
//       fill_buffer_x(6*workers_mpi.layer_size, workers_mpi.end[X]-2, workers_mpi.end[X]);
//   
//       MPI_Send(buffer_boundary_x + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_a, workers_mpi.right_node, LTAG,   MPI_COMM_WORLD);
//       
//       source  = workers_mpi.right_node;
//       msgtype = RTAG;
//       
//       MPI_Recv(buffer_boundary_x + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_a, source, msgtype, MPI_COMM_WORLD, &status);
//       
//       fill_gridinfo_x(9*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
//     }
//   } else {
//     if (workers_mpi.lastx == 0) {
//       source  = workers_mpi.right_node;
//       msgtype = RTAG;
//       MPI_Recv(buffer_boundary_x + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_a, source, msgtype, MPI_COMM_WORLD, &status);
//       
//       fill_gridinfo_x(9*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
//       
//       fill_buffer_x(6*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
//       
//       MPI_Send(buffer_boundary_x + 6*workers_mpi.rows_y*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.right_node,    LTAG, MPI_COMM_WORLD);
//     }
//     if(workers_mpi.firstx ==0) {
//       source  = workers_mpi.left_node;
//       msgtype = LTAG;
//       MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
//       
//       fill_gridinfo_x(0, 0, 2);
//       
//       fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
//       
//       MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.left_node, RTAG,    MPI_COMM_WORLD);
//     }
//   }
// }

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
void mpiexchange_top_bottom_stress(long taskid) {
 long i;
 if (workers_mpi.rank_y%2) {
    if (workers_mpi.firsty ==0) {
     
      fill_buffer_y_stress(3*9, 3, 5);
      
      MPI_Send(buffer_boundary_y_stress + 3*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
      source = workers_mpi.bottom_node;
      msgtype = BTAG;
    
      MPI_Recv(buffer_boundary_y_stress + 0*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y_stress(0, 0, 2);
      
    }
    if (workers_mpi.lasty == 0) {
      
      fill_buffer_y_stress(6*9, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
      
      MPI_Send(buffer_boundary_y_stress + 6*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
      
      source  = workers_mpi.top_node;
      msgtype = TTAG;
      
      MPI_Recv(buffer_boundary_y_stress + 9*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y_stress(9*9, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);
    }
  } else {
    if (workers_mpi.lasty == 0) {
      
      source = workers_mpi.top_node;
      msgtype = TTAG;

      MPI_Recv(buffer_boundary_y_stress + 9*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y_stress(9*9, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);
      
      fill_buffer_y_stress(6*9, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
      
      MPI_Send(buffer_boundary_y_stress + 6*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
      
    }
    if(workers_mpi.firsty ==0) {
      source = workers_mpi.bottom_node;
      msgtype = BTAG;
      
      MPI_Recv(buffer_boundary_y_stress + 0*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_y_stress(0, 0, 2);
      
      fill_buffer_y_stress(3*9, 3, 5);
      
      MPI_Send(buffer_boundary_y_stress + 3*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
    }
  }
}

void mpiexchange_front_back(long taskid) {
 long i;
 if (workers_mpi.rank_z%2) {
    if (workers_mpi.firstz ==0) {
     
      fill_buffer_z(3*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_z + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.back_node,   TTAG, MPI_COMM_WORLD);
      source = workers_mpi.back_node;
      msgtype = BTAG;
    
      MPI_Recv(buffer_boundary_z + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z(0, 0, 2);
      
    }
    if (workers_mpi.lastz == 0) {
      
      fill_buffer_z(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
      
      MPI_Send(buffer_boundary_z + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.front_node,   BTAG, MPI_COMM_WORLD);
      
      source  = workers_mpi.front_node;
      msgtype = TTAG;
      
      MPI_Recv(buffer_boundary_z + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);
    }
  } else {
    if (workers_mpi.lastz == 0) {
      
      source = workers_mpi.front_node;
      msgtype = TTAG;

      MPI_Recv(buffer_boundary_z + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);
      
      fill_buffer_z(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
      
      MPI_Send(buffer_boundary_z + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.front_node,   BTAG, MPI_COMM_WORLD);
      
    }
    if(workers_mpi.firstz ==0) {
      source = workers_mpi.back_node;
      msgtype = BTAG;
      
      MPI_Recv(buffer_boundary_z + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z(0, 0, 2);
      
      fill_buffer_z(3*SIZE_STRUCT_FIELDS, 3, 5);
      
      MPI_Send(buffer_boundary_z + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.back_node,   TTAG, MPI_COMM_WORLD);
    }
  }
}
void mpiexchange_front_back_stress(long taskid) {
 long i;
 if (workers_mpi.rank_z%2) {
    if (workers_mpi.firstz ==0) {
     
      fill_buffer_z_stress(3*9, 3, 5);
      
      MPI_Send(buffer_boundary_z_stress + 3*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.back_node,   TTAG, MPI_COMM_WORLD);
      source = workers_mpi.back_node;
      msgtype = BTAG;
    
      MPI_Recv(buffer_boundary_z_stress + 0*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z_stress(0, 0, 2);
      
    }
    if (workers_mpi.lastz == 0) {
      
      fill_buffer_z_stress(6*9, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
      
      MPI_Send(buffer_boundary_z_stress + 6*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.front_node,   BTAG, MPI_COMM_WORLD);
      
      source  = workers_mpi.front_node;
      msgtype = TTAG;
      
      MPI_Recv(buffer_boundary_z_stress + 9*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z_stress(9*9, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);
    }
  } else {
    if (workers_mpi.lastz == 0) {
      
      source = workers_mpi.front_node;
      msgtype = TTAG;

      MPI_Recv(buffer_boundary_z_stress + 9*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z_stress(9*9, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);
      
      fill_buffer_z_stress(6*9, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
      
      MPI_Send(buffer_boundary_z_stress + 6*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.front_node,   BTAG, MPI_COMM_WORLD);
      
    }
    if(workers_mpi.firstz ==0) {
      source = workers_mpi.back_node;
      msgtype = BTAG;
      
      MPI_Recv(buffer_boundary_z_stress + 0*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);
      
      fill_gridinfo_z_stress(0, 0, 2);
      
      fill_buffer_z_stress(3*9, 3, 5);
      
      MPI_Send(buffer_boundary_z_stress + 3*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.back_node,   TTAG, MPI_COMM_WORLD);
    }
  }
}
void mpiboundary_left_right(long taskid) {
  if (workers_mpi.firstx ==1) {
//     printf("rank_x=%d, rank=%ld, first_x=%d, last_x=%d, workers_mpi.left_node=%d, workers_mpi.layer_size=%ld\n", workers_mpi.rank_x, taskid, workers_mpi.firstx, workers_mpi.lastx, workers_mpi.left_node, workers_mpi.layer_size);
    fill_buffer_x(3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3, 5);
    
    MPI_Send(buffer_boundary_x + 3*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);
    source  = workers_mpi.left_node;
    msgtype = LTAG;
    MPI_Recv(buffer_boundary_x + 0*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x(0, 0, 2);
  }
  if(workers_mpi.lastx ==1) {
//      printf("rank_x=%d, rank=%ld, first_x=%d, last_x=%d, workers_mpi.right_node=%d, workers_mpi.layer_size=%ld\n", workers_mpi.rank_x, taskid, workers_mpi.firstx, workers_mpi.lastx, workers_mpi.right_node, workers_mpi.layer_size);
    source  = workers_mpi.right_node;
    msgtype = RTAG;
    MPI_Recv(buffer_boundary_x + 9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x(9*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    
    fill_buffer_x(6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, workers_mpi.end[X]-2, workers_mpi.end[X]);
    
    MPI_Send(buffer_boundary_x + 6*workers_mpi.layer_size*SIZE_STRUCT_FIELDS, 3*workers_mpi.layer_size, MPI_gridinfo, workers_mpi.right_node, LTAG, MPI_COMM_WORLD);
  }
}
void mpiboundary_left_right_stress(long taskid) {
  if (workers_mpi.firstx ==1) {
//     printf("rank_x=%d, rank=%ld, first_x=%d, last_x=%d, workers_mpi.left_node=%d, workers_mpi.layer_size=%ld\n", workers_mpi.rank_x, taskid, workers_mpi.firstx, workers_mpi.lastx, workers_mpi.left_node, workers_mpi.layer_size);
    fill_buffer_x_stress(3*workers_mpi.layer_size*9, 3, 5);
    
    MPI_Send(buffer_boundary_x_stress + 3*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.left_node, RTAG, MPI_COMM_WORLD);
    source  = workers_mpi.left_node;
    msgtype = LTAG;
    MPI_Recv(buffer_boundary_x_stress + 0*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x_stress(0, 0, 2);
  }
  if(workers_mpi.lastx ==1) {
//      printf("rank_x=%d, rank=%ld, first_x=%d, last_x=%d, workers_mpi.right_node=%d, workers_mpi.layer_size=%ld\n", workers_mpi.rank_x, taskid, workers_mpi.firstx, workers_mpi.lastx, workers_mpi.right_node, workers_mpi.layer_size);
    source  = workers_mpi.right_node;
    msgtype = RTAG;
    MPI_Recv(buffer_boundary_x_stress + 9*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_x_stress(9*workers_mpi.layer_size*9, workers_mpi.end[X]+1, workers_mpi.end[X]+3);
    
    fill_buffer_x_stress(6*workers_mpi.layer_size*9, workers_mpi.end[X]-2, workers_mpi.end[X]);
    
    MPI_Send(buffer_boundary_x_stress + 6*workers_mpi.layer_size*9, 3*workers_mpi.layer_size, MPI_iter_gridinfo, workers_mpi.right_node, LTAG, MPI_COMM_WORLD);
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
void mpiboundary_top_bottom_stress(long taskid) {
  long i;
  if (workers_mpi.firsty ==1) {
    fill_buffer_y_stress(3*9, 3, 5);
        
    MPI_Send(buffer_boundary_y_stress + 3*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.bottom_node,   TTAG, MPI_COMM_WORLD);
    source = workers_mpi.bottom_node;
    msgtype = BTAG;
  
    MPI_Recv(buffer_boundary_y_stress + 0*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_y_stress(0, 0, 2);
    
  }
  if(workers_mpi.lasty ==1) {
    source = workers_mpi.top_node;
    msgtype = TTAG;
    MPI_Recv(buffer_boundary_y_stress + 9*9, 1, MPI_gridinfo_vector_b_stress, source, msgtype, MPI_COMM_WORLD, &status);

    fill_gridinfo_y_stress(9*9, workers_mpi.end[Y]+1, workers_mpi.end[Y]+3);  
          
    fill_buffer_y_stress(6*9, workers_mpi.end[Y]-2, workers_mpi.end[Y]);
    
    MPI_Send(buffer_boundary_y_stress + 6*9, 1, MPI_gridinfo_vector_b_stress, workers_mpi.top_node,   BTAG, MPI_COMM_WORLD);
  }
}

void mpiboundary_front_back(long taskid) {
  long i;
  if (workers_mpi.firstz ==1) {
    fill_buffer_z(3*SIZE_STRUCT_FIELDS, 3, 5);
        
    MPI_Send(buffer_boundary_z + 3*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.back_node,   FTAG, MPI_COMM_WORLD);
    source = workers_mpi.back_node;
    msgtype = BATAG;
  
    MPI_Recv(buffer_boundary_z + 0*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_z(0, 0, 2);
    
  }
  if(workers_mpi.lastz ==1) {
    source = workers_mpi.front_node;
    msgtype = FTAG;
    MPI_Recv(buffer_boundary_z + 9*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, source, msgtype, MPI_COMM_WORLD, &status);

    fill_gridinfo_z(9*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);  
          
    fill_buffer_z(6*SIZE_STRUCT_FIELDS, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
    
    MPI_Send(buffer_boundary_z + 6*SIZE_STRUCT_FIELDS, 1, MPI_gridinfo_vector_c, workers_mpi.front_node, BATAG, MPI_COMM_WORLD);
  }
}
void mpiboundary_front_back_stress(long taskid) {
  long i;
  if (workers_mpi.firstz ==1) {
    fill_buffer_z_stress(3*9, 3, 5);
        
    MPI_Send(buffer_boundary_z_stress + 3*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.back_node,   FTAG, MPI_COMM_WORLD);
    source = workers_mpi.back_node;
    msgtype = BATAG;
  
    MPI_Recv(buffer_boundary_z_stress + 0*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);
    
    fill_gridinfo_z_stress(0, 0, 2);
    
  }
  if(workers_mpi.lastz ==1) {
    source = workers_mpi.front_node;
    msgtype = FTAG;
    MPI_Recv(buffer_boundary_z_stress + 9*9, 1, MPI_gridinfo_vector_c_stress, source, msgtype, MPI_COMM_WORLD, &status);

    fill_gridinfo_z_stress(9*9, workers_mpi.end[Z]+1, workers_mpi.end[Z]+3);  
          
    fill_buffer_z_stress(6*9, workers_mpi.end[Z]-2, workers_mpi.end[Z]);
    
    MPI_Send(buffer_boundary_z_stress + 6*9, 1, MPI_gridinfo_vector_c_stress, workers_mpi.front_node, BATAG, MPI_COMM_WORLD);
  }
}

void fill_buffer_x(long start_j, long x_start, long x_end) {
  long j, x, y, z;
  long index;
  j=start_j;
  for (x=x_start; x <= x_end; x++) {
    for (z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y < workers_mpi.rows_y; y++) {
        index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
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
}
void fill_buffer_x_stress(long start_j, long x_start, long x_end) {
  long j, x, y, z;
  long index;
  int dim, k;
  j=start_j;
  for (x=x_start; x <= x_end; x++) {
    for (z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y < workers_mpi.rows_y; y++) {
        index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
        for (dim=0; dim < 3; dim++) {
          for (k=0; k < 3; k++) {
            buffer_boundary_x_stress[j] = iter_gridinfo_w[index].disp[dim][k];
            j++;
          }
        }
      }
    }
  }
}
void fill_buffer_y(long start_j, long y_start, long y_end) {
  long j, x, y, z;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (z=0; z < workers_mpi.rows_z; z++) {
    for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
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
}
void fill_buffer_y_stress(long start_j, long y_start, long y_end) {
  long j, x, y, z;
  long index;
  long k;
  int dim;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (z=0; z < workers_mpi.rows_z; z++) {
    for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
      for (dim=0; dim < 3; dim++) {
        for (k=0; k < 3; k++) {
          buffer_boundary_y_stress[j] = iter_gridinfo_w[index].disp[dim][k];
          j++;
         }
       }
     }
     j += 9*9;
    }
  }
}
void fill_buffer_z(long start_j, long z_start, long z_end) {
  long j, x, y, z;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (y=0; y < workers_mpi.rows_y; y++) {
     for (z=z_start; z <= z_end; z++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
      for (a=0; a < NUMPHASES; a++) {
        buffer_boundary_z[j] = gridinfo_w[index].phia[a];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_z[j] = gridinfo_w[index].compi[k];
        j++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_boundary_z[j] = gridinfo_w[index].composition[k];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        buffer_boundary_z[j] = gridinfo_w[index].deltaphi[a];
        j++;
      }
      buffer_boundary_z[j] = gridinfo_w[index].temperature;
      j++;
     }
    j += 9*SIZE_STRUCT_FIELDS;
   }
  }
}
void fill_buffer_z_stress(long start_j, long z_start, long z_end) {
  long j, x, y, z;
  long index;
  long k;
  int dim;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (y=0; y < workers_mpi.rows_y; y++) {
     for (z=z_start; z <= z_end; z++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
       for (dim=0; dim < 3; dim++) {
        for (k=0; k < 3; k++) {
          buffer_boundary_z_stress[j] = iter_gridinfo_w[index].disp[dim][k];
          j++;
         }
       }
     }
    j += 9*9;
   }
  }
}

void fill_gridinfo_x(long start_j, long x_start, long x_end) {
  long j, x, y, z;
  long index;
  j = start_j;
  for (x=x_start; x <= x_end; x++) {
    for (z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y < workers_mpi.rows_y; y++) {
        index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
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
}
void fill_gridinfo_x_stress(long start_j, long x_start, long x_end) {
  long j, x, y, z;
  long index;
  long k;
  int dim;
  j = start_j;
  for (x=x_start; x <= x_end; x++) {
    for (z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y < workers_mpi.rows_y; y++) {
        index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
        for (dim=0; dim < 3; dim++) {
          for (k=0; k < 3; k++) {
            iter_gridinfo_w[index].disp[dim][k] = buffer_boundary_x_stress[j];
            j++;
          }
        }
      }
    }
  }
}
void fill_gridinfo_y(long start_j, long y_start, long y_end) {
  long j, x, y, z;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (z=0; z < workers_mpi.rows_z; z++) {
     for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
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
}
void fill_gridinfo_y_stress(long start_j, long y_start, long y_end) {
  long j, x, y, z;
  long index;
  long k;
  int dim;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (z=0; z < workers_mpi.rows_z; z++) {
     for (y=y_start; y <= y_end; y++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
      for (dim=0; dim < 3; dim++) {
        for (k=0; k < 3; k++) {
          iter_gridinfo_w[index].disp[dim][k] = buffer_boundary_y_stress[j] ;
          j++;
         }
       }
      }
     j += 9*9;
    }
  }
}
void fill_gridinfo_z(long start_j, long z_start, long z_end) {
  long j, x, y, z;
  long index;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (y=0; y < workers_mpi.rows_y; y++) {
     for (z=z_start; z <= z_end; z++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].phia[a] = buffer_boundary_z[j];
        j++;
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        gridinfo_w[index].compi[k] = buffer_boundary_z[j];
        j++;
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        gridinfo_w[index].composition[k] = buffer_boundary_z[j];
        j++;
      }
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_w[index].deltaphi[a] = buffer_boundary_z[j];
        j++;
      }
      gridinfo_w[index].temperature = buffer_boundary_z[j];
      j++;
     }
     j += 9*SIZE_STRUCT_FIELDS;
    }
  }
}
void fill_gridinfo_z_stress(long start_j, long z_start, long z_end) {
  long j, x, y, z;
  long index;
  long k;
  int dim;
  j = start_j;
  for (x=0; x < workers_mpi.rows_x; x++) {
   for (y=0; y < workers_mpi.rows_y; y++) {
     for (z=z_start; z <= z_end; z++) {
      index = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
      for (dim=0; dim < 3; dim++) {
        for (k=0; k < 3; k++) {
          iter_gridinfo_w[index].disp[dim][k] = buffer_boundary_z_stress[j] ;
          j++;
         }
       }
      }
     j += 9*9;
    }
  }
}
#endif

