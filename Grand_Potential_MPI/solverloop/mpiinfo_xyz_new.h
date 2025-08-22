#ifndef MPIINFO_NEW_XYZ_H_
#define MPIINFO_NEW_XYZ_H_

#include <mpi.h>

void ms_mpi_Init_numworkers(int argc, char * argv[])
{
  int mpi_size[3] = {0,0,0};

  if (argc == 4) {
    fprintf(stderr, "mpi - creating numworkers");
    ms_mpi_try(
      MPI_Dims_create(numtasks, (int)DIMENSION, mpi_size)
    );
  }
  else {
    assert(((DIMENSION == 2) && (argc >= 6)) || (DIMENSION == 3) && (argc == 7));
    mpi_size[X] = atoi(argv[4]);
    mpi_size[Y] = atoi(argv[5]);
    mpi_size[Z] = (DIMENSION == 2) ? 1 : atoi(argv[6]);
  }

  numworkers_x = mpi_size[X];
  numworkers_y = mpi_size[Y];
  numworkers_z = (DIMENSION == 2) ? 1 : mpi_size[Z];  
  
  if (!taskid)
  {
    fprintf(stderr, "numworkers = %d %d %d\n", numworkers_x, numworkers_y, numworkers_z);
  }

  if(numtasks != numworkers_x*numworkers_y*numworkers_z)
  {
    fprintf(stderr, "**ERROR** in numworkers.\n- Run as:\n  mpirun -np <num procs> ./microsim_gp InpFile FillFile OutName numworkers_x numworkers_y numworkers_z\n- Eg for 2D : mpirun -np 4 ./microsim_gp Input.in Fill.in out-data 2 2\n- Eg for 3D : mpirun -np 8 ./microsim_gp Input.in Fill.in out-data 2 2 2\n\n");
    ms_Abort();
  }
}

int ms_mpi_GetRankRelative(int * coord, int x, int y, int z, MPI_Comm COMM)
{
  int newrank, newcoord[3];
  newcoord[X] = coord[X] + x;
  newcoord[Y] = coord[Y] + y;
  newcoord[Z] = coord[Z] + z;
  MPI_Cart_rank(COMM, newcoord, &newrank);
  return newrank;
}

void ms_mpi_CreateUniverse() {

  if (!taskid)fprintf(stderr, "\nUSE_NEW_MPIINFO = True\n");

  int j, x, y, z, index, index_w;

  /**
   * Decompose uniformly
   */

  rows_x     = MESH_X + 6;
  rows_y     = MESH_Y + 6;
  rows_z     = MESH_Z + 6;

  averow[X] = (rows_x)/numworkers_x;
  averow[Y] = (rows_y)/numworkers_y;
  averow[Z] = (rows_z)/numworkers_z;
    
  extra[X]  = (rows_x)%numworkers_x;
  extra[Y]  = (rows_y)%numworkers_y;
  extra[Z]  = (rows_z)%numworkers_z;

  if (DIMENSION == 2)
  {
    averow[Z]     = 1;
    extra[Z]      = 0;
    numworkers_z  = 1;
  }
  
  /**
   * Declare topology
   */

  int mpi_topo_periods[3]  = {1,1,1};
  int mpi_self_coord[3]    = {0,0,0};
  int mpi_topo_dims[3]     = {numworkers_x, numworkers_y, numworkers_z};
  
  MPI_Cart_create (MPI_COMM_WORLD, 3, mpi_topo_dims, mpi_topo_periods, true, &mpi_topo_comm);

  MPI_Comm_rank   (mpi_topo_comm, &taskid);       // resetting taskid
  MPI_Comm_size   (mpi_topo_comm, &numtasks);     // resetting numtasks just in case
  
  if (taskid == MASTER)
  {
    workers_mpi_all = msMalloc(struct workers, numtasks);
  }

  for(int rank = 0; rank < numtasks; rank++)
  {
    if( rank == taskid || taskid == MASTER )
    {
      MPI_Cart_coords (mpi_topo_comm, rank, 3, mpi_self_coord);
      
      workers_mpi.rank   = rank;
      workers_mpi.rank_x = mpi_self_coord[X];
      workers_mpi.rank_y = mpi_self_coord[Y];
      workers_mpi.rank_z = mpi_self_coord[Z];

      workers_mpi.left_node   = ms_mpi_GetRankRelative(mpi_self_coord, -1,  0,  0, mpi_topo_comm);
      workers_mpi.right_node  = ms_mpi_GetRankRelative(mpi_self_coord,  1,  0,  0, mpi_topo_comm);

      workers_mpi.bottom_node = ms_mpi_GetRankRelative(mpi_self_coord,  0, -1,  0, mpi_topo_comm);
      workers_mpi.top_node    = ms_mpi_GetRankRelative(mpi_self_coord,  0,  1,  0, mpi_topo_comm);

      workers_mpi.front_node  = ms_mpi_GetRankRelative(mpi_self_coord,  0,  0,  1, mpi_topo_comm);
      workers_mpi.back_node   = ms_mpi_GetRankRelative(mpi_self_coord,  0,  0, -1, mpi_topo_comm);

      workers_mpi.firstx  = mpi_self_coord[X] == 0;
      workers_mpi.lastx   = mpi_self_coord[X] == (mpi_topo_dims[X]-1);
      workers_mpi.firsty  = mpi_self_coord[Y] == 0;
      workers_mpi.lasty   = mpi_self_coord[Y] == (mpi_topo_dims[Y]-1);
      workers_mpi.firstz  = mpi_self_coord[Z] == 0;
      workers_mpi.lastz   = mpi_self_coord[Z] == (mpi_topo_dims[Z]-1);
      
      workers_mpi.rows[X] = workers_mpi.rank_x < extra[X] ? averow[X] + 1 : averow[X];
      workers_mpi.rows[Y] = workers_mpi.rank_y < extra[Y] ? averow[Y] + 1 : averow[Y];
      workers_mpi.rows[Z] = workers_mpi.rank_z < extra[Z] ? averow[Z] + 1 : averow[Z];

      workers_mpi.rows_x  =  (workers_mpi.firstx || workers_mpi.lastx) ? workers_mpi.rows[X]+3 : workers_mpi.rows[X]+6;
      workers_mpi.rows_y  =  (workers_mpi.firsty || workers_mpi.lasty) ? workers_mpi.rows[Y]+3 : workers_mpi.rows[Y]+6;
      workers_mpi.rows_z  =  (workers_mpi.firstz || workers_mpi.lastz) ? workers_mpi.rows[Z]+3 : workers_mpi.rows[Z]+6;

      workers_mpi.start[X]    = 3;
      workers_mpi.start[Y]    = 3;
      workers_mpi.start[Z]    = 3;

      workers_mpi.offset[X] = workers_mpi.rank_x < extra[X] ? workers_mpi.rank_x*(averow[X] + 1) : workers_mpi.rank_x*averow[X] + extra[X];
      workers_mpi.offset[Y] = workers_mpi.rank_y < extra[Y] ? workers_mpi.rank_y*(averow[Y] + 1) : workers_mpi.rank_y*averow[Y] + extra[Y];
      workers_mpi.offset[Z] = workers_mpi.rank_z < extra[Z] ? workers_mpi.rank_z*(averow[Z] + 1) : workers_mpi.rank_z*averow[Z] + extra[Z];

      workers_mpi.offset_x =  workers_mpi.firstx ?  0 : 3;
      workers_mpi.offset_y =  workers_mpi.firsty ?  0 : 3;
      workers_mpi.offset_z =  workers_mpi.firstz ?  0 : 3;

      if (workers_mpi.firstx && workers_mpi.lastx)
      {
        workers_mpi.offset_x   = 0;
        workers_mpi.rows_x     = workers_mpi.rows[X];
      }

      if (workers_mpi.firsty && workers_mpi.lasty)
      {
        workers_mpi.offset_y = 0;
        workers_mpi.rows_y   = workers_mpi.rows[Y];
      }

      if (workers_mpi.firstz && workers_mpi.lastz)
      {
        workers_mpi.offset_z = 0;
        workers_mpi.rows_z   = workers_mpi.rows[Z];
      }

      workers_mpi.end[X]	= workers_mpi.rows_x - 4;
      workers_mpi.end[Y]	= workers_mpi.rows_y - 4;
      workers_mpi.end[Z]	= workers_mpi.rows_z - 4;

      workers_mpi.layer_size  = workers_mpi.rows_y * workers_mpi.rows_z;
      workers_mpi.index_count = workers_mpi.layer_size * workers_mpi.rows_x;
    }

    if (taskid == MASTER)
    {
      workers_mpi_all[rank] = workers_mpi;
    }
  }

  if(taskid == MASTER)
  {
    workers_mpi = workers_mpi_all[taskid];
  }

  if (!taskid)fprintf(stderr, "\nCREATED workers_mpi\n");

  /**
   * print info
   */
//   print_workers(taskid, &workers_mpi);

//   if(taskid == MASTER){
//     // printf("\n- MASTER printing workers_mpi_all");
//     for(int i  = 0; i < numtasks; i++){
//       print_workers(taskid, &workers_mpi_all[i]);
//     }
//   }

  /**
   * Set boundary options
   */

  boundary_worker = workers_mpi.firstx || workers_mpi.firsty || workers_mpi.lastx || workers_mpi.lasty || workers_mpi.firstz || workers_mpi.lastz;

  if(boundary_worker)
  {
    assign_boundary_points_mpi();
  }

  /**
   * Allocate Fields
   */
	GLOBAL_FIELDS_SIZE = 2*(NUMPHASES) + 2*(NUMCOMPONENTS-1);
	LBM_FIELDS_SIZE    = (DIMENSION == 2) ? 9 : 27;

	gridinfo_w				= msCalloc(struct fields, 	workers_mpi.index_count);
	GLOBAL_FIELDS_ALL	= msCalloc(double,					workers_mpi.index_count*GLOBAL_FIELDS_SIZE);

	for (long i = 0; i < workers_mpi.index_count; i++)
	{
		double *tmp = GLOBAL_FIELDS_ALL + i*GLOBAL_FIELDS_SIZE;
		gridinfo_w[i].phia        = tmp;
		gridinfo_w[i].compi       = tmp + NUMPHASES;
		gridinfo_w[i].composition = tmp + NUMPHASES   + (NUMCOMPONENTS-1);
		gridinfo_w[i].deltaphi    = tmp + NUMPHASES + 2*(NUMCOMPONENTS-1);
	}	
  
	if (ELASTICITY)
	{
		iter_gridinfo_w = msCalloc(struct iter_variables, workers_mpi.index_count);
	}

	if (LBM)
	{
		lbm_gridinfo_w 	= msCalloc(struct lbm_fields, workers_mpi.index_count);
		LBM_FIELDS_ALL 	= msCalloc(double,   	  			workers_mpi.index_count*(LBM_FIELDS_SIZE));
		for (long i = 0; i < workers_mpi.index_count; i++)
		{
			lbm_gridinfo_w[i].fdis = LBM_FIELDS_ALL + i*(LBM_FIELDS_SIZE);
		}

		gradc = ms_Malloc_3D(workers_mpi.layer_size, DIMENSION, NUMCOMPONENTS-1);
	}
 

  if (!taskid)fprintf(stderr, "\nAllocated fields\n");

  for (int layer=0; layer < 4; layer++)
  {
    gradient1[layer] = (struct gradlayer *)malloc((workers_mpi.layer_size)*(sizeof(*gradient1[layer])));
    for (int index=0; index < workers_mpi.layer_size; index++)
    {
      allocate_memory_gradlayer(&gradient1[layer][index]);
    }
  }

  if (!taskid)fprintf(stderr, "\nAllocated gradlayer\n");

  gradient = gradient1+1;

  /**
   * Allocate boundary layers
   */

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
  
  if (LBM) {
    buffer_boundary_x_lbm = (double *)malloc((workers_mpi.rows_z*workers_mpi.rows_y*12)*SIZE_LBM_FIELDS*sizeof(double));
    buffer_boundary_y_lbm = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_z*12)*SIZE_LBM_FIELDS*sizeof(double));
    buffer_boundary_z_lbm = (double *)malloc((workers_mpi.rows_x*workers_mpi.rows_y*12)*SIZE_LBM_FIELDS*sizeof(double));
    
    MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_z, 3*SIZE_LBM_FIELDS, 12*SIZE_LBM_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_b_lbm);
    MPI_Type_commit(&MPI_gridinfo_vector_b_lbm);
    
    MPI_Type_vector(workers_mpi.rows_x*workers_mpi.rows_y, 3*SIZE_LBM_FIELDS, 12*SIZE_LBM_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_c_lbm);
    MPI_Type_commit(&MPI_gridinfo_vector_c_lbm);
  }

  if (!taskid)fprintf(stderr, "\nAllocated boundary buffers\n");

  /**
   * Allocate max min things
   */

  workers_max_min.phi_max           = (double*)malloc(NUMPHASES*sizeof(double));
  workers_max_min.phi_min           = (double*)malloc(NUMPHASES*sizeof(double));
  workers_max_min.mu_max            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  workers_max_min.mu_min            = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  workers_max_min.rel_change_phi    = (double*)malloc((NUMPHASES)*sizeof(double));
  workers_max_min.rel_change_mu     = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));

  workers_max_min.rel_change_phi_max    = (double*)malloc((NUMPHASES)*sizeof(double));
  workers_max_min.rel_change_mu_max     = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));


  // printf("\n- [%d] Field alloc done", taskid);   
  /**
   * Master neeeds to send fields to workers
   */

  if(taskid == MASTER && RESTART == 0)
  {
	// For all ranks
    for(int rank = 0; rank < numtasks; rank ++)
    { 
      // Do not send master to master
      if(rank == MASTER){
        continue;
      }
      dest = rank;

      struct workers rankInfo = workers_mpi_all[rank];

      index_count = rankInfo.rows[Z] * rankInfo.rows[Y] * rankInfo.rows[X];

      // allocate the buffer;         
      buffer = (double *)malloc((index_count)*(SIZE_STRUCT_FIELDS)*sizeof(double));

      // fill the buffer
      j = 0;
      for (x=0; x < rankInfo.rows[X]; x++) {
        for (z=0; z < rankInfo.rows[Z]; z++) {
          for (y=0; y < rankInfo.rows[Y]; y++) {
            
            index = (x + rankInfo.offset[X])*rows_z*rows_y + (z+rankInfo.offset[Z])*rows_y + (y+rankInfo.offset[Y]);
            
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

      // printf("\n- [%d] in master sending to rank %d", taskid, dest);  

      // Send buffer
      MPI_Datatype MPI_gridinfo_vector;
      MPI_Type_vector(index_count, SIZE_STRUCT_FIELDS, SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector);
      MPI_Type_commit(&MPI_gridinfo_vector);
      MPI_Send(buffer, 1, MPI_gridinfo_vector, dest, BEGIN, mpi_topo_comm);
      MPI_Type_free(&MPI_gridinfo_vector); 

      // free buffer
      free(buffer);
      // printf("\n- [%d] in master sending to rank ss %d", taskid, dest); 

      if(LBM) {
        buffer_lbm = (double *)malloc((index_count)*(SIZE_LBM_FIELDS)*sizeof(double));
        j=0;
        for (x=0; x < rankInfo.rows[X]; x++) {
          for (z=0; z < rankInfo.rows[Z]; z++) {
            for (y=0; y < rankInfo.rows[Y]; y++) {
              index = (x + rankInfo.offset[X])*layer_size + (z+rankInfo.offset[Z])*rows_y + (y+rankInfo.offset[Y]);
              for (k=0; k < (SIZE_LBM_FIELDS-4); k++) {
                buffer_lbm[j] = lbm_gridinfo[index].fdis[k];
                j++;
              }
              buffer_lbm[j] = lbm_gridinfo[index].rho;
              j++;
              for(dim=0; dim < 3; dim++) {
                buffer_lbm[j] = lbm_gridinfo[index].u[dim];
                j++;
              }
            }
          }
        }
        // printf("\n- [%d] in master sending to rank %d lbm", taskid, dest); 
        MPI_Datatype MPI_gridinfo_vector_lbm;
        MPI_Type_vector(index_count, SIZE_LBM_FIELDS, SIZE_LBM_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_lbm);
        MPI_Type_commit(&MPI_gridinfo_vector_lbm);
        MPI_Send(buffer_lbm, 1, MPI_gridinfo_vector_lbm, dest, BEGIN, mpi_topo_comm);
        MPI_Type_free(&MPI_gridinfo_vector_lbm);
        
        free(buffer_lbm);

        // printf("\n- [%d] in master sending to rank lbm done %d", taskid, dest); 
      }
      // printf("\n- [0] master has sent to [%d] ", rank); 
    } // rank for loop ends
  } // taskid == MASTER) ends

  if (taskid != MASTER && RESTART == 0)
  {
    /**
     * Workers need to recieve from master and fill
     */
    source = MASTER;
    msgtype = BEGIN;
    /**
     * Recieve the buffer
     */

    buffer = (double *)malloc((workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z])*(SIZE_STRUCT_FIELDS)*sizeof(double));

    MPI_Datatype MPI_gridinfo_vector;
    MPI_Type_vector(workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z], SIZE_STRUCT_FIELDS, SIZE_STRUCT_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector);
    MPI_Type_commit(&MPI_gridinfo_vector);
    MPI_Recv(buffer, 1, MPI_gridinfo_vector, source, msgtype, mpi_topo_comm, &status);
    MPI_Type_free(&MPI_gridinfo_vector);
    // printf("\n- [%d] recieved fields ", taskid); 

    if(LBM) {
      // printf("\n- [%d] recieved fields lbm", taskid); 
      buffer_lbm = (double *)malloc((workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z])*(SIZE_LBM_FIELDS)*sizeof(double));

      MPI_Datatype MPI_gridinfo_vector_lbm;
      MPI_Type_vector(workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z], SIZE_LBM_FIELDS, SIZE_LBM_FIELDS, MPI_DOUBLE, &MPI_gridinfo_vector_lbm);
      MPI_Type_commit(&MPI_gridinfo_vector_lbm);
      MPI_Recv(buffer_lbm, 1, MPI_gridinfo_vector_lbm, source, msgtype, mpi_topo_comm, &status);
      MPI_Type_free(&MPI_gridinfo_vector_lbm);
      // printf("\n- [%d] recieved fields lbm done ", taskid); 
    }
    // printf("\n- [%d] Fields recieved", taskid);  

    // Copy fields from buffer
    j = 0;
    for (int i=0; i < workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z]; i++) {
      index = (i/(workers_mpi.rows[Y]*workers_mpi.rows[Z]) + workers_mpi.offset_x)*workers_mpi.layer_size + 
        ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z])/workers_mpi.rows[Y]) + workers_mpi.offset_z)*workers_mpi.rows_y + 
        ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z]))%workers_mpi.rows[Y] + workers_mpi.offset_y);

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
    free(buffer);
    
    if(LBM) {
      j=0;
      for (int i=0; i < workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z]; i++) {
        index = (i/(workers_mpi.rows[Y]*workers_mpi.rows[Z]) + workers_mpi.offset_x)*workers_mpi.layer_size + 
          ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z])/workers_mpi.rows[Y]) + workers_mpi.offset_z)*workers_mpi.rows_y + 
          ((i%(workers_mpi.rows[Y]*workers_mpi.rows[Z]))%workers_mpi.rows[Y] + workers_mpi.offset_y);
        for (int k=0; k<(SIZE_LBM_FIELDS-4); k++) {
          lbm_gridinfo_w[index].fdis[k] = buffer_lbm[j];
          j++;
        }
        lbm_gridinfo_w[index].rho = buffer_lbm[j];
        j++;
        for (int dim=0; dim<3; dim++) {
          lbm_gridinfo_w[index].u[dim] = buffer_lbm[j];
          j++;
        }
      }
      free(buffer_lbm);
    }
  }// taskid != MASTER ends

  /**
   * Master needs to fill from itself
   */
  
  if(taskid == MASTER && RESTART == 0)
  {
    // printf("\n- [%d] MASTER Field filling", taskid);  

    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index_w             = (x+workers_mpi.offset_x)*workers_mpi.layer_size + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y+workers_mpi.offset_y);
          index               = (x+workers_mpi.offset[X])*rows_y*rows_z + (z + workers_mpi.offset[Z])*rows_y + (y+workers_mpi.offset[Y]);
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

    if(LBM) {
      for (x=0; x < workers_mpi.rows[X]; x++) {
        for (z=0; z < workers_mpi.rows[Z]; z++) {
          for (y=0; y < workers_mpi.rows[Y]; y++) {
            index_w             =  (x+workers_mpi.offset_x)*workers_mpi.layer_size + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y+workers_mpi.offset_y);
            index               =  (x+workers_mpi.offset[X])*rows_y*rows_z + (z + workers_mpi.offset[Z])*rows_y + (y+workers_mpi.offset[Y]);
            for (k=0; k<(SIZE_LBM_FIELDS-4); k++) {
              lbm_gridinfo_w[index_w].fdis[k]     = lbm_gridinfo[index].fdis[k];
            }
            lbm_gridinfo_w[index_w].rho = lbm_gridinfo[index].rho;
            for(dim=0; dim < 3; dim++) {
              lbm_gridinfo_w[index_w].u[dim] = lbm_gridinfo[index].u[dim];
            }
          }
        }
      }
    }

    // printf("\n- [%d] MASTER Field filled done", taskid);  
  }// taskid == MASTER ends


  /**
   * Fill other stuff
   */

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
  if(LBM){
    workers_max_min.rho_max = -10;
    workers_max_min.rho_min = 10;
  }

  if (taskid == MASTER){
    for (a=0; a<NUMPHASES; a++) {
      global_max_min.phi_max[a] = workers_max_min.phi_max[a];
      global_max_min.phi_min[a] = workers_max_min.phi_min[a];
    }
    for (k=0; k<NUMCOMPONENTS-1; k++) {
      global_max_min.mu_max[k] = workers_max_min.mu_max[k];
      global_max_min.mu_min[k] = workers_max_min.mu_min[k];
    }
    if(LBM){
      global_max_min.rho_max = workers_max_min.rho_max;
      global_max_min.rho_min = workers_max_min.rho_min;
    }
  }
  printf("\n- [%d] Function done", taskid);  

}


#endif