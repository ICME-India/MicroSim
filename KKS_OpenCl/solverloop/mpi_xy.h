void mpi_xy(int rank, char *argv[]) {
  
  
  int cart_grid_dims[cart_grid_ndim];
  int cart_grid_periods[cart_grid_ndim];
  int cart_grid_reorder; 
  int cart_grid_coordns[cart_grid_ndim];
  int cart_grid_directn[cart_grid_ndim];
  int cart_grid_displac;
  //int dim0_rank_p, dim0_rank_m;

  int count_x, remainder_x, global_startx, global_endx;
  int x, y, z, gx, gy, gz, index, idx, b_idx, global_index, idxa, idxb, i1, i2;

  
  cart_grid_dims[0] = numtasks;
  cart_grid_periods[0] = 1;
  cart_grid_directn[0] = 0;
  cart_grid_reorder = 0;
  cart_grid_displac = 1;
  
  MPI_Cart_create(comm, cart_grid_ndim, cart_grid_dims, cart_grid_periods, cart_grid_reorder, &COMM_NEW);

  MPI_Cart_coords(COMM_NEW, rank, cart_grid_ndim, cart_grid_coordns);

  MPI_Cart_rank(COMM_NEW, cart_grid_coordns, &rank);

  MPI_Cart_shift(COMM_NEW, cart_grid_directn[0], cart_grid_displac, &dim0_rank_m, &dim0_rank_p);

  printf("rank-, rank, rank+: %d,%d,%d\n", dim0_rank_m, rank, dim0_rank_p);

  count_x = (int) rows_x/numtasks; 
  
  remainder_x = rows_x % numtasks;

  //printf("Rank = %d, count_x = %d, remainder_x = %d\n", rank, count_x, remainder_x);

  if ( rank < remainder_x ) { 
    // The first 'remainder' ranks get 'count + 1' tasks each
    global_startx = rank * (count_x + 1); 
    global_endx   = global_startx + count_x + 1;

    mpiparam.rows_x = ( global_endx - global_startx ) + 2;

  }
  else {
    // The remaining 'size - remainder' ranks get 'count' task each 
    global_startx = rank * count_x + remainder_x; 
    global_endx   = global_startx + (count_x - 1) + 1;

    mpiparam.rows_x = ( global_endx - global_startx ) + 2;

  }
  MPI_Barrier(COMM_NEW);
  
  mpiparam.rows_y = rows_y + 2;
  mpiparam.rows_z = rows_z + 2;

  printf("Size in local: Rank = %d, y, x, z: %ld, %ld, %ld, global_startx = %d, global_endx = %d \n", rank, mpiparam.rows_y, mpiparam.rows_x, mpiparam.rows_z, global_startx, global_endx);


  nx = mpiparam.rows_x;
  ny = mpiparam.rows_y;
  nz = mpiparam.rows_z;

  mpiparam.startx = 0;
  mpiparam.endx = mpiparam.rows_x; 

  mpiparam.starty = 0; 
  mpiparam.endy = mpiparam.rows_y;
  
  mpiparam.startz = 0; 
  mpiparam.endz = mpiparam.rows_z;

  sizempixyz = mpiparam.rows_x*mpiparam.rows_y*mpiparam.rows_z;
  sizebufyz = mpiparam.rows_y*mpiparam.rows_z;

  gridinfomN = (struct fields* )malloc((sizempixyz)*sizeof(*gridinfomN));
  gridinfomO = (struct fields* )malloc((sizempixyz)*sizeof(*gridinfomO));

//   if (ELASTICITY) {
    iter_gridinfom = (struct iter_variables*)malloc((sizempixyz)*(sizeof(*iter_gridinfom)));
//   }
//   else {
//     iter_gridinfom = (struct iter_variables*)malloc((1)*(sizeof(*iter_gridinfom)));
//   }



  if ((STARTTIME == 0) && (RESTART == 0)) {
    printf("Starting a new simulation\n");

    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          index = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          for ( i1 = 0; i1 < 3; i1++ ) {
            for ( i2 = 0; i2 < 3; i2++ ) {
              iter_gridinfom[index].disp[i1][i2] = 0.0;
            }
          }
        }
      }
    }

    //printf("%d, %d, %d, %d, %le, %le, %le\n", 2, 2, 2, index, iter_gridinfom[62].disp[X][2], iter_gridinfom[62].disp[Y][2], iter_gridinfom[62].disp[Z][2]);
    gx = global_startx;
    for ( x = 1; x < mpiparam.rows_x-1; x++ ) {
      gz = 0;
      for ( z = 1; z < mpiparam.rows_z-1; z++ ) {
        gy = 0;
        for ( y = 1; y < mpiparam.rows_y-1; y++ ) {
          index = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          global_index = gy + rows_y*(gz + gx*rows_z);
          gridinfomN[index] = gridinfo[global_index];
          gridinfomO[index] = gridinfomN[index];
          gy++;
          if ( gy > rows_y ) {
            printf("gy > rows_y: %d > %ld\n", gy, rows_y);
            exit(1);
          }
        }
        gz++;
        if ( gz > rows_z ) {
          printf("gz > rows_z: %d > %ld\n", gz, rows_z);
          exit(1);
        }
      }
      gx++;
      if ( gx > rows_x ) {
        printf("gx > rows_x: %d > %ld\n", gx, rows_x);
        exit(1);
      }
    }
  }
  else {
    printf("Restarting a simulation\n");
    if (ASCII) {
      readfromfile_serialmpi(gridinfomN, argv, STARTTIME);
    } else {
      readfromfile_serialmpi_binary(gridinfomN, argv, STARTTIME);
    }
  }
  MPI_Barrier(COMM_NEW);

  gridinfo_buf1Dss = (struct fields* )malloc((sizebufyz)*sizeof(*gridinfo_buf1Dss));
  gridinfo_buf1Dse = (struct fields* )malloc((sizebufyz)*sizeof(*gridinfo_buf1Dse));
  gridinfo_buf1Drs = (struct fields* )malloc((sizebufyz)*sizeof(*gridinfo_buf1Drs));
  gridinfo_buf1Dre = (struct fields* )malloc((sizebufyz)*sizeof(*gridinfo_buf1Dre));

  iter_gridinfo_buf2Dss = (struct iter_variables* )malloc((2*sizebufyz)*sizeof(*iter_gridinfo_buf2Dss));
  iter_gridinfo_buf2Dse = (struct iter_variables* )malloc((2*sizebufyz)*sizeof(*iter_gridinfo_buf2Dse));
  iter_gridinfo_buf2Drs = (struct iter_variables* )malloc((2*sizebufyz)*sizeof(*iter_gridinfo_buf2Drs));
  iter_gridinfo_buf2Dre = (struct iter_variables* )malloc((2*sizebufyz)*sizeof(*iter_gridinfo_buf2Dre));

  cscl_buf1Dss     = (struct csle*)malloc(sizebufyz*sizeof(struct csle));
  cscl_buf1Dse     = (struct csle*)malloc(sizebufyz*sizeof(struct csle));
  cscl_buf1Dss     = (struct csle*)malloc(sizebufyz*sizeof(struct csle));
  cscl_buf1Dse     = (struct csle*)malloc(sizebufyz*sizeof(struct csle));

  size_buffer = sizebufyz*SIZE_STRUCT_FIELDS;
    
  buffer_sxs = (double *)malloc((size_buffer)*sizeof(double));
  buffer_sxe = (double *)malloc((size_buffer)*sizeof(double));
  buffer_rxs = (double *)malloc((size_buffer)*sizeof(double));
  buffer_rxe = (double *)malloc((size_buffer)*sizeof(double));

  iter_size_buffer = 2*sizebufyz*9;
  iter_buffer_sxs = (double *)malloc((iter_size_buffer)*sizeof(double));
  iter_buffer_sxe = (double *)malloc((iter_size_buffer)*sizeof(double));
  iter_buffer_rxs = (double *)malloc((iter_size_buffer)*sizeof(double));
  iter_buffer_rxe = (double *)malloc((iter_size_buffer)*sizeof(double));

  x = 1;
  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
      for (a=0; a<NUMPHASES; a++) {
        buffer_sxs[b_idx] = gridinfomN[idx].phia[a];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxs[b_idx] = gridinfomN[idx].compi[k];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxs[b_idx] = gridinfomN[idx].composition[k];
        b_idx++;
      }
      buffer_sxs[b_idx] = gridinfomN[idx].temperature;
      b_idx++;
    }
  }

  x = mpiparam.rows_x-2;
  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
      for (a=0; a<NUMPHASES; a++) {
        buffer_sxe[b_idx] = gridinfomN[idx].phia[a];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxe[b_idx] = gridinfomN[idx].compi[k];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxe[b_idx] = gridinfomN[idx].composition[k];
        b_idx++;
      }
      buffer_sxe[b_idx] = gridinfomN[idx].temperature;
      b_idx++;
    }
  }

  MPI_Sendrecv(&buffer_sxs[0], size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_m, &buffer_rxe[0], size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_m, COMM_NEW, &status1);
  MPI_Sendrecv(&buffer_sxe[0], size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_p, &buffer_rxs[0], size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_p, COMM_NEW, &status2);
  MPI_Barrier(COMM_NEW);
  
  x = 0;
  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
      for (a=0; a<NUMPHASES; a++) {
        gridinfomN[idx].phia[a] = buffer_rxs[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfomN[idx].compi[k] = buffer_rxs[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfomN[idx].composition[k] = buffer_rxs[b_idx];
        b_idx++;
      }
      gridinfomN[idx].temperature = buffer_rxs[b_idx];
      b_idx++;
    }
  }
  
  x = mpiparam.rows_x-1;
  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
      for (a=0; a<NUMPHASES; a++) {
        gridinfomN[idx].phia[a] = buffer_rxe[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfomN[idx].compi[k] = buffer_rxe[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfomN[idx].composition[k] = buffer_rxe[b_idx];
        b_idx++;
      }
      gridinfomN[idx].temperature = buffer_rxe[b_idx];
      b_idx++;
    }
  }
  MPI_Barrier(COMM_NEW);
  
  // If not PERIODIC 
  if (rank == 0) { 
    x = 0;
    if ( boundary[1][0].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(1));
          for (a=0; a<NUMPHASES; a++) {
            gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
          }
        }
      }
    }
    if ( boundary[1][1].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(1));
          for (k=0; k<NUMCOMPONENTS-1; k++) {
            gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
            gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
          }
        }
      }
    }
    if ( boundary[1][2].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(1));
          gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
        }
      }
    }
  }

  if (rank == numtasks-1) { 
    x = mpiparam.rows_x-1;
    if ( boundary[0][0].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(nx-2));
          for (a=0; a<NUMPHASES; a++) {
            gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
          }
        }
      }
    }
    if ( boundary[0][1].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(nx-2));
          for (k=0; k<NUMCOMPONENTS-1; k++) {
            gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
            gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
          }
        }
      }
    }
    if ( boundary[0][2].type == 1 ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          idxb = y + mpiparam.rows_y*(z + mpiparam.rows_z*(nx-2));
          gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
        }
      }
    }
  }

  // Boundary condition to accommodate boundary layers along y (row_y + 2)
  y = 0; // at y- phi
  if ( boundary[3][0].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (1) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }
  else if ( boundary[3][0].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }
  
  y = mpiparam.rows_y-1; // at y+ phi 
  if ( boundary[2][0].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }
  else if ( boundary[2][0].type == 3 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = 1 + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }

  y = 0; // at y- composition and mu
  if ( boundary[3][1].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (1) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }
  else if ( boundary[3][1].type == 3 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }
  
  y = mpiparam.rows_y-1; // at y+ composition and mu
  if ( boundary[2][1].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }
  else if ( boundary[2][1].type == 3 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = 1 + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }


  y = 0; // at y- temperature
  if ( boundary[3][2].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (1) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  else if ( boundary[3][2].type == 3 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  
  y = mpiparam.rows_y-1; // at y+ temperature 
  if ( boundary[2][2].type == 1 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = (ny-2) + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  else if ( boundary[2][2].type == 3 ) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = 1 + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  



  // Boundary condition to accommodate boundary layers along z (row_z + 2)
  z = 0; // at z- phi
  if ( boundary[5][0].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((1) + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }
  else if ( boundary[5][0].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }

  z = mpiparam.rows_z-1; // at z+ phi
  if ( boundary[4][0].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }
  else if ( boundary[4][0].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*(1 + mpiparam.rows_z*x);
        for (a=0; a<NUMPHASES; a++) {
          gridinfomN[idxa].phia[a] = gridinfomN[idxb].phia[a];
        }
      }
    }
  }

  z = 0; // at z- composition and mu
  if ( boundary[5][1].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((1) + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }
  else if ( boundary[5][1].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }

  z = mpiparam.rows_z-1; // at z+ composition and mu
  if ( boundary[4][1].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }
  else if ( boundary[4][1].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*(1 + mpiparam.rows_z*x);
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          gridinfomN[idxa].compi[k]       = gridinfomN[idxb].compi[k];
          gridinfomN[idxa].composition[k] = gridinfomN[idxb].composition[k];
        }
      }
    }
  }

  z = 0; // at z- temperature
  if ( boundary[5][2].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((1) + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  else if ( boundary[5][2].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }

  z = mpiparam.rows_z-1; // at z+ temperature
  if ( boundary[4][2].type == 1 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*((nz-2) + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }
  else if ( boundary[4][2].type == 3 ) {
    for ( x = 0; x < mpiparam.rows_x; x++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        idxa = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        idxb = y + mpiparam.rows_y*(1 + mpiparam.rows_z*x);
        gridinfomN[idxa].temperature = gridinfomN[idxb].temperature;
      }
    }
  }


  if (ISOTHERMAL) { 
    for ( x = 0; x < mpiparam.rows_x; x++ ) { 
      for ( z = 0; z < mpiparam.rows_z; z++ ) {
        for ( y = 0; y < mpiparam.rows_y; y++ ) {
          index = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
          gridinfomN[index].temperature = T;
        }
      }
    }
  }

  for ( x = 0; x < mpiparam.rows_x; x++ ) {
    for ( z = 0; z < mpiparam.rows_z; z++ ) {
      for ( y = 0; y < mpiparam.rows_y; y++ ) {
        index = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);
        gridinfomO[index] = gridinfomN[index];
      }
    }
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

void mpi_exchange_dim0(int rank) {

  int y, b_idx, idx, z;

  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*z;
      for (a=0; a<NUMPHASES; a++) {
        buffer_sxs[b_idx] = gridinfo_buf1Dss[idx].phia[a];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxs[b_idx] = gridinfo_buf1Dss[idx].compi[k];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxs[b_idx] = gridinfo_buf1Dss[idx].composition[k];
        b_idx++;
      }
      buffer_sxs[b_idx] = gridinfo_buf1Dss[idx].temperature;
      b_idx++;
    }
  }

  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*z;
      for (a=0; a<NUMPHASES; a++) {
        buffer_sxe[b_idx] = gridinfo_buf1Dse[idx].phia[a];
        // if (rank==2)
        // printf("%lf %lf\t", buffer_sxe[b_idx], gridinfo_buf1Dse[y].phia[a]);
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxe[b_idx] = gridinfo_buf1Dse[idx].compi[k];
        // if (rank==2)
        // printf("%lf %lf\t", buffer_sxe[b_idx], gridinfo_buf1Dse[y].compi[k]);
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        buffer_sxe[b_idx] = gridinfo_buf1Dse[idx].composition[k];
        // if (rank==2)
        // printf("%lf %lf\t", buffer_sxe[b_idx], gridinfo_buf1Dse[y].composition[k]);
        b_idx++;
      }
      buffer_sxe[b_idx] = gridinfo_buf1Dse[idx].temperature;
      // if (rank==2)
      //   printf("%lf %lf\t", buffer_sxe[b_idx], gridinfo_buf1Dse[y].temperature);
      b_idx++;
    }
  }

  // if (rank==2) {
  //   for (i=0; i<size_buffer; i++) {
  //     printf("%lf \t", buffer_sxs[i]);
  //   }
  // }
  // printf("\n");

  //printf("##### %d \n ", size_buffer);

  MPI_Sendrecv(&buffer_sxs[0], size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_m, &buffer_rxe[0], size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_m, COMM_NEW, &status1);
  MPI_Sendrecv(&buffer_sxe[0], size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_p, &buffer_rxs[0], size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_p, COMM_NEW, &status2);
  MPI_Barrier(COMM_NEW);

  // if (rank==1) {
  //   for (i=0; i<size_buffer; i++) {
  //     printf("%lf \t", buffer_rxe[i]);
  //   }
  // }

  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*z;
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_buf1Drs[idx].phia[a] = buffer_rxs[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_buf1Drs[idx].compi[k] = buffer_rxs[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_buf1Drs[idx].composition[k] = buffer_rxs[b_idx];
        b_idx++;
      }
      gridinfo_buf1Drs[idx].temperature = buffer_rxs[b_idx];
      b_idx++;
    }
  }

  b_idx = 0;
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*z;
      for (a=0; a<NUMPHASES; a++) {
        gridinfo_buf1Dre[idx].phia[a] = buffer_rxe[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_buf1Dre[idx].compi[k] = buffer_rxe[b_idx];
        b_idx++;
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        gridinfo_buf1Dre[idx].composition[k] = buffer_rxe[b_idx];
        b_idx++;
      }
      gridinfo_buf1Dre[idx].temperature = buffer_rxe[b_idx];
      b_idx++;
    }
  }
  MPI_Barrier(COMM_NEW);

}



void mpi_exchange_dim0_iter(int rank) {

  int y, b_idx, idx, z, i1, i2, x;

  b_idx = 0;
  for ( x = 0; x < 2; x++ ) {
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);

      for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
          iter_buffer_sxs[b_idx] = iter_gridinfo_buf2Dss[idx].disp[i1][i2];
          b_idx++;
        }
      }

    }
  }
  }

  b_idx = 0;
  for ( x = 0; x < 2; x++ ) {
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);

      for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
          iter_buffer_sxe[b_idx] = iter_gridinfo_buf2Dse[idx].disp[i1][i2];
        }
      }

    }
  }
  }

  // if (rank==2) {
  //   for (i=0; i<size_buffer; i++) {
  //     printf("%lf \t", buffer_sxs[i]);
  //   }
  // }
  // printf("\n");

  //printf("##### %d \n ", size_buffer);

  MPI_Sendrecv(&iter_buffer_sxs[0], iter_size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_m, &iter_buffer_rxe[0], iter_size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_m, COMM_NEW, &status1);
  MPI_Sendrecv(&iter_buffer_sxe[0], iter_size_buffer, MPI_DOUBLE, dim0_rank_p, tag_dim0_p, &iter_buffer_rxs[0], iter_size_buffer, MPI_DOUBLE, dim0_rank_m, tag_dim0_p, COMM_NEW, &status2);
  MPI_Barrier(COMM_NEW);

  // if (rank==1) {
  //   for (i=0; i<size_buffer; i++) {
  //     printf("%lf \t", buffer_rxe[i]);
  //   }
  // }

  b_idx = 0;
  for ( x = 0; x < 2; x++ ) {
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);

      for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
          iter_gridinfo_buf2Drs[idx].disp[i1][i2] = iter_buffer_rxs[b_idx];
          b_idx++;
        }
      }

    }
  }
  }

  b_idx = 0;
  for ( x = 0; x < 2; x++ ) {
  for ( z = 0; z < mpiparam.rows_z; z++ ) {
    for ( y = 0; y < mpiparam.rows_y; y++ ) {
      idx = y + mpiparam.rows_y*(z + mpiparam.rows_z*x);

      for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
          iter_gridinfo_buf2Dre[idx].disp[i1][i2] = iter_buffer_rxe[b_idx];
          b_idx++;
        }
      }

    }
  }
  }
  MPI_Barrier(COMM_NEW);

}
