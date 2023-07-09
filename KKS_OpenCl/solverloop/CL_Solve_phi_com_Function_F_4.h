void CL_Solve_phi_com_Function_F_4() {
  
  //printf("In CL_Solve_phi_com_Function_F_4\n");
  
  tstep[0] = t;
  ret  = clEnqueueWriteBuffer(cmdQ, d_tstep, CL_TRUE, 0, sizeof(long), tstep, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_tstep %d\n", ret);
    exit(1);
  }

  if ( (t > tNoiseStart) && NOISE_PHASEFIELD ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_addNoise, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("ker_addNoise enq problem  %d\n",ret);
      exit(1);
    }
  }

  ret = clEnqueueNDRangeKernel(cmdQ, ker_copy_New_To_Old, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("ker_copy_New_To_Old enq problem  %d\n",ret);
    exit(1);
  }

  if ( tstep[0] > nsmooth ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverCsClEq_F4, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("ker_SolverCsClEq_F4 enq problem  %d in rank = %d \n",ret, rank);
      exit(1);
    }
  }

  if ( tstep[0] <= nsmooth ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverPhi_F4_smooth, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverPhi_F4, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F4 enq problem  %d at %ld \n",ret, tstep[0]);
    exit(1);
  }

  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomN, CL_TRUE,      1*nynz*sizeof(struct fields), nynz*sizeof(struct fields), gridinfo_buf1Dss, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomN at 1*nynz for mpi communication\n%d\n", ret);
    exit(1);
  }

  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomN, CL_TRUE, (nx-2)*nynz*sizeof(struct fields), nynz*sizeof(struct fields), gridinfo_buf1Dse, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomN at (nx-2)*nynz for mpi communication\n%d\n", ret);
    exit(1);
  }

  // MPI communication for ghost buffer update
  mpi_exchange_dim0(rank);

  //printf("Done phi exchanged rank = %d\n", rank);

  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, 0, 1*nynz*sizeof(struct fields), gridinfo_buf1Drs, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomN at 0 after mpi communication%d\n", ret);
    exit(1);
  }

  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, (nx-1)*nynz*sizeof(struct fields), 1*nynz*sizeof(struct fields), gridinfo_buf1Dre, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomN at (nx-1)*nynz after mpi communication%d\n", ret);
    exit(1);
  }
  
  
  //BC for x+ i.e. Nx
  if (rank == numtasks-1) {
  if ( boundary[0][0].type == 1 ) { 
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_xn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[0][0].type == 3 ) {
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_xn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_xn enq problem  %d\n",ret);
    exit(1);
  }
  }

  //BC for x- i.e. 0
  if (rank == 0) {
  if ( boundary[1][0].type == 1 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_x0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[1][0].type == 3 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_x0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_x0 enq problem  %d\n",ret);
    exit(1);
  }
  }
  
  //BC for y+ i.e. Ny
  if ( boundary[2][0].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_yn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[2][0].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_yn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_yn enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for y- i.e. 0
  if ( boundary[3][0].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_y0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[3][0].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_y0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_y0 enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for z+ i.e. Nz
  if ( boundary[4][0].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_zn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[4][0].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_zn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  z+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_zn enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for z- i.e. 0
  if ( boundary[5][0].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_z0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[5][0].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_phi_z0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for z- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_phi_z0 enq problem  %d\n",ret);
    exit(1);
  }

  if ( tstep[0] <= nsmooth ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverCatr_F4_smooth, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverCatr_F4, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F4 enq problem  %d at %ld \n",ret, tstep[0]);
    exit(1);
  }

  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomN, CL_TRUE,      1*nynz*sizeof(struct fields), nynz*sizeof(struct fields), gridinfo_buf1Dss, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomN at 1*nynz for mpi communication\n%d\n", ret);
    exit(1);
  }


  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomN, CL_TRUE, (nx-2)*nynz*sizeof(struct fields), nynz*sizeof(struct fields), gridinfo_buf1Dse, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomN at (nx-2)*nynz for mpi communication\n%d\n", ret);
    exit(1);
  }

  // MPI communication for ghost buffer update
  mpi_exchange_dim0(rank);

  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, 0, 1*nynz*sizeof(struct fields), gridinfo_buf1Drs, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomN at 0 after mpi communication%d\n", ret);
    exit(1);
  }

  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, (nx-1)*nynz*sizeof(struct fields), 1*nynz*sizeof(struct fields), gridinfo_buf1Dre, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomN at (nx-1)*nynz after mpi communication%d\n", ret);
    exit(1);
  }
  
  
  //BC for x+ i.e. Nx
  if (rank == numtasks-1) {
  if ( boundary[0][1].type == 1 ) { 
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_xn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[0][1].type == 3 ) {
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_xn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_xn enq problem  %d\n",ret);
    exit(1);
  }
  }

  //BC for x- i.e. 0
  if (rank == 0) {
  if ( boundary[1][1].type == 1 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_x0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[1][1].type == 3 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_x0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_x0 enq problem  %d\n",ret);
    exit(1);
  }
  }
  
  //BC for y+ i.e. Ny
  if ( boundary[2][1].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_yn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[2][1].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_yn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_yn enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for y- i.e. 0
  if ( boundary[3][1].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_y0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[3][1].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_y0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_y0 enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for z+ i.e. Nz
  if ( boundary[4][1].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_zn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[4][1].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_zn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  z+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_zn enq problem  %d\n",ret);
    exit(1);
  }
  
  //BC for z- i.e. 0
  if ( boundary[5][1].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_z0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[5][1].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_com_z0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for z- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_com_z0 enq problem  %d\n",ret);
    exit(1);
  }

  if ( ELASTICITY && ( tstep[0] > nsmooth ) ) {
    for(iter=1; iter <= MAX_ITERATIONS; iter++) {

      if ( DIMENSION == 3 ) {
      ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverStress_iterative, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
      if (ret!=CL_SUCCESS) {
        printf("ker_SolverStress_iterative enq problem  %d at time step %ld \n",ret, tstep[0]);
        exit(1);
      }
      }
      else if ( DIMENSION == 2 ) {
      ret = clEnqueueNDRangeKernel(cmdQ, ker_SolverStress_iterative_2D, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
      if (ret!=CL_SUCCESS) {
        printf("ker_SolverStress_iterative_2D enq problem  %d at time step %ld \n",ret, tstep[0]);
        exit(1);
      }
      }


      //printf("%d, %d, %d\n", 2*nynz, 2*nynz, ny*nx*nz);
      ret = clEnqueueReadBuffer(cmdQ, d_iter_gridinfom, CL_TRUE,      2*nynz*sizeof(struct iter_variables), 2*nynz*sizeof(struct iter_variables), iter_gridinfo_buf2Dss, 0, NULL, NULL);
      if (ret != CL_SUCCESS) {
        printf("Error: Failed to read d_iter_gridinfo_buf2Dss at 2*nynz for mpi communication\n%d\n", ret);
        exit(1);
      }

      ret = clEnqueueReadBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, (nx-4)*nynz*sizeof(struct iter_variables), 2*nynz*sizeof(struct iter_variables), iter_gridinfo_buf2Dse, 0, NULL, NULL);
      if (ret != CL_SUCCESS) {
        printf("Error: Failed to read d_iter_gridinfo_buf2Dse at (nx-4)*nynz for mpi communication\n%d\n", ret);
        exit(1);
      }

      // MPI communication for ghost buffer update
      mpi_exchange_dim0_iter(rank);

      ret  = clEnqueueWriteBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, 0, 2*nynz*sizeof(struct iter_variables), iter_gridinfo_buf2Drs, 0, NULL, NULL);
      if (ret!=CL_SUCCESS) {
        printf("enq buffer write error d_iter_gridinfom at 0 after mpi communication%d\n", ret);
        exit(1);
      }

      ret  = clEnqueueWriteBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, (nx-2)*nynz*sizeof(struct iter_variables), 2*nynz*sizeof(struct iter_variables), iter_gridinfo_buf2Dre, 0, NULL, NULL);
      if (ret!=CL_SUCCESS) {
        printf("enq buffer write error d_iter_gridinfom at (nx-2)*nynz after mpi communication%d\n", ret);
        exit(1);
      }

        //BC for x+ i.e. Nx
  if (rank == numtasks-1) {
  if ( boundary[0][3].type == 1 ) {
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_xn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[0][3].type == 3 ) {
    ret  = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_xn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_xn enq problem  %d\n",ret);
    exit(1);
  }
  }

  //BC for x- i.e. 0
  if (rank == 0) {
  if ( boundary[1][3].type == 1 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_x0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[1][3].type == 3 ) {
    ret |= clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_x0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  x- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_x0 enq problem  %d\n",ret);
    exit(1);
  }
  }

  //BC for y+ i.e. Ny
  if ( boundary[2][3].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_yn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[2][3].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_yn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_yn enq problem  %d\n",ret);
    exit(1);
  }

  //BC for y- i.e. 0
  if ( boundary[3][3].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_y0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[3][3].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_y0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  y- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_y0 enq problem  %d\n",ret);
    exit(1);
  }

  //BC for z+ i.e. Nz
  if ( boundary[4][3].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_zn_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[4][3].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_zn_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for  z+ \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_zn enq problem  %d\n",ret);
    exit(1);
  }

  //BC for z- i.e. 0
  if ( boundary[5][3].type == 1 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_z0_noflux, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else if ( boundary[5][3].type == 3 ) {
    ret = clEnqueueNDRangeKernel(cmdQ, ker_apply_BC_ela_z0_periodic, work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
  }
  else {
    printf("No Boundary condition for z- \n");
    exit(1);
  }
  if (ret!=CL_SUCCESS) {
    printf("ker_apply_BC_ela_z0 enq problem  %d\n",ret);
    exit(1);
  }



    }
  }


  
}
