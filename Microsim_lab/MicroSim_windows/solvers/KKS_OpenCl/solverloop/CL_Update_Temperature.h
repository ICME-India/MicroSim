void CL_Update_Temperature(long t) {
  
  if ( !ISOTHERMAL ) {
    tstep[0] = t;
    //printf("Passed t = %ld\n", t);
    ret  = clEnqueueWriteBuffer(cmdQ, d_tstep, CL_TRUE, 0, sizeof(long), tstep, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("enq buffer write error d_tstep %d\n", ret);
      exit(1);
    }

    ret = clEnqueueNDRangeKernel(cmdQ, kernel10[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel10 enq problem  %d\n",ret);
      exit(1);
    }

    if ( boundary[2][2].type == 1 ) {
      ret = clEnqueueNDRangeKernel(cmdQ, kernel11[0], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    }
    else if ( boundary[2][2].type == 3 ) {
      //printf("Only no flux boundary is defined for Temperature\n");
      //printf("because of the imposition of 1D temperature gradient implementation\n");
      //printf("Using No Flux boundary condition for Temperature and proceeding\n");
      ret = clEnqueueNDRangeKernel(cmdQ, kernel11[0], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    }
    if (ret!=CL_SUCCESS) {
      printf("kernel11[0] enq problem  %d\n",ret);
      exit(1);
    }

    if ( boundary[3][2].type == 1 ) {
      ret = clEnqueueNDRangeKernel(cmdQ, kernel11[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    }
    else if ( boundary[3][2].type == 3 ) {
      //printf("Only no flux boundary is defined for Temperature\n");
      //printf("because of the imposition of 1D temperature gradient implementation\n");
      //printf("Using No Flux boundary condition for Temperature and proceeding\n");
      ret = clEnqueueNDRangeKernel(cmdQ, kernel11[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    }
    if (ret!=CL_SUCCESS) {
      printf("kernel11[1] enq problem  %d\n",ret);
      exit(1);
    }
  }
}