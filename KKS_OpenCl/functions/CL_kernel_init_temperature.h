void CL_kernel_init_temperature() {

  if ( ISOTHERMAL ) { 
    ret = clEnqueueNDRangeKernel(cmdQ, kernel10[0], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel10 enq problem  %d\n",ret);
      exit(1);
    }
  }
  else if ( !ISOTHERMAL ) {
    tstep[0] = tstart;
    ret  = clEnqueueWriteBuffer(cmdQ, d_tstep, CL_TRUE, 0, sizeof(int), tstep, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("enq buffer write error d_tstep %d\n", ret);
      exit(1);
    }

    ret = clEnqueueNDRangeKernel(cmdQ, kernel10[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel10 enq problem  %d\n",ret);
      exit(1);
    }

    ret = clEnqueueNDRangeKernel(cmdQ, kernel11[0], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel13 enq problem  %d\n",ret);
      exit(1);
    }

    ret = clEnqueueNDRangeKernel(cmdQ, kernel11[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel13 enq problem  %d\n",ret);
      exit(1);
    }
  }
  else {
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
    printf("o        Temperature condition is not choosen correctly         o\n");
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
    exit(1);
  }

}
