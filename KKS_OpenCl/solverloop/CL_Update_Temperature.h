void CL_Update_Temperature() {
  
  if ( !ISOTHERMAL ) {
    tstep[0] = t;
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
      printf("kernel11 enq problem  %d\n",ret);
      exit(1);
    }

    ret = clEnqueueNDRangeKernel(cmdQ, kernel11[1], work_dim, NULL, globaldim, NULL, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("kernel11 enq problem  %d\n",ret);
      exit(1);
    }
  }
}