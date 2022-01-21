void CL_buffer_allocation() {
  
  /* Create buffer object */
  d_gridOld = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxny*sizeof(struct grid), gridOld, &ret);
  if ( !d_gridOld ) {
    printf("falied to allocate device memory d_gridOld in rank = %d, Error code = %d \n", rank, ret);
    exit(1);
  }
  d_gridNew = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxny*sizeof(struct grid), gridNew, &ret);
  if ( !d_gridNew ) {
    printf("falied to allocate device memory d_gridNew %d\n", ret);
    exit(1);
  }
  d_cscl = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxny*sizeof(struct csle), cscl, &ret);
  if ( !d_cscl ) {
    printf("falied to allocate device memory d_csl %d\n", ret);
    exit(1);
  }
  d_pfmdat = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(struct pfmval), &pfmdat, &ret);
  if ( !d_pfmdat ) {
    printf("falied to allocate device memory d_pfmdat %d\n", ret);
    exit(1);
  }
  d_pfmvar = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(struct pfmpar), &pfmvar, &ret);
  if ( !d_pfmvar ) {
    printf("falied to allocate device memory d_pfmvar %d\n", ret);
    exit(1);
  }
  d_temp = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nx*sizeof(double), temp, &ret);//Changed to nx, According to MESH_X
  if ( !d_temp ) {
    printf("falied to allocate device memory d_temp %d\n", ret);
    exit(1);
  }
  d_tstep = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(long), tstep, &ret);
  if ( !d_tstep ) {
    printf("falied to allocate device memory d_tstep %d\n", ret);
    exit(1);
  }
  printf("Device buffers created. \n");

  /* Copy input data to memory buffer */
  ret  = clEnqueueWriteBuffer(cmdQ, d_gridOld, CL_TRUE, 0, nxny*sizeof(struct grid), gridOld, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridOld in rank = %d. Error code = %d\n", rank, ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_gridNew, CL_TRUE, 0, nxny*sizeof(struct grid), gridNew, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridNew %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_cscl, CL_TRUE, 0, nxny*sizeof(struct csle), cscl, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_cscl %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_pfmdat, CL_TRUE, 0, sizeof(struct pfmval), &pfmdat, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_pfmdat %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_pfmvar, CL_TRUE, 0, sizeof(struct pfmpar), &pfmvar, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_pfmvar %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_temp, CL_TRUE, 0, nx*sizeof(double), temp, 0, NULL, NULL);//Changed to nx, According to MESH_X
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_temp %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_tstep, CL_TRUE, 0, sizeof(long), tstep, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_tstep %d\n", ret);
    exit(1);
  }
  printf("Write to device completed. \n");

}