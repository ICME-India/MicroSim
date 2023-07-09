void CL_buffer_allocation() {

  /* Create buffer object */
//  if (ELASTICITY) {
    d_iter_gridinfom = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxnynz*sizeof(struct iter_variables), iter_gridinfom, &ret);
    if ( !d_iter_gridinfom ) {
      printf("falied to allocate device memory d_iter_gridinfom in rank = %d, Error code = %d \n", rank, ret);
      exit(1);
    }
//  }
//  else {
//    d_iter_gridinfom = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1*sizeof(struct iter_variables), iter_gridinfom, &ret);
//    if ( !d_iter_gridinfom ) {
//      printf("falied to allocate device memory d_iter_gridinfom in rank = %d, Error code = %d \n", rank, ret);
//      exit(1);
//    }
//  }
    d_eigen_strain_phase =  clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, NUMPHASES*sizeof(struct symmetric_tensor), eigen_strain_phase, &ret);
    if ( !d_eigen_strain_phase ) {
      printf("falied to allocate device memory d_eigen_strain_phase in rank = %d, Error code = %d \n", rank, ret);
      exit(1);
    }
    d_stiffness_phase = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, NUMPHASES*sizeof(struct Stiffness_cubic), stiffness_phase, &ret);
    if ( !d_stiffness_phase ) {
      printf("falied to allocate device memory d_stiffness_phase in rank = %d, Error code = %d \n", rank, ret);
      exit(1);
    }
    d_stiffness_phase_n = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, NUMPHASES*sizeof(struct Stiffness_cubic), stiffness_phase_n, &ret);
    if ( !d_stiffness_phase_n ) {
      printf("falied to allocate device memory d_stiffness_phase_n in rank = %d, Error code = %d \n", rank, ret);
      exit(1);
    }

  d_gridinfomO = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxnynz*sizeof(struct fields), gridinfomO, &ret);
  if ( !d_gridinfomO ) {
    printf("falied to allocate device memory d_gridinfomO in rank = %d, Error code = %d \n", rank, ret);
    exit(1);
  }
  d_gridinfomN = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxnynz*sizeof(struct fields), gridinfomN, &ret);
  if ( !d_gridinfomN ) {
    printf("falied to allocate device memory d_gridinfomN %d\n", ret);
    exit(1);
  }
  d_cscl = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nxnynz*sizeof(struct csle), cscl, &ret);
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
  //d_temp = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nx*sizeof(double), temp, &ret);//Changed to nx, According to MESH_X
  //if ( !d_temp ) {
  //  printf("falied to allocate device memory d_temp %d\n", ret);
  //  exit(1);
  //}
  d_tstep = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(long), tstep, &ret);
  if ( !d_tstep ) {
    printf("falied to allocate device memory d_tstep %d\n", ret);
    exit(1);
  }
  d_propf3 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(struct propmatf3), &propf3, &ret);
  if ( !d_propf3 ) {
    printf("falied to allocate device memory d_propf3 %d\n", ret);
    exit(1);
  }
  d_propf4 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(struct propmatf4), &propf4, &ret);
  if ( !d_propf4 ) {
    printf("falied to allocate device memory d_propf4 %d\n", ret);
    exit(1);
  }
  d_propf4spline = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nx*sizeof(struct propmatf4spline), propf4spline, &ret);
  if ( !d_propf4spline ) {
    printf("falied to allocate device memory d_propf4spline %d\n", ret);
    exit(1);
  }
  d_propf4spline1 = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nx*sizeof(struct propmatf4spline), propf4spline1, &ret);
  if ( !d_propf4spline1 ) {
    printf("falied to allocate device memory d_propf4spline1 %d\n", ret);
    exit(1);
  }
  printf("Device buffers created. \n");

  /* Copy input data to memory buffer */
  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomO, CL_TRUE, 0, nxnynz*sizeof(struct fields), gridinfomO, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomO in rank = %d. Error code = %d\n", rank, ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_gridinfomN, CL_TRUE, 0, nxnynz*sizeof(struct fields), gridinfomN, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_gridinfomN %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_cscl, CL_TRUE, 0, nxnynz*sizeof(struct csle), cscl, 0, NULL, NULL);
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
  //ret  = clEnqueueWriteBuffer(cmdQ, d_temp, CL_TRUE, 0, nx*sizeof(double), temp, 0, NULL, NULL);//Changed to nx, According to MESH_X
  //if (ret!=CL_SUCCESS) {
  //  printf("enq buffer write error d_temp %d\n", ret);
  //  exit(1);
  //}
  ret  = clEnqueueWriteBuffer(cmdQ, d_tstep, CL_TRUE, 0, sizeof(long), tstep, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_tstep %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_propf3, CL_TRUE, 0, sizeof(struct propmatf3), &propf3, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_propf3 %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_propf4, CL_TRUE, 0, sizeof(struct propmatf4), &propf4, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_propf4 %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_propf4spline, CL_TRUE, 0, nx*sizeof(struct propmatf4spline), propf4spline, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_propf4spline %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_propf4spline1, CL_TRUE, 0, nx*sizeof(struct propmatf4spline), propf4spline1, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_propf4spline1 %d\n", ret);
    exit(1);
  }

//   if (ELASTICITY) {
    ret  = clEnqueueWriteBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, 0, nxnynz*sizeof(struct iter_variables), iter_gridinfom, 0, NULL, NULL);
    if (ret!=CL_SUCCESS) {
      printf("enq buffer write error d_iter_gridinfom in rank = %d. Error code = %d\n", rank, ret);
      exit(1);
    }
//   }
//   else {
//     ret  = clEnqueueWriteBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, 0, 1*sizeof(struct iter_variables), iter_gridinfom, 0, NULL, NULL);
//     if (ret!=CL_SUCCESS) {
//       printf("enq buffer write error d_iter_gridinfom in rank = %d. Error code = %d\n", rank, ret);
//       exit(1);
//     }
//   }

  ret  = clEnqueueWriteBuffer(cmdQ, d_eigen_strain_phase, CL_TRUE, 0, NUMPHASES*sizeof(struct symmetric_tensor), eigen_strain_phase, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_eigen_strain_phase in rank = %d. Error code = %d\n", rank, ret);
    exit(1);
  }

  ret  = clEnqueueWriteBuffer(cmdQ, d_stiffness_phase, CL_TRUE, 0, NUMPHASES*sizeof(struct Stiffness_cubic), stiffness_phase, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_stiffness_phase in rank = %d. Error code = %d\n", rank, ret);
    exit(1);
  }

  ret  = clEnqueueWriteBuffer(cmdQ, d_stiffness_phase_n, CL_TRUE, 0, NUMPHASES*sizeof(struct Stiffness_cubic), stiffness_phase_n, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_stiffness_phase_n in rank = %d. Error code = %d\n", rank, ret);
    exit(1);
  }

  printf("Write to device completed. \n");

}
