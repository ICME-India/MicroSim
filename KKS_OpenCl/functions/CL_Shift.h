void CL_Shift() {

  int x, y, z;
  long index;
  long INTERFACE_POS;
  
  CL_DeviceToHost(); 

  ret = clEnqueueReadBuffer(cmdQ, d_cscl, CL_TRUE, 0, nxny*sizeof(struct csle), cscl, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read cscl in mvframe \n%d\n", ret);
    exit(1);
  }

  for(x=start[X]-2; x<=end[X]+2; x++) {
    if (x > 1) {
      if (SHIFT) {
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > MAX_INTERFACE_POS) {
          MAX_INTERFACE_POS = INTERFACE_POS;
        }
      }
    }
  }

  if(SHIFT) {
    if(MAX_INTERFACE_POS > shiftj) {
      shift_ON = 1;
    }
    if (shift_ON) {
      apply_shiftY(gridinfo, MAX_INTERFACE_POS); 
      shift_OFFSET += (MAX_INTERFACE_POS - shiftj);
      shift_ON=0;
    }
  }

  for (x=0; x<rows_x; x++) {
    for(z=0; z<rows_z; z++) {
      for (y=0; y<rows_y; y++) {

        index = x*layer_size + z*rows_y + y;

        gridNew[index].c1 = gridinfo[index].composition[0];

        gridNew[index].mu[0] = gridinfo[index].compi[0];

        gridNew[index].phi = gridinfo[index].phia[0];

        temp[y] = gridinfo[index].temperature;
      }
    }
  }

  pfmdat.shift_OFFSET = shift_OFFSET;

  ret  = clEnqueueWriteBuffer(cmdQ, d_pfmdat, CL_TRUE, 0, sizeof(struct pfmval), &pfmdat, 0, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_pfmdat %d\n", ret);
    exit(1);
  }
  ret  = clEnqueueWriteBuffer(cmdQ, d_temp, CL_TRUE, 0, nx*sizeof(double), temp, 0, NULL, NULL);//Changed to nx, According to MESH_X
  if (ret!=CL_SUCCESS) {
    printf("enq buffer write error d_temp %d\n", ret);
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

}


