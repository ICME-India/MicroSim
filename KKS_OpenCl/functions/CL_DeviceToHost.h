void CL_DeviceToHost() {

  int x, y, z;
  long index;

  ret = clEnqueueReadBuffer(cmdQ, d_gridNew, CL_TRUE, 0, nxny*sizeof(struct grid), gridNew, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridNew \n%d\n", ret);
    exit(1);
  }
  ret = clEnqueueReadBuffer(cmdQ, d_temp, CL_TRUE, 0, ny*sizeof(double), temp, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_temp \n%d\n", ret);
    exit(1);
  }

  /*printf("Checking for NAN\n");
  for(i=istart; i < iend; i++) {
    for(j=0; j < nx; j++) {
      index = (i*nx + j);
      if ( isnan(gridNew[index].phi) ) {
        printf("1NAN while writing at timestep=%d, i=%d, j=%d, rank=%d\n", t, i, j, rank);
        exit(1);
      }
    }
  }*/

  for (x=0; x<rows_x; x++) {
    for(z=0; z<rows_z; z++) {
      for (y=0; y<rows_y; y++) {

        index = x*layer_size + z*rows_y + y;

        gridinfo[index].compi[0] = gridNew[index].c1;

        gridinfo[index].phia[0] = gridNew[index].phi;

      }
    }
  }


}

void CL_Global_Max_Min() { 

  int x, y, z;
  long index;

  CL_DeviceToHost();

  ret = clEnqueueReadBuffer(cmdQ, d_gridOld, CL_TRUE, 0, nxny*sizeof(struct grid), gridOld, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridOld \n%d\n", ret);
    exit(1);
  }

  global_max_min.rel_change_phi[0] = 0.0;
  global_max_min.rel_change_com[0] = 0.0;

  for (x=0; x<rows_x; x++) {
    for(z=0; z<rows_z; z++) {
      for (y=0; y<rows_y; y++) {

        index = x*layer_size + z*rows_y + y;

        gridinfoO[index].compi[0] = gridOld[index].c1;

        gridinfoO[index].phia[0] = gridOld[index].phi;

        global_max_min.rel_change_phi[0] += (gridinfo[index].phia[0]-gridinfoO[index].phia[0])*(gridinfo[index].phia[0]-gridinfoO[index].phia[0]);

        if (gridinfo[index].phia[0] > global_max_min.phi_max[0]) { 
          global_max_min.phi_max[0] = gridinfo[index].phia[0];
        }
        if (gridinfo[index].phia[0] < global_max_min.phi_min[0]) { 
          global_max_min.phi_min[0] = gridinfo[index].phia[0];
        }

        global_max_min.rel_change_com[0] += (gridinfo[index].compi[0]-gridinfoO[index].compi[0])*(gridinfo[index].compi[0]-gridinfoO[index].compi[0]);

        if (gridinfo[index].compi[0] > global_max_min.com_max[0]) { 
          global_max_min.com_max[0] = gridinfo[index].compi[0];
        }
        if (gridinfo[index].compi[0] < global_max_min.com_min[0]) { 
          global_max_min.com_min[0] = gridinfo[index].compi[0];
        }

      }
    }
  }

}

/*
void savetimestep(struct grid *gridNew, struct pfmval *pfmdat, struct pfmpar *pfmvar, double *temp, int t) {

  int i, j, k;
  int index;
  int nx;
  int ny;
  int nz;
  int rank;
  
  FILE *fp0;
  FILE *fp1;
  FILE *fp2;
  FILE *fp3;
  FILE *fp4;
  
  char fname0[100];
  char fname1[100];
  char fname2[100];
  char fname3[100];
  char fname4[100];
  
  char command0[150];
  char command1[150];
  char command2[150];
  char command3[150];
  char command4[150];

  int istart;
  int iend;
  int totny;
  int nxtotny;
  int indexa;
  int ia;
  
  int ig;
  double tmp0;
  
  ig=0;
  
  rank = pfmdat->myrank;
  
  nx = pfmdat->jNx;
  ny = pfmdat->iNy;
  
  istart = 0;
  iend = ny;
  
  sprintf(fname2, "c_1_%011d_%04d.dat", t, rank);

  ////Writing phi
  sprintf(fname1, "phi_%011d_%04d.dat", t, rank);
  fp1 = fopen(fname1,"w");
  for(i = istart+1; i < iend-1; i++) {  //Boundary removed
    for(j=1; j < nx-1; j++) {           //Boundary removed
      index = (i*nx + j);
      fprintf(fp1,"%le\t", gridNew[index].phi);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  sprintf(command1, "gzip %s", fname1);
  system(command1);

  fp2 = fopen(fname2,"w");
  for(i = istart+1; i < iend-1; i++) {  //Boundary removed
    for(j=1; j < nx-1; j++) {           //Boundary removed
      index = (i*nx + j);
      fprintf(fp2,"%le\t", gridNew[index].c1);
    }
    fprintf(fp2,"\n");
  }
  fclose(fp2);

  sprintf(command2, "gzip %s", fname2);
  system(command2);

  if ( !ISOTHERMAL ) {
    sprintf(fname4, "tem_%011d_%04d.dat", t, rank);
    fp4 = fopen(fname4,"w");
    for(i=1; i < ny-1; i++) {           //Boundary removed
      fprintf(fp4,"%le\n", temp[i]);  
    }
    fclose(fp4);
    
    sprintf(command4, "gzip %s", fname4);
    system(command4);

  }
  
}
*/
