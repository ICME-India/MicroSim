void CL_DeviceToHost() {

  int x, y, z;
  long index;
  int ip, is;
  double tmpmx, tmpmn;

  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomN, CL_TRUE, 0, nxnynz*sizeof(struct fields), gridinfomN, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomN \n%d\n", ret);
    exit(1);
  }

  ret = clEnqueueReadBuffer(cmdQ, d_iter_gridinfom, CL_TRUE, 0, nxnynz*sizeof(struct iter_variables), iter_gridinfom, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_iter_gridinfom \n%d\n", ret);
    exit(1);
  }
  
}

void CL_Global_Max_Min() { 

  int x, y, z, ip, is;
  long index;

  CL_DeviceToHost();

  ret = clEnqueueReadBuffer(cmdQ, d_gridinfomO, CL_TRUE, 0, nxnynz*sizeof(struct fields), gridinfomO, 0, NULL, NULL);
  if (ret != CL_SUCCESS) {
    printf("Error: Failed to read d_gridinfomO \n%d\n", ret);
    exit(1);
  }
  
  for ( ip = 0; ip < npha; ip++ ) { 
    global_max_min.rel_change_phi[ip] = 0.0;
  }
  for ( is = 0; is < nsol; is++ ) { 
    global_max_min.rel_change_com[is] = 0.0;
  }

  for (x=0; x<nx; x++) {
    for(z=0; z<nz; z++) {
      for (y=0; y<ny; y++) {

        index = x*layer_size + z*ny + y;
        
        for ( ip = 0; ip < npha; ip++ ) { 
        global_max_min.rel_change_phi[ip] += (gridinfomN[index].phia[ip]-gridinfomO[index].phia[ip])*(gridinfomN[index].phia[ip]-gridinfomO[index].phia[ip]);
        }
        
        //for ( ip = 0; ip < npha; ip++ ) { 
        //
        //if (gridinfo[index].phia[ip] > global_max_min.phi_max[ip]) { 
        //  global_max_min.phi_max[ip] = gridinfo[index].phia[ip];
        //}
        //if (gridinfo[index].phia[ip] < global_max_min.phi_min[ip]) { 
        //  global_max_min.phi_min[ip] = gridinfo[index].phia[ip];
        //}
        //
        //}
        
        for ( is = 0; is < nsol; is++ ) { 

        global_max_min.rel_change_com[is] += (gridinfomN[index].compi[is]-gridinfomO[index].compi[is])*(gridinfomN[index].compi[is]-gridinfomO[index].compi[is]);
        
        }
        
        for ( is = 0; is < nsol; is++ ) { 

        global_max_min.rel_change_composition[is] += (gridinfomN[index].composition[is]-gridinfomO[index].composition[is])*(gridinfomN[index].composition[is]-gridinfomO[index].composition[is]);
        
        }
        
        //for ( is = 0; is < nsol; is++ ) { 
        //
        //if (gridinfo[index].compi[is] > global_max_min.com_max[is]) { 
        //  global_max_min.com_max[is] = gridinfo[index].compi[is];
        //}
        //if (gridinfo[index].compi[is] < global_max_min.com_min[is]) { 
        //  global_max_min.com_min[is] = gridinfo[index].compi[is];
        //}
        //
        //}

      }
    }
  }
  
  
  for ( is = 0; is < nsol; is++ ) { 
    
    global_max_min.com_max[is] = gridinfomN[0].compi[is];
    global_max_min.com_min[is] = gridinfomN[0].compi[is];
    
    global_max_min.composition_max[is] = gridinfomN[0].composition[is];
    global_max_min.composition_min[is] = gridinfomN[0].composition[is];

    for (x=0; x<nx; x++) {
      for(z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
          
          index = x*layer_size + z*ny + y;
          
          if (gridinfomN[index].compi[is] > global_max_min.com_max[is]) {
            global_max_min.com_max[is] = gridinfomN[index].compi[is];
          }
          if (gridinfomN[index].compi[is] < global_max_min.com_min[is]) {
            global_max_min.com_min[is] = gridinfomN[index].compi[is];
          }
          
          if (gridinfomN[index].composition[is] > global_max_min.composition_max[is]) {
            global_max_min.composition_max[is] = gridinfomN[index].composition[is];
          }
          if (gridinfomN[index].composition[is] < global_max_min.composition_min[is]) {
            global_max_min.composition_min[is] = gridinfomN[index].composition[is];
          }
        
        }
      }
    }
  }
  
  for ( ip = 0; ip < npha; ip++ ) { 
    
    global_max_min.phi_max[ip] = gridinfomN[0].phia[ip];
    global_max_min.phi_min[ip] = gridinfomN[0].phia[ip];

    for (x=0; x<nx; x++) {
      for(z=0; z<nz; z++) {
        for (y=0; y<ny; y++) {
          
          index = x*layer_size + z*ny + y;
          
          if (gridinfomN[index].phia[ip] > global_max_min.phi_max[ip]) {
            global_max_min.phi_max[ip] = gridinfomN[index].phia[ip];
          }
          if (gridinfomN[index].phia[ip] < global_max_min.phi_min[ip]) {
            global_max_min.phi_min[ip] = gridinfomN[index].phia[ip];
          }
        
        }
      }
    }
  }

}


void savetimestepcsle(struct csle *cscl, struct pfmval *pfmdat, int t) {

  int i, j, k;
  int index;
  int nx;
  int ny;
  int nz;
  int rank;
  
  FILE *fp0;
  FILE *fp1;
  FILE *fp2;
  
  char fname0[100];
  char fname1[100];
  char fname2[100];
  
  char command0[150];
  char command1[150];
  char command2[150];
  
  int istart;
  int iend;
  int totny;
  int nxtotny;
  int indexa;
  int ia;
  
  int ig;
  float tmp0;
  int is, js, ip;
  
  ig=0;
  
  rank = pfmdat->myrank;
  
  nx = pfmdat->Nx;
  ny = pfmdat->Ny;
  
  istart = 0;
  iend = ny;

  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is = 0; is < (NUMCOMPONENTS-1); is++ ) { 
      sprintf(fname2, "cl_p%01d_s%01d_%011d_%04d.csl", ip+1, is+1, t, rank);
      fp2 = fopen(fname2,"w");
      for(i = istart+1; i < iend-1; i++) {
        for(j=1; j < nx-1; j++) {
          index = (i*nx + j);
          fprintf(fp2,"%f\t", cscl[index].comie[ip][is]);
        }
        fprintf(fp2,"\n");
      }
      fclose(fp2);
      sprintf(command2, "gzip -f %s", fname2);
      //system(command2);
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

