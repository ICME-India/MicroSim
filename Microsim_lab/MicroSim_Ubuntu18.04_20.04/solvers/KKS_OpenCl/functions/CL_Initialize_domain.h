void CL_Initialize_domain() {

  //int restart=0;
  int x, y, z;
  long index;

  for (x=0; x<rows_x; x++) {
    for(z=0; z<rows_z; z++) {
      for (y=0; y<rows_y; y++) {

        index = x*layer_size + z*rows_y + y;

        cscl[index].c1l = ( pfmdat.c1l_1stguess );
        cscl[index].c1s = ( pfmdat.c1s_1stguess );

        gridOld[index].c1  = gridinfo[index].compi[0];
        gridNew[index].c1  = gridinfo[index].compi[0];

        gridOld[index].phi = gridinfo[index].phia[0];
        gridNew[index].phi = gridinfo[index].phia[0];

      }
    }
  }

  int i0;
  int rank=0;
  int itimestep=STARTTIME;  

  if ( !ISOTHERMAL ) {
    for ( i = 0; i < nx; i++ ) { //Changed to nx, According to MESH_X
      //i0 = i + (rank*(ny-2))-1;
      //temp[i] = pfmdat.Toffset + pfmdat.TG*((i0-pfmdat.PosOffset)*pfmvar.dx - pfmdat.Vp*itimestep*pfmvar.dt);

      i0 = i + (rank*(nx-2));//Changed to nx, According to MESH_X
      temp[i] = pfmdat.Toffset + pfmdat.TGRADIENT*((i0-pfmdat.TPosOffset+pfmdat.shift_OFFSET)*pfmvar.deltax-(pfmdat.velocity*itimestep*pfmvar.deltat));
      //printf("CL:%le\n", temp[i]);

    }
    if ( rank == 0 ) {
      i = 0;
      i0 = 0;
      temp[i] = temp[i+1];
    }
    else if ( rank == pfmdat.nproc-1 ) { 
      i = (ny-1);
      i0 = (rank+1)*(ny-2)-1;
      temp[i] = temp[i-1];
    }
  }
  else if ( ISOTHERMAL ) {
    for ( i = 0; i < ny; i++ ) {
      temp[i] = pfmdat.T0;
    }
  }
  else {
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
    printf("o      Temperature conditions are not chosen correctly          o\n");
    printf("ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
    exit(1);
  }

  
  //x=1;
  //z=1;
  //for (y=0; y<rows_y; y++) {

    //index = x*layer_size + z*rows_y + y;

    //temp[y] = gridinfo[index].temperature;
    //printf("MS:%le\n", temp[y]);

  //}


}