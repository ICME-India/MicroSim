void print_input_parameters(char *argv[]){
  int i,j;
  char tempbuff[1000];
  
  char tmpstr1[1000];
  char tmpstr2[1000];
  long length;
  
  FILE *fr;
//   char **tmp;
//   
//   bool decision;
//   
//   char *str1, *str2, *token, *subtoken;
//   char *saveptr1, *saveptr2;
//   
//   long k, j;
//   long index;
//   long length;
//   long phase;
  
  char outfile[2000];
  
  strcpy(tmpstr2, argv[1]);
//   sscanf(argv[1], "%100s.%100[^;]", tmpstr1, tmpstr2);
  strcpy(tmpstr1,strtok(tmpstr2, "."));
  
  sprintf(outfile, "%s.out", tmpstr1);
  
  fr = fopen(outfile, "w");
  
  char key[1000];
  
  strcpy(key, "MESH_X");
  
  PRINT_LONG(key, MESH_X, fr);
  
  strcpy(key, "MESH_Y");
  
  PRINT_LONG(key, MESH_Y, fr);
  
  strcpy(key, "MESH_Z");
  
  PRINT_LONG(key, MESH_Z, fr);
  
  strcpy(key,"DELTA_X");
  
  PRINT_DOUBLE(key, deltax, fr);
  
  strcpy(key,"DELTA_Y");
  
  PRINT_DOUBLE(key, deltay, fr);
  
  strcpy(key,"DELTA_Z");
  
  PRINT_DOUBLE(key, deltaz, fr);
  
  strcpy(key, "NUMPHASES");
  
  PRINT_INT(key, NUMPHASES, fr);

  strcpy(key, "NUMCOMPONENTS");
  
  PRINT_INT(key, NUMCOMPONENTS, fr);
  
  strcpy(key, "DIMENSION");
  
  PRINT_INT(key, DIMENSION, fr);
  
  strcpy(key, "NTIMESTEPS");
  
  PRINT_LONG(key, ntimesteps, fr);
  
  strcpy(key, "SAVET");
  
  PRINT_LONG(key, saveT, fr);
  
  strcpy(key, "NSMOOTH");
  
  PRINT_LONG(key, nsmooth, fr);

  strcpy(key,"R");
  
  PRINT_DOUBLE(key, R, fr);
  
  strcpy(key,"V");
  
  PRINT_DOUBLE(key, V, fr);
  
  strcpy(key,"TEMPERATURE_SCALE");
  
  PRINT_DOUBLE(key, TEMPERATURE_SCALE, fr);
  
  strcpy(key,"TEMPERATURE");
  
  PRINT_DOUBLE(key, T, fr);
  
  
  strcpy(key,"Equilibrium_temperature");
  
  PRINT_DOUBLE(key, Teq, fr);
  
  strcpy(key,"Filling_temperature");
  
  PRINT_DOUBLE(key, Tfill, fr);
  
  strcpy(key, "GAMMA");
  
  PRINT_MATRIX(key, Gamma, NUMPHASES, NUMPHASES, fr);
  
  strcpy(key, "dab");
  
  PRINT_MATRIX(key, dab, NUMPHASES, NUMPHASES, fr);
  
  strcpy(key, "Tau");
  
  PRINT_MATRIX(key, tau_ab, NUMPHASES, NUMPHASES, fr);
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "DIFFUSIVITY[%s]",Phases[i]);
    PRINT_MATRIX(key, Diffusivity[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "Slopes[Solidus,%s]",Phases[i]);
    PRINT_VECTOR(key, slopes[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "Slopes[Liquidus,%s]",Phases[i]);
    PRINT_VECTOR(key, slopes[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "A[%s]",Phases[i]);
    PRINT_MATRIX(key, A[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    for (j=0; j<NUMPHASES; j++) {
      if(i!=j) {
        sprintf(key, "Rotation_matrix[%s][%s]",Phases[i],Phases[j]);
        PRINT_MATRIX(key, Rotation_matrix[i][j], 3, 3, fr);
      }
    }
  }
  
  for (i=0; i<NUMPHASES; i++) {
    for (j=0; j<NUMPHASES; j++) {
      if(i!=j) {
        sprintf(key, "Inv_Rotation_matrix[%s][%s]",Phases[i],Phases[j]);
        PRINT_MATRIX(key, Inv_Rotation_matrix[i][j], 3, 3, fr);
      }
    }
  }
  
  strcpy(key, "Function_F");
  
  PRINT_INT(key, FUNCTION_F, fr);
  
  strcpy(key, "Function_W");
  
  PRINT_INT(key, FUNCTION_W, fr);
  
  strcpy(key, "Function_anisotropy");
  
  PRINT_INT(key, FUNCTION_ANISOTROPY, fr);
  
  strcpy(key, "Anisotropy_type");
  
  PRINT_INT(key, FOLD, fr);
  
  strcpy(key, "COMPONENTS");
  
  PRINT_STRING_ARRAY(key, Components, NUMCOMPONENTS, fr);
  
  strcpy(key, "PHASES");
  
  PRINT_STRING_ARRAY(key, Phases, NUMPHASES, fr);
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "ceq[Solidus,%s]",Phases[i]);
    PRINT_VECTOR(key, ceq[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "ceq[Liquidus,%s]",Phases[i]);
    PRINT_VECTOR(key, ceq[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "cfill[Liquidus,%s]",Phases[i]);
    PRINT_VECTOR(key, cfill[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "EIGEN_STRAIN[%s]",Phases[i]);
    PRINT_SYMMETRIC_TENSOR(key, eigen_strain_phase[i], fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "Stiffness_cubic[%s]",Phases[i]);
    PRINT_VOIGT_CUBIC(key, stiffness_phase[i], fr);
  }
  
  
  
  max_length = strlen(Phases[0]); 
  for (a=1; a<NUMPHASES; a++) {
    length = strlen(Phases[a]);
    if (length > max_length) {
      max_length = length;
    }
  }
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    length = strlen(Components[k]);
    if (length > max_length) {
      max_length = length;
    }
  }
  
//   PRINT_BOUNDARY_CONDITIONS(fr);
  
  fclose(fr);
}
