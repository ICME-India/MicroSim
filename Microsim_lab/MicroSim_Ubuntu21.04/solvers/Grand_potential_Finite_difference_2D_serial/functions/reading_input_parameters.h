void reading_input_parameters(char *argv[]) {
  FILE * fr = fopen(argv[1], "rt");
  
  if(fr == NULL) {
    printf("file %s not found", argv[1]);
  }
  int i;
  char tempbuff[1000];
  
  char tmpstr1[100];
  char tmpstr2[100];
  char **tmp;
  
  bool decision;
  
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  
  long k, j;
  long index;
  long length;
  long phase;
  
  while(fgets(tempbuff,1000,fr)) {
    sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);
//     printf("%s\n",  tmpstr1);
//     printf("%s\n",  tmpstr2);
    
    if(tmpstr1[0] != '#') {
      if (strcmp(tmpstr1,"DIMENSION")==0) {
        DIMENSION = atoi(tmpstr2);
        start   = (long*)malloc(3*sizeof(*start));
        end     = (long*)malloc(3*sizeof(*end));
        averow  = (long*)malloc(3*sizeof(*averow));
        rows    = (long*)malloc(3*sizeof(*rows));
        offset  = (long*)malloc(3*sizeof(*offset));
        extra   = (long*)malloc(3*sizeof(*extra));
      }
      else if (strcmp(tmpstr1,"MESH_X")==0) {
         MESH_X = atol(tmpstr2);
      } 
      else if (strcmp(tmpstr1,"MESH_Y")==0) {
         MESH_Y = atol(tmpstr2);
      }
      else if(strcmp(tmpstr1,"MESH_Z")==0) {
        MESH_Z = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"DELTA_X")==0) {
        deltax = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"DELTA_Y")==0) {
        deltay = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"DELTA_Z")==0) {
        deltaz = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"DELTA_t")==0) {
        deltat = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"NUMPHASES")==0) {
        NUMPHASES = atoi(tmpstr2);
        Phases = (char**)malloc(sizeof(char*)*NUMPHASES);
        for (i = 0; i < NUMPHASES; ++i) {
          Phases[i] = (char*)malloc(sizeof(char)*51);
        }
        filling_type_phase = (struct filling_type* )malloc((NUMPHASES)*sizeof(*filling_type_phase));
      }
      else if (strcmp(tmpstr1,"NUMCOMPONENTS")==0) {
        NUMCOMPONENTS = atoi(tmpstr2);
        Components = (char**)malloc(sizeof(char*)*NUMCOMPONENTS);
        for (i = 0; i < NUMCOMPONENTS; ++i) {
          Components[i] = (char*)malloc(sizeof(char)*51);
        }
        if (((NUMCOMPONENTS-1)>0) && (NUMPHASES>0)) {
          Diffusivity         = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
          ceq                 = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
          cfill               = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
          ceq_coeffs          = Malloc3M(NUMPHASES, NUMCOMPONENTS-1,               4);
          slopes              = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
          dcbdT               = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
          A                   = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
          
          DELTA_T             = MallocM(NUMPHASES,  NUMPHASES);
          DELTA_C             = MallocM(NUMPHASES,  NUMCOMPONENTS-1);
          dcbdT_phase         = MallocM(NUMPHASES,  NUMCOMPONENTS-1);
          B                   = MallocM(NUMPHASES,  NUMCOMPONENTS-1);
          Beq                 = MallocM(NUMPHASES,  NUMCOMPONENTS-1);
          dBbdT               = MallocM(NUMPHASES,  NUMCOMPONENTS-1);
          C                   = (double *)malloc(NUMPHASES*sizeof(double));
          
          cmu                 = Malloc3M(NUMPHASES, NUMCOMPONENTS-1,   NUMCOMPONENTS-1);
	  muc                 = Malloc3M(NUMPHASES, NUMCOMPONENTS-1,   NUMCOMPONENTS-1);
	  Rotation_matrix     = Malloc4M(NUMPHASES,       NUMPHASES,   3,             3);
	  Inv_Rotation_matrix = Malloc4M(NUMPHASES,       NUMPHASES,   3,             3);
	  Rotated_qab         = MallocV(3);
          
          
          eigen_strain    = (struct symmetric_tensor *)malloc(NUMPHASES*sizeof(*eigen_strain));
          Stiffness_c     = (struct Stiffness_cubic *)malloc(NUMPHASES*sizeof(*Stiffness_c));
          Stiffness_t     = (struct Stiffness_tetragonal *)malloc(NUMPHASES*sizeof(*Stiffness_t));
          for (i=0; i < 6; i++) {
            boundary[i] = (struct bc_scalars*)malloc(3*sizeof(*boundary[i])); //3=Number of scalar fields
          }
        }
      }
      else if (strcmp(tmpstr1,"epsilon")==0) {
        epsilon = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"tau")==0) {
        tau = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"NTIMESTEPS")==0) {
        ntimesteps = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"SAVET")==0) {
        saveT = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"NSMOOTH")==0) {
        nsmooth = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"STARTTIME")==0) {
        STARTTIME = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"RESTART")==0) {
        RESTART = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"R")==0) {
        R = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"V")==0) {
        V = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"WRITEFORMAT")==0) {
        if(strcmp(tmpstr2,"ASCII") == 0) {
          ASCII = 1;
        }
        else if (strcmp(tmpstr2,"BINARY") == 0) {
          ASCII = 0;
        }
      }
      else if (strcmp(tmpstr1,"TRACK_PROGRESS")==0) {
        time_output = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"TEMPERATURE_SCALE")==0) {
        TEMPERATURE_SCALE = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"Equilibrium_temperature")==0) {
        Teq = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"Filling_temperature")==0) {
        Tfill = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"T")==0) {
        T = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"tilt_angle")==0) {
        tilt_angle = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"Rtheta")==0) {
        Rtheta = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"COMPONENTS")==0) {
        populate_string_array(Components, tmpstr2, NUMCOMPONENTS);
      }
      else if (strcmp(tmpstr1,"PHASES")==0) {
        populate_string_array(Phases, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "Function_anisotropy") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        FUNCTION_ANISOTROPY = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Anisotropy_type") == 0) && (FUNCTION_ANISOTROPY !=0)) {
         FOLD = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Rotation_matrix") == 0) && (FUNCTION_ANISOTROPY !=0)) {
         populate_rotation_matrix(Rotation_matrix, Inv_Rotation_matrix, tmpstr2);
      }
      else if (strcmp(tmpstr1, "OBSTACLE") == 0) {
        OBSTACLE = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Function_W") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)){
        FUNCTION_W = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "ISOTHERMAL") == 0) {
        ISOTHERMAL = atoi(tmpstr2);
        if (!ISOTHERMAL) {
          TEMPGRADY = 1;
        }
      }
      else if ((strcmp(tmpstr1, "Tempgrady") == 0) && (TEMPGRADY)) {
        tmp = (char**)malloc(sizeof(char*)*5);
        for (i = 0; i < 5; ++i) {
          tmp[i] = (char*)malloc(sizeof(char)*10);
        }
        for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL) {
          token = strtok_r(str1, "{,}", &saveptr1);
          if (token == NULL)
              break;
          strcpy(tmp[i],token);
        }
        temperature_gradientY.base_temp          = atof(tmp[0]);
        temperature_gradientY.DeltaT             = atof(tmp[1]);
        temperature_gradientY.Distance           = atof(tmp[2]);
        temperature_gradientY.gradient_OFFSET    = atof(tmp[3]);
        temperature_gradientY.velocity           = atof(tmp[4]);
        for (i = 0; i < 5; ++i) {
          free(tmp[i]);
        }
        free(tmp);
      }
      else if (strcmp(tmpstr1, "Shift") == 0) {
        SHIFT = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "Shiftj") == 0) {
        if (SHIFT) {
         shiftj = atol(tmpstr2);
        }
      }
      else if (strcmp(tmpstr1, "BINARY") == 0) {
        BINARY = atoi(tmpstr2);
        TERNARY = 0;
        DILUTE = 0;
      }
      else if (strcmp(tmpstr1, "TERNARY") == 0) {
        TERNARY = atoi(tmpstr2);
        BINARY = 0;
        DILUTE = 0;
      }
      else if (strcmp(tmpstr1, "DILUTE") == 0) {
        DILUTE = atoi(tmpstr2);
        BINARY = 0;
        TERNARY = 0;
      }
      else if (strcmp(tmpstr1, "Writecomposition") ==0) {
        WRITECOMPOSITION = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "Noise_phasefield") ==0) {
        NOISE_PHASEFIELD = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "Amp_Noise_Phase") ==0) {
        if (NOISE_PHASEFIELD) {
          AMP_NOISE_PHASE = atof(tmpstr2);
        }
      }
      else if ((strcmp(tmpstr1, "GAMMA") == 0) && (NUMPHASES > 0)) {
        Gamma = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(Gamma, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "Tau") == 0) && (NUMPHASES > 0)) {
        tau_ab = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(tau_ab, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "dab") == 0) && (NUMPHASES > 0)) {
        dab = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(dab, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "ec") == 0) && (NUMPHASES > 0)) {
        ec = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(ec, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "e2") == 0) && (NUMPHASES > 0)) {
        e2 = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(e2, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "e4") == 0) && (NUMPHASES > 0)) {
        e4 = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(e4, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "beta") == 0) && (NUMPHASES > 0)) {
        beta = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(beta, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "Gamma_abc") == 0) && (NUMPHASES > 0)) {
        Gamma_abc = Malloc3M(NUMPHASES, NUMPHASES, NUMPHASES);
        populate_matrix3M(Gamma_abc, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "DIFFUSIVITY") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_diffusivity_matrix(Diffusivity, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "ceq") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(ceq, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "cfill") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(cfill, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "Function_F") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        FUNCTION_F = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "A") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_A_matrix(A, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "slopes") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(slopes, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "Function_W") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        FUNCTION_W = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "EIGEN_STRAIN") == 0) && (NUMPHASES > 0)) {
        populate_symmetric_tensor(eigen_strain, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "VOIGT_ISOTROPIC") == 0) && (NUMPHASES > 0)) {
        populate_cubic_stiffness(Stiffness_c, tmpstr2);
      }
      else {
        printf("Unrecongized parameter : \"%s\"\n", tmpstr1);
      }
    }
  }
  fclose(fr);
  
  char outfile[100];
  
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
  
  strcpy(key, "STARTTIME");
  
  PRINT_LONG(key, STARTTIME, fr);
  
  strcpy(key, "RESTART");
  
  PRINT_LONG(key, RESTART, fr);

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
  
  PRINT_STRING_ARRAY(key, Phases, NUMCOMPONENTS, fr);
  
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
    PRINT_SYMMETRIC_TENSOR(key, eigen_strain[i], fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "Stiffness_cubic[%s]",Phases[i]);
    PRINT_VOIGT_CUBIC(key, Stiffness_c[i], fr);
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