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


          for (i=0; i < 6; i++) {
            boundary[i] = (struct bc_scalars*)malloc(3*sizeof(*boundary[i])); //3=Number of scalar fields
          }
          
        }
      }
      else if (strcmp(tmpstr1,"epsilon")==0) {
        epsilon = atof(tmpstr2);
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
      else if (strcmp(tmpstr1,"T")==0) {
        T = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"tilt_angle")==0) {
        tilt_angle = atof(tmpstr2);
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
      else if ((strcmp(tmpstr1, "Rotation_matrix") == 0) && (FUNCTION_ANISOTROPY !=0)) {
        RotAngles = MallocV(3);
        Read_Rotation_Angles(RotAngles, tmpstr2);
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
      else if (strcmp(tmpstr1, "BINARY") == 0) {
        BINARY = atoi(tmpstr2);
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
      else if ((strcmp(tmpstr1, "dab") == 0) && (NUMPHASES > 0)) {
        dab = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(dab, tmpstr2, NUMPHASES);
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
      else if (strcmp(tmpstr1,"tNoiseStart")==0) {
        tNoiseStart = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"DirectionalSolidification")==0) {
        DirectionalSolidification = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1,"TLiquidus")==0) {
        TLiquidus = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"TemperatureGradient")==0) {
        TemperatureGradient = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"PullingVelocity")==0) {
        PullingVelocity = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"PositionOffset")==0) {
        PositionOffset = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1,"T_Offset")==0) {
        T_Offset = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1, "atr") == 0) {
        atr = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "CLplatformID") == 0) {
        platformnumber = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "CLdeviceID") == 0) {
        devicenumber = atoi(tmpstr2);
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
  
  strcpy(key, "GAMMA");
  
  PRINT_MATRIX(key, Gamma, NUMPHASES, NUMPHASES, fr);
  
  strcpy(key, "dab");
  
  PRINT_MATRIX(key, dab, NUMPHASES, NUMPHASES, fr);
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "DIFFUSIVITY[%s]",Phases[i]);
    PRINT_MATRIX(key, Diffusivity[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
  }
  
  strcpy(key, "Function_anisotropy");
  
  PRINT_INT(key, FUNCTION_ANISOTROPY, fr);
  
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
    sprintf(key, "cfill[Solidus,%s]",Phases[i]);
    PRINT_VECTOR(key, cfill[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "cfill[Liquidus,%s]",Phases[i]);
    PRINT_VECTOR(key, cfill[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
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
