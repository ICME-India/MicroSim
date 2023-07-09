void reading_input_parameters(char *argv[]) {
  FILE * fr = fopen(argv[1], "rt");
  
  if(fr == NULL) {
    printf("file %s not found", argv[1]);
  }
  printf("Processor:%d--Reading Infile\n",rank);
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
        phase_map = (char**)malloc(sizeof(char*)*NUMPHASES);
        for (i = 0; i < NUMPHASES; ++i) {
          phase_map[i] = (char*)malloc(sizeof(char)*51);
        }
        filling_type_phase = (struct filling_type* )malloc((NUMPHASES)*sizeof(*filling_type_phase));
      }
      else if (strcmp(tmpstr1,"num_thermo_phases")==0) {
        NUM_THERMO_PHASES = atoi(tmpstr2);
        Phases_tdb = (char**)malloc(sizeof(char*)*NUM_THERMO_PHASES);
        for (i = 0; i < NUM_THERMO_PHASES; ++i) {
          Phases_tdb[i] = (char*)malloc(sizeof(char)*51);
        }
      }
      
      else if (strcmp(tmpstr1,"NUMCOMPONENTS")==0) {
        NUMCOMPONENTS = atoi(tmpstr2);
        Components = (char**)malloc(sizeof(char*)*NUMCOMPONENTS);
        for (i = 0; i < NUMCOMPONENTS; ++i) {
          Components[i] = (char*)malloc(sizeof(char)*51);
        }
        if (((NUMCOMPONENTS-1)>0) && (NUMPHASES>0)) {
          Diffusivity         = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
          DiffusivityInv      = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
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

          eigen_strain_phase = (struct symmetric_tensor *)malloc(NUMPHASES*sizeof(*eigen_strain_phase));
          stiffness_phase    = (struct Stiffness_cubic *)malloc(NUMPHASES*sizeof(*stiffness_phase));
          stiffness_phase_n  = (struct Stiffness_cubic *)malloc(NUMPHASES*sizeof(*stiffness_phase_n));
          stiffness_t_phase  = (struct Stiffness_tetragonal *)malloc(NUMPHASES*sizeof(*stiffness_t_phase));


          for (i=0; i < 6; i++) {
            boundary[i] = (struct bc_scalars*)malloc(4*sizeof(*boundary[i])); //3=Number of scalar fields
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
      else if (strcmp(tmpstr1,"STARTTIME")==0) {
        STARTTIME = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"RESTART")==0) {
        RESTART = atol(tmpstr2);
      }
      else if ((strcmp(tmpstr1,"numworkers")==0) && RESTART) {
        numworkers = atol(tmpstr2);
        printf("%d\n", numworkers);
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
      else if (strcmp(tmpstr1,"WRITEHDF5")==0) {
        WRITEHDF5 = atoi(tmpstr2);
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
      else if (strcmp(tmpstr1,"temperature")==0) {
        temperature = atof(tmpstr2);
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
      else if (strcmp(tmpstr1,"tdb_phases")==0) {
        populate_string_array(Phases_tdb, tmpstr2, NUM_THERMO_PHASES);
      }
      else if (strcmp(tmpstr1,"phase_map")==0) {
        populate_string_array(phase_map, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "Function_anisotropy") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        FUNCTION_ANISOTROPY = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Anisotropy_type") == 0) && (FUNCTION_ANISOTROPY !=0)) {
         FOLD = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Rotation_matrix") == 0) && (FUNCTION_ANISOTROPY !=0)) {
        //RotAngles = MallocV(3);
        RotAngles = Malloc3M(NUMPHASES, NUMPHASES, 3);
        //Read_Rotation_Angles(RotAngles, tmpstr2);
        populate_rotation_matrix(Rotation_matrix, Inv_Rotation_matrix, tmpstr2);
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
      else if ((strcmp(tmpstr1, "Gamma_abc") == 0) && (NUMPHASES > 0)) {
        Gamma_abc = Malloc3M(NUMPHASES, NUMPHASES, NUMPHASES);
        populate_matrix3M(Gamma_abc, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "DIFFUSIVITY") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_diffusivity_matrix(Diffusivity, tmpstr2, NUMCOMPONENTS);
        for ( a = 0; a < NUMPHASES; a++ ) { 
          matinvnew(Diffusivity[a], DiffusivityInv[a], NUMCOMPONENTS-1);
        }
      }
      else if ((strcmp(tmpstr1, "ceq") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(ceq, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "cfill") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(cfill, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "Function_F") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        FUNCTION_F = atoi(tmpstr2);
        // if (FUNCTION_F == 2) {
          c_guess = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
        // }
      }
      else if ((strcmp(tmpstr1, "A") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_A_matrix(A, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "slopes") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(slopes, tmpstr2, NUMCOMPONENTS);
      }
      else if((strcmp(tmpstr1, "c_guess") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)){
        populate_thermodynamic_matrix(c_guess, tmpstr2, NUMCOMPONENTS);
      }
      else if (strcmp(tmpstr1,"tNoiseStart")==0) {
        tNoiseStart = atol(tmpstr2);
      }
      /*else if (strcmp(tmpstr1,"DirectionalSolidification")==0) {
        DirectionalSolidification = atoi(tmpstr2);
      }*/
      else if (strcmp(tmpstr1,"TLiquidus")==0) {
        TLiquidus = atof(tmpstr2);
      }
      /*else if (strcmp(tmpstr1,"TemperatureGradient")==0) {
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
      }*/
      else if (strcmp(tmpstr1, "atr") == 0) {
        atr = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "CLplatformID") == 0) {
        platformnumber = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1, "CLdeviceID") == 0) {
        devicenumber = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1,"tdbfname")==0) {
        strcpy(tdbfname,tmpstr2);
      }
      else if ((strcmp(tmpstr1, "ELASTICITY") == 0)) {
        ELASTICITY = atoi(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "EIGEN_STRAIN") == 0) && (NUMPHASES > 0) && ELASTICITY) {
        populate_symmetric_tensor(eigen_strain_phase, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "VOIGT_ISOTROPIC") == 0) && (NUMPHASES > 0) && ELASTICITY) {
        populate_cubic_stiffness(stiffness_phase, tmpstr2);
      }
      else if ((strcmp(tmpstr1, "VOIGT_CUBIC") == 0) && (NUMPHASES > 0) && ELASTICITY) {
        populate_cubic_stiffness(stiffness_phase, tmpstr2);
      }
      else if ((strcmp(tmpstr1, "rho") == 0) && ELASTICITY) {
        rho = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "damping_factor") == 0) && ELASTICITY) {
        damping_factor = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "tolerance") == 0) && ELASTICITY) {
        tolerance = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "max_iterations") == 0) && ELASTICITY) {
        MAX_ITERATIONS = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "deltat_e") == 0) && ELASTICITY) {
        deltat_e = atof(tmpstr2);
      }
//       else {
//         printf("Unrecongized parameter : \"%s\"\n", tmpstr1);
//       }
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
  
  strcpy(key,"DELTA_t");
  
  PRINT_DOUBLE(key, deltat, fr);
  
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

  strcpy(key,"T");
  
  PRINT_DOUBLE(key, T, fr);

  strcpy(key,"epsilon");
  
  PRINT_DOUBLE(key, epsilon, fr);
  
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
    sprintf(key, "cfill[Solidus,%s]",Phases[i]);
    PRINT_VECTOR(key, cfill[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "cfill[Liquidus,%s]",Phases[i]);
    PRINT_VECTOR(key, cfill[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
  }

  /*******/

  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "ceq[%le,%s]",T,Phases[i]);
    PRINT_VECTOR(key, ceq[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "ceq[%le,%s]",T,Phases[i]);
    PRINT_VECTOR(key, ceq[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "cfill[%le,%s]",TLiquidus, Phases[i]);
    PRINT_VECTOR(key, cfill[i][i], NUMCOMPONENTS-1, fr);
  }
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "cfill[%le,%s]",TLiquidus, Phases[i]);
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

  strcpy(key, "tNoiseStart");
  
  PRINT_INT(key, tNoiseStart, fr);

  strcpy(key,"TLiquidus");
  
  PRINT_DOUBLE(key, TLiquidus, fr);
  
  strcpy(key, "atr");
  
  PRINT_INT(key, atr, fr);
  
  strcpy(key, "CLplatformID");
  
  PRINT_INT(key, platformnumber, fr);
  
  strcpy(key, "CLdeviceID");
  
  PRINT_INT(key, devicenumber, fr);

  fprintf(fr, "tdbfname = %s\n\n",tdbfname);

  //fprintf(fr, "Time step (Dimensional) = %le s\n", deltat*(Gamma[0][1]*V/(TLiquidus*R))*(Gamma[0][1]*V/(TLiquidus*R))/Diffusivity[1][0][0]); 

  //fprintf(fr, "Grid size (x) (Dimensional) = %le m\n", deltax*Gamma[0][1]*V/(TLiquidus*R));

  //fprintf(fr, "Grid size (y) (Dimensional) = %le m\n", deltay*Gamma[0][1]*V/(TLiquidus*R));


//   PRINT_BOUNDARY_CONDITIONS(fr);
  
  fclose(fr);
}
