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
          //Diffusivity         = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
          AtomicMobility      = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
          ceq                 = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
          cfill               = Malloc3M(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);

        }
      }
      /*else if (strcmp(tmpstr1,"epsilon")==0) {
        epsilon = atof(tmpstr2);
      }*/
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
      else if (strcmp(tmpstr1,"COMPONENTS")==0) {
        populate_string_array(Components, tmpstr2, NUMCOMPONENTS);
      }
      else if (strcmp(tmpstr1,"PHASES")==0) {
        populate_string_array(Phases, tmpstr2, NUMCOMPONENTS);
      }
      //else if ((strcmp(tmpstr1, "GAMMA") == 0) && (NUMPHASES > 0)) {
        //Gamma = MallocM(NUMPHASES, NUMPHASES);
        //populate_matrix(Gamma, tmpstr2, NUMPHASES);
      //}
      //else if ((strcmp(tmpstr1, "DIFFUSIVITY") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        //populate_diffusivity_matrix(Diffusivity, tmpstr2, NUMCOMPONENTS);
      //}
      else if ((strcmp(tmpstr1, "ceq") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(ceq, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "cfill") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_thermodynamic_matrix(cfill, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "AtomicMobility") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        populate_diffusivity_matrix(AtomicMobility, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "Kappa_phi") == 0) && (NUMPHASES > 0)) {
        Kappa_phi = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(Kappa_phi, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "L_phi") == 0) && (NUMPHASES > 0)) {
        L_phi = MallocM(NUMPHASES, NUMPHASES);
        populate_matrix(L_phi, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "Kappa_c") == 0) && ((NUMCOMPONENTS-1) >0)) {
        Kappa_c = MallocM(NUMCOMPONENTS, NUMCOMPONENTS);
        populate_matrix(Kappa_c, tmpstr2, NUMCOMPONENTS);
      }
      else if ((strcmp(tmpstr1, "B_fp") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
        //B_fp = atof(tmpstr2);
        //printf("B_fp=%le\n", B_fp);
        B_fp = MallocV(NUMPHASES-1);
        populate_vector(B_fp, tmpstr2, (NUMPHASES-1));
      }
      else if ((strcmp(tmpstr1, "A_fm") == 0) && (NUMPHASES > 0)) {
        //A_fm = atof(tmpstr2);
        //printf("A_fm=%le\n", A_fm);
        A_fm = MallocV(1);
        populate_vector(A_fm, tmpstr2, 1);
      }
      else if (strcmp(tmpstr1,"spinodal")==0) {
        SPINODAL = atoi(tmpstr2);
      }
      else if (strcmp(tmpstr1,"tdbflag")==0) {
        tdbflag = atoi(tmpstr2);
      }
      else if (tdbflag && (strcmp(tmpstr1,"tdbfname")==0)) {
        strcpy(tdbfname,tmpstr2);
        printf("%s\n",tdbfname);
      }
      else if (strcmp(tmpstr1,"temperature")==0) {
        temperature = atof(tmpstr2);
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
  
  //strcpy(key, "GAMMA");
  
  //PRINT_MATRIX(key, Gamma, NUMPHASES, NUMPHASES, fr);

  //for (i=0; i<NUMPHASES; i++) {
    //sprintf(key, "DIFFUSIVITY[%s]",Phases[i]);
    //PRINT_MATRIX(key, Diffusivity[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
  //}
  
  
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
  
  
  strcpy(key, "Kappa_phi");
  
  PRINT_MATRIX(key, Kappa_phi, NUMPHASES, NUMPHASES, fr);
  
  strcpy(key, "L_phi");
  
  PRINT_MATRIX(key, L_phi, NUMPHASES, NUMPHASES, fr);
  
  strcpy(key, "Kappa_c");
  
  PRINT_MATRIX(key, Kappa_c, NUMCOMPONENTS, NUMCOMPONENTS, fr);
  
  for (i=0; i<NUMPHASES-1; i++) {
    sprintf(key, "B_fp[%s]",Phases[i]);
    PRINT_VECTOR(key, B_fp, NUMPHASES-1, fr);
  }
  
  for (i=0; i<1; i++) {
    sprintf(key, "A_fm[%s]",Phases[NUMPHASES-1]);
    PRINT_VECTOR(key, A_fm, NUMPHASES-1, fr);
  }
  
  for (i=0; i<NUMPHASES; i++) {
    sprintf(key, "AtomicMobility[%s]",Phases[i]);
    PRINT_MATRIX(key, AtomicMobility[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
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