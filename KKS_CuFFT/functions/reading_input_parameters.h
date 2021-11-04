void reading_input_parameters(char *argv[]){

    FILE *fr;
    if (fr = fopen(argv[1], "rt"))
        printf("\nReading input parameters from %s\n", argv[1]);
    else
        printf("\nFile %s not found\n", argv[1]);

    int i;
    char tempbuff[1000];
    char tmpstr1[100];
    char tmpstr2[100];

    while (fgets(tempbuff,1000,fr)) {
        sscanf(tempbuff,"%100s = %100[^;];", tmpstr1, tmpstr2);

        if (tmpstr1[0] != '#') {
            if (strcmp(tmpstr1, "DIMENSION") == 0) {
                DIMENSION = atoi(tmpstr2);
                start   = (long*)malloc(3*sizeof(*start));
                end     = (long*)malloc(3*sizeof(*end));
            }
            else if (strcmp(tmpstr1, "MESH_X") == 0) {
                MESH_X = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "MESH_Y") == 0) {
                MESH_Y = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "MESH_Z") == 0) {
                MESH_Z = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_X") == 0) {
                DELTA_X = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_Y") == 0) {
                DELTA_Y = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_Z") == 0) {
                DELTA_Z = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "DELTA_t") == 0) {
                DELTA_t = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "NUMSTEPS") == 0) {
                numsteps = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "NSMOOTH") == 0) {
                nsmooth = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "SAVET") == 0) {
                saveT = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "lambda") == 0) {
                Ln = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "relax_coeff") == 0) {
                relax_coeff = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "Noise_phasefield") == 0) {
                Noise_phasefield = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "Amp_Noise_Phase") == 0) {
                Amp_Noise_Phase = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "T") == 0) {
                T = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "ELAST_INT") == 0) {
                elast_int = atoi(tmpstr2);
            }
            else if (strcmp(tmpstr1, "Writecomposition") ==0) {
                WRITECOMPOSITION = atoi(tmpstr2);
            }
            else if (strcmp(tmpstr1, "EIGEN_STRAIN") == 0) {
                populate_symmetric_tensor(eigen_strain, tmpstr2, NUMPHASES);
            }
            else if (strcmp(tmpstr1, "NUMPHASES") == 0) {
                NUMPHASES = atoi(tmpstr2);

                if (NUMPHASES > 2) {
                    printf("MultiPhase code with more than two phases"
                    " will be released in future."
                    " Reverting to two phases\n");
                    NUMPHASES = 2;
                }

                PHASES = (char**)malloc(sizeof(char*)*NUMPHASES);
                for(i = 0; i < NUMPHASES; i++)
                    PHASES[i] = (char*)malloc(sizeof(char)*51);
            }
            else if (strcmp(tmpstr1, "NUMCOMPONENTS") == 0) {
                NUMCOMPONENTS = atoi(tmpstr2);
                if (NUMCOMPONENTS > 2) {
                    printf("Multicomponent code will be released in future."
                    " Reverting to two components\n");
                    NUMCOMPONENTS = 2;
                }
                COMPONENTS = (char**)malloc(sizeof(char*)*NUMCOMPONENTS);

                for(i = 0; i < NUMCOMPONENTS; i++)
                    COMPONENTS[i] = (char*)malloc(sizeof(char)*51);

                if (((NUMCOMPONENTS-1) > 0) && (NUMPHASES > 0)) {
                    Diffusivity  = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
                    F0  = Malloc3M(NUMCOMPONENTS, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
                    ceq  = Malloc3M(NUMPHASES, NUMPHASES, NUMCOMPONENTS-1);
                    eigen_strain = (struct symmetric_tensor*)malloc(NUMPHASES*sizeof(*eigen_strain));
                    Stiffness_c = (struct Stiffness_cubic*)malloc(NUMPHASES*sizeof(*Stiffness_c));
                }
            }
            else if (strcmp(tmpstr1, "V") == 0) {
                Vm = atof(tmpstr2);
            }
            else if(strcmp(tmpstr2, "ASCII") == 0) {
                ASCII = 0;
            }
            else if (strcmp(tmpstr1,"TRACK_PROGRESS")==0) {
                time_output = atol(tmpstr2);
            }
            else if (strcmp(tmpstr1, "alpha") == 0) {
                alpha = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "c0") == 0) {
                c0 = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "PPT_RADIUS") == 0) {
                ppt_radius = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "VF") == 0) {
                vf = atof(tmpstr2);
            }
            else if (strcmp(tmpstr1, "INITFLAG") == 0){
                initflag = atoi(tmpstr2);
            }
            else if (strcmp(tmpstr1, "INITCOUNT") == 0){
                initcount = atoi(tmpstr2);
            }
            else if((strcmp(tmpstr1, "GAMMA") == 0) && (NUMPHASES > 0)) {
                Gamma = MallocM(NUMPHASES, NUMPHASES);
                populate_matrix(Gamma, tmpstr2, NUMPHASES);
            }
            else if ((strcmp(tmpstr1, "f0") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) > 0)) {
                populate_A_matrix(F0, tmpstr2, NUMCOMPONENTS);
            }
            else if ((strcmp(tmpstr1, "DIFFUSIVITY") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
                populate_diffusivity_matrix(Diffusivity, tmpstr2, NUMCOMPONENTS);
            }
            else if ((strcmp(tmpstr1, "ceq") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) >0)) {
                populate_thermodynamic_matrix(ceq, tmpstr2, NUMCOMPONENTS);
            }
            else if (strcmp(tmpstr1, "COMPONENTS")==0) {
                populate_string_array(COMPONENTS, tmpstr2, NUMCOMPONENTS);
            }
            else if (strcmp(tmpstr1, "PHASES")==0) {
                populate_string_array(PHASES, tmpstr2, NUMCOMPONENTS);
            }
            else if ((strcmp(tmpstr1, "VOIGT_CUBIC") == 0) && (NUMPHASES > 0)) {
                populate_cubic_stiffness(Stiffness_c, tmpstr2);
            }
            else {
                printf("Unrecognized parameter : \"%s\"\n", tmpstr1);
            }
        }
    }

    fclose(fr);

    char outfile[1000];
    strcpy(tmpstr2, argv[1]);

    strcpy(tmpstr1, strtok(tmpstr2, "."));

    sprintf(outfile, "DATA/%s.out", tmpstr1);

    fr = fopen(outfile, "w");

    char key[1000];

    strcpy(key, "DIMENSION");
    PRINT_INT(key, DIMENSION, fr);

    strcpy(key, "MESH_X");
    PRINT_LONG(key, MESH_X, fr);

    strcpy(key, "MESH_Y");
    PRINT_LONG(key, MESH_Y, fr);

    strcpy(key, "MESH_Z");
    PRINT_LONG(key, MESH_Z, fr);

    strcpy(key, "DELTA_X");
    PRINT_DOUBLE(key, DELTA_X, fr);

    strcpy(key, "DELTA_Y");
    PRINT_DOUBLE(key, DELTA_Y, fr);

    strcpy(key, "DELTA_Z");
    PRINT_DOUBLE(key, DELTA_Z, fr);

    strcpy(key, "DELTA_t");
    PRINT_DOUBLE(key, DELTA_t, fr);

    strcpy(key, "NUMPHASES");
    PRINT_INT(key, NUMPHASES, fr);

    strcpy(key, "NUMCOMPONENTS");
    PRINT_INT(key, NUMCOMPONENTS, fr);

    strcpy(key, "NUMSTEPS");
    PRINT_LONG(key, numsteps, fr);

    strcpy(key, "NSMOOTH");
    PRINT_LONG(key, nsmooth, fr);

    strcpy(key, "SAVET");
    PRINT_LONG(key, saveT, fr);

    strcpy(key, "COMPONENTS");
    PRINT_STRING_ARRAY(key, COMPONENTS, NUMCOMPONENTS, fr);

    strcpy(key, "PHASES");
    PRINT_STRING_ARRAY(key, PHASES, NUMCOMPONENTS, fr);

    strcpy(key, "Gamma");
    PRINT_MATRIX(key, Gamma, NUMPHASES, NUMPHASES, fr);

    for (i = 0; i < NUMPHASES; i++) {
        sprintf(key, "Diffusivity[%s]",PHASES[i]);
        PRINT_MATRIX(key, Diffusivity[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
    }

    strcpy(key, "V");
    PRINT_DOUBLE(key, Vm, fr);

    strcpy(key, "relax_coeff");
    PRINT_DOUBLE(key, relax_coeff, fr);

    strcpy(key, "ELAST_INT");
    PRINT_INT(key, elast_int, fr);

    for (i = 0; i < NUMPHASES; i++) {
        sprintf(key, "Stiffness_cubic[%s]", PHASES[i]);
        PRINT_VOIGT_CUBIC(key, Stiffness_c[i], fr);
    }

    for (i = 0; i < NUMPHASES; i++) {
        sprintf(key, "EIGEN_STRAIN[%s]", PHASES[i]);
        PRINT_SYMMETRIC_TENSOR(key, eigen_strain[i], fr);
    }

    strcpy(key, "TRACK_PROGRESS");
    PRINT_LONG(key, time_output, fr);

    strcpy(key, "Noise_phasefield");
    PRINT_INT(key, Noise_phasefield, fr);

    strcpy(key, "Amp_Noise_Phase");
    PRINT_INT(key, Amp_Noise_Phase, fr);

    for (i = 0; i < NUMPHASES; i++) {
        sprintf(key, "ceq[%s]", PHASES[i]);
        PRINT_VECTOR(key, ceq[i][i], NUMCOMPONENTS-1, fr);
    }

    strcpy(key, "c0");
    PRINT_DOUBLE(key, c0, fr);

    for (i = 0; i < NUMCOMPONENTS; i++) {
        sprintf(key, "f0[%s]", COMPONENTS[i]);
        PRINT_MATRIX(key, F0[i], NUMCOMPONENTS-1, NUMCOMPONENTS-1, fr);
    }
    for (i = 0; i < NUMPHASES; i++) {
        sprintf(key, "ceq[%s]", PHASES[i]);
        PRINT_VECTOR(key, ceq[i][NUMPHASES-1], NUMCOMPONENTS-1, fr);
    }

    strcpy(key, "alpha");
    PRINT_DOUBLE(key, alpha, fr);

    strcpy(key, "lambda");
    PRINT_DOUBLE(key, Ln, fr);

    fclose(fr);
}
