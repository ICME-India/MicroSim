#ifndef MICROSIM_READ_INP_PARAMS
#define MICROSIM_READ_INP_PARAMS
#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <sys/stat.h>

#include "global_vars.h"
#include "utility_functions_new.h"

#define msExpectedSizeIs(SIZE)\
    ms_Throw_SizeExpectationError((size_t)(SIZE), Name, Tokens, NumTokens);

void ms_GlobalInit_PhaseFieldMatrices()
{
    start   = msMalloc(int, 3);
    end     = msMalloc(int, 3);
    averow  = msMalloc(long,3);
    rows    = msMalloc(long,3);
    offset  = msMalloc(long,3);
    extra   = msMalloc(long,3);   

    Diffusivity         = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);
    ceq                 = ms_Malloc_3D(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
    cfill               = ms_Malloc_3D(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
    c_guess             = ms_Malloc_3D(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
    ceq_coeffs          = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1,               4);
    slopes              = ms_Malloc_3D(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
    dcbdT               = ms_Malloc_3D(NUMPHASES, NUMPHASES,       NUMCOMPONENTS-1);
    A                   = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);

    DELTA_T             = ms_Malloc_2D(NUMPHASES,  NUMPHASES);
    DELTA_C             = ms_Malloc_2D(NUMPHASES,  NUMCOMPONENTS-1);
    dcbdT_phase         = ms_Malloc_2D(NUMPHASES,  NUMCOMPONENTS-1);
    B                   = ms_Malloc_2D(NUMPHASES,  NUMCOMPONENTS-1);
    Beq                 = ms_Malloc_2D(NUMPHASES,  NUMCOMPONENTS-1);
    dBbdT               = ms_Malloc_2D(NUMPHASES,  NUMCOMPONENTS-1);
    C                   = ms_Malloc_1D(NUMPHASES);

    cmu                 = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1,   NUMCOMPONENTS-1);
    muc                 = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1,   NUMCOMPONENTS-1);
    
    Rotation_matrix     = ms_Malloc_4D(NUMPHASES,       NUMPHASES,   3,             3);
    Inv_Rotation_matrix = ms_Malloc_4D(NUMPHASES,       NUMPHASES,   3,             3);
    Rotated_qab         = ms_Malloc_1D(3);

    // ms_Init_RotationMatrix(Rotation_matrix, Inv_Rotation_matrix);
    
    eigen_strain_phase  = msMalloc(struct symmetric_tensor,      NUMPHASES);
    stiffness_phase     = msMalloc(struct Stiffness_cubic,       NUMPHASES);
    stiffness_phase_n   = msMalloc(struct Stiffness_cubic,       NUMPHASES);
    stiffness_t_phase   = msMalloc(struct Stiffness_tetragonal,  NUMPHASES);
    
    for (i = 0; i < 6; i++)
    {
        /// 5 fields : phi, mu, T, u and F
        boundary[i] = msMalloc(struct bc_scalars, 5);
    }

    ms_BoundaryConditions_PopulateDefault();
    ms_BoundaryValues_PopulateDefault();

    filling_type_phase = msMalloc(struct filling_type, 6);

    Gamma       = ms_Malloc_2D(NUMPHASES, NUMPHASES);
    tau_ab      = ms_Malloc_2D(NUMPHASES, NUMPHASES);
    dab         = ms_Malloc_2D(NUMPHASES, NUMPHASES);
    fab         = ms_Malloc_2D(NUMPHASES, NUMPHASES);

    Gamma_abc   = ms_Malloc_3D(NUMPHASES, NUMPHASES, NUMPHASES);
}

void ms_GlobalFree_PhaseFieldMatrices()
{
    free(start);
    free(end);
    free(averow);
    free(rows);
    free(offset);
    free(extra);

    ms_Free_3D(Diffusivity);
    ms_Free_3D(ceq);
    ms_Free_3D(cfill);
    ms_Free_3D(c_guess);
    ms_Free_3D(ceq_coeffs);
    ms_Free_3D(slopes);
    ms_Free_3D(dcbdT);
    ms_Free_3D(A);

    ms_Free_2D(DELTA_T);
    ms_Free_2D(DELTA_C);
    ms_Free_2D(dcbdT_phase);
    ms_Free_2D(B);
    ms_Free_2D(Beq);
    ms_Free_2D(dBbdT);
    ms_Free_1D(C);

    ms_Free_3D(cmu);
    ms_Free_3D(muc);

    ms_Free_4D(Rotation_matrix);
    ms_Free_4D(Inv_Rotation_matrix);
    ms_Free_1D(Rotated_qab);

    ms_Free_2D(Gamma);
    ms_Free_2D(tau_ab);
    ms_Free_2D(dab);
    ms_Free_2D(fab);

    ms_Free_3D(Gamma_abc);

    ms_ArrayString_Free(Components, NUMCOMPONENTS);
    ms_ArrayString_Free(Phases,     NUMPHASES);
    ms_ArrayString_Free(Phases_tdb, NUMPHASES);
    ms_ArrayString_Free(phase_map,  NUMPHASES);
}

// long ms_readvar_l(char *file_path, char * var_name)
// {    
//     long var;
//     char line_buffer[1024];
//     char name_buffer[64];
//     char value_buffer[64];
//     FILE *FPtr = ms_FILE_Open(file_path, "rt");
//     while(fgets(line_buffer,1024,FPtr))
//     {
//         sscanf(line_buffer, "%64s = %64[^;];", name_buffer, value_buffer);
//         if(name_buffer[0] != '#')
//         {
//             if (!strcmp(name_buffer,var_name))
//             {
//                 var = atol(value_buffer);
//             }
//         }
//     }
//     fclose(FPtr);
//     return var;
// }

void ms_ReadInputParameters(char *file_path)
{
    static bool     IS_FIRST_READING = true;
    static time_t   LAST_READ_TIME   = 0;

    struct stat file_stat;
    if (stat(file_path, &file_stat) != 0)
    {
        perror("ms_ReadInputParameters::stat::error");
        exit(EXIT_FAILURE);
    }

    if (LAST_READ_TIME == file_stat.st_mtime)
    {
        return;
    }

    LAST_READ_TIME = file_stat.st_mtime;

    if (!IS_FIRST_READING && (taskid == MASTER))
    {
        fprintf(stderr, "re-reading input file\n");
    }
    
    FILE *FPtr = ms_FILE_Open(file_path, "rt");

    /// static because they are large
    static char Name [128];
    static char Value[4098];

    char    *Line    = NULL,
            **Tokens = NULL;

    size_t  LineLen = 0,
            NumTokens = 0;

    ssize_t read;

    bool    lineIsVar = false,
            lineIsArr = false;

    regex_t reg_var,
            reg_arr;

    regmatch_t    reg_match[3];

    static const char pattern_var[] = "^([a-zA-Z][a-zA-Z0-9_]*)=([^;{}]+);";
    static const char pattern_arr[] = "^([a-zA-Z][a-zA-Z0-9_]*)=\\{([^}]*)\\};";

    if (regcomp(&reg_var, pattern_var, REG_EXTENDED)!=0 || regcomp(&reg_arr, pattern_arr, REG_EXTENDED)!=0)
    {
        fprintf(stderr, "could not compile regex");
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&Line, &LineLen, FPtr)) != -1)
    {
        if (LineLen < 5 || (strchr("#`-*|/", Line[0]) != NULL))
        {
            continue;
        }

        ms_String_RemoveWhitespaces(Line);

        lineIsVar = (0 == regexec(&reg_var, Line, 3, reg_match, 0));
        
        if (!lineIsVar)
        {
            lineIsArr = (0 == regexec(&reg_arr, Line, 3, reg_match, 0));    
        }

        if (!lineIsArr && !lineIsVar)
        {
            continue;
        }

        ms_Regex_ExtractMatch(Line, reg_match[1], Name,  sizeof(Name));
        ms_Regex_ExtractMatch(Line, reg_match[2], Value, sizeof(Value));

        if (Tokens != NULL)
        {
            ms_ArrayString_Free(Tokens, NumTokens);
            Tokens    = NULL;
            NumTokens = 0;
        }       

        if (lineIsArr)
        {   
            NumTokens = ms_String_CountChar(Value, ',') + 1;
            Tokens    = ms_ArrayString_Create(NumTokens, 64);
            ms_ArrayString_Populate(Tokens, Value, NumTokens);
        }

        if (IS_FIRST_READING && taskid == MASTER)
        {
            if (lineIsVar)
            {
                fprintf(stderr, "%s = %s\n", Name, Value);
            }
            else
            {
                ms_ArrayString_Print(stderr, Name, Tokens, NumTokens);
            }
        }

        /// the Barrier only exists to bind the print statements together
        MPI_Barrier(MPI_COMM_WORLD);

        /**
         * Input variables to be read one on starting the sim
         */

        if (IS_FIRST_READING)
        {
            if (!strcmp(Name, "DIMENSION"))
            {
                DIMENSION = atol(Value);
                assert(DIMENSION == 2 || DIMENSION == 3);
            }
            else if (!strcmp(Name, "MESH_X"))    { MESH_X = atol(Value); } 
            else if (!strcmp(Name, "MESH_Y"))    { MESH_Y = atol(Value); }
            else if (!strcmp(Name, "MESH_Z"))    { MESH_Z = (DIMENSION == 3) ? atol(Value) : 1; }
            else if (!strcmp(Name, "DELTA_X"))   { deltax = atof(Value); }
            else if (!strcmp(Name, "DELTA_Y"))   { deltay = atof(Value); }
            else if (!strcmp(Name, "DELTA_Z"))   { deltaz = (DIMENSION == 3) ? atof(Value) : 1; }

            else if (!strcmp(Name, "NUMPHASES"))
            {
                NUMPHASES = atoi(Value);
                assert(NUMPHASES > 1);
            }
            else if (!strcmp(Name, "NUMCOMPONENTS"))
            {
                NUMCOMPONENTS = atoi(Value);
                assert(NUMCOMPONENTS > 1);

                BINARY  = NUMCOMPONENTS == 2;
                TERNARY = NUMCOMPONENTS == 3;

                ms_GlobalInit_PhaseFieldMatrices();
            }        

            else if (!strcmp(Name, "NSMOOTH"))      { nsmooth   = atol(Value); }
            else if (!strcmp(Name, "STARTTIME"))    { STARTTIME = atol(Value); }
            else if (!strcmp(Name, "RESTART"))      { RESTART   = atol(Value); }
            else if (!strcmp(Name, "INVERT_GSL"))   { USE_GSL_MATINV = atol(Value) == 1;}

            else if (!strcmp(Name, "WRITEFORMAT"))
            {
                ASCII =  (strcmp(Value,"ASCII")  == 0);
                ASCII = !(strcmp(Value,"BINARY") == 0);
            }

            else if (!strcmp(Name, "WRITEHDF5"))        { WRITEHDF5        = atoi(Value) == 1; }
            else if (!strcmp(Name, "Noise_phasefield")) { NOISE_PHASEFIELD = atoi(Value) == 1; }

            else if (!strcmp(Name, "COMPONENTS"))
            {
                msExpectedSizeIs(NUMCOMPONENTS);
                Components = Tokens;
                Tokens = NULL;
            }
            else if (!strcmp(Name, "PHASES"))
            {
                msExpectedSizeIs(NUMPHASES);
                Phases = Tokens;
                Tokens = NULL;
            }
            else if (!strcmp(Name, "num_thermo_phases"))
            {
                NUM_THERMO_PHASES = atoi(Value);
            }
            else if (!strcmp(Name, "tdb_phases"))
            {
                msExpectedSizeIs(NUM_THERMO_PHASES);
                Phases_tdb = Tokens;
                Tokens = NULL;
            }
            else if (!strcmp(Name, "phase_map"))
            {
                msExpectedSizeIs(NUMPHASES);
                phase_map = Tokens;
                Tokens = NULL;
            }

            else if (!strcmp(Name, "BOUNDARY"))
            {
                msExpectedSizeIs(7);
                ms_BoundaryConditions_Populate(Tokens);
            }

            else if (!strcmp(Name, "BOUNDARY_VALUE"))
            {
                msExpectedSizeIs(7);
                ms_BoundaryValues_Populate(Tokens);
            }

            else if (!strcmp(Name, "epsilon"))   { epsilon    = atof(Value); }
            else if (!strcmp(Name, "tau"))       { tau        = atof(Value); }
            else if (!strcmp(Name, "R"))         { R          = atof(Value); }
            else if (!strcmp(Name, "V"))         { V          = atof(Value); }
            else if (!strcmp(Name, "OBSTACLE"))  { OBSTACLE   = atoi(Value); }
            else if (!strcmp(Name, "Function_W")){ FUNCTION_W = atoi(Value); }

            else if (!strcmp(Name,  "DILUTE"))
            {
                DILUTE = atoi(Value) == 1;
            }
            else if (!strcmp(Name, "Tau"))
            {
                msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
                ms_PopulateMatrix_Symmetric_ab(tau_ab, NUMPHASES, Tokens);
            }
            else if (!strcmp(Name, "ceq"))
            {
                msExpectedSizeIs(2 + (NUMCOMPONENTS-1));
                ms_PopulateMatrix_Thermodynamic_abi(ceq, Tokens);
            }
            else if (!strcmp(Name, "cfill"))
            {
                msExpectedSizeIs(2 + (NUMCOMPONENTS-1));
                ms_PopulateMatrix_Thermodynamic_abi(cfill, Tokens);
            }

            else if (!strcmp(Name, "Function_F")) { FUNCTION_F = atoi(Value); }
            else if (!strcmp(Name,  "ISOTHERMAL"))
            {
                ISOTHERMAL = atoi(Value) == 1;
                TEMPGRADY  = !ISOTHERMAL;
            }
            
            else if (!strcmp(Name, "Equilibrium_temperature")) { Teq   = atof(Value); }
            else if (!strcmp(Name, "Filling_temperature"))     { Tfill = atof(Value); }

            else if ((FUNCTION_F == 5) && !strcmp(Name, "Latent_heat")) { Lf = atof(Value);}
            else if ((FUNCTION_F == 5) && !strcmp(Name, "Thermal_conductivity")) { therm_cond = atof(Value);}

            else if (!strcmp(Name, "A"))
            {
                msExpectedSizeIs(1 + (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1));
                ms_PopulateMatrix_A_aij(A, Tokens);
            }
            else if (!strcmp(Name, "slopes"))
            {
                msExpectedSizeIs(2 + (NUMCOMPONENTS-1));
                ms_PopulateMatrix_Thermodynamic_abi(slopes, Tokens);
            }
            else if(!strcmp(Name, "c_guess"))
            {
                msExpectedSizeIs(2 + (NUMCOMPONENTS-1));
                ms_PopulateMatrix_Thermodynamic_abi(c_guess, Tokens);
            }

            /**
             * Elastacity variables
             */

            else if (!strcmp(Name, "ELASTICITY")) 
            { 
                ELASTICITY = atoi(Value) == 1;
                if (ELASTICITY)
                {
                    ec   = ms_Malloc_2D(NUMPHASES, NUMPHASES);
                    e2   = ms_Malloc_2D(NUMPHASES, NUMPHASES);
                    e4   = ms_Malloc_2D(NUMPHASES, NUMPHASES);
                    // beta = ms_Malloc_2D(NUMPHASES, NUMPHASES);
                }
            }
            else if (ELASTICITY && !strcmp(Name, "ec"))
            {
                msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
                ms_PopulateMatrix_Symmetric_ab(ec, NUMPHASES, Tokens);
            }
            else if (ELASTICITY && !strcmp(Name, "e2"))
            {
                msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
                ms_PopulateMatrix_Symmetric_ab(e2, NUMPHASES, Tokens);
            }
            else if (ELASTICITY && !strcmp(Name, "e4"))
            {
                msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
                ms_PopulateMatrix_Symmetric_ab(e4, NUMPHASES, Tokens);
            }
            // else if (ELASTICITY && !strcmp(Name, "beta"))
            // {
            //     msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
            //     ms_PopulateMatrix_Symmetric_ab(beta, NUMPHASES, Tokens);
            // }
            else if (!strcmp(Name, "EIGEN_STRAIN")) 
            {
                msExpectedSizeIs(7);
                ms_PopulateSymmetricTensor(eigen_strain_phase, Tokens);
            }
            else if (!strcmp(Name, "VOIGT_ISOTROPIC"))
            {
                msExpectedSizeIs(4);
                ms_PopulateCubicStiffness(stiffness_phase, Tokens);
            }
            else if (!strcmp(Name, "VOIGT_CUBIC"))
            {
                msExpectedSizeIs(4);
                ms_PopulateCubicStiffness(stiffness_phase, Tokens);
            }

            else if (ELASTICITY && !strcmp(Name, "rho")) {
                rho = atof(Value);
            }
            else if (ELASTICITY && !strcmp(Name, "damping_factor")) {
                damping_factor = atof(Value);
            }
            else if (ELASTICITY && !strcmp(Name, "tolerance")) {
                tolerance = atof(Value);
            }
            else if (ELASTICITY && !strcmp(Name, "max_iterations")) {
                MAX_ITERATIONS = atof(Value);
            }
            else if (ELASTICITY && !strcmp(Name, "deltat_e")) {
                deltat_e = atof(Value);
            }
            else if (!strcmp(Name, "GRAIN_GROWTH")) { 
                GRAIN_GROWTH = atoi(Value);
            }
            else if (!strcmp(Name, "CONSTRAINED"))  {
                CONSTRAINED  = atoi(Value);
            }
            else if (CONSTRAINED && !strcmp(Name, "lambda")) {
                lambda = atof(Value);
            }

            /**
             * LBM variables
             */
            else if (!strcmp(Name, "LBM"))
            {
                LBM = atoi(Value) == 1;
                if (LBM)
                {
                    rho_LBM = msMalloc(double, NUMPHASES);
                    beta_c  = msMalloc(double, NUMCOMPONENTS-1);   
                }
            }

            else if (LBM && !strcmp(Name, "LBM_RESTART"))    { LBM_RESTART = atoi(Value); }
            else if (LBM && !strcmp(Name, "NU_LBM"))         { nu_lbm      = atof(Value); }

            else if (LBM && !strcmp(Name, "nu")) { nu   = atof(Value); }
            else if (LBM && !strcmp(Name, "W0")) { W_0  = atof(Value); }
            else if (LBM && !strcmp(Name, "rho_LBM"))
            {
                msExpectedSizeIs(NUMPHASES);
                for (int i=0; i < NUMPHASES; i++)
                {
                    rho_LBM[i] = atof(Tokens[i]); 
                }
            }
            else if (LBM && !strcmp(Name, "beta_c"))
            {
                msExpectedSizeIs(NUMCOMPONENTS-1);
                for (int i=0; i < NUMCOMPONENTS-1; i++)
                {
                    beta_c[i] = atof(Tokens[i]); 
                }
            }
            else if (!strcmp(Name, "eidt_mode"))
            {
                EIDT.mode           = atol(Value);
                EIDT_SWITCH         = EIDT.mode == 2 || EIDT.mode == 1;
                if (EIDT_SWITCH)
                {
                    EIDT.comp_far_field = msMalloc(double, (NUMCOMPONENTS-1));
                    EIDT.comp_rate      = msMalloc(double, (NUMCOMPONENTS-1));
                    EIDT.comp_center    = msMalloc(double, (NUMCOMPONENTS-1));
                    for (int i = 0; i < NUMCOMPONENTS-1; i++) {
                        EIDT.comp_center[i] = ceq[NUMPHASES-1][NUMPHASES-1][i];
                    }
                }

                EIDT.radius_coeff = 0;
                EIDT.step_coeff = 0;
                EIDT.step0 = 0;
            }
            else if (EIDT_SWITCH && !strcmp(Name, "comp_ff"))
            {
                for (int i = 0; i < NUMCOMPONENTS-1; i++)
                {
                    EIDT.comp_far_field[i]  = atof(Tokens[i]);
                }
            }
        }

        /**
         * Input variables to be read dynamically everytime the function is called
         */

        if (!strcmp(Name, "DELTA_t"))
        { 
            deltat = atof(Value); 
        }

        else if (!strcmp(Name, "NTIMESTEPS"))        { ntimesteps  = atol(Value); }
        else if (!strcmp(Name, "SAVET"))             { saveT       = atol(Value); }
        else if (!strcmp(Name, "TRACK_PROGRESS"))    { time_output = atol(Value); }
        
        else if (!strcmp(Name, "DIFFUSIVITY"))
        {
            // if (atoi(Tokens[0]) == 1) {
            //     msExpectedSizeIs(2 + (NUMCOMPONENTS-1));
            // }
            // else {
            //     msExpectedSizeIs(2 + (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1));
            // }
            ms_PopulateMatrix_Diffusivity_aij(Diffusivity, Tokens);
        }
        else if (!strcmp(Name, "GAMMA"))
        {
            msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
            ms_PopulateMatrix_Symmetric_ab(Gamma, NUMPHASES, Tokens);
        }
        else if (!strcmp(Name, "Gamma_abc") && (NUMPHASES > 2))
        {
            msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)*(NUMPHASES-2)/6);
            ms_PopulateMatrix_Symmetric_abc(Gamma_abc, NUMPHASES, Tokens);
        }
        else if (!strcmp(Name, "Function_anisotropy"))
        {
            FUNCTION_ANISOTROPY = atoi(Value);
            ANISOTROPY = FUNCTION_ANISOTROPY != 0;
        }
        else if (ANISOTROPY && !strcmp(Name, "Anisotropy_type"))
        {
            FOLD = atoi(Value);
        }
        else if (ANISOTROPY && !strcmp(Name, "dab"))
        {
            msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
            ms_PopulateMatrix_Symmetric_ab(dab, NUMPHASES, Tokens);
        }
        else if (ANISOTROPY && !strcmp(Name, "fab"))
        {
            msExpectedSizeIs(NUMPHASES*(NUMPHASES-1)/2);
            ms_PopulateMatrix_Symmetric_ab(fab, NUMPHASES, Tokens);
            USING_FAB = true;
        }
        else if (ANISOTROPY && !strcmp(Name, "Rotation_matrix"))
        {
            msExpectedSizeIs(5);
            ms_PopulateMatrix_Rotation_abqr(Rotation_matrix, Inv_Rotation_matrix, Tokens);
        }

        else if (!strcmp(Name, "T")) { T = atof(Value); }
        
        else if (TEMPGRADY && !strcmp(Name, "Tempgrady"))
        {
            msExpectedSizeIs(5);
            if (IS_FIRST_READING)
            {
                temperature_gradientY.base_temp         = atof(Tokens[0]);
                temperature_gradientY.Distance          = atof(Tokens[2]);
                temperature_gradientY.gradient_OFFSET   = atof(Tokens[3]);                
            }
            temperature_gradientY.DeltaT   = atof(Tokens[1]);
            temperature_gradientY.velocity = atof(Tokens[4]);
        }
        
        else if (!strcmp(Name, "Shift"))             { SHIFT  = atoi(Value) == 1;}
        else if (SHIFT && !strcmp(Name, "Shiftj"))   { shiftj = atol(Value); assert(shiftj>3);}
        
        else if (NOISE_PHASEFIELD && !strcmp(Name, "Amp_Noise_Phase")) 
        {
            AMP_NOISE_PHASE = atof(Value);
        }

        else if (LBM && !strcmp(Name, "LBM_SAVE_FREQ"))  { lbmSaveFreq = atoi(Value); }
        else if (LBM && !strcmp(Name, "gy")) { g_y  = atof(Value); }
        else if (LBM && !strcmp(Name, "dt")) { dt   = atof(Value); }

        else if ((EIDT.mode == 1) && !strcmp(Name, "comp_ff_rate"))
        {
            for (int i = 0; i < NUMCOMPONENTS-1; i++)
            {
                EIDT.comp_rate[i]  = atof(Tokens[i]);
            }
        }
        else if ((EIDT.mode == 2) && !strcmp(Name, "rad_coeff"))
        {
            EIDT.radius_coeff = atof(Value);
        }
        else if ((EIDT.mode == 2) && !strcmp(Name, "step_coeff"))
        {
            EIDT.step_coeff = atof(Value);   
        }
        else if ((EIDT.mode == 2) && !strcmp(Name, "step_zero"))
        {
            EIDT.step0 = atol(Value);
        }
    }

    if (IS_FIRST_READING)
    {
        IS_FIRST_READING = false;
    }
    
    regfree(&reg_var);
    regfree(&reg_arr);
    free(Line);
    fclose(FPtr);
}


#ifdef __cplusplus
}
#endif
#endif