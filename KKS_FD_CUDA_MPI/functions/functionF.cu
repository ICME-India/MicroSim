#include "functionF.cuh"

void function_F_04_init_propertymatrices(domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    //Initialize property matrices
    FILE *fp;
    long a;
    long i, j, k;
    char filename[1000];
    long numlines, lines;
    char *file_contents;

    struct stat sb;

    for (i = 0; i < simDomain->numThermoPhases-1; i++)
    {
        sprintf(filename, "tdbs_encrypted/Composition_%s.csv", simDomain->phases_tdb[i+1]);
        fp = fopen(filename, "r");
        if (stat(filename, &sb) == -1)
        {
            perror("filename");
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < simDomain->numThermoPhases; i++)
    {
        sprintf(filename, "tdbs_encrypted/HSN_%s.csv", simDomain->phases_tdb[i]);
        fp = fopen(filename, "r");
        if (stat(filename, &sb) == -1)
        {
            perror("stat");
            exit(EXIT_FAILURE);
        }
    }

    file_contents = (char*)malloc(sb.st_size);

    numlines = 0;
    while (fscanf(fp, "%[^\n] ", file_contents) != EOF)
    {
        numlines++;
    }
    fclose(fp);

    double ****comp_ES = malloc4M(simDomain->numThermoPhases, simDomain->numComponents-1, 2,                          numlines);
    double ****ThF     = malloc4M(simDomain->numThermoPhases, simDomain->numComponents-1, simDomain->numComponents-1, numlines);
    double   **T_ES    = malloc2M(simDomain->numThermoPhases, numlines);
    double   **T_ThF   = malloc2M(simDomain->numThermoPhases, numlines);

    for (a = 0; a < simDomain->numThermoPhases-1; a++)
    {
        sprintf(filename, "tdbs_encrypted/Composition_%s.csv", simDomain->phases_tdb[a+1]);
        fp = fopen(filename, "r");
        fscanf(fp, "%*[^\n]\n");

        for (lines = 0; lines < numlines-1; lines++)
        {
            fscanf(fp, "%le,", &T_ES[a][lines]);
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                fscanf(fp, "%le,", &comp_ES[a][k][0][lines]);
            }
            for (k=0; k < simDomain->numComponents-1; k++)
            {
                fscanf(fp, "%le,", &comp_ES[a][k][1][lines]);
            }
        }
        fclose(fp);
    }

    for (a = 0; a < simDomain->numThermoPhases; a++)
    {
        sprintf(filename,"tdbs_encrypted/HSN_%s.csv", simDomain->phases_tdb[a]);
        fp = fopen(filename, "r");
        fscanf(fp, "%*[^\n]\n");

        for (lines = 0; lines < numlines-1; lines++)
        {
            fscanf(fp, "%le,", &T_ThF[a][lines]);
            for (j = 0; j < simDomain->numComponents-1; j++)
            {
                for (k = 0; k < simDomain->numComponents-1; k++)
                {
                    if (j <= k)
                    {
                        fscanf(fp, "%le,", &ThF[a][j][k][lines]);
                    }
                    else
                    {
                        ThF[a][j][k][lines] = ThF[a][k][j][lines];
                    }
                }
            }
        }
        fclose(fp);
    }

    gsl_interp_accel ****acc_ES     = (gsl_interp_accel****)malloc((simDomain->numThermoPhases-1)*sizeof(gsl_interp_accel***));
    gsl_spline       ****spline_ES  = (gsl_spline ****)malloc((simDomain->numThermoPhases-1)*sizeof(gsl_spline***));
    gsl_interp_accel ****acc_ThF    = (gsl_interp_accel****)malloc(simDomain->numThermoPhases*sizeof(gsl_interp_accel***));
    gsl_spline       ****spline_ThF = (gsl_spline****)malloc(simDomain->numThermoPhases*sizeof(gsl_spline***));

    for (a = 0; a < simDomain->numThermoPhases; a++)
    {
        acc_ThF[a]    = (gsl_interp_accel***)malloc((simDomain->numComponents-1)*sizeof(gsl_interp_accel**));
        spline_ThF[a] = (gsl_spline***)malloc((simDomain->numComponents-1)*sizeof(gsl_spline**));

        for (k = 0; k < simDomain->numComponents-1; k++)
        {
            acc_ThF[a][k]    = (gsl_interp_accel**)malloc((simDomain->numComponents-1)*sizeof(gsl_interp_accel*));
            spline_ThF[a][k] = (gsl_spline**)malloc((simDomain->numComponents-1)*sizeof(gsl_spline*));

            for (j = 0; j < simDomain->numComponents-1; j++)
            {
                acc_ThF[a][k][j]    = gsl_interp_accel_alloc ();
                spline_ThF[a][k][j] = gsl_spline_alloc (gsl_interp_cspline, numlines-1);
                gsl_spline_init (spline_ThF[a][k][j], T_ThF[a], ThF[a][k][j], numlines-1);
            }
        }
    }

    for (a = 0; a < simDomain->numThermoPhases-1; a++)
    {
        acc_ES[a]     = (gsl_interp_accel***)malloc((simDomain->numComponents-1)*sizeof(gsl_interp_accel**));
        spline_ES[a]  = (gsl_spline***)malloc((simDomain->numComponents-1)*sizeof(gsl_spline**));

        for (k = 0; k < simDomain->numComponents-1; k++)
        {
            acc_ES[a][k]     = (gsl_interp_accel**)malloc((simDomain->numComponents-1)*sizeof(gsl_interp_accel*));
            spline_ES[a][k]  = (gsl_spline**)malloc((simDomain->numComponents-1)*sizeof(gsl_spline*));

            for (i = 0; i < 2; i++)
            {
                acc_ES[a][k][i]    = gsl_interp_accel_alloc();
                spline_ES[a][k][i] = gsl_spline_alloc (gsl_interp_cspline, numlines-1);
                gsl_spline_init(spline_ES[a][k][i],  T_ES[a], comp_ES[a][k][i], numlines-1);
            }
        }
    }

    free4M(comp_ES, simDomain->numThermoPhases, simDomain->numComponents-1, 2);
    free4M(ThF, simDomain->numThermoPhases, simDomain->numComponents-1, simDomain->numComponents-1);
    free2M(T_ES, simDomain->numThermoPhases);
    free2M(T_ThF, simDomain->numThermoPhases);

    for (a = 0; a < simDomain->numPhases; a++)
    {
        for(i = 0; i<simDomain->numComponents-1; i++)
        {
            for(j = 0; j<simDomain->numComponents-1; j++)
            {
                if (i == j)
                {
                    simParams->F0_A_host[a][i][j] = 0.5*gsl_spline_eval(spline_ThF[simDomain->thermo_phase_host[a]][i][j], simParams->T, acc_ThF[simDomain->thermo_phase_host[a]][i][j]) / (simParams->molarVolume * 1.602*1e8);
                }
                else
                {
                    simParams->F0_A_host[a][i][j] = gsl_spline_eval(spline_ThF[simDomain->thermo_phase_host[a]][i][j], simParams->T, acc_ThF[simDomain->thermo_phase_host[a]][i][j]) / (simParams->molarVolume * 1.602*1e8);
                }
            }
        }
    }

    double c_mat[simDomain->numComponents-1];
    double c_ppt[simDomain->numComponents-1];

    double sum_c = 0.0;
    simParams->F0_B_host[0][0] = 0.0;
    simParams->F0_C_host[0] = 0.0;

    for (a = 1; a < simDomain->numPhases; a++)
    {
        for (i = 0; i < simDomain->numComponents-1; i++)
        {
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                c_mat[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][1], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][1]);
                c_ppt[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][0], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][0]);

                printf("%lf\t%lf\n", c_mat[k], c_ppt[k]);

                if (k != i)
                    sum_c += simParams->F0_A_host[0][k][i]*c_mat[k] - simParams->F0_A_host[a][k][i]*c_ppt[k];
            }
                printf("%lf\t%lf\n", simParams->F0_A_host[0][i][i], simParams->F0_A_host[a][i][i]);

            simParams->F0_B_host[a][i] = 2.0*(simParams->F0_A_host[0][i][i]*c_mat[i] - simParams->F0_A_host[a][i][i]*c_ppt[i]) + sum_c;
        }

        sum_c = 0.0;

        for (k = 0; k < simDomain->numComponents-1; k++)
        {
            c_mat[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][1], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][1]);
            c_ppt[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][0], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][0]);

        }

        for (i = 0; i < simDomain->numComponents-1; i++)
        {
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                if (i <= k)
                    sum_c += simParams->F0_A_host[a][i][k]*c_ppt[i]*c_ppt[k] - simParams->F0_A_host[0][i][k]*c_mat[i]*c_mat[k];
            }
        }
        simParams->F0_C_host[a] = sum_c;

        sum_c = 0.0;
    }
}

void calcFreeEnergyCoeffs(domainInfo *simDomain, controls *simControls, simParameters *simParams)
{

    if (simControls->FUNCTION_F == 1)
    {
        simParams->F0_B_host[0][0] = 0.0;
        simParams->F0_C_host[0]    = 0.0;

        for (int i = 1; i < simDomain->numPhases; i++)
        {
            simParams->F0_C_host[i] = 0.0;

            for (int j = 0; j < simDomain->numComponents-1; j++)
            {
                simParams->F0_B_host[i][j] = 2.0*(simParams->F0_A_host[0][j][j]*simParams->ceq_host[0][0][j] - simParams->F0_A_host[i][j][j]*simParams->ceq_host[i][i][j]);

                for (int k = 0; k < simDomain->numComponents-1; k++)
                {
                    if (j == k)
                        continue;
                    simParams->F0_B_host[i][j] += simParams->F0_A_host[i][j][k]*simParams->ceq_host[i][i][k] - simParams->F0_A_host[0][j][k]*simParams->ceq_host[0][0][k];
                }

                for (int k = 0; k <= j; k++)
                {
                    simParams->F0_C_host[i] += simParams->F0_A_host[i][j][k]*simParams->ceq_host[i][i][j]*simParams->ceq_host[i][i][k] - simParams->F0_A_host[0][j][k]*simParams->ceq_host[0][0][j]*simParams->ceq_host[0][0][k];
                }
            }
        }
    }
    else if (simControls->FUNCTION_F == 2)
    {
        for (int a = 0; a < simDomain->numPhases; a++)
        {
            for (int b = 0; b < simDomain->numThermoPhases; b++)
            {
                if (strcmp(simDomain->phase_map[a], simDomain->phases_tdb[b]) == 0)
                {
                    simDomain->thermo_phase_host[a] = b;
                }
            }
        }
    }
    else if (simControls->FUNCTION_F == 3)
    {
        // Mapping thermodynamic phases to the simulated phases
        for (int a = 0; a < simDomain->numPhases; a++)
        {
            for (int b = 0; b < simDomain->numThermoPhases; b++)
            {
                if (strcmp(simDomain->phase_map[a], simDomain->phases_tdb[b]) == 0)
                {
                    simDomain->thermo_phase_host[a] = b;
                }
            }
        }

        // Calculating curvatures for each phase
        for (int a = 0; a < simDomain->numPhases; a++)
            for (int b = 0; b < simDomain->numComponents-1; b++)
                simParams->F0_A_host[a][b][b] = 0.5*evalFunc(dmudc_tdb[simDomain->thermo_phase_host[a]], simParams->cguess_host[a][a][b], simParams->T)/simParams->molarVolume;

        // Calculating linear and constant part of the parabolic free energy function for each phase
        // Set matrix values of the above to 0
        simParams->F0_B_host[0][0] = 0.0;
        simParams->F0_C_host[0]    = 0.0;

        for (int i = 1; i < simDomain->numPhases; i++)
        {
            simParams->F0_C_host[i] = 0.0;

            for (int j = 0; j < simDomain->numComponents-1; j++)
            {
                simParams->F0_B_host[i][j] = 2.0*(simParams->F0_A_host[0][j][j]*simParams->cguess_host[0][0][j] - simParams->F0_A_host[i][j][j]*simParams->cguess_host[i][i][j]);

                for (int k = 0; k < simDomain->numComponents-1; k++)
                {
                    if (j == k)
                        continue;
                    simParams->F0_B_host[i][j] += simParams->F0_A_host[i][j][k]*simParams->cguess_host[i][i][k] - simParams->F0_A_host[0][j][k]*simParams->cguess_host[0][0][k];
                }

                for (int k = 0; k <= j; k++)
                {
                    simParams->F0_C_host[i] += simParams->F0_A_host[i][j][k]*simParams->cguess_host[i][i][j]*simParams->cguess_host[i][i][k] - simParams->F0_A_host[0][j][k]*simParams->cguess_host[0][0][j]*simParams->cguess_host[0][0][k];
                }
            }
        }
    }
    else if (simControls->FUNCTION_F == 4)
    {
        if (simControls->count == simControls->startTime)
        {
            for (int a = 0; a < simDomain->numPhases; a++)
            {
                for (int b = 0; b < simDomain->numThermoPhases; b++)
                {
                    if (strcmp(simDomain->phase_map[a], simDomain->phases_tdb[b]) == 0)
                    {
                        simDomain->thermo_phase_host[a] = b;
                    }
                }
            }
        }

        function_F_04_init_propertymatrices(simDomain, simControls, simParams);
    }

    FILE *fp = fopen("DATA/FreeEnergy.out", "w");

    fprintf(fp, "Count = %d\n", simControls->count);
    for (int i = 0; i < simDomain->numPhases; i++)
    {
        for (int j = 0; j < simDomain->numComponents-1; j++)
        {
            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                fprintf(fp, "A %d\t%d\t%d\t%lf\n", i, j, k, simParams->F0_A_host[i][j][k]*simParams->molarVolume * 1.602*1e8);
            }
            fprintf(fp, "B %d\t%d\t%lf\n", i, j, simParams->F0_B_host[i][j]);
        }
        fprintf(fp, "C %d\t%lf\n", i, simParams->F0_C_host[i]);
    }

    fclose(fp);
}
