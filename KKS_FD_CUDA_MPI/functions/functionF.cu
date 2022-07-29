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
            perror(filename);
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
                    simParams->F0_A_host[a][i][j] = 0.5*gsl_spline_eval(spline_ThF[simDomain->thermo_phase_host[a]][i][j], simParams->T, acc_ThF[simDomain->thermo_phase_host[a]][i][j]);
                }
                else
                {
                    simParams->F0_A_host[a][i][j] = gsl_spline_eval(spline_ThF[simDomain->thermo_phase_host[a]][i][j], simParams->T, acc_ThF[simDomain->thermo_phase_host[a]][i][j]);
                }
            }
        }
    }

    double c_mat[simDomain->numComponents-1];
    double c_ppt[simDomain->numComponents-1];

    double sum_c = 0.0;


    for (a = 0; a < simDomain->numPhases; a++)
    {
        if (a == 0)
        {
            for (i = 0; i < simDomain->numComponents-1; i++)
                simParams->F0_B_host[0][i] = 0.0;
            simParams->F0_C_host[0] = 0.0;
        }
        else
        {
            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                for (k = 0; k < simDomain->numComponents-1; k++)
                {
                    c_mat[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][1], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][1]);
                    c_ppt[k] = gsl_spline_eval(spline_ES[simDomain->thermo_phase_host[a-1]][k][0], simParams->T, acc_ES[simDomain->thermo_phase_host[a-1]][k][0]);

                    if (k != i)
                        sum_c += simParams->F0_A_host[0][k][i]*c_mat[k] - simParams->F0_A_host[a][k][i]*c_ppt[k];
                }

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

        double *dmudc = (double*)malloc(sizeof(double)*(simDomain->numComponents-1)*(simDomain->numComponents-1));
        long i, j, index;
        long a,k;
        double **dmudc_matrix;
        double sum = 0.0;

        double *y = (double*)malloc(simDomain->numComponents);
        dmudc_matrix = malloc2M((simDomain->numComponents-1),(simDomain->numComponents-1));

        for (a = 0; a < simDomain->numPhases; a++)
        {
            sum = 0.0;
            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                y[k] = simParams->ceq_host[a][a][k];
                sum += y[k];
            }

            y[simDomain->numComponents-1] = 1.0 - sum;

            (*dmudc_tdb[simDomain->thermo_phase_host[a]])((simParams->T + simParams->Teq)/2.0, y, dmudc);

            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                for (j = 0; j < simDomain->numComponents-1; j++)
                {
                    index = i*(simDomain->numComponents-1) + j;
                    dmudc_matrix[i][j] = dmudc[index];
                }
            }

            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                for (j = 0; j < simDomain->numComponents-1; j++)
                {
                    if (i == j)
                        simParams->F0_A_host[a][i][j] = 0.5*dmudc_matrix[i][j];
                    else
                        simParams->F0_A_host[a][i][j] = dmudc_matrix[i][j];

                }
            }
        }

        for (a = simDomain->numPhases-1; a > 0; a--)
        {
            // long   index;
            //double **dmudc_matrix;
            double sum_cs[simDomain->numComponents-1];
            double sum_cl[simDomain->numComponents-1];
            double mu_s[simDomain->numComponents-1];
            double dmu_s[simDomain->numComponents-1];
            double mu_l[simDomain->numComponents-1];
            double dmu_l[simDomain->numComponents-1];
            double f_s, f_l;
            double df_s, df_l;
            double DT = 1;
            double dS;
            double dmuS, dmuL;

            for (k = 0; k<simDomain->numComponents-1; k++)
            {
                y[k] = simParams->ceq_host[a][0][k];
                sum += y[k];
                sum_cl[k] = 0.0;
            }

            y[simDomain->numComponents-1] = 1.0 - sum;

            (*free_energy_tdb[simDomain->thermo_phase_host[0]])((simParams->T + simParams->Teq)/2.0 - DT, y, &f_l);
            (*free_energy_tdb[simDomain->thermo_phase_host[0]])((simParams->T + simParams->Teq)/2.0 + DT, y, &df_l);

            (*Mu_tdb[simDomain->thermo_phase_host[0]])((simParams->T + simParams->Teq)/2.0-DT, y, mu_l);
            (*Mu_tdb[simDomain->thermo_phase_host[0]])((simParams->T + simParams->Teq)/2.0+DT, y, dmu_l);


            sum = 0.0;
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                y[k] = simParams->ceq_host[a][a][k];
                sum += y[k];
                sum_cs[k] = 0.0;
            }

            y[simDomain->numComponents-1] = 1.0 - sum;

            (*free_energy_tdb[simDomain->thermo_phase_host[a]])((simParams->T + simParams->Teq)/2.0 - DT, y, &f_s);
            (*free_energy_tdb[simDomain->thermo_phase_host[a]])((simParams->T + simParams->Teq)/2.0 + DT, y, &df_s);

            (*Mu_tdb[simDomain->thermo_phase_host[a]])((simParams->T + simParams->Teq)/2.0-DT, y, mu_s);
            (*Mu_tdb[simDomain->thermo_phase_host[a]])((simParams->T + simParams->Teq)/2.0+DT, y, dmu_s);

            dS   = 0.5*(df_s - f_s)/DT - 0.5*(df_l - f_l)/DT;

            dmuS = 0.0;
            dmuL = 0.0;

            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                dmuS += 0.5*(dmu_s[k] - mu_s[k])*(simParams->ceq_host[a][a][k] - simParams->ceq_host[a][0][k])/DT;
                dmuL += 0.5*(dmu_l[k] - mu_l[k])*(simParams->ceq_host[a][a][k] - simParams->ceq_host[a][0][k])/DT;
            }

            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                for (int l = 0; l < simDomain->numComponents-1; l++)
                {
                    if (k == l)
                    {
                        sum_cs[k] += 2.0*(simParams->ceq_host[a][a][l] - simParams->ceq_host[a][0][l])*simParams->F0_A_host[a][l][k];
                        sum_cl[k] += 2.0*(simParams->ceq_host[a][a][l] - simParams->ceq_host[a][0][l])*simParams->F0_A_host[0][l][k];
                    }
                    else
                    {
                        sum_cs[k] += (simParams->ceq_host[a][a][l] - simParams->ceq_host[a][0][l])*simParams->F0_A_host[a][l][k];
                        sum_cl[k] += (simParams->ceq_host[a][a][l] - simParams->ceq_host[a][0][l])*simParams->F0_A_host[0][l][k];
                    }

                }

                simParams->slopes[a][a][k] = sum_cs[k]/(dS - dmuS);
                simParams->slopes[a][0][k] = sum_cl[k]/(dS - dmuL);
                simParams->slopes[0][a][k] = simParams->slopes[a][0][k];

            }
        }

        for (a = 0; a < simDomain->numPhases; a++)
        {
            sum = 0.0;
            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                y[k] = simParams->cguess_host[a][a][k];
                sum += y[k];
            }

            y[simDomain->numComponents-1] = 1.0 - sum;

            (*dmudc_tdb[simDomain->thermo_phase_host[a]])(simParams->T, y, dmudc);

            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                for (j = 0; j < simDomain->numComponents-1; j++)
                {
                    index = i*(simDomain->numComponents-1) + j;
                    dmudc_matrix[i][j] = dmudc[index];
                }
            }

            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                for (j = 0; j < simDomain->numComponents-1; j++)
                {
                    if (i == j)
                        simParams->F0_A_host[a][i][j] = 0.5*dmudc_matrix[i][j];
                    else
                        simParams->F0_A_host[a][i][j] = dmudc_matrix[i][j];
                }
            }
        }


        for (a = 1; a < simDomain->numPhases; a++)
        {
            simParams->DELTA_T[a][0] = 0.0;
            simParams->DELTA_T[a][a] = 0.0;

            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                simParams->DELTA_T[a][a] += simParams->slopes[a][a][k]*(simParams->ceq_host[a][0][k] - simParams->ceq_host[a][a][k]);
                simParams->DELTA_T[a][0] += simParams->slopes[a][0][k]*(simParams->ceq_host[a][0][k] - simParams->ceq_host[a][a][k]);
                simParams->DELTA_C[a][k]  = simParams->ceq_host[a][0][k] - simParams->ceq_host[a][a][k];
            }
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                simParams->dcbdT[a][a][k] = simParams->DELTA_C[a][k]/simParams->DELTA_T[a][a];
                simParams->dcbdT[a][0][k] = simParams->DELTA_C[a][k]/simParams->DELTA_T[a][0];
            }
            for (k = 0; k < simDomain->numComponents-1; k++)
            {
                simParams->dBbdT[a][k] = 2.0*(simParams->F0_A_host[0][k][k]*simParams->dcbdT[a][0][k] - simParams->F0_A_host[a][k][k]*simParams->dcbdT[a][a][k]);
                for (i = 0; i < simDomain->numComponents-1; i++)
                {
                    if (k != i)
                    {
                        simParams->dBbdT[a][k] += (simParams->F0_A_host[0][k][i]*simParams->dcbdT[a][0][i] - simParams->F0_A_host[a][k][i]*simParams->dcbdT[a][a][i]);
                    }
                }
            }
        }
        // Got here

        for (k = 0; k < simDomain->numComponents-1; k++)
        {
            simParams->dBbdT[0][k] = 0.0;
        }

        for (a = 0; a < simDomain->numPhases; a++)
        {
            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                double c_Mat[simDomain->numComponents-1];
                double c_PPT[simDomain->numComponents-1];

                double sum_c = 0.0;

                if (a != 0)
                {
                    for (k = 0; k < simDomain->numComponents-1; k++)
                    {
                        c_Mat[k] = simParams->ceq_host[a][0][k];
                        c_PPT[k] = simParams->ceq_host[a][a][k];

                        if (k != i)
                        {
                            sum_c += simParams->F0_A_host[0][k][i]*c_Mat[k] - simParams->F0_A_host[a][k][i]*c_PPT[k];
                        }
                    }
                    simParams->F0_Beq_host[a][i] = (2.0*(simParams->F0_A_host[0][i][i]*c_Mat[i] - simParams->F0_A_host[a][i][i]*c_PPT[i]) + sum_c);
                }
            }
        }

        for (a = 0; a < simDomain->numPhases; a++)
        {
            for (i = 0; i < simDomain->numComponents-1; i++)
            {
                double c_Mat[simDomain->numComponents-1];
                double c_PPT[simDomain->numComponents-1];

                double sum_c = 0.0;

                if (a != 0)
                {
                    for (k = 0; k < simDomain->numComponents-1; k++)
                    {
                        c_Mat[k] = simParams->ceq_host[a][0][k] - (simParams->DELTA_C[a][k])*(simParams->Teq-simParams->T)/(simParams->DELTA_T[a][0]);
                        c_PPT[k] = simParams->ceq_host[a][a][k] - (simParams->DELTA_C[a][k])*(simParams->Teq-simParams->T)/(simParams->DELTA_T[a][a]);
                        if (k != i)
                        {
                            sum_c += simParams->F0_A_host[0][k][i]*c_Mat[k] - simParams->F0_A_host[a][k][i]*c_PPT[k];
                        }
                    }
                    simParams->F0_B_host[a][i] = (2.0*(simParams->F0_A_host[0][i][i]*c_Mat[i] - simParams->F0_A_host[a][i][i]*c_PPT[i]) + sum_c);
                }
            }
        }

        for (a = 0; a < simDomain->numPhases; a++)
        {
            double c_Mat[simDomain->numComponents-1];
            double c_PPT[simDomain->numComponents-1];

            double sum_c = 0.0;

            if (a != 0)
            {
                for (k = 0; k < simDomain->numComponents-1; k++)
                {
                    c_Mat[k] = simParams->ceq_host[a][0][k] - (simParams->DELTA_C[a][k])*(simParams->Teq-simParams->T)/(simParams->DELTA_T[a][0]);
                    c_PPT[k] = simParams->ceq_host[a][a][k] - (simParams->DELTA_C[a][k])*(simParams->Teq-simParams->T)/(simParams->DELTA_T[a][a]);
                }

                for (i = 0; i < simDomain->numComponents-1; i++) {
                    for (j = 0; j < simDomain->numComponents-1; j++) {
                        if (i <= j) {
                            sum_c += (simParams->F0_A_host[a][i][j]*c_PPT[i]*c_PPT[j] - simParams->F0_A_host[0][i][j]*c_Mat[i]*c_Mat[j]);
                        }
                    }
                }
            }

            simParams->F0_C_host[a] = sum_c;
        }


        free2M(dmudc_matrix, simDomain->numComponents-1);
        free(dmudc);
        free(y);
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
}
