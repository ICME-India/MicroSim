#include "utilityFunctions.h"

#if THERMO == 1
void testThermoFuncs(domainInfo simDomain, simParameters simParams)
{
    FILE *fp[6];

    long N = 3000;

    double y[2], ans[2];

    y[0] = -1.0, y[1] = 2.0;

    fp[0] = fopen("DATA/GE0.bin", "w");
    fp[1] = fopen("DATA/GE1.bin", "w");
    fp[2] = fopen("DATA/mu0.bin", "w");
    fp[3] = fopen("DATA/mu1.bin", "w");
    fp[4] = fopen("DATA/dMu0.bin", "w");
    fp[5] = fopen("DATA/dMu1.bin", "w");

    for (int i = 0; i < N; i++)
    {
        y[0] += (double)3.0/N;
        y[1] -= (double)3.0/N;

        (*free_energy_tdb[simDomain.thermo_phase_host[0]])(simParams.Teq, y, ans);
        fprintf(fp[0], "%le\t%le\n", y[0], ans[0]);

        (*Mu_tdb[simDomain.thermo_phase_host[0]])(simParams.Teq, y, ans);
        fprintf(fp[2], "%le\t%le\n", y[0], ans[0]);

        (*dmudc_tdb[simDomain.thermo_phase_host[0]])(simParams.Teq, y, ans);
        fprintf(fp[4], "%le\t%le\n", y[0], ans[0]);

        (*free_energy_tdb[simDomain.thermo_phase_host[1]])(simParams.Teq, y, ans);
        fprintf(fp[1], "%le\t%le\n", y[0], ans[0]);

        (*Mu_tdb[simDomain.thermo_phase_host[1]])(simParams.Teq, y, ans);
        fprintf(fp[3], "%le\t%le\n", y[0], ans[0]);

        (*dmudc_tdb[simDomain.thermo_phase_host[1]])(simParams.Teq, y, ans);
        fprintf(fp[5], "%le\t%le\n", y[0], ans[0]);
    }

    for (int i = 0; i < 6; i++)
        fclose(fp[i]);
}
#endif

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    long i, j, k;

    tmp = (char**)malloc(sizeof(char*)*NUMPHASES*(NUMPHASES-1)*0.5);

    for (i = 0; i < NUMPHASES*(NUMPHASES-1)*0.5; ++i)
    {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i], token);
    }
    k = 0;
    for(i = 0; i < NUMPHASES; i++)
    {
        for (j = i+1; j < NUMPHASES; j++)
        {
            Mat[i][i] = 0.0;
            Mat[i][j] = atof(tmp[k]);
            Mat[j][i] = Mat[i][j];
            k++;
        }
    }

    for (i = 0; i < NUMPHASES*(NUMPHASES-1)*0.5; ++i)
    {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void populate_matrix3M(double ***Mat, char *tmpstr, long NUMPHASES) {
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    long i,j,k,l;
    long len = NUMPHASES*(NUMPHASES-1)*(NUMPHASES-2)/6;

    if (len != 0)
    {
        tmp = (char**)malloc(sizeof(char*)*len);

        for (i = 0; i < len; ++i) {
            tmp[i] = (char*)malloc(sizeof(char)*10);
        }
        for (i = 0, str1 = tmpstr; ; i++, str1 = NULL) {
            token = strtok_r(str1, "{,}", &saveptr1);
            if (token == NULL)
                break;
            strcpy(tmp[i],token);
        }
        l=0;
        for(i=0; i < NUMPHASES; i++) {
            for (j=i+1; j < NUMPHASES; j++) {
                for (k=j+1; k < NUMPHASES; k++) {
                    Mat[i][i][i] = 0.0;
                    Mat[i][j][j] = 0.0;
                    Mat[i][k][k] = 0.0;

                    Mat[i][j][k] = atof(tmp[l]);
                    Mat[i][k][j] = Mat[i][j][k];
                    Mat[j][i][k] = Mat[i][j][k];
                    Mat[j][k][i] = Mat[i][j][k];
                    Mat[k][i][j] = Mat[i][j][k];
                    Mat[k][j][i] = Mat[i][j][k];

                    l++;
                }
            }
        }
        for (i = 0; i < len; ++i) {
            free(tmp[i]);
        }
        free(tmp);
        tmp = NULL;
    }
}

void populate_thetaij_matrix(double **Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *token;
    char *saveptr1;

    long i, j, k;

    tmp = (char**)malloc(sizeof(char*) * NUMPHASES*(NUMPHASES-1)*0.5);

    for (i = 0; i < NUMPHASES*(NUMPHASES-1)*0.5; i++)
        tmp[i] = (char*)malloc(sizeof(char)*10);

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    k = 0;

    for (i = 0; i < NUMPHASES-1; i++)
    {
        Mat[i][i] = 0.0;

        for (j = i+1; j < NUMPHASES; j++)
        {
            Mat[i][j] = atof(tmp[k]);
            Mat[j][i] = Mat[i][j];
            k++;
        }
    }

    Mat[NUMPHASES-1][NUMPHASES-1] = 0.0;

    for (i = 0; i < NUMPHASES*(NUMPHASES-1)*0.5; i++)
        free(tmp[i]);

    free(tmp);
    tmp = NULL;
}

void populate_thetai_matrix(double *Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *token;
    char *saveptr1;

    long i;
    long len = 2;
    long phase;

    tmp = (char**)malloc(sizeof(char*)*len);

    for (i = 0; i < len; i++)
        tmp[i] = (char*)malloc(sizeof(char)*10);

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i], token);
    }

    phase = atol(tmp[0]);
    if (phase < NUMPHASES)
        Mat[phase] = atof(tmp[1]);

    for (i = 0; i < len; ++i)
        free(tmp[i]);

    free(tmp);
    tmp = NULL;
}

void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    long i,j,k,l;
    long len = (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1) +2;
    long phase;

    //   length = (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1) + 2;
    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i) {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }
    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL) {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }
    if (strcmp(tmp[0],"1")==0) {
        //     printf("The matrix is diagonal\n");
        //Read only the diagonal components of the diffusivity matrix
        l=1;
        for(i=0; i < NUMCOMPONENTS-1; i++) {
            phase = atol(tmp[1]);
            Mat[phase][i][i] = atof(tmp[1+l]);
            l++;
        }
        for(i=0; i < NUMCOMPONENTS-1; i++) {
            for (j=i+1; j < NUMCOMPONENTS-1; j++) {
                Mat[phase][i][j] = 0.0;
                Mat[phase][j][i] = 0.0;
                l++;
            }
        }
    } else {
        l=1;
        for(i=0; i < NUMCOMPONENTS-1; i++) {
            phase = atol(tmp[1]);
            Mat[phase][i][i] = atof(tmp[1+l]);
            l++;
        }
        for(i=0; i < NUMCOMPONENTS-1; i++) {
            for (j=0; j < NUMCOMPONENTS-1; j++) {
                Mat[phase][i][j] = atof(tmp[1+l]);
                l++;
            }
        }
    }
    for (i = 0; i < len; ++i) {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void populate_A_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    long i,j,k,l;
    long len = (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1) +1;
    long phase;

    //   length = (NUMCOMPONENTS-1)*(NUMCOMPONENTS-1) + 2;
    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i) {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }
    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL) {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }
    l=1;
    for(i=0; i < NUMCOMPONENTS-1; i++) {
        phase = atol(tmp[0]);
        Mat[phase][i][i] = atof(tmp[l]);
        l++;
    }
    for(i=0; i < NUMCOMPONENTS-1; i++) {
        for (j=i+1; j < NUMCOMPONENTS-1; j++) {
            Mat[phase][i][j] = atof(tmp[l]);
            Mat[phase][j][i] = Mat[phase][i][j];
            l++;
        }
    }
    for (i = 0; i < len; ++i) {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}


void populate_string_array(char **string, char *tmpstr, long size)
{
    char *str1, *token;
    char *saveptr1;

    long i;

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(string[i],token);

        if (string[i][0] == ' ')
            memmove(string[i], string[i]+1, strlen(string[i]));

        if (string[i][strlen(string[i])-1] == ' ')
            string[i][strlen(string[i])-1] = '\0';
    }
}

void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    long i,j,k,l;
    long len = (NUMCOMPONENTS-1) + 2;
    long phase1, phase2;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i)
    {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }
    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }
    phase1 = atoi(tmp[0]);
    phase2 = atoi(tmp[1]);

    l=1;
    for (i=0; i < NUMCOMPONENTS-1; i++)
    {
        Mat[phase1][phase2][i] = atof(tmp[l+1]);
        l++;
    }

    for (i = 0; i < len; ++i)
    {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void get_Rotation_Matrix(double **R, double theta, int axis) {
    double costheta, sintheta;

    costheta = cos(theta*M_PI/180.0);
    sintheta = sin(theta*M_PI/180.0);

    if (axis == 0) {
        R[0][0] = 1.0;
        R[0][1] = 0.0;
        R[0][2] = 0.0;
        R[1][0] = 0.0;
        R[2][0] = 0.0;

        R[1][1] = costheta;
        R[1][2] = -sintheta;
        R[2][1] = sintheta;
        R[2][2] = costheta;
    }

    if (axis == 1) {
        R[0][0] = costheta;
        R[0][1] = 0.0;
        R[0][2] = sintheta;
        R[1][0] = 0.0;
        R[2][0] = -sintheta;

        R[1][1] = 1.0;
        R[1][2] = 0.0;
        R[2][1] = 0.0;
        R[2][2] = costheta;
    }

    if (axis == 2) {
        R[0][0] = costheta;
        R[0][1] = -sintheta;
        R[0][2] = 0.0;
        R[1][0] = sintheta;
        R[2][0] = 0.0;

        R[1][1] = costheta;
        R[1][2] = 0.0;
        R[2][1] = 0.0;
        R[2][2] = 1.0;
    }
}

void populate_rotation_matrix(double ****Mat, double ****Mat_Inv, char *tmpstr)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;
    int i,l,j;

    long len = 5;
    long phase1;
    long phase2;
    double thetax, thetay, thetaz;
    double **Rx;
    double **Ry;
    double **Rz;
    double **mult;
 //   double **LUTemp;

    Rx     = MallocM(3,3);
    Ry     = MallocM(3,3);
    Rz     = MallocM(3,3);
    mult   = MallocM(3,3);
  //  LUTemp = MallocM(3,3);

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i) {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL) {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    phase1 = atol(tmp[0]);
    phase2 = atol(tmp[1]);

    thetax = atof(tmp[2]);
    thetay = atof(tmp[3]);
    thetaz = atof(tmp[4]);

    //   printf("phase1=%ld, phase2=%ld\n", phase1, phase2);

    get_Rotation_Matrix(Rx, thetax, 0);
    get_Rotation_Matrix(Ry, thetay, 1);
    get_Rotation_Matrix(Rz, thetaz, 2);


    matrixMultiply(Rx, Ry, mult, 3);

    matrixMultiply(mult, Rz, Mat[phase1][phase2], 3);

    int P[4];

    // for (i=0; i<3; i++) {
    //     for (j=0; j<3; j++) {
    //         LUTemp[i][j] = Mat[phase1][phase2][i][j];
    //     }
    // }

    // LUPDecompose(LUTemp, 3, 0.000001, P);
    // LUPInvert(LUTemp, P, 3, Mat_Inv[phase1][phase2]);
    matinvnew(Mat[phase1][phase2], Mat_Inv[phase1][phase2], 3);

    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            Mat[phase2][phase1][i][j]     = Mat[phase1][phase2][i][j];
            Mat_Inv[phase2][phase1][i][j] = Mat_Inv[phase1][phase2][i][j];
                   //printf("%le ",Mat[phase1][phase2][i][j]);
        }
            // printf("\n");
    }

    //   exit(0);

    FreeM(Rx, 3);
    FreeM(Ry, 3);
    FreeM(Rz, 3);
    FreeM(mult, 3);
  //  FreeM(LUTemp, 3);

    for (i = 0; i < len; ++i) {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void populate_symmetric_tensor(symmetric_tensor *Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    int i, j, k, l;
    long len = 7;
    long phase;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i)
    {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }
    phase = atoi(tmp[0]);

    Mat[phase].xx = atof(tmp[1]);
    Mat[phase].yy = atof(tmp[2]);
    Mat[phase].zz = atof(tmp[3]);
    Mat[phase].yz = atof(tmp[4]);
    Mat[phase].xz = atof(tmp[5]);
    Mat[phase].xy = atof(tmp[6]);


    for (i = 0; i < len; ++i)
    {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

void populate_cubic_stiffness(Stiffness_cubic *Mat, char *tmpstr)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;
    int i, l;

    long len = 4;
    long phase;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; ++i)
    {
        tmp[i] = (char*)malloc(sizeof(char)*10);
    }

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    phase = atoi(tmp[0]);

    Mat[phase].C11 = atof(tmp[1]);
    Mat[phase].C12 = atof(tmp[2]);
    Mat[phase].C44 = atof(tmp[3]);

    for (i = 0; i < len; ++i)
    {
        free(tmp[i]);
    }
    free(tmp);
    tmp = NULL;
}

double* MallocV(long m)
{
    double *Mat;
    Mat = (double*)malloc(m*sizeof(*Mat));
    return Mat;
}

double** MallocM(long a, long b)
{
    long i, j;
    double** Mat;

    Mat = (double**)malloc(a*sizeof(**Mat));

    for (i = 0; i < a; i++)
    {
        Mat[i] = (double*)malloc(b*sizeof(*Mat[i]));
        for (j = 0; j < b; j++)
            Mat[i][j] = 0.0;
    }

    return Mat;
}

double*** Malloc3M(long a, long b, long c)
{
    long i, j, k;
    double*** Mat;

    Mat = (double***)malloc(a*sizeof(***Mat));

    for (i = 0; i < a; i++)
    {
        Mat[i] = (double**)malloc(b*sizeof(**Mat[i]));
        for (j = 0; j < b; j++)
        {
            Mat[i][j] = (double *)malloc(c*sizeof(*Mat[i][j]));
            for (k = 0; k < c; k++)
                Mat[i][j][k] = 0.0;
        }
    }

    return Mat;
}

double**** Malloc4M(long a, long b, long c, long d)
{
    long i, j, p;
    double**** Mat;

    Mat = (double****)malloc(a*sizeof(****Mat));
    for (i = 0; i < a; i++)
    {
        Mat[i] = (double***)malloc(b*sizeof(***Mat[i]));
        for (j = 0; j < b; j++)
        {
            Mat[i][j] = (double**)malloc(c*sizeof(**Mat[i][j]));
            for (p = 0; p < c; p++)
            {
                Mat[i][j][p] = (double*)malloc(d*sizeof(*Mat[i][j][p]));
                for (int q = 0; q < d; q++)
                    Mat[i][j][p][q] = 0.0;
            }
        }
    }

    return Mat;
}

void FreeM(double **Mat, long a)
{
    long i;

    for (i = 0; i < a; i++)
        free(Mat[i]);

    free(Mat);
    Mat = NULL;
}

void Free3M(double ***Mat, long a, long b)
{
    long i, j;

    for (i = 0; i < a; i++)
    {
        for (j = 0; j < b; j++)
            free(Mat[i][j]);
        free(Mat[i]);
    }

    free(Mat);
    Mat = NULL;
}

void Free4M(double ****Mat, long a, long b, long c)
{
    long i, j, l;

    for(i = 0; i < a; i++)
    {
        for(j = 0; j < b; j++)
        {
            for(l = 0; l < c; l++)
            {
                free(Mat[i][j][l]);
            }
            free(Mat[i][j]);
        }
        free(Mat[i]);
    }

    free(Mat);
    Mat = NULL;
}

void allocOnDev(double **arr, double ***arr2, long N, long stride)
{
    cudaMalloc((void**)arr, sizeof(double)*N*stride);

    *arr2 = (double**)malloc(sizeof(double*) * N);
    for (long i = 0; i < N; i++)
        *arr2[i] = *arr + stride;
}

void freeOnDev(double **arr, double ***arr2)
{
    cudaFree(*arr);
    free(*arr2);
}

void freeVars(domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    Free3M(simParams->slopes, simDomain->numPhases, simDomain->numPhases);
    FreeM(simParams->DELTA_T, simDomain->numPhases);

    FreeM(simParams->gamma_host, simDomain->numPhases);
    cudaFree(simParams->gamma_dev);

    FreeM(simParams->kappaPhi_host, simDomain->numPhases);
    cudaFree(simParams->kappaPhi_dev);

    FreeM(simParams->relax_coeff_host, simDomain->numPhases);
    cudaFree(simParams->relax_coeff_dev);

    FreeM(simParams->Tau_host, simDomain->numPhases);

    Free3M(simParams->diffusivity_host, simDomain->numPhases, simDomain->numComponents-1);
    cudaFree(simParams->diffusivity_dev);

    Free3M(simParams->mobility_host, simDomain->numPhases, simDomain->numComponents-1);
    cudaFree(simParams->mobility_dev);

    Free3M(simParams->F0_A_host, simDomain->numPhases, simDomain->numComponents-1);
    cudaFree(simParams->F0_A_dev);

    FreeM(simParams->F0_B_host, simDomain->numPhases);
    cudaFree(simParams->F0_B_dev);

    free(simParams->F0_C_host);
    cudaFree(simParams->F0_C_dev);

    Free3M(simParams->ceq_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->ceq_dev);

    Free3M(simParams->cfill_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->cfill_dev);

    Free3M(simParams->cguess_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->cguess_dev);

    Free3M(simParams->theta_ijk_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->theta_ijk_dev);

    FreeM(simParams->theta_ij_host, simDomain->numPhases);
    cudaFree(simParams->theta_ij_dev);

    free(simParams->theta_i_host);
    cudaFree(simParams->theta_i_dev);

    Free4M(simParams->Rotation_matrix_host, simDomain->numPhases, simDomain->numPhases, 3);
    cudaFree(simParams->Rotation_matrix_dev);

    Free4M(simParams->Inv_Rotation_matrix_host, simDomain->numPhases, simDomain->numPhases, 3);
    cudaFree(simParams->Inv_Rotation_matrix_dev);

    FreeM(simParams->dab_host, simDomain->numPhases);
    cudaFree(simParams->dab_dev);

    free(simParams->eigen_strain);
    free(simParams->Stiffness_c);

    for (long i = 0; i < simDomain->numThermoPhases; i++)
        free(simDomain->phases_tdb[i]);

    free(simDomain->phases_tdb);

    for (long i = 0; i < simDomain->numPhases; i++)
    {
        free(simDomain->phaseNames[i]);
        free(simDomain->phase_map[i]);
    }

    for (long i = 0; i < simDomain->numComponents; i++)
        free(simDomain->componentNames[i]);

    free(simDomain->phaseNames);
    free(simDomain->componentNames);
    free(simDomain->phase_map);
    free(simDomain->thermo_phase_host);
    cudaFree(simDomain->thermo_phase_dev);

    free(simControls->eigenSwitch);

    for (long i = 0; i < 6; i++)
    {
        free(simControls->boundary[i]);
    }
}
