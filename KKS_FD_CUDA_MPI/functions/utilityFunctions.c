#include "utilityFunctions.h"

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *str2, *token;
    char *saveptr1, *saveptr2;

    int i, j, k;

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

  int i,j,k,l;
  long len = NUMPHASES*(NUMPHASES-1)*(NUMPHASES-2)/6;

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

void populate_thetaij_matrix(double **Mat, char *tmpstr, long NUMPHASES)
{
    char **tmp;
    char *str1, *token;
    char *saveptr1;

    int i, j, k;

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

    int i;
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

    int i,j,k,l;
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

    int i,j,k,l;
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

    int i;

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

    int i,j,k,l;
    int len = (NUMCOMPONENTS-1) + 2;
    int phase1, phase2;

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

double** malloc2M(int a, int b)
{
    int i;
    double** Mat;

    Mat = (double**)malloc(a*sizeof(**Mat));

    for (i = 0; i < a; i++)
        Mat[i] = (double*)malloc(b*sizeof(*Mat[i]));

    return Mat;
}

double*** malloc3M(int a, int b, int c)
{
    int i, j;
    double*** Mat;

    Mat = (double***)malloc(a*sizeof(***Mat));

    for (i = 0; i < a; i++)
    {
        Mat[i] = (double**)malloc(b*sizeof(**Mat[i]));
        for (j = 0; j < b; j++)
            Mat[i][j] = (double *)malloc(c*sizeof(*Mat[i][j]));
    }

    return Mat;
}

double**** malloc4M(int a, int b, int c, int d)
{
    int i, j, p;
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
            }
        }
    }

    return Mat;
}

void free2M(double **Mat, int a)
{
    int i;

    for (i = 0; i < a; i++)
        free(Mat[i]);

    free(Mat);
    Mat = NULL;
}

void free3M(double ***Mat, int a, int b)
{
    int i, j;

    for (i = 0; i < a; i++)
    {
        for (j = 0; j < b; j++)
            free(Mat[i][j]);
        free(Mat[i]);
    }

    free(Mat);
    Mat = NULL;
}

void free4M(double ****Mat, int a, int b, int c)
{
    int i, j, l;

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

void allocOnDev(double **arr, double ***arr2, int N, int stride)
{
    cudaMalloc((void**)arr, sizeof(double)*N*stride);

    *arr2 = (double**)malloc(sizeof(double*) * N);
    for (int i = 0; i < N; i++)
        *arr2[i] = *arr + stride;
}

void freeOnDev(double **arr, double ***arr2)
{
    cudaFree(*arr);
    free(*arr2);
}

void freeVars(domainInfo *simDomain, simParameters *simParams)
{
    free2M(simParams->gamma_host, simDomain->numPhases);
    cudaFree(simParams->gamma_dev);

    free2M(simParams->kappaPhi_host, simDomain->numPhases);
    cudaFree(simParams->kappaPhi_dev);

    free2M(simParams->relax_coeff_host, simDomain->numPhases);
    cudaFree(simParams->relax_coeff_dev);

    free2M(simParams->Tau_host, simDomain->numPhases);

    free3M(simParams->diffusivity_host, simDomain->numPhases, simDomain->numComponents-1);
    cudaFree(simParams->diffusivity_dev);

    free3M(simParams->F0_A_host, simDomain->numPhases, simDomain->numComponents-1);
    cudaFree(simParams->F0_A_dev);

    free2M(simParams->F0_B_host, simDomain->numPhases);
    cudaFree(simParams->F0_B_dev);

    free(simParams->F0_C_host);
    cudaFree(simParams->F0_C_dev);

    free3M(simParams->ceq_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->ceq_dev);

    free3M(simParams->cfill_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->cfill_dev);

    free3M(simParams->cguess_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->cguess_dev);

    free3M(simParams->theta_ijk_host, simDomain->numPhases, simDomain->numPhases);
    cudaFree(simParams->theta_ijk_dev);

    free2M(simParams->theta_ij_host, simDomain->numPhases);
    cudaFree(simParams->theta_ij_dev);

    free(simParams->theta_i_host);
    cudaFree(simParams->theta_i_dev);

    for (int i = 0; i < simDomain->numThermoPhases; i++)
        free(simDomain->phases_tdb[i]);

    free(simDomain->phases_tdb);

    for (int i = 0; i < simDomain->numPhases; i++)
    {
        free(simDomain->phaseNames[i]);
        free(simDomain->phase_map[i]);
    }

    for (int i = 0; i < simDomain->numComponents; i++)
        free(simDomain->componentNames[i]);

    free(simDomain->phaseNames);
    free(simDomain->componentNames);
    free(simDomain->phase_map);
    free(simDomain->thermo_phase_host);
    cudaFree(simDomain->thermo_phase_dev);
}
