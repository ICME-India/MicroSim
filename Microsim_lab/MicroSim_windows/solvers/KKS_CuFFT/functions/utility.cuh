#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double swap_bytes(double value)
{
    double  src_num = value;
    int64_t tmp_num = htobe64(le64toh(*(int64_t*)&src_num));
    double  dst_num = *(double*)&tmp_num;
    return  dst_num;
}

double ran(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

__global__ void checkOverlap_cuda_sph(long MESH_X, long MESH_Y, long MESH_Z, double centx, double centy, double centz,
                             double rad_shield, int *occupancy, int *overlap)
{
    long i = threadIdx.x + blockDim.x*blockIdx.x;
    long j = threadIdx.y + blockDim.y*blockIdx.y;
    long k = threadIdx.z + blockDim.z*blockIdx.z;

    if (((double)((i - centx)*(i - centx)) + (double)((j - centy)*(j - centy)) + (double)((k - centz)*(k - centz)))
        <= rad_shield)
        if (occupancy[k+j*MESH_Z+i*MESH_Y*MESH_Z] == 1)
            *overlap = 1;
    __syncthreads();
}

__global__ void checkOverlap_cuda_cyl(long MESH_X, long MESH_Y, long MESH_Z, double centx, double centy,
                             double rad_shield, int *occupancy, int *overlap)
{
    long i = threadIdx.x + blockDim.x*blockIdx.x;
    long j = threadIdx.y + blockDim.y*blockIdx.y;
    long k = threadIdx.z + blockDim.z*blockIdx.z;

    if (((double)((i - centx)*(i - centx)) + (double)((j - centy)*(j - centy)))
        <= rad_shield)
        if (occupancy[k+j*MESH_Z+i*MESH_Y*MESH_Z] == 1)
            *overlap = 1;
    __syncthreads();
}

int checkOverlap(long MESH_X, long MESH_Y, long MESH_Z, long centx, long centy, long centz, double rad, long shield_dist, int *occupancy)
{
    if (!shield_dist)
        return 0;

    for (int i = 0; i < MESH_X; i++)
        for (int j = 0; j < MESH_Y; j++)
            for (int k = 0; k < MESH_Z; k++)
                if (occupancy[k + j*MESH_Z + i*MESH_Y*MESH_Z] == 1)
                    if ((i - centx)*(i - centx) + (j - centy)*(j - centy)
                        + (k - centz)*(k - centz) <= (shield_dist*rad)*(shield_dist*rad))
                    return 1;
    return 0;
}

double** MallocM(long m, long n)
{
    int i;
    double **Mat;
    Mat = (double **)malloc(m*sizeof(**Mat));
    for (i = 0; i < m; i++)
        Mat[i] = (double *)malloc(n*sizeof(*Mat[i]));
    return Mat;
}

double*** Malloc3M(long m, long n, long k)
{
    int i, j;
    double*** Mat;
    Mat = (double ***)malloc(m*sizeof(***Mat));
    for (i = 0; i < m; i++)
    {
        Mat[i] = (double **)malloc(n*sizeof(**Mat[i]));
        for (j = 0; j < n; j++)
            Mat[i][j] = (double *)malloc(k*sizeof(*Mat[i][j]));
    }
    return Mat;
}

void FreeM(double **Mat, long m)
{
    int i;
    for(i = 0; i < m; i++)
        free(Mat[i]);
    free(Mat);
    Mat = NULL;
}

void Free3M(double ***Mat, long m, long n)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            free(Mat[i][j]);
        free(Mat[i]);
    }
    free(Mat);
    Mat = NULL;
}


/*
 * Populate functions
 */

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES)
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

    for (i = 0; i < NUMPHASES; i++)
    {
        for (j = i+1; j < NUMPHASES; j++)
        {
            Mat[i][i] = 0.0;
            Mat[i][j] = atof(tmp[k]);
            Mat[j][i] = Mat[i][j];
            k++;
        }
    }

    for (i = 0; i < NUMPHASES*(NUMPHASES-1)*0.5; i++)
        free(tmp[i]);

    free(tmp);
    tmp = NULL;
}

void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1, *token;
    char *saveptr1;

    int i, j, l;
    int len = (NUMCOMPONENTS - 1) * (NUMCOMPONENTS - 1) + 2;
    int phase;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; i++)
        tmp[i] = (char*)malloc(sizeof(char)*10);

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    if (strcmp(tmp[0], "1") == 0)
    {
        l = 1;
        for(i = 0; i < NUMCOMPONENTS - 1; i++)
        {
            phase = atol(tmp[1]);
            Mat[phase][i][i] = atof(tmp[1+l]);
            l++;
        }

        for(i = 0; i < NUMCOMPONENTS - 1; i++)
        {
            for (j = i + 1; j < NUMCOMPONENTS - 1; j++)
            {
                Mat[phase][i][j] = 0.0;
                Mat[phase][j][i] = 0.0;
                l++;
            }
        }
    }
    else
    {
        l = 1;
        for(i = 0; i < NUMCOMPONENTS - 1; i++)
        {
            phase = atol(tmp[1]);
            Mat[phase][i][i] = atof(tmp[1+l]);
            l++;
        }
        for(i = 0; i < NUMCOMPONENTS - 1; i++)
        {
            for (j = 0; j < NUMCOMPONENTS - 1; j++)
            {
                Mat[phase][i][j] = atof(tmp[1+l]);
                l++;
            }
        }
    }
    for (i = 0; i < len; i++)
        free(tmp[i]);
    free(tmp);
    tmp = NULL;
}

void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1,  *token;
    char *saveptr1;

    int i, l;
    long len = (NUMCOMPONENTS - 1) + 2;
    long phase1, phase2;

    tmp = (char**)malloc(sizeof(char*)*len);
    for (i = 0; i < len; i++)
        tmp[i] = (char*)malloc(sizeof(char)*10);

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    phase1 = atoi(tmp[0]);
    phase2 = atoi(tmp[1]);

    l = 1;
    for (i = 0; i < NUMCOMPONENTS - 1; i++)
    {
        Mat[phase1][phase2][i] = atof(tmp[l + 1]);
        l++;
    }

    for (i = 0; i < len; i++)
        free(tmp[i]);
    free(tmp);
    tmp = NULL;
}

void populate_A_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS)
{
    char **tmp;
    char *str1, *token;
    char *saveptr1;

    int i, j, l;
    long len = (NUMCOMPONENTS - 1) * (NUMCOMPONENTS - 1) + 1;
    long phase;

    tmp = (char**)malloc(sizeof(char*)*len);

    for (i = 0; i < len; ++i)
        tmp[i] = (char*)malloc(sizeof(char)*10);

    for (i = 0, str1 = tmpstr; ; i++, str1 = NULL)
    {
        token = strtok_r(str1, "{,}", &saveptr1);
        if (token == NULL)
            break;
        strcpy(tmp[i],token);
    }

    l = 1;

    for(i = 0; i < NUMCOMPONENTS - 1; i++)
    {
        phase = atol(tmp[0]);
        Mat[phase][i][i] = atof(tmp[l]);
        l++;
    }

    for(i = 0; i < NUMCOMPONENTS - 1; i++)
    {
        for (j = i + 1; j < NUMCOMPONENTS - 1; j++)
        {
            Mat[phase][i][j] = atof(tmp[l]);
            Mat[phase][j][i] = Mat[phase][i][j];
            l++;
        }
    }

    for (i = 0; i < len; ++i)
        free(tmp[i]);

    free(tmp);
    tmp = NULL;
}

void populate_cubic_stiffness(struct Stiffness_cubic *Mat, char *tmpstr)
{
  char **tmp;
  char *str1, *token;
  char *saveptr1;
  int i;

  long len = 4;
  long phase;

  tmp = (char**)malloc(sizeof(char*)*len);
  for (i = 0; i < len; ++i)
    tmp[i] = (char*)malloc(sizeof(char)*10);

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
    free(tmp[i]);

  free(tmp);
  tmp = NULL;
}

void populate_string_array(char**string, char *tmpstr, long size)
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
    }
}

void populate_symmetric_tensor(struct symmetric_tensor *Mat, char *tmpstr, long NUMPHASES)
{
  char **tmp;
  char *str1, *token;
  char *saveptr1;

  int i;
  long len = 7;
  long phase;

  tmp = (char**)malloc(sizeof(char*)*len);
  for (i = 0; i < len; ++i)
    tmp[i] = (char*)malloc(sizeof(char)*10);

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
    free(tmp[i]);

  free(tmp);
  tmp = NULL;
}

/*
 * Print functions
 */

void PRINT_INT(char *key, int value, FILE *fp)
{
    fprintf(fp,"%s = %d\n", key, value);
    fprintf(fp,"\n");
}

void PRINT_DOUBLE(char *key, double value, FILE *fp)
{
    fprintf(fp,"%s = %le\n", key, value);
    fprintf(fp,"\n");
}

void PRINT_LONG(char *key, long value, FILE *fp)
{
    fprintf(fp, "%s = %ld\n", key, value);
    fprintf(fp,"\n");
}

void PRINT_VECTOR(char *key, double *Mat, long m, FILE *fp)
{
    long i;
    fprintf(fp, "%s = [", key);
    for (i = 0; i < m; i++)
        if (i < (m - 1))
            fprintf(fp, "%le,", Mat[i]);
        else
            fprintf(fp, "%le]\n", Mat[i]);
        fprintf(fp,"\n");
}

void PRINT_MATRIX(char *key, double **Mat, long m, long n, FILE *fp)
{
    int i, j;
    fprintf(fp, "%s = [", key);
    strcat(key, "= [ ");
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            if (j < (n-1))
                fprintf(fp, "%le,", Mat[i][j]);
            else
                fprintf(fp, "%le", Mat[i][j]);

            if (i < (m-1))
                fprintf(fp, "\n%*s", (int)strlen(key), "");
            else
                fprintf(fp, "]\n");
    }
    fprintf(fp, "\n");
}

void PRINT_STRING_ARRAY(char *key, char **str, long m, FILE *fp)
{
    long i;
    fprintf(fp, "%s = [", key);

    for (i = 0; i < m; i++)
        if (i < (m-1))
            fprintf(fp, "%s,", str[i]);
        else
            fprintf(fp, "%s]\n", str[i]);

        fprintf(fp,"\n");
}

void PRINT_SYMMETRIC_TENSOR(char *key, struct symmetric_tensor e, FILE *fp)
{
    char tmp[100];
    sprintf(tmp, "%s.xx",key);
    fprintf(fp, "%s = %le\n", tmp, e.xx);
    sprintf(tmp, "%s.yy",key);
    fprintf(fp, "%s = %le\n", tmp, e.yy);
    sprintf(tmp, "%s.zz",key);
    fprintf(fp, "%s = %le\n", tmp, e.zz);
    sprintf(tmp, "%s.yz",key);
    fprintf(fp, "%s = %le\n", tmp, e.yz);
    sprintf(tmp, "%s.xz",key);
    fprintf(fp, "%s = %le\n", tmp, e.xz);
    sprintf(tmp, "%s.xy",key);
    fprintf(fp, "%s = %le\n", tmp, e.xy);
    fprintf(fp, "\n");
}

void PRINT_VOIGT_CUBIC(char* key, struct Stiffness_cubic Stiffness, FILE *fp) {
  char tmp[100];
  sprintf(tmp, "%s.C11",key);
  fprintf(fp, "%s = %le\n", tmp, Stiffness.C11);
  sprintf(tmp, "%s.C12",key);
  fprintf(fp, "%s = %le\n", tmp, Stiffness.C12);
  sprintf(tmp, "%s.C44",key);
  fprintf(fp, "%s = %le\n", tmp, Stiffness.C44);
  fprintf(fp,"\n");
}

void gridinfo_to_gpu(struct fields *gridinfo, cufftDoubleComplex **comp, cufftDoubleComplex **phi)
{
    long a, b;
    long x, y, z;
    long index;

    cufftDoubleComplex *hostWriteTemp = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z);

    for (a = 0; a < NUMCOMPONENTS-1; a++)
    {
        for (x = 0; x < rows_x; x++)
        {
            for (y = 0; y < rows_y; y++)
            {
                for (z = 0; z < rows_z; z++)
                {
                    index = x*layer_size + y*rows_z + z;
                    hostWriteTemp[index].x = gridinfo[index].compi[a];
                    hostWriteTemp[index].y = 0.0;
                }
            }
        }
        cudaMemcpy(comp[a], hostWriteTemp, sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyHostToDevice);
    }

    for (b = 1; b < NUMPHASES; b++)
    {
        for (x = 0; x < rows_x; x++)
        {
            for (y = 0; y < rows_y; y++)
            {
                for (z = 0; z < rows_z; z++)
                {
                    index = x*layer_size + y*rows_z + z;
                    hostWriteTemp[index].x = gridinfo[index].phia[b];
                    hostWriteTemp[index].y = 0.0;
                }
            }
        }
        cudaMemcpy(phi[b], hostWriteTemp, sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyHostToDevice);
        //cudaMemcpy(dfdphiHost[b], hostWriteTemp, sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyHostToDevice);
    }
    free(hostWriteTemp);
}

void gpu_to_gridinfo(struct fields *gridinfo, cufftDoubleComplex **comp, cufftDoubleComplex **phi)
{
    long a, b;
    long x, y, z;
    long index;

    cufftDoubleComplex *hostPrintTemp = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z);

    for (x = 0; x < rows_x; x++)
    {
        for (y = 0; y < rows_y; y++)
        {
            for (z = 0; z < rows_z; z++)
            {
                index = x*layer_size + y*rows_z + z;
                gridinfo[index].phia[0] = 1.0;
            }
        }
    }

    for (b = 1; b < NUMPHASES; b++)
    {
        cudaMemcpy(hostPrintTemp, phi[b], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyDeviceToHost);
        for (x = 0; x < rows_x; x++)
        {
            for (y = 0; y < rows_y; y++)
            {
                for (z = 0; z < rows_z; z++)
                {
                    index = x*layer_size + y*rows_z + z;
                    gridinfo[index].phia[b] = hostPrintTemp[index].x;
                    gridinfo[index].phia[0] -= gridinfo[index].phia[b];
                }
            }
        }
    }

    for (a = 0; a < NUMCOMPONENTS-1; a++)
    {
        cudaMemcpy(hostPrintTemp, comp[a], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyDeviceToHost);
        for (x = 0; x < rows_x; x++)
        {
            for (y = 0; y < rows_y; y++)
            {
                for (z = 0; z < rows_z; z++)
                {
                    index = x*layer_size + y*rows_z + z;
                    gridinfo[index].compi[a] = hostPrintTemp[index].x;
                }
            }
        }
    }
    free(hostPrintTemp);
}
