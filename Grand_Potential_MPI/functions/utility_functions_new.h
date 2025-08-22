#ifndef MICROSIM_UTILITY_FUNCTIONS_HEADER
#define MICROSIM_UTILITY_FUNCTIONS_HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <ctype.h>
#include <assert.h>
#include <strings.h>

#include <mpi.h>

#include "global_vars.h"

#define msSelf(var) #var
#define msMalloc(T, arr_size) ((T*)malloc(sizeof(T)*(arr_size)))
#define msCalloc(T, arr_size) ((T*)calloc((arr_size),sizeof(T)))

void        ms_Abort();
void        ms_Assert_Terminate(int err);

double*     ms_Malloc_1D(size_t m);
double**    ms_Malloc_2D(size_t m, size_t n);
double***   ms_Malloc_3D(size_t m, size_t n, size_t k);
double****  ms_Malloc_4D(size_t m, size_t n, size_t k, size_t l);

void        ms_Free_1D(double*     ptr);
void        ms_Free_2D(double**    ptr);
void        ms_Free_3D(double***   ptr);
void        ms_Free_4D(double****  ptr);

FILE*   ms_FILE_Open(const char * path, const char * mode);
bool    ms_FILE_IsRegular(const char * path);

long    ms_readvar_l(char *file_path, char * var_name);

size_t  ms_String_CountChar         (const char * str, const char c);
void    ms_String_RemoveWhitespaces (char * Line);

void    ms_Regex_ExtractMatch(const char * Line, regmatch_t match, char * dest, size_t dest_size);

void    ms_Throw_SizeExpectationError(size_t expected_size, const char *Name, char **Tokens, size_t NumTok);

char**  ms_ArrayString_Create   (size_t arr_size, size_t string_size);
void    ms_ArrayString_Free     (char** tmp, size_t arr_size);
void    ms_ArrayString_Populate (char **arr, char *values, size_t arr_size);
void    ms_ArrayString_Print    (FILE * OUT, char * NAME, char** tmp, size_t arr_size);

void    ms_PopulateMatrix_Rotation_abqr     (double**** Mat, double**** Inv, char **tokens);
void    ms_PopulateMatrix_Symmetric_ab      (double**   Mat, size_t arr_size, char **tokens);
void    ms_PopulateMatrix_Symmetric_abc     (double***  Mat, size_t arr_size, char **tokens);
void    ms_PopulateMatrix_Diffusivity_aij   (double***  Mat, char **tokens);
void    ms_PopulateMatrix_Thermodynamic_abi (double***  Mat, char **tokens);
void    ms_PopulateMatrix_A_aij             (double***  Mat, char **tokens);
void    ms_PopulateSymmetricTensor          (struct symmetric_tensor *Mat, char **tokens);
void    ms_PopulateCubicStiffness           (struct Stiffness_cubic *Mat, char **tokens);

void    ms_BoundaryConditions_Populate(char **tokens);
void    ms_BoundaryValues_Populate(char **tokens);


void ms_Assert_Terminate(int err)
{
  if (err)
  {
    ms_Abort();
  }
}

double *ms_Malloc_1D(size_t m)
{
    double *Mat = msMalloc(double, m);
    for (size_t i = 0; i < m; i++)
    {
        Mat[i] = 0;
    }
    return Mat;
}

double** ms_Malloc_2D(size_t m, size_t n)
{
    double **Mat = NULL;
    Mat    = msMalloc(double*, m);
    Mat[0] = ms_Malloc_1D(m*n);
    for(size_t i = 1; i < m; i++)
    {
        Mat[i] = Mat[0] + i*n;
    }   
    return Mat;
}

double*** ms_Malloc_3D(size_t m, size_t n, size_t k)
{
    double ***Mat = NULL;
    Mat       = msMalloc(double**,  m);
    Mat[0]    = msMalloc(double*,   m*n);
    Mat[0][0] = ms_Malloc_1D(m*n*k);
    
    size_t i, j;
    for(i = 0; i < m; i++)
    {
        Mat[i] = Mat[0] + i*n;

        for (j = 0; j < n; j++)
        {
            Mat[i][j] = Mat[0][0] + i*(n*k) + j*(k);
        }
    }
    return Mat;
}

double**** ms_Malloc_4D(size_t m, size_t n, size_t k, size_t l)
{
    double ****Mat = NULL;
    Mat          = msMalloc(double***,  m); 
    Mat[0]       = msMalloc(double**,   m * n);
    Mat[0][0]    = msMalloc(double*,    m * n * k);
    Mat[0][0][0] = ms_Malloc_1D(m * n * k * l);

    size_t i, j, p;
    for(i = 0; i < m; i++)
    {
        Mat[i] = Mat[0] + i * n;

        for (j = 0; j < n; j++)
        {
            Mat[i][j] = Mat[0][0] + i * (n * k) + j * k;

            for (p = 0; p < k; p++)
            {
                Mat[i][j][p] = Mat[0][0][0] + i * (n * k * l) + j * (k * l) + p * l;
            }
        }
    }
    return Mat;
}

void ms_Free_1D(double *ptr)
{
    free(ptr);
    ptr = NULL;
}

void ms_Free_2D(double **ptr)
{
    free(ptr[0]);
    free(ptr);
    ptr = NULL;
}

void ms_Free_3D(double ***ptr)
{
    free(ptr[0][0]);
    free(ptr[0]);
    free(ptr);
    ptr = NULL;
}

void ms_Free_4D(double ****ptr)
{
    free(ptr[0][0][0]);
    free(ptr[0][0]);
    free(ptr[0]);
    free(ptr);
    ptr = NULL;
}

FILE *ms_FILE_Open(const char * path, const char * mode)
{
  FILE *fr = fopen(path, mode);
  if(fr == NULL)
  {
    fprintf(stderr, "cannot open file at path %s", path);
    ms_Assert_Terminate(true);
  }
  return fr;
}

bool ms_FILE_IsRegular(const char * path)
{
  struct stat _file_stat;
  if (stat(path, &_file_stat))
  {
    perror("stat");
    fprintf(stderr, "file %s does not exist (or some internal error)", path);
    return false;
  }
  if (S_ISREG(_file_stat.st_mode))
  {
    return true;
  } else {
    fprintf(stderr, "file %s is not a regular file", path);
    return false;
  }
}

size_t ms_String_CountChar(const char * str, const char c)
{
    size_t count = 0 ;
    for (size_t i = 0; str[i] != '\0'; i++)
    {
        if (str[i] == c)
        {
            count ++;   
        }
    }
    return count;
}

void ms_String_RemoveWhitespaces(char * Line)
{
    char    *read  = Line,
            *write = Line;
    while (*read)
    {
        if (!isspace((unsigned char)*read))
        {
            *write++ = *read;
        }
        read++;
    }
    *write = '\0';
}

void ms_Regex_ExtractMatch(const char * Line, regmatch_t match, char * dest, size_t dest_size)
{
    int LineLen = match.rm_eo - match.rm_so;
    if (LineLen >= dest_size)
    {
        LineLen = dest_size -1;
    }
    strncpy(dest, Line + match.rm_so, LineLen);
    dest[LineLen] = '\0';
}

void ms_Abort()
{
    fprintf(stderr, "\nAborting Execution!\nPlease consult our documentation!\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(EXIT_FAILURE);
}

void ms_Throw_SizeExpectationError(size_t expected, const char * Name, char **Tokens, size_t numtok)
{
    if (numtok != expected)
    {
        fprintf(stderr, "\n**ERROR** Expected %ld values in array %s\n", expected, Name);
        fprintf(stderr, "You entered %ld values\n", numtok);
        for (size_t i = 0; i < numtok; i++)
        {
            fprintf(stderr, "%s[%ld] = %s\n", Name, i, Tokens[i]);
        }
        ms_Abort();
    }
}

void ms_mpi_try(int mpi_error)
{
    static char msg[MPI_MAX_ERROR_STRING];
    msg[0]  = '\0';
    int msg_len = 0;
    if (mpi_error != MPI_SUCCESS)
    {
        fprintf(stderr, "\n**ERROR** MPI Error!\n");
        MPI_Error_string(mpi_error, msg, &msg_len);
        fprintf(stderr, "%s[%d]", msg, msg_len);
        ms_Abort();
    }    
}

char **ms_ArrayString_Create(size_t arr_size, size_t string_size)
{
    char **M = msMalloc(char*, arr_size);
    for (size_t i = 0; i < arr_size; ++i)
    {
        M[i] = msMalloc(char, string_size+1);
    }
    return M;
}

void ms_ArrayString_Free(char** tmp, size_t arr_size)
{
    if (tmp)
    {
        for (size_t i = 0; i < arr_size; ++i) {
            free(tmp[i]);
        }
        free(tmp);
        tmp = NULL;        
    }
}

void ms_ArrayString_Populate(char **arr, char *values, size_t arr_size)
{
    if (arr == NULL)
    {
        arr = ms_ArrayString_Create(arr_size, 64);
    }

    size_t count = 0;
    char *token = strtok(values, ",");

    while ((token != NULL) && (count < arr_size))
    {
        strcpy(arr[count++], token);
        token = strtok(NULL, ",");
    }
}

void ms_ArrayString_Print(FILE *OUT, char *NAME, char **tmp, size_t arr_size)
{
    fprintf(OUT, "%s = {", NAME);
    for (size_t i = 0; i < arr_size; i++)
    {
        fprintf(OUT,"%s", tmp[i]);
        if (i != arr_size-1) {fprintf(OUT, ", ");}
    }
    fprintf(OUT, "};\n");
}

void ms_PopulateMatrix_Rotation_abqr(double ****Mat, double ****Inv, char **tmp)
{
    double **Rx     = ms_Malloc_2D(3,3);
    double **Ry     = ms_Malloc_2D(3,3);
    double **Rz     = ms_Malloc_2D(3,3);
    double **mult   = ms_Malloc_2D(3,3);

    long    phase1 = atol(tmp[0]),
            phase2 = atol(tmp[1]);
    
    double  thetax = atof(tmp[2]),
            thetay = atof(tmp[3]),
            thetaz = atof(tmp[4]);

    get_Rotation_Matrix(Rx, thetax, 0);
    get_Rotation_Matrix(Ry, thetay, 1);
    get_Rotation_Matrix(Rz, thetaz, 2);

    multiply2d(Rx,   Ry, mult,                3);
    multiply2d(mult, Rz, Mat[phase1][phase2], 3);

    if (USE_GSL_MATINV)
    {
        matinv_gsl_t rotinv;
        matinv_gsl_Malloc(&rotinv, 3);
        matinv_gsl_Copy_to_Mat(rotinv.MAT, Mat[phase1][phase2]);

        matinv_gsl_Invert(&rotinv);

        matinv_gsl_Copy_from_Mat(Inv[phase1][phase2], rotinv.INV);
        matinv_gsl_Free(&rotinv);
    }
    else
    {
        matinvnew(Mat[phase1][phase2], Inv[phase1][phase2], 3);
    }
    
    int i, j;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Mat[phase2][phase1][i][j] = Mat[phase1][phase2][i][j];
            Inv[phase2][phase1][i][j] = Inv[phase1][phase2][i][j];
        }
    }
}

void ms_PopulateMatrix_Symmetric_ab(double **Mat, size_t size, char **tokens)
{
    assert(Mat != NULL);
    assert(size > 1);

    size_t i, j, k = 0;
    for (i = 0; i < size; i++)
    {
        for (j = i+1; j < size; j++)
        {
            Mat[i][i] = 0.0;
            Mat[i][j] = atof(tokens[k++]);
            Mat[j][i] = Mat[i][j];
        }
    }
}

void ms_PopulateMatrix_Symmetric_abc(double ***Mat, size_t size, char **tokens)
{
    assert(Mat != NULL);
    assert(size > 2);

    size_t i, j, k, l = 0;
    for(i=0; i < size; i++)
    {
        for (j=i+1; j < size; j++)
        {
            for (k=j+1; k < size; k++)
            {
                Mat[i][i][i] = 0.0;
                Mat[i][j][j] = 0.0;
                Mat[i][k][k] = 0.0;

                Mat[i][j][k] = atof(tokens[l++]);
                
                Mat[i][k][j] = Mat[i][j][k];
                Mat[j][i][k] = Mat[i][j][k];
                Mat[j][k][i] = Mat[i][j][k];
                Mat[k][i][j] = Mat[i][j][k];
                Mat[k][j][i] = Mat[i][j][k];
            }
        }
    }
}

void ms_PopulateMatrix_Diffusivity_aij(double *** Mat, char **tokens)
{
    int  i, j, l = 0;

    int diagonal = atoi(tokens[l++]) == 1,
        phase    = atoi(tokens[l++]);

    for (i = 0; i < NUMCOMPONENTS-1; i++)
    {
        Mat[phase][i][i] = atof(tokens[l++]);
    }

    if (!diagonal)
    {
        for (i=0; i < NUMCOMPONENTS-1; i++)
        {
            for (j=0; j < NUMCOMPONENTS-1; j++)
            {
                if (i != j)
                {
                    Mat[phase][i][j] = atof(tokens[l++]);    
                }
            }
        }
    }
}

void ms_PopulateMatrix_Thermodynamic_abi(double ***Mat, char **tokens)
{
    int  i, j = 0;
    int phase1 = atoi(tokens[j++]);
    int phase2 = atoi(tokens[j++]);
    for (i=0; i < NUMCOMPONENTS-1; i++)
    {
        double value = atof(tokens[j++]);
        Mat[phase1][phase2][i] = value;
        Mat[phase2][phase1][i] = value;
    }
}

void ms_PopulateMatrix_A_aij(double *** Mat, char **tokens)
{
    int  i, j, k = 0;

    int phase = atoi(tokens[k++]);

    for(i=0; i < NUMCOMPONENTS-1; i++)
    {
        Mat[phase][i][i] = atof(tokens[k++]);
    }

    for(i=0; i < NUMCOMPONENTS-1; i++)
    {
        for (j=i+1; j < NUMCOMPONENTS-1; j++)
        {
            Mat[phase][i][j] = atof(tokens[k++]);
            Mat[phase][j][i] = Mat[phase][i][j];
        }
    }
}

void ms_PopulateSymmetricTensor(struct symmetric_tensor *Mat, char **tokens)
{
    int phase = atoi(tokens[0]);
    Mat[phase].xx = atof(tokens[1]);
    Mat[phase].yy = atof(tokens[2]);
    Mat[phase].zz = atof(tokens[3]);
    Mat[phase].yz = atof(tokens[4]);
    Mat[phase].xz = atof(tokens[5]);
    Mat[phase].xy = atof(tokens[6]);
}

void ms_PopulateCubicStiffness(struct Stiffness_cubic *Mat, char **tokens)
{
    int phase = atoi(tokens[0]);  
    Mat[phase].C11 = atof(tokens[1]);
    Mat[phase].C12 = atof(tokens[2]);
    Mat[phase].C44 = atof(tokens[3]);
}

void ms_BoundaryConditions_Populate(char ** tmp)
{
    BOUNDARY_LEFT   = atoi(tmp[2]);
    BOUNDARY_RIGHT  = atoi(tmp[1]);
    BOUNDARY_FRONT  = atoi(tmp[3]);
    BOUNDARY_BACK   = atoi(tmp[4]);
    BOUNDARY_TOP    = atoi(tmp[5]);
    BOUNDARY_BOTTOM = atoi(tmp[6]);

    if (strcmp(tmp[0], "phi") == 0)
    {
        assign_buffer_points_conditions(0, BOUNDARY_LEFT, BOUNDARY_RIGHT,
            BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM);
    }
    if (strcmp(tmp[0], "mu") == 0)
    {
        assign_buffer_points_conditions(1, BOUNDARY_LEFT, BOUNDARY_RIGHT,
            BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM);
    }
    if (strcmp(tmp[0], "T") == 0)
    {
        assign_buffer_points_conditions(2, BOUNDARY_LEFT, BOUNDARY_RIGHT,
            BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM);
    }
    if (strcmp(tmp[0], "u") == 0)
    {
        assign_buffer_points_conditions(3, BOUNDARY_LEFT, BOUNDARY_RIGHT,
            BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM);
    }
    if (strcmp(tmp[0], "F") == 0)
    {
        assign_buffer_points_conditions(4, BOUNDARY_LEFT, BOUNDARY_RIGHT,
            BOUNDARY_FRONT, BOUNDARY_BACK, BOUNDARY_TOP, BOUNDARY_BOTTOM);
    }
}

void ms_BoundaryConditions_PopulateDefault()
{
    /// phi
    assign_buffer_points_conditions(0, 3, 3, 3, 3, 3, 3);
    /// mu
    assign_buffer_points_conditions(1, 3, 3, 3, 3, 3, 3);
    /// T
    assign_buffer_points_conditions(2, 3, 3, 3, 3, 3, 3);
    /// u
    assign_buffer_points_conditions(3, 3, 3, 3, 3, 3, 3);
    /// F
    assign_buffer_points_conditions(4, 3, 3, 3, 3, 3, 3);
}

void ms_BoundaryValues_Populate(char **tmp)
{
    double val[6];

    val[0] = atof(tmp[1]);
    val[1] = atof(tmp[2]);
    val[2] = atof(tmp[3]);
    val[3] = atof(tmp[4]);
    val[4] = atof(tmp[5]);
    val[5] = atof(tmp[6]);
  
    if (strcmp(tmp[0], "phi") == 0)
    {
        assign_boundary_points_values(0, val[1], val[0], val[2], val[3], val[4], val[5]);
    }
    else if (strcmp(tmp[0], "mu") == 0)
    {
        assign_boundary_points_values(1, val[1], val[0], val[2], val[3], val[4], val[5]);
    }
    else if (strcmp(tmp[0], "T") == 0)
    {
        assign_boundary_points_values(2, val[1], val[0], val[2], val[3], val[4], val[5]);
    }
    else if (strcmp(tmp[0], "u") == 0)
    {
        assign_boundary_points_values(3, val[1], val[0], val[2], val[3], val[4], val[5]);
    }
    else if (strcmp(tmp[0], "F") == 0)
    {
        assign_boundary_points_values(4, val[1], val[0], val[2], val[3], val[4], val[5]);
    }
}

void ms_BoundaryValues_PopulateDefault()
{
    /// phi
    assign_boundary_points_values(0, 0, 0, 0, 0, 0, 0);
    /// mu
    assign_boundary_points_values(1, 0, 0, 0, 0, 0, 0);
    /// T
    assign_boundary_points_values(2, 0, 0, 0, 0, 0, 0);
    /// u
    assign_boundary_points_values(3, 0, 0, 0, 0, 0, 0);
    /// F
    assign_boundary_points_values(4, 0, 0, 0, 0, 0, 0);
}

#ifdef __cplusplus
}
#endif

#endif