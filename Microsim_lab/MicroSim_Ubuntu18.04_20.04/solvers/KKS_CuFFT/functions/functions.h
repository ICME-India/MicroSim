#ifndef FUNCTIONS_DEFS_H_
#define FUNCTIONS_DEFS_H_

double** MallocM(long m, long n);
double*** Malloc3M(long m, long n, long k);

void FreeM(double **a, long m);
void Free3M(double ***a, long m, long n);

void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_F0_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_symmetric_tensor(struct symmetric_tensor *Mat, char *tmpstr, long NUMPHASES);

void PRINT_INT(char *key, int value, FILE *fp);
void PRINT_LONG(char *key, long value, FILE *fp);
void PRINT_DOUBLE(char *key, double value, FILE *fp);
void PRINT_MATRIX(char *key, double **Mat, long m, long n, FILE *fp);
void PRINT_VECTOR(char *key, double *Mat, long m, FILE *fp);

void populate_string_array(char**string, char *tmpstr, long size);

void PRINT_STRING_ARRAY(char *key, char **str, long m, FILE *fp);

void devPropQuery();

double ran(long *idum);
int checkOverlap(long MESH_X, long MESH_Y, double centx, double centy, double rad, int *occupancy);

void writetofile_serial2D(struct fields *gridinfo, char *argv[], long t);
void writetofile_serial2D_binary(struct fields *gridinfo, char *argv[], long t);
void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo);
void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo);

void gpu_to_gridinfo(struct fields *gridinfo, cufftDoubleComplex **comp, cufftDoubleComplex **phi);
void gridinfo_to_gpu(struct fields *gridinfo, cufftDoubleComplex **comp, cufftDoubleComplex **phi);

void calc_strn(double eigen_strn[3][3], double eps11, double eps22, double eps33, double eps12, double eps13, double eps23);
void eval_lambda(double C[6][6], double lambda[3][3][3][3]);
void calc_inv_omega(double kx, double ky, double kz, double lambda[3][3][3][3], double inv_omega[3][3]);
void calc_inverse(double a[3][3], double ainv[3][3]);
void calc_stress(double eigen_strn[3][3], double lambda[3][3][3][3], double stress[3][3]);
double elast(double kx, double ky, double kz, double omega[3][3], double eigen_strn1[3][3], double eigen_strn2[3][3], double stress1[3][3], double stress2[3][3], double lambda[3][3][3][3]);
void calculate_Bn(double *B, double eigen_strn1[3][3],double eigen_strn2[3][3], double stress1[3][3], double stress2[3][3]);

void evolve_FFT(char *argv[]);
void evolve_FD(char *argv[]);

#endif
