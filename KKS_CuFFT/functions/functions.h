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

#endif
