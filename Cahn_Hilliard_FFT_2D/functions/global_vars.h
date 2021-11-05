#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_
#include<complex.h>
#include <gsl/gsl_rng.h>
#include "fftw3.h"
//Global variables
long MESH_X; 
long MESH_Y;
long MESH_Z;
long layer_size;
int DIMENSION;
double deltax;
double deltay;
double deltaz;
double deltat;

int NUMPHASES;
int NUMCOMPONENTS;
int t;

long ntimesteps;
long saveT;
long nsmooth;

double R; 
double V;

char **Components;
char **Phases;

// double ceq[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***ceq;
// double cfill[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***cfill;

// double Diffusivity[NUMPHASES][NUMCOMPONENTS-1][NUMCOMPONENTS-1];
//double ***Diffusivity;

// double Gamma[NUMPHASES][NUMPHASES];
double **Gamma;

long rows_x, rows_y, rows_z;
long *start, *end;
long *averow, *rows, *offset, *extra;


int t, to, n, starttime;
long a,b;
int ASCII=0;
long time_output;
int max_length;

struct max_min {
  double *phi_max;
  double *phi_min;
  double *com_max;
  double *com_min;
  double *rel_change_phi;
  double *rel_change_com;
};

struct max_min global_max_min;

struct fill_cube {
  long x_start;
  long x_end;
  long y_start;
  long y_end;
  long z_start; 
  long z_end;
};
struct fill_cube fill_cube_parameters;

struct fill_cylinder {
  long x_center;
  long y_center;
  long z_start;
  long z_end;
  double radius;
};
struct fill_cylinder fill_cylinder_parameters;

struct fill_ellipse {
  long x_center;
  long y_center;
  long z_center;
  double major_axis;
  double eccentricity;
  double rot_angle;
};

struct fill_ellipse fill_ellipse_parameters;

struct fill_sphere {
  long x_center;
  long y_center;
  long z_center;
  double radius;
};

struct fill_sphere fill_sphere_parameters;
  
struct filling_type {
  long NUMCUBES;
  long NUMCIRCLES;
  long NUMTRIANGLES;
  double volume_fraction;
  long length;
};
struct filling_type *filling_type_phase;

struct fields {
  double *phia;
  double *compi;
};

struct fields *gridinfo;

fftw_complex **phi; 
fftw_complex **dfdphi;
fftw_complex **com;
fftw_complex **dfdc;

fftw_plan planF;
fftw_plan planB;

long i, j, k;
long x, y, z;
long half_MESH_X;
long half_MESH_Y;
long index_count;
double sum1;
double tmp1, tmp2;
double k2;
double inv_denom;
double delta_kx;
double delta_ky;
double dWdphi;
double *W;
double *kx;
double *ky;
double *P;

double **Kappa_phi, **Kappa_c, **L_phi;
double ***AtomicMobility;
double *A_fm, *B_fp;

#define X 0
#define Y 1
#define Z 2
#define TRUE 1


#endif
