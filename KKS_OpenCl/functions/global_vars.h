#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_

//Global variables
long MESH_X; 
long MESH_Y;
long MESH_Z;
long layer_size;
int DIMENSION;
double deltax=1;
double deltay=1;
double deltaz=1;
double deltat=1;

int NUMPHASES;
int NUMCOMPONENTS;
int t;


long ntimesteps;
long saveT;
long nsmooth;

double Teq;
double Tfill;
double T;
int TEMPGRADY=0;
int ISOTHERMAL=1;
double TEMPERATURE_SCALE;

double epsilon;
double tau;
double R; 
double V;

double tilt_angle;
double Rtheta;
int FOLD;
int ANISOTROPY = 0;
int ANISOTROPY_GRADIENT = 0;
double ****Rotation_matrix;
double *Rotated_qab;
double ****Inv_Rotation_matrix;
double Rth_phase;


int FUNCTION_W = 1;
int FUNCTION_ANISOTROPY = 1;

int OBSTACLE = 0;
int WELL = 0;
int SHIFT =0;
int shiftj=0;
int BINARY=0;
int TERNARY=0;
int DILUTE=0;
int WRITECOMPOSITION=0;
int NOISE_PHASEFIELD=0;
double AMP_NOISE_PHASE=0.0;
char **Components;
char **Phases;

double t_smooth;
// double C[NUMPHASES];
double *C;
// double ceq[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***ceq;
// double cfill[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***cfill;


// double c_old[NUMCOMPONENTS-1],c_new[NUMCOMPONENTS-1],c[NUMCOMPONENTS-1];
double *c_old, *c_new, *c;
// double Diffusivity[NUMPHASES][NUMCOMPONENTS-1][NUMCOMPONENTS-1];
double ***Diffusivity;
long bulk_phase;
// double Gamma[NUMPHASES][NUMPHASES];
double **Gamma;
// double dab[NUMPHASES][NUMPHASES];
double **dab;

int BOUNDARY_LEFT;
int BOUNDARY_RIGHT;
int BOUNDARY_FRONT;
int BOUNDARY_BACK;
int BOUNDARY_TOP;
int BOUNDARY_BOTTOM;


struct max_min {
  double *phi_max;
  double *phi_min;
  double *mu_max;
  double *mu_min;
  double *rel_change_phi;
  double *rel_change_mu;
};

struct max_min global_max_min;


struct bc_scalars {
  int type;
  long points[3];
  long proxy[3];
  double value[3];
};

struct bc_scalars *boundary[6];


long rows_x, rows_y, rows_z;
long *start, *end;
long *averow, *rows, *offset, *extra;
long MAX_INTERFACE_POS=0;


int t, to, n, starttime;
long a,b;
double sum_lambdaphi, sum_dhphi;
long active_phases, count_phases;
double Deltaphi;
long interface;
double scalprod;
//long k;
double mu;
int ASCII=0;
long time_output;
int max_length;

struct Tempgrad {
 double base_temp;
 double DeltaT;
 double Distance;
 double gradient_OFFSET;
 double velocity;
 double GRADIENT;
};
struct Tempgrad temperature_gradientY;



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


long shift_OFFSET=0;
long shift_OFFSET_GLOBAL=0;
  
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
  double temperature;
};

struct fields *gridinfo;

long tNoiseStart;
double TLiquidus;
int DirectionalSolidification=0;
double TemperatureGradient=1e6;
double PullingVelocity=0.05;
int PositionOffset=20;
double T_Offset=1696.0;
double *RotAngles;
int atr=1;

#define X 0
#define Y 1
#define Z 2
#define TRUE 1

// #define PHI 0
// #define MU  1
// #define T   2

char *Scalars[] = {"PHI", "MU", "T"}; 


#endif
