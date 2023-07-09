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

double deltat_e = 0.2;

int NUMPHASES;
int NUM_THERMO_PHASES;
int NUMCOMPONENTS;
long t;


long ntimesteps;
long saveT;
long nsmooth;
long STARTTIME=0;
long RESTART=0;
long WRITEHDF5=0;

double Teq;
double Tfill;
double T=-1.0;
double temperature=-1.0;
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
char **Phases_tdb;
char **phase_map;
long *thermo_phase;

double t_smooth;

double ***cmu,***muc;
// double A[NUMPHASES][NUMCOMPONENTS-1][NUMCOMPONENTS-1];
double ***A;
// double DELTA_T[NUMPHASES][NUMPHASES];
double **DELTA_T;
// double DELTA_C[NUMPHASES][NUMCOMPONENTS-1];
double **DELTA_C;
// double dcbdT[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***dcbdT;
// double dcbdT_phase[NUMPHASES][NUMCOMPONENTS-1];
double **dcbdT_phase;
// double B[NUMPHASES][NUMCOMPONENTS-1];
double **B;
// double Beq[NUMPHASES][NUMCOMPONENTS-1];
double **Beq;
// double dBbdT[NUMPHASES][NUMCOMPONENTS-1];
double **dBbdT;
// double C[NUMPHASES];
double *C;
// double ceq[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***ceq;
// double cfill[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***cfill;
// double ceq_coeffs[NUMPHASES][NUMCOMPONENTS-1][4];
double ***ceq_coeffs;
// double slopes[NUMPHASES][NUMPHASES][NUMCOMPONENTS-1];
double ***slopes;
double ***c_guess;
double ****comp_ES;
double ****ThF;
double **T_ES;
double **T_ThF;
gsl_interp_accel ****acc_ES;
gsl_spline ****spline_ES;
gsl_interp_accel ****acc_ThF;
gsl_spline ****spline_ThF;

// double c_old[NUMCOMPONENTS-1],c_new[NUMCOMPONENTS-1],c[NUMCOMPONENTS-1];
double *c_old, *c_new, *c;
// double Diffusivity[NUMPHASES][NUMCOMPONENTS-1][NUMCOMPONENTS-1];
double ***Diffusivity;
double ***DiffusivityInv;
double **dcdmu, **inv_dcdmu, *deltamu, *deltac, *sum, ***dcdmu_phase, **Ddcdmu;
long bulk_phase;
// double Gamma[NUMPHASES][NUMPHASES];
double **Gamma;
// double dab[NUMPHASES][NUMPHASES];
double **dab;
// double Gamma_abc[NUMPHASES][NUMPHASES][NUMPHASES];
double ***Gamma_abc;

int BOUNDARY_LEFT;
int BOUNDARY_RIGHT;
int BOUNDARY_FRONT;
int BOUNDARY_BACK;
int BOUNDARY_TOP;
int BOUNDARY_BOTTOM;


struct max_min {
  double *phi_max;
  double *phi_min;
  double *com_max;
  double *com_min;
  double *composition_max;
  double *composition_min;
  double *rel_change_phi;
  double *rel_change_com;
  double *rel_change_composition;
};

struct max_min global_max_min;
struct max_min workers_max_min;
struct max_min global_max_min1;


struct bc_scalars {
  int type;
  long points[3];
  long proxy[3];
  double value[3];
};

struct bc_scalars *boundary[6];



struct mpiparameters {
  long rows_x;
  long rows_y;
  long rows_z;
  long startx;
  long starty;
  long startz;
  long endx;
  long endy;
  long endz;
};

struct mpiparameters mpiparam;
#define cart_grid_ndim (1)
#define MASTER      (0)
double *buffer_sxs;
double *buffer_sxe;
double *buffer_rxs;
double *buffer_rxe;
double *buffer_sxs2;
double *buffer_sxe2;
double *buffer_rxs2;
double *buffer_rxe2;
double *iter_buffer_sxs;
double *iter_buffer_sxe;
double *iter_buffer_rxs;
double *iter_buffer_rxe;

MPI_Comm comm=MPI_COMM_WORLD, COMM_NEW;
MPI_Status status1, status2;
MPI_Datatype MPI_gridinfo;
int tag_dim0_m, tag_dim0_p;
int numtasks;
int rank, numworkers;
int sizempixyz, size_buffer, sizebufyz, size_buffer2, iter_size_buffer;
int dim0_rank_p, dim0_rank_m;


long rows_x, rows_y, rows_z;
long *start, *end;
long *averow, *rows, *offset, *extra;


int to, n, starttime;
long a,b;
double sum_lambdaphi, sum_dhphi;
long active_phases, count_phases;
double Deltaphi;
long interface;
double scalprod;
//long k;
long iter;
double mu;
double global_error=0.0;
double error;
double tolerance=1e-12;
long MAX_ITERATIONS=5;
int ASCII=0;
//herr_t status_h;
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
double BASE_POS=0;
double GRADIENT;
double temp_bottom;

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
int shift_ON=0;
long shift_position=0;
long position;
long time_file;
long file_iter;

long INTERFACE_POS_GLOBAL=0;
long MAX_INTERFACE_POS=0;

  
struct filling_type {
  long NUMCUBES;
  long NUMCIRCLES;
  long NUMTRIANGLES;
  double volume_fraction;
  long length;
};
struct filling_type *filling_type_phase;

// struct tensor {
//   
// }
// 
struct symmetric_tensor {
  double xx;
  double yy;
  double zz;
  double yz;
  double xz;
  double xy;
};

struct symmetric_tensor *eigen_strain_phase;
struct symmetric_tensor strain[3];
struct symmetric_tensor eigen_strain[3];


struct elast_properties {
  double Poisson_ratio;
  double Y_modulus;
  double B_modulus;
  double S_modulus;
  char symmetry[10];
};

struct Stiffness_cubic {
  double C11;
  double C12;
  double C44;
};

struct Stiffness_cubic *stiffness_phase;
struct Stiffness_cubic *stiffness_phase_n;
struct Stiffness_cubic  stiffness_c[3];

struct Stiffness_tetragonal {
   double C11;
   double C12;
   double C13;
   double C33;
   double C44;
   double C66;
};

struct Stiffness_tetragonal *stiffness_t_phase;

struct fields {
  double phia[npha];
  double compi[nsol];
  double composition[nsol];
  double temperature;
};
double rho=10;
double damping_factor=0.8;

struct iter_variables {
 double disp[3][3];
};

// struct iter_variables *iter_gridinfo;
struct iter_variables *iter_gridinfo;
struct iter_variables *iter_gridinfom;
struct iter_variables *iter_gridinfo_buf2Dss;
struct iter_variables *iter_gridinfo_buf2Dse;
struct iter_variables *iter_gridinfo_buf2Drs;
struct iter_variables *iter_gridinfo_buf2Dre;


double *buffer;
double *buffer_boundary_x;
double *buffer_boundary_y;
double *buffer_boundary_z;
double *buffer_boundary_x_stress;
double *buffer_boundary_y_stress;
double *buffer_boundary_z_stress;


struct fields *gridinfo;
struct fields *gridinfomN;
struct fields *gridinfomO;
struct fields *gridinfo_buf1Dss;
struct fields *gridinfo_buf1Dse;
struct fields *gridinfo_buf1Drs;
struct fields *gridinfo_buf1Dre;
struct fields *gridinfo_buf2Dss;
struct fields *gridinfo_buf2Dse;
struct fields *gridinfo_buf2Drs;
struct fields *gridinfo_buf2Dre;

#define SIZE_STRUCT_FIELDS (NUMPHASES+2*(NUMCOMPONENTS-1)+1)

struct gradlayer {
 double **gradphi;
 double **gradphi_c;
 double **phistagg;
 double **gradchempot;
 double ***Dmid;
 double **jat;
 double *deltaphi;
 double **phase_comp;
 double **dcbdT_phase;
 double ***dcdmu_phase;
 int interface;
 int bulk_phase;
//  struct symmetric_tensor strain[3];
//  struct symmetric_tensor eigen_strain[3];
//  struct Stiffness_cubic  stiffness_c[3];
};




int ELASTICITY=0;
int GRAIN_GROWTH=0;

struct gradlayer **gradient;
struct gradlayer *gradient1[4];
struct gradlayer *tmp;
struct gradlayer *grad, *grad_right, *grad_left, *grad_back, *grad_front, *grad_boundary;
struct gradlayer test;

long tNoiseStart;
double TLiquidus;
int DirectionalSolidification=0;
double TemperatureGradient=1e6;
double PullingVelocity=0.05;
int PositionOffset=20;
double T_Offset=1696.0;
double ***RotAngles;
int atr=1;
char tdbfname[100];

#define X 0
#define Y 1
#define Z 2
#define TRUE 1

int FUNCTION_F = 3;

// #define PHI 0
// #define MU  1
// #define T   2




double (*free_energy)(double *c, double T, long a);
// double (*Mu)(double *c, double T, long a, long i);
void (*Mu)(double *c, double T, long a, double *Mu);
// double (*c_mu)(double *mu, double T, long a, long i);
void (*c_mu)(double *mu, double *c, double T, long a, double *c_guess);
// double (*dc_dmu)(double *mu, double T, long a, long i, long j);
void (*dc_dmu)(double *mu, double *phase_comp, double T, long a, double **dcdmu);
// double (*dpsi)(double *mu, double T, double *phi, long a);
double (*dpsi)(double *mu, double **phase_comp, double T, double *phi, long a);
void (*function_A)(double T, double ***c);
double (*function_B)(double T, long i, long a);
double (*function_C)(double T, long a);
void (*init_propertymatrices)(double T);

void CL_define_npha_com();
void CL_thermo_function_calls();

char *Scalars[] = {"PHI", "MU", "T", "U"};
char dirname[1000];
char **coordNames;
long size_fields;
FILE *fp;
long position,file_iter;
long time_file;

#endif
