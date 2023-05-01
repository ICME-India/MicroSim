#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_

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

double deltat_e = 0.2;

int NUMPHASES;
int NUM_THERMO_PHASES;
int NUMCOMPONENTS;
int t;


long ntimesteps;
long saveT;
long nsmooth;
long STARTTIME=0;
long RESTART=0;

double Teq;
double Tfill;
double T;
int TEMPGRADY=0;
int ISOTHERMAL=0;
double TEMPERATURE_SCALE;
double DELTAT;

double epsilon;
double tau;
double R; 
double V;
double Lf; //Latent-heat; functionF=5
double therm_cond;

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
// double divphi[NUMPHASES],lambda_phi[NUMPHASES];
double *divphi, *lambda_phi;
// double divflux[NUMCOMPONENTS-1], divjat[NUMCOMPONENTS-1];
double *divflux, *divjat;
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


double DET;

// double c_old[NUMCOMPONENTS-1],c_new[NUMCOMPONENTS-1],c[NUMCOMPONENTS-1];
double *c_old, *c_new, *c, *c_tdt;
// double Diffusivity[NUMPHASES][NUMCOMPONENTS-1][NUMCOMPONENTS-1];
double ***Diffusivity;
double **dcdmu, **inv_dcdmu, *deltamu, *deltac, *sum, ***dcdmu_phase, **Ddcdmu;
long bulk_phase;
// double Gamma[NUMPHASES][NUMPHASES];
double **Gamma;
// double tau_ab[NUMPHASES][NUMPHASES];
double **tau_ab;
// double dab[NUMPHASES][NUMPHASES];
double **dab;
// double ec[NUMPHASES][NUMPHASES];
double **ec;
// double e2[NUMPHASES][NUMPHASES];
double **e2;
// double e4[NUMPHASES][NUMPHASES];
double **e4;
// double beta[NUMPHASES][NUMPHASES];
double **beta;
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
  double *mu_max;
  double *mu_min;
  double *rel_change_phi;
  double *rel_change_mu;
  long INTERFACE_POS_MAX;
  long INTERFACE_POS_MIN;
};

struct max_min global_max_min;
struct max_min workers_max_min;


struct bc_scalars {
  int type;
  long points[3];
  long proxy[3];
  double value[3];
};

struct bc_scalars *boundary[6];


long rows_x, rows_y, rows_z;
int *start, *end;
long *averow, *rows, *offset, *extra;

//Variables for mpi........................................................
int     taskid,                                                                  /* this task's unique id */
        numworkers,                                                              /* number of worker processes */
        numworkers_x,
        numworkers_y,
        numworkers_z,
        numtasks,                                                                /* number of tasks */
//         averow[DIMENSION],rows[DIMENSION],offset[DIMENSION],extra[DIMENSION],    /* for sending rows of data */
//         rows_x, rows_y,
        dest, source,                                                            /* to - from for message send-receive */
        left_node,right_node,                                                    /* neighbor tasks */
        top_node, bottom_node,
        front_node, back_node,
        msgtype,                                                                 /* for message types */
        rc;
// long    start[DIMENSION],end[DIMENSION];                                         /* misc */
        MPI_Status status;
        MPI_Datatype MPI_gridinfo;                                               //New Datatype to send structure
        MPI_Datatype MPI_gridinfo_vector_b;
        MPI_Datatype MPI_gridinfo_vector_c;
        MPI_Datatype MPI_iter_gridinfo;                                               //New Datatype to send structure
        MPI_Datatype MPI_gridinfo_vector_b_stress;
        MPI_Datatype MPI_gridinfo_vector_c_stress;
long    offset_x, offset_y, offset_z;
int boundary_worker=0;

struct workers {
  int lastx;
  int firstx;
  int lasty;
  int firsty;
  int lastz;
  int firstz;
  int rank_x;
  int rank_y;
  int rank_z;
  long offset[3];
  long offset_x;
  long offset_y;
  long offset_z;
  long rows[3];
  long rows_x;
  long rows_y;
  long rows_z;
  int left_node;
  int right_node;
  int top_node;
  int bottom_node;
  int front_node;
  int back_node;
  long start[3];
  long end[3];
  long layer_size;
};
struct workers workers_mpi;
//Variables for mpi........................................................

//Variables for mpi........................................................

double gradphixf, gradphixb, gradphiyb, gradphiyf;
double gradphix, gradphiy, gradphiz;
double gradphix_l, gradphiy_l, normgradphix_l, normgradphiy_l, normgradphiz_l;
double gradphix_l_y, gradphix_l_z;
double gradphiy_l_x, gradphiy_l_z;
double gradphiz_l_x, gradphiz_l_y;
double normgradphi;
int t, to, n, starttime;
long gidy, center, center1, front, back, left, right, top, bottom;
long index_, index_count;
double s_phi_center,s_phi_right,s_phi_front, s_phi_top, s_phi_bottom;
long a,b;
double sum_lambdaphi, sum_dhphi;
long active_phases, count_phases;
double Deltaphi;
long interface;
double scalprod;
long i,j,k;
long x;
long iter;
double mu;
double global_error=0.0;
double error;
double tolerance=1e-12;
long MAX_ITERATIONS=5;
int ASCII=0;
int WRITEHDF5;
herr_t status_h;
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
// struct symmetric_tensor *eigen_strain_field;
// 
// eigen_strain_field = (struct symmetric_tensor *)malloc(MESH_X*MESH_Y*sizeof(*eigen_strain_field));

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
  double *phia;
  double *compi;
  double *composition;
  double *deltaphi;
  double temperature;
};
double rho=10;
double damping_factor=0.8;

struct iter_variables {
 double disp[3][3];
};

// struct iter_variables *iter_gridinfo;
struct iter_variables *iter_gridinfo_w;
struct iter_variables *iter_gridinfo_w_instance;

double *buffer;
double *buffer_boundary_x;
double *buffer_boundary_y;
double *buffer_boundary_z;
double *buffer_boundary_x_stress;
double *buffer_boundary_y_stress;
double *buffer_boundary_z_stress;
#define SIZE_STRUCT_FIELDS (2*NUMPHASES+2*(NUMCOMPONENTS-1)+1)

// struct fields {
//   double phia[2];
//   double compi[1];
//   double deltaphi[2];
//   double temperature;
// };

struct fields *gridinfo;
struct fields *gridinfo_w;
struct fields *gridinfo_instance;


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
//  double ***strain;
//  double **stiffness;
//  double ***eigen_strain;
 struct symmetric_tensor strain[3];
 struct symmetric_tensor eigen_strain[3];
 struct Stiffness_cubic  stiffness_c[3];
};
 //double strain[DIMENSION][DIMENSION][DIMENSION];
 //double stiffness[DIMENSION][3];			//2 to 3
 //double eigen_strain[DIMENSION][DIMENSION][DIMENSION];
// struct symmetric_tensor *eigen_strain_phase;
// struct Stiffness_cubic  *stiffness_phase;

int ELASTICITY=0;
int GRAIN_GROWTH=0;

struct gradlayer **gradient;
struct gradlayer *gradient1[4];
struct gradlayer *tmp;
struct gradlayer *grad, *grad_right, *grad_left, *grad_back, *grad_front, *grad_boundary, *grad_top, *grad_bottom;
struct gradlayer test;

// #define MAXWORKER   8                  /* maximum number of worker tasks */
// #define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define BTAG        555                /* message tag */
#define TTAG        666                /* message tag */
#define FTAG        777                /* message tag */
#define BATAG       888                /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */
#define SHIFT_SIGNAL 5
#define SHIFT_POS    6

#define X 0
#define Y 1
#define Z 2
#define TRUE 1

int FUNCTION_F = 1;

// void(*function_F_dpsi[])(double *mu, double T, double *phi, long a)      = {function_f_01_dpsi};
// void(*function_F_free_energy[])(double *c, double T, long a)             = {function_f_01_free_energy};
// void(*function_F_Mu[])(double *c, double T, long a, long i)              = {function_f_01_free_energy};
// void(*function_F_c_mu[])(double *mu, double T, long a, long i)           = {function_f_01_c_mu};
// void(*function_F_dc_dmu[])(double *mu, double T, long a, long i, long j) = {function_f_01_dc_dmu};
// void(*function_F_function_A[])(double T1, long i, long j, long a)        = {function_f_01_function_A};
// void(*function_F_function_B[])(double T, long i, long a)                 = {function_f_01_function_B};
// void(*function_F_function_C[])(double T, long a)                         = {function_f_01_function_C};
// void(*function_F_compute_chemicalpotential[])(double T, long a)          = {function_f_01_compute_chemicalpotential};


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
// void (*compute_chemicalpotential)(struct fields* gridinfo);
// double (*free_energy)(double *c, double T, long a);

double (*dwdphi)(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a);
double (*dwdphi_smooth)(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a);

double (*dAdphi)(double *phi, struct gradlayer **gradient, long gidy, long a);
double (*dAdphi_smooth)(double *phi, struct gradlayer **gradient, long gidy, long a);
double (*divdAdgradphi)(struct gradlayer **gradient, long index, long gidy, long a);
double (*divdAdgradphi_smooth)(struct gradlayer **gradient, long index, long gidy, long a);


void (*dAdq)(double *qab, double* dadq, long a, long b);
double (*function_ac)(double *qab, long a, long b);

void (*calculate_gradients)(long x, struct gradlayer **gradient);
// void (*calculate_gradients_phasefield)(long x, struct gradlayer **gradient);
void (*calculate_gradients_phasefield)(long x, struct gradlayer **gradient, int CALCULATE_COMPOSITION);
void (*calculate_gradients_concentration)(long x, struct gradlayer **gradient);
void (*calculate_fluxes_concentration)(long x, struct gradlayer **gradient);
void (*calculate_divergence_concentration)(long x, struct gradlayer **gradient);
void (*calculate_divergence_concentration_smooth)(long x, struct gradlayer **gradient);
void (*calculate_divergence_concentration_smooth_concentration)(long x, struct gradlayer **gradient);
void (*calculate_divergence_phasefield)(long x, struct gradlayer **gradient);
void (*calculate_divergence_phasefield_smooth)(long x, struct gradlayer **gradient);
void (*calculate_divergence_stress)(long x, struct gradlayer **gradient);
void (*calculate_gradients_stress)(long x, struct gradlayer **gradient);
void (*compute_error)(long x, double *error);
void (*solverloop_phasefield)(long *start, long *end);
void (*solverloop_concentration)(long *start, long *end);
void (*writetofile_mpi)(struct fields* gridinfo, char *argv[], long t);
void (*writetofile_mpi_binary)(struct fields* gridinfo, char *argv[], long t);
void (*writetofile_mpi_hdf5)(struct fields* gridinfo, char *argv[], long t);
void (*readfromfile_mpi_hdf5)(struct fields* gridinfo, char *argv[], long numworkers, long t);
void (*readfromfile_mpi)(struct fields* gridinfo, char *argv[], long t);
void (*readfromfile_mpi_binary)(struct fields* gridinfo, char *argv[], long t);
double (*df_elast)(struct gradlayer *gradient, struct symmetric_tensor sigma, double *phi, long a);

// #define PHI 0
// #define MU  1
// #define T   2

char *Scalars[] = {"PHI", "MU", "T", "U"}; 
char dirname[1000];
char **coordNames;
long size_fields;
FILE *fp;
long position,file_iter;
long time_file;
#endif
