#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_

#define PI acos(-1.0) 
#define Tolerance 1.0e-06
#define COMPERR 1.0e-8

#define X 0
#define Y 1
#define Z 2

int NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z;
int Blocks_X, Blocks_Y, Blocks_Z;

int DIMENSION;

// Variables : composition field, derivative of bulk free energy w.r.t. composition
cufftDoubleComplex **dfdcDev, **compDev;
cufftDoubleComplex **dfdcHost, **compHost;

// Variables : structural order parameter field, derivative of bulk free energy w.r.t. order parameter
cufftDoubleComplex **dfdphiDev, **phiDev;
cufftDoubleComplex **dfdphiHost, **phiHost;

//FFT Handle
cufftHandle plan;

//total number of simulation steps
int NUMCOMPONENTS;
int NUMPHASES;
long numsteps;
long time_output;
long nsmooth;
long saveT;
char **COMPONENTS;
char **PHASES;

int WRITECOMPOSITION = 0;
int ASCII = 0;
//Configuration to be initialized or to be read
int  restart, initcount = 0;

// Alloy composition, amplitude of white noise to be added to the system
double  alloycomp, Noise_phasefield, Amp_Noise_Phase; 

//Grid size along x, timestep
double  DELTA_X, DELTA_Y, DELTA_Z, DELTA_t;

//Total simulation time (nondimensional)
double     sim_time = 0.0, total_time;

//System dimensions along x and y
long        MESH_X, MESH_Y, MESH_Z, nx_half, ny_half, nz_half;

//Bulk free energy coefficients
double     f0A, f0B, Vm, sigma, width, alpha, w, vf;

//Gradient energy coefficients associated with composition and structural
//order parameter fields
double     kappa_phi;
//Mobility of solute required for CH equation (mobility)
//Relaxation coefficient for CA equation (relax_coeff)
double     ***Diffusivity, relax_coeff, **Gamma, interface_width;
//Excess solute concentration in matrix phase
double     c0, alloy_comp;
double     ***F0;
double     ***ceq;
double     c_alpha_eq, c_beta_eq;
double     wn, T, Ln, En, Tn;
long       SEED;
FILE       *fpout;
double     dkx, dky, dkz;
double     *kx, *ky, *kz;
double     *kx_h, *ky_h, *kz_h;
double     sizescale;
double     temperature;

//Elasticity-related variables
double     *B;
double     mubar = 800, nubar = 0.3333333, Aniso = 3.0;
double     eigen_strn1[3][3], stress1[3][3];
int        elast_int = 0;

struct symmetric_tensor
{
  double xx;
  double yy;
  double zz;
  double yz;
  double xz;
  double xy;
};

struct symmetric_tensor *eigen_strain;

struct Stiffness_cubic
{
  double C11;
  double C12;
  double C44;
};

struct Stiffness_cubic *Stiffness_c;

int calc_interface_energy;

int Num_of_blocks;

dim3 Gridsize, Blocksize;

long layer_size;

long rows_x, rows_y, rows_z;
long *start, *end;

struct fields
{
    double *phia;
    double *compi;
    double *deltaphi;
    double temperature;
};

struct fields *gridinfo;

struct fill_cylinder
{
    long x_center;
    long y_center;
    long z_start;
    long z_end;
    double radius;
};

struct fill_cylinder fill_cylinder_parameters;

struct fill_sphere
{
    long x_center;
    long y_center;
    long z_center;
    double radius;
};

struct fill_sphere fill_sphere_parameters;

size_t double_size;
size_t complex_size;

// Function pointers
short solvertype = 0;
void (*evolve)(char *argv[]);

// Thermo-writer vars
short thermo_writer;

// mGPU vars
int nGPUs, *whichGPUs;
int GPU_ID, GPU_width;

#endif
