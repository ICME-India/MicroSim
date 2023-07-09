#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
//#include "solverloop/defines1.h"

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif


int devicenumber=0;
int platformnumber=0;

// struct grid {
//   double phi[npha];
//   double mu[nsol];
//   double com[nsol];
// };
struct csle {
  double comie[npha][nsol];
};
struct pfmpar {
  double surfTen[npha*npha];
  double ee[npha*npha];
  double w[npha*npha];
  double wh[npha*npha];
  double eesqrt[npha*npha];
  double w_abc[npha*npha*npha];
  double deltax;
  double deltay;
  double deltaz;
  double deltat;
  double Er;
  double E0;
};
struct pfmval {
  double Rotation_matrix[npha][npha][3][3];
  double Inv_Rotation_matrix[npha][npha][3][3];
  double D[npha][nsol*nsol];
  double DInv[npha][nsol*nsol];
  double cguess[npha*npha][nsol];
  double c_eq[npha*npha][nsol];
  double c_fill[npha*npha][nsol];
  double angle[npha*npha*3];
  double gamma[npha*npha];
  double epsc[npha*npha];
  double Tr;
  double Vm;
  double phisolid;
  double philiquid;
  double Rg;
  double T0;
  double Teq;
  double Tfill;
  double lrep;
  double com_1_Initial;
  double com_2_Initial;
  double com_1_1stguess;
  double com_2_1stguess;
  double a2;
  double rad;
  double intwidth;
  double dxIntfcPoints;
  double dtParam;
  double RefD;
  double TLiquidus;
  double Toffset;
  double TPosOffset;
  double TGRADIENT;
  double velocity;
  double NoiseFac;
  double interfaceUplimit;
  double interfaceDownlimit;
  double deltat_e;
  double rho;
  double damping_factor;
  long shift_OFFSET;
  int thermophase[npha];
  int   nproc;
  int   Nx;
  int   Ny;
  int   Nz;
  int   ntimesteps;
  int   savetime;
  int   myrank;
  int   ISOTHERMAL;
  int atr;
  int Noise_phasefield;
  int Function_anisotropy;
  int anisotropy_type;
  int ELASTICITY;
  int DIMENSION;
};
struct propmatf3 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
  double A[npha][nsol][nsol];
  double DELTA_T[npha][npha];
  double DELTA_C[npha][nsol];
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
};
struct propmatf4 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
  double A[npha][nsol][nsol];
  double DELTA_T[npha][npha];
  double DELTA_C[npha][nsol];
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
};
struct propmatf4spline {
  double A[npha][nsol][nsol];
  double B[npha][nsol];
  double C[npha];
};


//long t;

struct grid *gridNew;
struct grid *gridOld;
struct csle *cscl;
struct pfmval pfmdat;
struct pfmpar pfmvar;
struct propmatf3 propf3;
struct propmatf4 propf4;
struct propmatf4spline *propf4spline;
struct propmatf4spline *propf4spline1;

struct csle *cscl_buf1Dss;
struct csle *cscl_buf1Dse;
struct csle *cscl_buf1Dss;
struct csle *cscl_buf1Dse;

double *temp;
int timestep;
long *tstep;

int device_run;

int tstart=1;

int nx;
int ny;
int nz;

int nynz;
int nxnynz;
int ntNoNoise;

int rank;

int globaldim0;
int globaldim1;

int np;

int i, j, k, i1;

//cl_platform_id platform_id = NULL;
//cl_device_id device_id = NULL;
cl_context context = NULL;
cl_command_queue cmdQ = NULL;
cl_program program;
cl_kernel ker_SolverCsClEq_F2;
cl_kernel ker_SolverCsClEq_F3;
cl_kernel ker_SolverCsClEq_F4;
cl_kernel ker_SolverPhi_F2_smooth;
cl_kernel ker_SolverPhi_F3_smooth;
cl_kernel ker_SolverPhi_F4_smooth;
cl_kernel ker_SolverPhi_F2;
cl_kernel ker_SolverPhi_F3;
cl_kernel ker_SolverPhi_F4;
cl_kernel ker_SolverCatr_F2_smooth;
cl_kernel ker_SolverCatr_F3_smooth;
cl_kernel ker_SolverCatr_F4_smooth;
cl_kernel ker_SolverCatr_F2;
cl_kernel ker_SolverCatr_F3;
cl_kernel ker_SolverCatr_F4;
cl_kernel ker_SolverStress_iterative;
cl_kernel ker_SolverStress_iterative_2D;
cl_kernel ker_apply_BC_phi_y0_noflux;
cl_kernel ker_apply_BC_phi_yn_noflux;
cl_kernel ker_apply_BC_phi_y0_periodic;
cl_kernel ker_apply_BC_phi_yn_periodic;
cl_kernel ker_apply_BC_phi_z0_noflux;
cl_kernel ker_apply_BC_phi_zn_noflux;
cl_kernel ker_apply_BC_phi_z0_periodic;
cl_kernel ker_apply_BC_phi_zn_periodic;
cl_kernel ker_apply_BC_phi_x0_noflux;
cl_kernel ker_apply_BC_phi_xn_noflux;
cl_kernel ker_apply_BC_phi_x0_periodic;
cl_kernel ker_apply_BC_phi_xn_periodic;
cl_kernel ker_apply_BC_com_y0_noflux;
cl_kernel ker_apply_BC_com_yn_noflux;
cl_kernel ker_apply_BC_com_y0_periodic;
cl_kernel ker_apply_BC_com_yn_periodic;
cl_kernel ker_apply_BC_com_z0_noflux;
cl_kernel ker_apply_BC_com_zn_noflux;
cl_kernel ker_apply_BC_com_z0_periodic;
cl_kernel ker_apply_BC_com_zn_periodic;
cl_kernel ker_apply_BC_com_x0_noflux;
cl_kernel ker_apply_BC_com_xn_noflux;
cl_kernel ker_apply_BC_com_x0_periodic;
cl_kernel ker_apply_BC_com_xn_periodic;
cl_kernel ker_apply_BC_ela_y0_noflux;
cl_kernel ker_apply_BC_ela_yn_noflux;
cl_kernel ker_apply_BC_ela_y0_periodic;
cl_kernel ker_apply_BC_ela_yn_periodic;
cl_kernel ker_apply_BC_ela_z0_noflux;
cl_kernel ker_apply_BC_ela_zn_noflux;
cl_kernel ker_apply_BC_ela_z0_periodic;
cl_kernel ker_apply_BC_ela_zn_periodic;
cl_kernel ker_apply_BC_ela_x0_noflux;
cl_kernel ker_apply_BC_ela_xn_noflux;
cl_kernel ker_apply_BC_ela_x0_periodic;
cl_kernel ker_apply_BC_ela_xn_periodic;
cl_kernel kernel4;
cl_kernel kernel5[9];    
cl_kernel kernel6[9];  
cl_kernel kernel7[9];  
cl_kernel ker_addNoise;  
cl_kernel ker_copy_New_To_Old;
cl_kernel kernel10[3];
cl_kernel kernel11[2];
cl_kernel kernel12;
cl_uint ret_num_devices;
cl_uint ret_num_platforms;
cl_int ret;
cl_ulong local_size;
cl_int cl_local_size;
size_t kernel_code_size;
size_t local_size_size;
cl_event kerndone;
cl_char device_name[1024] = {0}; 
cl_char vendor_name[1024] = {0};
cl_device_type device_type; 
cl_char string[10240] = {0};
cl_uint work_dim=3;

size_t globaldim[3];

// cl_mem d_gridNew;
// cl_mem d_gridOld;
cl_mem d_iter_gridinfom;
cl_mem d_gridinfomN;
cl_mem d_gridinfomO;
cl_mem d_cscl;
cl_mem d_pfmdat;
cl_mem d_pfmvar;
cl_mem d_temp;
cl_mem d_tstep;
cl_mem d_propf3;
cl_mem d_propf4;
cl_mem d_propf4spline;
cl_mem d_propf4spline1;
cl_mem d_eigen_strain_phase;
cl_mem d_stiffness_phase;
cl_mem d_stiffness_phase_n;
//cl_mem d_eigen_strain;
//cl_mem d_strain;
//cl_mem d_stiffness_c;
//cl_mem d_sigma;

char *source_str;

void CL_build();
void CL_create_kernels();
void CL_device_selection();
void CL_initialize_variables();
void CL_memory_apis();
void CL_Initialize_domain();
void CL_kernel_init_temperature();
void CL_Update_Temperature(long t);
void CL_Solve_phi_com();
void CL_DeviceToHost();
void CL_Global_Max_Min();
void CL_Shift();
void (*CL_Solve_phi_com_Function)();

void initialize_domain(struct grid *, struct grid *, struct csle *, struct pfmval *, struct pfmpar *, double *);
void savetimestep(struct grid *, struct pfmval *, struct pfmpar *, double *, int);
void savetimestepcsle(struct csle *, struct pfmval *, int);
int RandRange(int , int );
void copy_New_To_Old(struct grid *, struct grid *, struct pfmval *);
void savetimestepcscl(struct csle *, struct pfmval *, struct pfmpar *, int);
void apply_BC(struct grid *, struct pfmval *);
void read_rank(struct grid *, struct csle *, double *temp, struct pfmval *, struct pfmpar *, int);
void apply_BC_temp_ib_noflux(double *, struct pfmval *);
void apply_BC_temp_it_noflux(double *, struct pfmval *);
void propf3Hostupdate(struct propmatf3 *);
void propf4Hostupdate(struct propmatf4 *);
void FunctionF_4_SplineCPU(struct propmatf4spline *propf4spline, double *temp, int );
void mpi_exchange_dim0(int rank);
void mpi_exchange_dim0_iter(int rank);
