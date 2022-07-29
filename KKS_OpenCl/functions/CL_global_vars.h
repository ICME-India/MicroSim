#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include "solverloop/defines1.h"

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__           
#include <OpenCL/opencl.h>           
#else          
#include <CL/cl.h>           
#endif

int devicenumber=0;
int platformnumber=0;

struct grid {
  double phi;
  double c1;
  double mu[1];
};
struct csle {
  double c1l;
  double c1s;
};
struct pfmpar {
  double E0;
  double surfTen;
  double ee;
  double w;
  double deltax;
  double deltay;
  double deltat;
  //double dx;
  //double dy;
  //double dt;
  double eesqrt;
  double Er;
  double IntMob;
  double IntMobInv;
};
struct pfmval {
  double Rotation_matrix[2][2][3][3];
  double Inv_Rotation_matrix[2][2][3][3];
  double Tr;
  double sigma;
  double Vm;
  double D11l;
  double D11s;
  double phisolid;
  double philiquid;
  double Rg;
  double T0;
  double Teq;
  double Tfill;
  double lrep;
  double c1l_Initial;
  double c1s_Initial;
  double c1l_1stguess;
  double c1s_1stguess;
  double a2;
  double rad;
  double epsc;
  double intwidth;
  double dxIntfcPoints;
  double dtParam;
  double epsm;
  double InterfaceMobility;
  double RefD;
  double angle;
  double TLiquidus;
  double Toffset;
  //double PosOffset;
  //double TG;
  //double Vp;
  double TPosOffset; 
  double TGRADIENT;
  double velocity;
  double NoiseFac;
  double interfaceUplimit;
  double interfaceDownlimit;
  long shift_OFFSET;
  int thermophase[npha];
  int   nproc;
  int   jNx;
  int   iNy;
  int   jDimX;
  int   iDimY;
  int   ntimesteps;
  int   savetime;
  int   myrank;
  int   ISOTHERMAL;
};
struct propmatf3 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double A[npha][nsol][nsol]; 
  double DELTA_T[npha][npha]; 
  double DELTA_C[npha][nsol]; 
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
};
struct propmatf4 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double A[npha][nsol][nsol]; 
  double DELTA_T[npha][npha]; 
  double DELTA_C[npha][nsol]; 
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
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

double *temp;
int timestep;
long *tstep;

int device_run;

int tstart=1;

int nx;
int ny;
int NX;
int NY;

int nxny;
int NXNY;
int ntNoNoise;

int rank;
int totny;
int nxtotny;
int totsliceny;
int nxtotsliceny;
int istart;
int iend;

int globaldim0;
int globaldim1;

int np;

int i;
int j;
int k;
int i1;

//cl_platform_id platform_id = NULL;
//cl_device_id device_id = NULL;
cl_context context = NULL;
cl_command_queue cmdQ = NULL;
cl_program program;
cl_kernel kernel1;
cl_kernel kernel1_2;
cl_kernel kernel1_3;
cl_kernel kernel1_4;
cl_kernel kernel2;
cl_kernel kernel3;
cl_kernel kernel2a;
cl_kernel kernel2a_2;
cl_kernel kernel2a_3;
cl_kernel kernel2a_4;
cl_kernel kernel2b;
cl_kernel kernel2b_2;
cl_kernel kernel2b_3;
cl_kernel kernel2b_4;
cl_kernel kernel3a;
cl_kernel kernel3a_2;
cl_kernel kernel3a_3;
cl_kernel kernel3a_4;
cl_kernel kernel3b;
cl_kernel kernel3b_2;
cl_kernel kernel3b_3;
cl_kernel kernel3b_4;
cl_kernel kernel4;
cl_kernel kernel5[9];    
cl_kernel kernel6[9];  
cl_kernel kernel7[9];  
cl_kernel kernel8;  
cl_kernel kernel9;
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
cl_uint work_dim=2;

size_t globaldim[2];

cl_mem d_gridNew;
cl_mem d_gridOld;
cl_mem d_cscl;
cl_mem d_pfmdat;
cl_mem d_pfmvar;
cl_mem d_temp;
cl_mem d_tstep;
cl_mem d_propf3;
cl_mem d_propf4;
cl_mem d_propf4spline;
cl_mem d_propf4spline1;

char *source_str;

void CL_build();
void CL_create_kernels();
void CL_device_selection();
void CL_initialize_variables();
void CL_memory_apis();
void CL_Initialize_domain();
void CL_kernel_init_temperature();
void CL_Update_Temperature(long t);
void CL_Solve_phi_com(long t);
void CL_DeviceToHost();
void CL_Global_Max_Min();
void CL_Shift();
void (*CL_Solve_phi_com_Function)(long t);

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
