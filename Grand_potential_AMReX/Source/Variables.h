#ifndef _VARIABLES_H_
#define _VARIABLES_H_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

using namespace amrex;

//#ifndef phasecount
#define phasecount 7
//#endif

//#ifndef compcount
#define compcount 3
//#endif

#ifndef X
#define X 0
#endif

#ifndef Y
#define Y 1
#endif

#ifndef Z
#define Z 2
#endif

#ifndef cent
// center
#define cent 0
#endif


#ifndef iph
// i + 1/2
#define iph 1
#endif

#ifndef imh
// i - 1/2
#define imh 2
#endif

#ifndef jph
// j + 1/2
#define jph 3
#endif

#ifndef jmh
// j - 1/2
#define jmh 4
#endif

#ifndef kph
// k + 1/2
#define kph 5
#endif

#ifndef kmh
// k - 1/2
#define kmh 6
#endif


//Geometrical dimensions
std::string s_dim;
int dim;
std::string s_ncellx;
int ncellx;
std::string s_ncelly;
int ncelly;
std::string s_ncellz;
int ncellz;
//Descretization
std::string s_dx;
Real dx;
std::string s_dy;
Real dy;
std::string s_dz;
Real dz;
std::string s_dt;
Real dt;
//Number of phases
std::string s_nump;
int nump;
std::string s_numcom;
int numcom;
//Running and saving
std::string s_nsteps;
int nsteps;
std::string s_nsmooth;
int nsmooth;
std::string s_savet;
int savet;
std::string s_startt;
int startt;
std::string s_restart;
int restart;
std::string s_numworkers;
int numworkers;
//Component and phase names
std::string s_comp;
Vector <std::string> comp;
std::string s_phase;
Vector <std::string> phase;
//Material properties
std::string s_gammaa;
Vector<Real> gammaa;
Vector<Vector<Real>> gam;
Vector <std::string> val;
Vector <Real> dbl;
Vector <std::string> s_diff;
Vector <Vector<double>> diffu;
Vector<Vector<Vector<double>>> diff;
//Gas constant and molar volume
std::string s_R;
Real R;
std::string s_Vm;
Real Vm;
//Elastic parameters
Vector <std::string> s_egstr;
Vector <Vector<double>> egstr;
Vector <std::string> s_voigiso;
Vector <Vector<double>> voigiso;
//Boundary
Vector <std::string> s_bound;
Vector <Vector<string>> bound;
//Boundary value
Vector <std::string> s_boundval;
Vector <Vector<string>> boundval;

//Type of simulation
std::string s_isothermal;
int isothermal;
std::string s_binary;
int binary;
//Ternary
std::string s_dilute;
int dilute;
std::string s_T;
Real T;
//Filewriting
std::string s_writeformat;
std::string writeformat;
std::string s_writehdf5;
int writehdf5;
std::string s_trackprog;
int trackprog;
//Model specific GP model
std::string s_eps;
Real eps;
std::string s_tau;
Real tau;
std::string s_Tau;
Real Tau;
Real tau_final; 
//Anisotropy
std::string s_funcANI;
int funcANI;
std::string s_ANItype;
int ANItype;
std::string s_dab;
Vector<Real> dab;
Vector<Vector<Real>> dabb;
//Rotation matrix
Vector<std::string> s_rotmat;
Vector<Vector<Real>> rotmat;
//Potential function
std::string s_funcW;
int funcW;
std::string gamma_abc;
Vector<Real> gammaa_abc;
Vector<Vector<Vector<Real>>> gam_abc;
//Shifting of domain
std::string s_shiftdom;
Real shiftdom;
std::string s_shiftj;
Real shiftj;
//Write composition and chemical ptential fields
std::string s_writecomp;
int writecomp;
//Noise
std::string s_noise_pf;
int noise_pf;
std::string s_amp_noise_phase;
Real amp_noise_phase;
//Temperature
std::string s_Teq;
Real Teq;
std::string s_Tfill;
Real Tfill;
//Temperature gradient
std::string s_tempgrady;
Vector<Real> tempgrady;
//Function_F
std::string s_funcf;
int funcf;
//A
std::vector <std::string> s_A;
Vector<Vector<Real>> A1;
Vector<Vector<Vector<Real>>> A{};
Vector<Vector<Vector<Real>>> Aeq{};
//B and D
Vector<Vector<Real>> B{};
Vector <Real> C{};
//Real BB;
//Real DD;

Vector<Vector<Vector<double>>> A_values;
Vector<Vector<double>> A_temp;
//Vector<std::array<int, 1024>> Aliq;

//ceq
std::vector <std::string> s_ceq;
Vector<Vector<Real>> ceq;
Vector <Vector<Vector<Real>>> c_eq;
//cfill
std::vector <std::string> s_cfill;
Vector<Vector<Real>> cfill;
Vector <Vector<Vector<Real>>> c_fill;
//cguess
std::vector <std::string> s_cguess;
Vector<Vector<Real>> cguess;
Vector <Vector<Vector<Real>>> c_guess;
//slopes
std::vector <std::string> s_slopes;
Vector<Vector<Real>> slopes;
//thermo phase
std::string s_ntp;
int ntp;
//tdbfname
std::string s_tdbname; 
std::string tdbname;
//tdb phase
std::string s_tdbphase;
Vector<std::string> tdb_phase;
//phase map
std::string s_phasemap;
Vector<std::string> phasemap;

//Filling 
Vector <std::string> s_cube;
Vector <std::string> s_sphere;
Vector <std::string> s_cylinder;
Vector <std::string> s_ellipse;
Vector <Vector<Real>> cube;
Vector <Vector<Real>> sphere;
Vector <Vector<Real>> cylinder;
Vector <Vector<Real>> ellipse;

//GSL
gsl_interp_accel *A_accel_ptr;
gsl_spline *A_spline_ptr;

//dcdmu
Vector<Vector<Vector<Real>>> dcdmu;
Vector<Vector<Vector<Real>>> dcdmu_eq;

//cmu
//Vector<Vector<Vector<Real>>> muc{};
//Vector<Vector<Vector<Real>>> cmu{};
Array3D<Real,0,phasecount,0,compcount-1,0,compcount-1> cmu{};
Array2D<Real,0,phasecount,0,compcount-1> c{};

//Vector<Vector<Vector<Real>>> muc_eq{};
Array3D<Real,0,phasecount,0,compcount-1,0,compcount-1> cmu_eq{};


//Dummy
// Real sum =0.0;
// Real sum1 =0.0;  

//Function pointers
void (*dc_dmu)();
//void (*c_mu)(int i, int j, int l, amrex::Array4<Real const> const& mu, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> &c, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA, int nc, int a);
void (*Mu)(MultiFab& ray);
void (*function_A)();
void (*function_B)();
void (*function_C)();
//void (*free_energy)(int i, int j, int k, int numphase, Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, Array1D <Real,0,phasecount-1> CC, Array1D<Real,0,phasecount-1> &fe , Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> &c, int numcomp, int a);
void (*dwdphi)(MultiFab& uhl, MultiFab& klk, Geometry const& mrt);
void (*aniso_term)(MultiFab& qtr, MultiFab& prt, Geometry const& trp);
void (*Chem_pot)(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& comp_new, MultiFab& comp_old, Geometry const& geom);
//void (*dAdq)(int i, int j , int k, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> daby, Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &r_qab,  Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &dadq,  Array1D<Real, 0, AMREX_SPACEDIM*2> &q2, int a, int b);
//void (*function_ac)(int i, int j , int k, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> daby, Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &r_qab,  Array1D<Real, 0, AMREX_SPACEDIM*2> &ac,  Array1D<Real, 0, AMREX_SPACEDIM*2> &q2, int a, int b);
void (*Initialise_phi)(MultiFab& qer);


//General variables
Real timee = 0.0;
std::string restart_chkfile = "";
std::string chk_file = "chk";
int stepnum = 0;
Vector<BoxArray> grids;
int n=0;
Vector <int> numb;
int val1=0;
int val2=0;
Vector<Real> fe;

double strt_time;
double stop_time;
double rst_time=0.0;

Vector<Vector<Real>> conc;
Vector<Vector<Real>> conceq;
Vector<Vector<Vector<double>>> conc_Sol;
Vector<Vector<Vector<double>>> conc_Liq;
Vector<Vector<double>> temprt;

Vector<Vector<Real>> tau_ab;


#endif
