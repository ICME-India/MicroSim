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
Real gammaa;
Vector <std::string> val;
Vector <Real> dbl;
Vector <std::string> s_diff;
Vector <Vector<double>> diffu;
Vector<double> diff;
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
Real dab;
//Rotation matrix
Vector<std::string> s_rotmat;
Vector<Vector<Real>> rotmat;
//Potential function
std::string s_funcW;
int funcW;
//std::string gamma_abc;
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
Vector<Vector<Vector<Real>>> A;
Vector<Vector<Vector<Real>>> Aeq;
//B and D
Vector <Real> B;
Vector <Real> D;
//Real BB;
//Real DD;

Vector<Vector<double>> A_values;
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
Vector <amrex::Real> dcdmu;
Vector <amrex::Real> dcdmu_eq;

//cmu
Vector <amrex::Real> cmu;

//Dummy
// Real sum =0.0;
// Real sum1 =0.0;  

//Function pointers
void (*dc_dmu)();
void (*c_mu)(MultiFab& sun);
void (*Mu)(MultiFab& ray);
void (*function_A)();
void (*function_B)();
void (*function_D)();
void (*free_energy)();
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

Vector<Vector<Vector<Real>>> conc;
Vector<Vector<Vector<Real>>> conceq;
Vector<Vector<double>> conc_Sol;
Vector<Vector<double>> conc_Liq;
Vector<Vector<double>> temprt;
#endif
