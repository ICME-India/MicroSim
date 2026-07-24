#ifndef AmrCoreGP_H_
#define AmrCoreGP_H_

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>

#include <fstream>
#include <string>
#include <limits>
#include <memory>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
using namespace amrex;

struct AmrCoreFill
{   AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*data*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
                    {		
                    	//do something for external Dirichlet (BCType::ext_dir)
                    }
};


class AmrCoreGP
    : public amrex::AmrCore
{
public:

AmrCoreGP ();
virtual ~AmrCoreGP();

virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,const amrex::DistributionMapping& dm) override;
virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,const amrex::DistributionMapping& dm) override;
virtual void ClearLevel (int lev) override;
virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,const amrex::DistributionMapping& dm) override;
virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

void InitData();
void Evolve();

private:
void AdvancePhiAtLevel(Real time,Real dt, int lev);
void AdvancePhiAllLevels(Real time,Real dt);
void AverageDown();

void Calculate_Tau();
void CalculateGrad(amrex::Real t, int lev);

void dmudt_2D(amrex::MultiFab& mu_new, amrex::MultiFab& mu_old, amrex::MultiFab& phi_new, amrex::MultiFab& phi_old, amrex::MultiFab& comp_new, amrex::MultiFab& comp_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx,amrex::Real dt);
void dmudt_3D(amrex::MultiFab& mu_new, amrex::MultiFab& mu_old, amrex::MultiFab& phi_new, amrex::MultiFab& phi_old, amrex::MultiFab& comp_new, amrex::MultiFab& comp_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx,amrex::Real dt);
void dpsi(amrex::MultiFab& mu_old, amrex::MultiFab& term3, amrex::MultiFab& phi_old, amrex::MultiFab& psi, amrex::GpuArray<Real,AMREX_SPACEDIM> dx);

double findC(int phase, double temp, int i);
double findA(int phase, double temp, int o, int r);
amrex::Real FunctionTau(int i, int j, int k, int numphase, Real Tauu, amrex::Array4<Real const> const& phi);

void function_F_04_function_A();
void function_F_04_function_B();
void function_F_04_function_C();
void function_F_04_dc_dmu();

void function_W_01_dwdphi(amrex::MultiFab& term2, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx);
void function_W_02_dwdphi(amrex::MultiFab& term2, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx);

void Function_Q(int i, int j, int k, amrex::Array4<Real const> const& phi,GpuArray<Real,AMREX_SPACEDIM> delta, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1, Order::C> &qab);
void function_A_00_iso_2D(amrex::MultiFab& term1, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM>& dx);
void function_A_00_iso_3D(amrex::MultiFab& term1, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM>& dx);
void function_A_01_ani_2D(amrex::MultiFab& term1, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM>& dx);
void function_A_01_ani_3D(amrex::MultiFab& term1, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM>& dx);

void FillPatch (std::string ref,int lev, Real time, MultiFab& mf, int icomp, int ncomp);
void FillCoarsePatch (std::string ref, int lev, Real time, MultiFab& mf, int icomp, int ncomp);

void GetData (std::string ref,int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime);
void getc();

void Init_liq(amrex::MultiFab& phi_new);
void Init_phi_cyl(amrex::MultiFab& phi_new);
void Init_phi_sph(MultiFab& phi_new);
void Init_mu(amrex::MultiFab& mu_new, amrex::MultiFab& phi_new);
void Init_comp(amrex::MultiFab& phi_new, amrex::MultiFab& comp_new);

void timeStepWithSubcycling(amrex::Real time, amrex::Real dt, int lev);
void timeStepNoSubcycling(amrex::Real time, amrex::Real dt);

void ReadInput();
void readA();
void readC();

void update_phi(amrex::MultiFab& phi_new, amrex::MultiFab& phi_old, amrex::MultiFab& term1, amrex::MultiFab& term2, amrex::MultiFab& term3,amrex::MultiFab& term4,amrex::MultiFab& lambad, amrex::GpuArray<Real,AMREX_SPACEDIM> dx,amrex::Real dt);

std::string PlotFileName (int lev) const;
amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;
amrex::Vector<std::string> PlotFileVarNames () const;
void WritePlotFile ();


// Required multifab
amrex::Vector<amrex::MultiFab>  phi_new;
amrex::Vector<amrex::MultiFab>  phi_old;
amrex::Vector<amrex::MultiFab>  Grad;
amrex::Vector<amrex::MultiFab>  mu_old;
amrex::Vector<amrex::MultiFab>  mu_new;
amrex::Vector<amrex::MultiFab>  comp_old;
amrex::Vector<amrex::MultiFab>  comp_new;
amrex::Vector<amrex::MultiFab>  term1;
amrex::Vector<amrex::MultiFab>  term2;
amrex::Vector<amrex::MultiFab>  term3;
amrex::Vector<amrex::MultiFab>  term4;
amrex::Vector<amrex::MultiFab>  psi;
amrex::Vector<amrex::MultiFab>  lambad;
amrex::Vector<amrex::MultiFab>  plot;

// Store boundary condition
amrex::Vector<amrex::BCRec> bcs;
amrex::Vector<amrex::BCRec> bcs_phi;
amrex::Vector<amrex::BCRec> bcs_mu; 
amrex::Vector<amrex::BCRec> bcs_comp; 
amrex::Vector<amrex::BCRec> bcs_disp; 
amrex::Vector<amrex::BCRec> bcs_str; 

//AMR Paraneter
amrex::Vector<amrex::Real> dtlev;
amrex::Vector<int> istep;
amrex::Vector<int> nsubsteps;    
amrex::Real dt;
amrex::Real time;
int regrid_int ;
int do_subcycle;
int max_step ;

// Plotting file
std::string plot_file {"plt"};
int plot_int ;

// General Parameter of Model
//#ifndef phasecount
#define phasecount 7
//#endif
//#ifndef compcount
#define compcount 4
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
#define cent 0
#endif


#ifndef iph
#define iph 1
#endif

#ifndef imh
#define imh 2
#endif

#ifndef jph
#define jph 3
#endif

#ifndef jmh
#define jmh 4
#endif

#ifndef kph
#define kph 5
#endif

#ifndef kmh
#define kmh 6
#endif

#ifndef xx
#define xx 0
#endif

#ifndef yy
#define yy 1
#endif

#ifndef zz
#define zz 2
#endif

#ifndef yz
#define yz 3
#endif

#ifndef xz
#define xz 4
#endif

#ifndef xy
#define xy 5
#endif

#ifndef C11
#define C11 0
#endif

#ifndef C12
#define C12 1
#endif

#ifndef C44
#define C44 2
#endif


//Component and phase names
int dim;
int nump;
int numcom;
amrex::Vector <std::string> comp;
amrex::Vector <std::string> phase;

//Material properties
amrex::Vector <amrex::Real> gammaa;
amrex::Vector <amrex::Vector<amrex::Real>> gam;
amrex::Vector <std::string> val;
amrex::Vector <amrex::Real> dbl;
amrex::Vector <amrex::Vector<double>> diffu;
amrex::Vector<amrex::Vector <amrex::Vector<double>>> diff;

//Gas constant and molar volume
Real R;
Real Vm;

//Elastic parameters
int ELASTICITY{0};
amrex::Vector <amrex::Vector<double>> egstr;
amrex::Vector <amrex::Vector<double>> voigiso;
amrex::Vector <amrex::Vector<double>> stiffness_n;
amrex::Real ro;
amrex::Real dampfac;
amrex::Real MAX_ITERATIONS;
amrex::Real deltae=0.2;

//Boundary
amrex::Vector <amrex::Vector<std::string>> bound;

//Boundary value
amrex::Vector <amrex::Vector<std::string>> boundval;

//Type of simulation
int isothermal;
int binary;
int ternary;

//Ternary
int dilute;
amrex::Real T;

//Filewriting
std::string s_writeformat;
std::string writeformat;
std::string s_writehdf5;
int writehdf5;
std::string s_trackprog;
int trackprog;

//Model specific GP model
amrex::Real eps;
amrex::Real tau;
amrex::Real Tau;
amrex::Real tau_final; 

//Anisotropy
int funcANI;
int ANItype;
amrex::Vector<amrex::Real> dab;
amrex::Vector<amrex::Vector<amrex::Real>> dabb;

//Rotation matrix
amrex::Vector<amrex::Vector<amrex::Real>> rotmat;

//Potential function
int funcW;
std::string gamma_abc;
amrex::Vector<amrex::Real> gammaa_abc;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> gam_abc;

//Shifting of domain
amrex::Real shiftdom;
amrex::Real shiftj;

//Write composition and chemical ptential fields
int writecomp;

//Noise
int noise_pf;
amrex::Real amp_noise_phase;

//Temperature
amrex::Real Teq;
amrex::Real Tfill;

//Temperature gradient
amrex::Vector<amrex::Real> tempgrady;

//Function_F
int funcf;

//A
amrex::Vector<amrex::Vector<amrex::Real>> A1;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> A{};
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> Aeq{};

//B and D
amrex::Vector<amrex::Vector <amrex::Real>> B;
amrex::Vector <amrex::Real> C;
//Real BB;
//Real DD;

amrex::Vector<amrex::Vector<amrex::Vector<double>>> A_values;
amrex::Vector<amrex::Vector<double>> A_temp;

//ceq
amrex::Vector<amrex::Vector<amrex::Real>> ceq;
amrex::Vector <amrex::Vector<amrex::Vector<amrex::Real>>> c_eq;

//cfill
amrex::Vector<amrex::Vector<amrex::Real>> cfill;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> c_fill;

//cguess
amrex::Vector<amrex::Vector<amrex::Real>> cguess;
amrex::Vector <amrex::Vector<amrex::Vector<amrex::Real>>> c_guess;

//slopes
amrex::Vector<amrex::Vector<amrex::Real>> slopes;

//thermo phase
int ntp;

//tdbfname
std::string tdbname;

//tdb phase
std::string s_tdbphase;
Vector<std::string> tdb_phase;

//phase map
Vector<std::string> phasemap;

//Filling 
amrex::Vector <amrex::Vector<amrex::Real>> cube;
amrex::Vector <amrex::Vector<amrex::Real>> sphere;
amrex::Vector <amrex::Vector<amrex::Real>> cylinder;
amrex::Vector <amrex::Vector<amrex::Real>> ellipse;
amrex::Vector <amrex::Vector<amrex::Real>> cylrand;
amrex::Vector <amrex::Vector<amrex::Real>> sphrand;
amrex::Vector<amrex::Real> cuberand;
amrex::Vector<amrex::Real> cubepat;

//GSL
gsl_interp_accel *A_accel_ptr;
gsl_spline *A_spline_ptr;

//dcdmu
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> dcdmu;
amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> dcdmu_eq;

//cmu
Array2D<Real,0,phasecount,0,compcount-1> c{};
Array3D<Real,0,phasecount,0,compcount-1,0,compcount-1> cmu_eq{};

// Genral variables
int n=0;
amrex::Vector <int> numb;
int val1=0;
int val2=0;
amrex::Vector<Real> fe;
double strt_time;
double stop_time;
double rst_time=0.0;

amrex::Vector<amrex::Vector<amrex::Real>> conc;
amrex::Vector<amrex::Vector<amrex::Real>> conceq;
amrex::Vector<amrex::Vector<amrex::Vector<double>>> conc_Sol;
amrex::Vector<amrex::Vector<amrex::Vector<double>>> conc_Liq;
amrex::Vector<amrex::Vector<double>> temprt;
amrex::Vector<amrex::Vector<amrex::Real>> tau_ab;

double maxt1, maxt2, maxt3, maxt4;
double mint1, mint2, mint3, mint4;
double te1;
double te2;
double te3;
double te4;
double aniso_max1, aniso_max2;
double aniso_min1, aniso_min2;

void mat_mul1D(amrex::Vector<amrex::Vector<amrex::Real>> &inv_dcdmu, amrex::Vector<amrex::Real> &deltac,amrex::Vector<amrex::Real> &deltamu, int sz);
void mat_mul2D(amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> &dc_dmu,amrex::Vector<amrex::Vector<amrex::Vector<double>>> &diffu,amrex::Vector<amrex::Vector<amrex::Real>> &prod, int sz);
void mat_inv(amrex::Array2D<amrex::Real,0,compcount-2,0,compcount-2, Order::C> &denom,amrex::Array2D<Real,0,compcount-2,0,compcount-2, Order::C> &inv_denom,int numcomp);
void Inv_vec_rot(int i, int j, int k, int a, int b,amrex::Array2D <amrex::Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec,amrex::Array2D <amrex::Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z);
void Vec_rot(int i, int j, int k, int a, int b,amrex::Array1D <amrex::Real, 0, AMREX_SPACEDIM-1> &vec,amrex::Array1D <amrex::Real, 0, AMREX_SPACEDIM-1> &r_vec,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z);
void Vec_rot(int i, int j, int k, int a, int b,amrex::Array2D<amrex::Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec,amrex::Array2D<amrex::Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,amrex::Array2D <amrex::Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z);
};

#endif




































                        
