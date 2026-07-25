#pragma once
#include <vector>
#include <string>


#define MAX_PHASES 5
#define MAX_COMP 15

#define Alpha 0
#define Beta 1

#define  PI    3.14159265358979323846  /* PI     */
#define INV_PI 0.31830988618379067154  /* 1 / PI */

#define XX 0
#define YY 1
#define ZZ 2
#define YZ 3
#define XZ 4
#define XY 5

#define IDX_CC 0
#define IDX_CE 1
#define IDX_WC 2
#define IDX_CN 3
#define IDX_SC 4
#define IDX_CF 5
#define IDX_BC 6

using std::vector, std::string;

struct Domain{
    int DIMENSION;
    int MESH_X;
    int MESH_Y;
    int MESH_Z;
    double DELTA_X;
    double DELTA_Y;
    double DELTA_Z;
    double DELTA_T;
};

struct Input_parameters{
    int NUM_PHASES;
    int NUM_COMPONENTS;

    int Nsmooth;
    int start_time;
    int ntime_steps;
    int saveT;
    int Num_thermo_phase;
    double epsilon;
    double tau;
    double R;
    double Volm;
    int Obstacle;
    int Function_W;
    int function_anisotropy;
    int FOLD;
    double Noise_phase_field;
    double T_eq;
    double Filling_temp;
    double T;
    
    // KKS parameters
    double refD;
    double Er;
    double epsc;
    double epsm;
    double W;
    double ee;
    double IntMobInv;
    double a2;
};

struct Settings{
    int periodic[3] = {1, 1, 1}; // default MPI wraps around the cart in all the three axis
    bool BINARY;
    bool TERNARY;
    int Function_F;
    bool restart = false;
    bool ASCII = false;
    bool HDF5 = true;
    bool elasticity = false;
    bool ISOTHERMAL =true;  
    bool Dilute = true;
    bool shift = false;
    bool Noise = false;
    bool anisotropy = true;
};

struct phases_components{
    vector<string> componenets;
    vector<string> phases;
    vector<string> phases_tdb;
    vector<string> phases_map;
};

struct Temp_gradient{
    double base_Temp;
    double Distance;
    double gradient_offset;
    double DeltaT;
    double offset;
    double velocity;
};


struct Phasefield_matraix{
    double *Diffusivity;
    double *Inv_Diffusivity;
    double *ceq;
    double *cfill;
    double *c_guess;
    double *ceq_coeffs;
    double *slopes;
    double *dcbdT;
    double *A;
    double *B;
    double *Beq;
    double *C;
    double *DELTA_T;
    double *DELTA_C;
    double *dcbdT_phase;
    double *dBbdT;
    double *cmu;
    double *muc;
    double *Rotation_matrix;
    double *Inv_Rotation_matrix;
    double *Rotated_qab;
    double *gamma;
    double *tau_ab;
    double *gamma_abc;
    double *dab;
    double *fab;
};

// -------------- END of Input parameters -----------------//

enum field{
    Order_Parameter,
    Chemical_potential,
    Composition,
    Temperature,
    Unknown
};

struct fields {
  // --- Shared Fields ---
  double *phia;        // Evolved: Allen-Cahn
  
  // --- Method Specifics ---
#ifdef __KKS
  // KKS: We evolve C, but need Cs/Cl for the solver history
  double *C;    // Evolved: Cahn-Hilliard
  double *C_phase;   // Auxiliary: For Newton-Raphson initial guess (stores concentrations for all phases)
  double *mu;

#elif defined(__FID)
  double *C;
  double *Cs;   // Evolved: Diffusion + Solute trapping
  double *Cl;   // Evolved: Diffusion + Solute trapping
  double *mus;
  double *mul;

#else  // GP - Default
  double *C;    
  double *mu;
  double *flux_x;
  double *flux_y;
  double *flux_z;
#endif
  double *deltaPhi;
  double *temperature;
};

struct OFFSETS{
    int NUM_FIELDS;

    int OFF_PHI;

#ifdef __KKS
    int OFF_COMP;
    int OFF_C_PHASE;
    int OFF_MU;

#elif defined(__FID)
    int OFF_COMP;
    int OFF_CS;
    int OFF_CL;
    int OFF_MUS;
    int OFF_MUL;

#else
    int OFF_COMP;
    int OFF_MU;
    int OFF_TERM1;
    int OFF_TERM2;
    int OFF_TERM3;
#endif
    int OFF_DELTAPHI;
    int OFF_TEMP;
};

struct FieldRange{
    int start_idx;
    int count;
};

// ---------------- Filling parameters --------------------//

enum fill {
    FILL_CUBE,
    FILL_CYLINDER,
    FILL_SPHERE,
    FILL_CYLINDER_RANDOM,
    FILL_SPHERE_RANDOM,
    FILL_CUBE_RANDOM,
    FILL_ELLIPSE,
    FILL_VORONOI2D,
    FILL_VORONOI3D
};

struct fillParameters{
    enum fill *fillType;
    int *xC, *yC, *zC;
    int *xS, *yS, *zS;
    int *xE, *yE, *zE;
    int *radius;
    int *phase;
    double *major_axis, *eccentricity, *rot_angle;
    int *seed;
    double *volFrac;
    int *shieldDist;
    int *shiftFrac;
    double *radVar;
    int countFill;
};

// ---------------- Boundary conditions ------------------//

enum BCType {
    //BC_PERIODIC = 0, // Handled by MPI topology
    //BC_NEUMANN  = 1, // Zero Flux (copy inner)
    //BC_DIRICHLET= 2, // Fixed Value
    //BC_FLUX     = 3  // Fixed Gradient (Inner + val)
    BC_FREE     = 0, // Don't know what it is but it was part of microsim
    BC_NEUMANN  = 1,
    BC_DIRICHLET= 2,
    BC_PERIODIC = 3,
    BC_FLUX     = 4
};

struct BCValue {
    BCType type;
    double value;    // Value for Dirichlet or Flux amount
};

// Container for BCs of ALL fields on a SINGLE face
struct FaceBCs {
    BCValue *phi;
    BCValue *comp;
    BCValue *mu;
    BCValue *temp;
};

// Container for the whole domain
struct DomainBCs {
    FaceBCs left, right;
    FaceBCs top, bottom;
    FaceBCs front, back;
};

struct DomainDecomp{
    int Dims[3]    = {0,0,0};
    int periods[3] = {0,0,0};
    int coords[3]  = {0,0,0};
};

// ------------------ MPI Decomposition ---------------//
struct workers {
    int left;
    int right;
    int up;
    int down;
    int front;
    int back;

    bool is_left_edge;
    bool is_right_edge;
    bool is_top_edge;
    bool is_bot_edge;
    bool is_front_edge;
    bool is_back_edge;

    int local_d;
    int local_h;
    int local_w;

    int ghost_d;
    int ghost_h;
    int ghost_w;

    size_t num_pts;
    size_t sz_yz;
    size_t sz_xz;
    size_t sz_xy;

    size_t off_l;
    size_t off_r;
    size_t off_t;
    size_t off_b;
    size_t off_f;
    size_t off_bk;

    size_t total_send_bytes;
};

struct SplineCoefficients {
    double *c0;  // constant term
    double *c1;  // linear term
    double *c2;  // quadratic term
    double *c3;  // cubic term
    double *x;   // breakpoints
    int n_points;
    
    // Default constructor
    SplineCoefficients() : c0(nullptr), c1(nullptr), c2(nullptr), 
                           c3(nullptr), x(nullptr), n_points(0) {}
};

