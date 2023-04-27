#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include "defines.h"

/*
 * Definition for a structure to hold information about the domain.
 */
typedef struct domainInfo
{
    // Domain size
    long MESH_X;
    long MESH_Y;
    long MESH_Z;

    long numCells;

    // Cell size
    double DELTA_X;
    double DELTA_Y;
    double DELTA_Z;

    long DIMENSION;

    long numPhases;
    long numComponents;
    char **phaseNames;
    char **componentNames;

    long numThermoPhases;
    char **phases_tdb;
    char **phase_map;
    long  *thermo_phase_dev, *thermo_phase_host;
} domainInfo;

/*
 * Definition for a structure to boundary condition instructions
 */
typedef struct bc_scalars
{
    int type;
    long points[3];
    long proxy[3];
    double value[3];
} bc_scalars;

/*
 * Definition for a structure to hold simulation controls
 */
typedef struct controls
{
    // MPI information
    long numworkers;

    // Simulation parameters
    double DELTA_t;
    long startTime;
    long count;

    // Number of iterations for all of the following
    long numSteps;
    long saveInterval;
    long trackProgress;
    long restart;
    long nsmooth;

    // Extent of padding at the block boundaries
    long padding;

    // Type of free-energy
    long FUNCTION_F;

    // MPMC or binary
    long multiphase;

    // Filewriting
    long writeFormat;
    long writeHDF5;

    // Anisotropy
    int FUNCTION_ANISOTROPY;
    int FOLD;
    int ANISOTROPY;
    int ANISOTROPY_GRADIENT;
    double Rth_phase;

    // Temperature behaviour
    long ISOTHERMAL;
    double dTdt;
    long T_update;

    // Boundary conditions
    bc_scalars *boundary[6];
} controls;

typedef struct simParameters
{
    // User input
    double **gamma_host, *gamma_dev;
    double ***diffusivity_host, *diffusivity_dev;
    double ***mobility_host, *mobility_dev;
    double **Tau_host;
    double **relax_coeff_host, *relax_coeff_dev;
    double ***F0_A_host, *F0_A_dev;
    double **F0_B_host, *F0_B_dev;
    double **F0_Beq_host, *F0_Beq_dev;
    double *F0_C_host, *F0_C_dev;

    long ISOTHERMAL;

    double **DELTA_T;
    double **DELTA_C;
    double ***dcbdT;
    double **dBbdT;

    double alpha;
    double epsilon;
    double ***ceq_host, *ceq_dev;
    double ***cfill_host, *cfill_dev;
    double ***cguess_host, *cguess_dev;

    double R;
    double molarVolume;
    double T, Teq, Tfill;

    double ***slopes;

    double *theta_i_host, *theta_i_dev;
    double **theta_ij_host, *theta_ij_dev;
    double ***theta_ijk_host, *theta_ijk_dev;

    long SEED;

    double **dab_host, *dab_dev;
    double tilt_angle;
    double Rtheta;
    double ****Rotation_matrix_host, *Rotation_matrix_dev;
    double ****Inv_Rotation_matrix_host, *Inv_Rotation_matrix_dev;

    // Calculated from user input
    double **kappaPhi_host, *kappaPhi_dev;
} simParameters;

/*
 * Definition for a structure to hold information about the domain that is
 * stored in any particular MPI process
 */
typedef struct subdomainInfo
{
    long xS, yS, zS;
    long xE, yE, zE;

    long sizeX, sizeY, sizeZ;

    long numCells;
    long numCompCells;
    long shiftPointer;
    long padding;
    long yStep, zStep;

    int rank;

    int nbLeft, nbRight;        //Along x-axis
    int nbUp, nbDown;           //Along y-axis
    int nbFront, nbBack;        //Along z-axis
} subdomainInfo;

/*
 * Definition for a structure to hold sphere filling information
 */
typedef struct sphere
{
    long xC, yC, zC;
    long radius;
    long phase;
} sphere;

/*
 * Definition for a structure to hold cylinder filling information
 */
typedef struct cylinder
{
    long xC, yC;
    long zS, zE;
    long radius;
    long phase;
} cylinder;

/*
 * Definition for a structure to hold cube filling information
 */
typedef struct cube
{
    long xS, yS, zS;
    long xE, yE, zE;
    long phase;
} cube;

/*
 * Definition for a structure to hold ellipse filling information
 */
typedef struct ellipse
{
    long xC, yC, zC;
    double major_axis;
    double eccentricity;
    double rot_angle;
    long phase;
} ellipse;

enum fill{FILLCYLINDER, FILLSPHERE, FILLCYLINDERRANDOM, FILLSPHERERANDOM, FILLCUBE, FILLELLIPSE};

/*
 * Definition for a structure to hold all domain filling information
 */
typedef struct fillParameters
{
    enum fill *fillType;
    long *xC, *yC, *zC;
    long *xS, *yS, *zS;
    long *xE, *yE, *zE;
    long *radius;
    long *phase;
    double *major_axis, *eccentricity, *rot_angle;
    long *seed;
    double *volFrac;
    long *shieldDist;
    double *radVar;

    long countFill;
} fillParameters;

#endif
