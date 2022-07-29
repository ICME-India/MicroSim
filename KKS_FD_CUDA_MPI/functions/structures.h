#ifndef STRUCTURES_H_
#define STRUCTURES_H_

/*
 * Definition for a structure to hold information about the domain.
 */
typedef struct domainInfo
{
// Domain size
    int MESH_X;
    int MESH_Y;
    int MESH_Z;

    int numCells;

// Cell size
    double DELTA_X;
    double DELTA_Y;
    double DELTA_Z;

    int DIMENSION;

    int numPhases;
    int numComponents;
    char **phaseNames;
    char **componentNames;

    int numThermoPhases;
    char **phases_tdb;
    char **phase_map;
    int  *thermo_phase_dev, *thermo_phase_host;
} domainInfo;

/*
 * Definition for a structure to hold simulation controls
 */
typedef struct controls
{
// Simulation parameters
    double DELTA_t;
    int startTime;
    int count;

// Number of iterations for all of the following
    int numSteps;
    int saveInterval;
    int trackProgress;
    int restart;

// Extent of padding at the block boundaries
    int padding;

// Type of free-energy
    int FUNCTION_F;

// MPMC or binary
    int multiphase;

// Filewriting
    int writeFormat;
    int writeHDF5;

// Temperature behaviour
    int ISOTHERMAL;
    double dTdt;
    int T_update;

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

// Calculated from user input
    double **kappaPhi_host, *kappaPhi_dev;
} simParameters;

/*
 * Definition for a structure to hold information about the domain that is
 * stored in any particular MPI process
 */
typedef struct subdomainInfo
{
    int xS, yS, zS;
    int xE, yE, zE;

    int sizeX, sizeY, sizeZ;

    int numCells;
    int numCompCells;
    int shiftPointer;
    int padding;

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
    int xC, yC, zC;
    int radius;
    int phase;
} sphere;

/*
 * Definition for a structure to hold cylinder filling information
 */
typedef struct cylinder
{
    int xC, yC;
    int zS, zE;
    int radius;
    int phase;
} cylinder;

enum fill{FILLCYLINDER, FILLSPHERE, FILLCYLINDERRANDOM, FILLSPHERERANDOM};

/*
 * Definition for a structure to hold all domain filling information
 */
typedef struct fillParameters
{
    enum fill *fillType;
    int *xC, *yC, *zC;
    int *zS, *zE;
    int *radius;
    int *phase;
    long *seed;
    double *volFrac;
    int *shieldDist;
    double *radVar;

    int countFill;
} fillParameters;

#endif
