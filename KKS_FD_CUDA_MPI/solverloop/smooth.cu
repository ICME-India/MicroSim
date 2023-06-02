#include "smooth.cuh"

#define phiAniso(phase, x, y, z) (phiAniso[(((phase)*3 + (x))*3 + (y))*3 + (z)])

__global__
void __smooth__(double **phi, double **phiNew,
                double *relaxCoeff, double *kappaPhi,
                double *dab, double *Rotation_matrix, double *Inv_rotation_matrix, int FUNCTION_ANISOTROPY,
                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                long xStep, long yStep, long padding,
                double DELTA_X, double DELTA_Y, double DELTA_Z,
                double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long x, y, z;

    long index[3][3][3];

    double phiAniso[MAX_NUM_PHASES*27];

    double aniso[MAX_NUM_PHASES] = {0.0};
    double dfdphiSum = 0.0;

    long phase, p;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        /*
         * index is a 3D matrix that stores location indices for the ...
         * 26 gridpoints around the gridpoint at which computation is being carried out.
         *
         * Eg:
         * index[0][1][2] is (i-1, j, k+1)
         *
         * for 2D simulations, the third index will be set to 1 (index[x][y][1])
         *
         */

        for (x = 0; x < 3; x++)
        {
            for (y = 0; y < 3; y++)
            {
                for (z = 0; z < 3; z++)
                {
                    index[x][y][z] = (k+z-1) + (j+y-1)*yStep + (i+x-1)*xStep;
                }
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            for (x = 0; x < 3; x++)
            {
                for (y = 0; y < 3; y++)
                {
                    if (DIMENSION == 3)
                    {
                        for (z = 0; z < 3; z++)
                        {
                            phiAniso(phase, x, y, z) = phi[phase][index[x][y][z]];
                        }
                    }
                    else
                    {
                        phiAniso(phase, x, y, 1) = phi[phase][index[x][y][1]];
                    }
                }
            }
        }


        if (FUNCTION_ANISOTROPY == 0)
        {
            for (phase = 0; phase < NUMPHASES; phase++)
            {
                if (DIMENSION == 1)
                {
                    aniso[phase] = (phiAniso(phase, 0, 1, 1) - 2.0*phiAniso(phase, 1, 1, 1) + phiAniso(phase, 2, 1, 1))/(DELTA_X*DELTA_X);
                }
                else if (DIMENSION == 2)
                {
                    //Centre
                    aniso[phase] = -3.0*phiAniso(phase, 1, 1, 1)/(DELTA_X*DELTA_Y);

                    // Nearest neighbours
                    aniso[phase] += 0.5*(phiAniso(phase, 0, 1, 1) + phiAniso(phase, 2, 1, 1))/(DELTA_X*DELTA_X);
                    aniso[phase] += 0.5*(phiAniso(phase, 1, 0, 1) + phiAniso(phase, 1, 2, 1))/(DELTA_Y*DELTA_Y);

                    // Second-nearest neighbours
                    aniso[phase] += 0.25*(phiAniso(phase, 0, 0, 1) + phiAniso(phase, 0, 2, 1) + phiAniso(phase, 2, 2, 1) + phiAniso(phase, 2, 0, 1))/(DELTA_X*DELTA_Y);
                }
                else if (DIMENSION == 3)
                {
                    // Centre
                    aniso[phase] = -4.0*phiAniso(phase, 1, 1, 1)/(DELTA_X*DELTA_X);

                    // Nearest neighbours
                    aniso[phase] += (phiAniso(phase, 0, 1, 1) + phiAniso(phase, 2, 1, 1))/(3.0*DELTA_X*DELTA_X);
                    aniso[phase] += (phiAniso(phase, 1, 0, 1) + phiAniso(phase, 1, 2, 1))/(3.0*DELTA_Y*DELTA_Y);
                    aniso[phase] += (phiAniso(phase, 1, 1, 0) + phiAniso(phase, 1, 1, 2))/(3.0*DELTA_Z*DELTA_Z);

                    // Second-nearest neighbours
                    aniso[phase] += (phiAniso(phase, 0, 0, 1) + phiAniso(phase, 0, 2, 1) + phiAniso(phase, 2, 2, 1) + phiAniso(phase, 2, 0, 1))/(6.0*DELTA_X*DELTA_Y);
                    aniso[phase] += (phiAniso(phase, 1, 0, 0) + phiAniso(phase, 1, 0, 2) + phiAniso(phase, 1, 2, 2) + phiAniso(phase, 1, 2, 0))/(6.0*DELTA_Y*DELTA_Z);
                    aniso[phase] += (phiAniso(phase, 0, 1, 0) + phiAniso(phase, 0, 1, 2) + phiAniso(phase, 2, 1, 2) + phiAniso(phase, 2, 1, 0))/(6.0*DELTA_Z*DELTA_X);
                }
            }
        }
        else if (FUNCTION_ANISOTROPY == 1 || FUNCTION_ANISOTROPY == 2)
        {
            for (phase = 0; phase < NUMPHASES; phase++)
            {
                aniso[phase] = calcAnisotropy_01(phiAniso, dab, kappaPhi, Rotation_matrix, Inv_rotation_matrix, phase, NUMPHASES, DIMENSION, DELTA_X, DELTA_Y, DELTA_Z);
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            dfdphiSum = 0.0;

            for (p = 0; p < NUMPHASES; p++)
            {
                if (p == phase)
                    continue;

                if (FUNCTION_ANISOTROPY == 0)
                {
                    dfdphiSum += 2.0*kappaPhi[phase*NUMPHASES + p]*(aniso[p] - aniso[phase]);
                }
                else if (FUNCTION_ANISOTROPY == 1 || FUNCTION_ANISOTROPY == 2)
                {
                    dfdphiSum += 2.0*(aniso[p] - aniso[phase]);
                }
            }

            phiNew[phase][index[1][1][1]] = phi[phase][index[1][1][1]] - DELTA_t*FunctionTau(phi, relaxCoeff, index[1][1][1], NUMPHASES)*dfdphiSum/(double)NUMPHASES;
        }
    }
}


void smooth(double **phi, double **phiNew,
            domainInfo* simDomain, controls* simControls,
            simParameters* simParams, subdomainInfo* subdomain,
            dim3 gridSize, dim3 blockSize)
{
    __smooth__<<<gridSize, blockSize>>>(phi, phiNew,
                                        simParams->relax_coeff_dev, simParams->kappaPhi_dev,
                                        simParams->dab_dev, simParams->Rotation_matrix_dev, simParams->Inv_Rotation_matrix_dev, 0,
                                        simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                        subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                        subdomain->xStep, subdomain->yStep, subdomain->padding,
                                        simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                        simControls->DELTA_t);
}
