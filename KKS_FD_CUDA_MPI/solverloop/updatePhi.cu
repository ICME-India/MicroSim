#include "updatePhi.cuh"

__global__
void __updatePhi__(double **phi, double **dfdphi, double **phiNew,
                   double *relaxCoeff, double *kappaPhi,
                   double *dab, double *Rotation_matrix, double *Inv_rotation_matrix, int FUNCTION_ANISOTROPY,
                   long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ,
                   long yStep, long zStep, long padding,
                   double DELTA_X, double DELTA_Y, double DELTA_Z,
                   double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long x, y, z;

    long index[3][3][3];

    double phiAniso[MAX_NUM_PHASES][3][3][3];

    double divphi = 0.0;

    double aniso[MAX_NUM_PHASES] = {0.0};
    double dfdphiSum = 0.0;

    long phase, p;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        for (x = 0; x < 3; x++)
        {
            for (y = 0; y < 3; y++)
            {
                for (z = 0; z < 3; z++)
                {
                    index[x][y][z] = (k+z-1)*zStep + (j+y-1)*yStep + (i+x-1);
                }
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            for (x = 0; x < 3; x++)
            {
                for (y = 0; y < 3; y++)
                {
                    for (z = 0; z < 3; z++)
                    {
                        phiAniso[phase][x][y][z] = phi[phase][index[x][y][z]];
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
                    aniso[phase] = (phiAniso[phase][0][1][1] - 2.0*phiAniso[phase][1][1][1] + phiAniso[phase][2][1][1])/(DELTA_X*DELTA_X);
                }
                else if (DIMENSION == 2)
                {
                    // Centre
                    aniso[phase] = -3.0*phiAniso[phase][1][1][1]/(DELTA_X*DELTA_Y);

                    // Nearest neighbours
                    aniso[phase] += 0.5*(phiAniso[phase][0][1][1] + phiAniso[phase][2][1][1])/(DELTA_X*DELTA_X);
                    aniso[phase] += 0.5*(phiAniso[phase][1][0][1] + phiAniso[phase][1][2][1])/(DELTA_Y*DELTA_Y);

                    // Second-nearest neighbours
                    aniso[phase] += 0.25*(phiAniso[phase][0][0][1] + phiAniso[phase][0][2][1] + phiAniso[phase][2][2][1] + phiAniso[phase][2][0][1])/(DELTA_X*DELTA_Y);
                }
                else if (DIMENSION == 3)
                {
                    aniso[phase] = 0.0;
                }
            }
        }
        else if (FUNCTION_ANISOTROPY == 1)
        {
            for (phase = 0; phase < NUMPHASES; phase++)
            {
                aniso[phase] = calcAnisotropy_01(phiAniso, dab, kappaPhi, Rotation_matrix, Inv_rotation_matrix, phase, NUMPHASES, DIMENSION, DELTA_X, DELTA_Y, DELTA_Z);
            }
        }
        else if (FUNCTION_ANISOTROPY == 2)
        {
            for (phase = 0; phase < NUMPHASES; phase++)
            {
                aniso[phase] = calcAnisotropy_02(phiAniso, dab, kappaPhi, Rotation_matrix, Inv_rotation_matrix, phase, NUMPHASES, DIMENSION, DELTA_X, DELTA_Y, DELTA_Z);
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            divphi = (phiAniso[phase][2][1][1] - 2.0*phiAniso[phase][1][1][1] + phiAniso[phase][0][1][1])/(DELTA_X*DELTA_X);
            if (DIMENSION >= 2)
                divphi += (phiAniso[phase][1][2][1] - 2.0*phiAniso[phase][1][1][1] + phiAniso[phase][1][0][1])/(DELTA_Y*DELTA_Y);
            if (DIMENSION == 3)
                divphi += (phiAniso[phase][1][1][2] - 2.0*phiAniso[phase][1][1][1] + phiAniso[phase][1][1][0])/(DELTA_Z*DELTA_Z);

            dfdphiSum = 0.0;

            for (p = 0; p < NUMPHASES; p++)
            {
                if (p == phase)
                    continue;

                dfdphiSum += (dfdphi[phase][index[1][1][1]] - dfdphi[p][index[1][1][1]]);

                if (FUNCTION_ANISOTROPY == 0)
                {
                    dfdphiSum += 2.0*kappaPhi[phase*NUMPHASES + p]*(aniso[p] - aniso[phase]);
                }
                else if (FUNCTION_ANISOTROPY == 1 || FUNCTION_ANISOTROPY == 2)
                {
                    dfdphiSum += (aniso[p] - aniso[phase]);
                }
            }

            if (fabs(divphi) > 0.0)
                phiNew[phase][index[1][1][1]] = phi[phase][index[1][1][1]] - DELTA_t*FunctionTau(phi, relaxCoeff, index[1][1][1], NUMPHASES)*dfdphiSum/(double)NUMPHASES;
        }
    }
}

void updatePhi(double **phi, double **dfdphi, double **phiNew, double **phaseComp,
               domainInfo* simDomain, controls* simControls,
               simParameters* simParams, subdomainInfo* subdomain,
               dim3 gridSize, dim3 blockSize)
{
    __updatePhi__<<<gridSize, blockSize>>>(phi, dfdphi, phiNew,
                                           simParams->relax_coeff_dev, simParams->kappaPhi_dev,
                                           simParams->dab_dev, simParams->Rotation_matrix_dev, simParams->Inv_Rotation_matrix_dev, simControls->FUNCTION_ANISOTROPY,
                                           simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                           subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                           subdomain->yStep, subdomain->zStep, subdomain->padding,
                                           simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                           simControls->DELTA_t);

    applyBoundaryCondition(phiNew, 0, simDomain->numPhases,
                           simDomain, simControls,
                           simParams, subdomain,
                           gridSize, blockSize);
}
