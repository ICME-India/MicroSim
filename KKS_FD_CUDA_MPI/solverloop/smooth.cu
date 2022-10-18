#include "smooth.cuh"

__global__
void __smooth__(double **phi, double **phiNew,
                double *relaxCoeff, double *kappaPhi,
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
    double aniso[MAX_NUM_PHASES] = {0.0};

    long phase, p;

    double dfdphiSum = 0.0;

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

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            dfdphiSum = 0.0;

            for (p = 0; p < NUMPHASES; p++)
            {
                if (p == phase)
                    continue;

                dfdphiSum += phiAniso[phase][1][1][1]*(4.0*phiAniso[phase][1][1][1]*phiAniso[phase][1][1][1] - 6.0*phiAniso[phase][1][1][1] + 2.0);
                dfdphiSum -= phiAniso[p][1][1][1]*(4.0*phiAniso[p][1][1][1]*phiAniso[p][1][1][1] - 6.0*phiAniso[p][1][1][1] + 2.0);

                dfdphiSum += kappaPhi[phase*NUMPHASES + p]*(aniso[p] - aniso[phase]);

            }

            phiNew[phase][index[1][1][1]] = phi[phase][index[1][1][1]] - DELTA_t*2.0*FunctionTau(phi, relaxCoeff, index[1][1][1], NUMPHASES)*dfdphiSum/(double)NUMPHASES;
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
                                        simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                        subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                        subdomain->yStep, subdomain->zStep, subdomain->padding,
                                        simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                        simControls->DELTA_t);
}
