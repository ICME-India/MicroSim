#include "updatePhi.cuh"

__global__
void __updatePhi__(double **phi, double **dfdphi, double **phiNew,
                   double *relaxCoeff,
                   int NUMPHASES, int NUMCOMPONENTS,
                   int sizeX, int sizeY, int sizeZ,
                   double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int phase, p;

    //if ((i < sizeX && j > 1 && j < sizeY-2 && sizeZ == 1) || (i < sizeX && j < sizeY && sizeZ > 1 && k > 1 && k < sizeZ-2))
    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        for (phase = 0; phase < NUMPHASES; phase++)
        {
            phiNew[phase][idx] = phi[phase][idx];
            for (p = 0; p < NUMPHASES; p++)
            {
                if (p == phase)
                    continue;

                phiNew[phase][idx] -= (DELTA_t*relaxCoeff[phase*NUMPHASES + p]/(double)NUMPHASES)*(dfdphi[phase][idx] - dfdphi[p][idx]);
            }
        }
    }
    __syncthreads();
}

__global__
void __updatePhiBinary__(double **phi, double **dfdphi, double **phiNew,
                         double *relaxCoeff, double kappaPhi,
                         int NUMPHASES, int NUMCOMPONENTS,
                         int sizeX, int sizeY, int sizeZ,
                         double DELTA_X, double DELTA_Y, double DELTA_Z,
                         double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp[2], xm[2], yp[2], ym[2], zp[2], zm[2];

    //if ((i < sizeX && j > 1 && j < sizeY-2 && sizeZ == 1) || (i < sizeX && j < sizeY && sizeZ > 1 && k > 1 && k < sizeZ-2))
    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp[0] = (j + k*sizeY)*sizeX + i+1;
        xp[1] = (j + k*sizeY)*sizeX + i+2;
        xm[0] = (j + k*sizeY)*sizeX + i-1;
        xm[1] = (j + k*sizeY)*sizeX + i-2;

        if (i == 0)
        {
            xm[0] = (j + k*sizeY)*sizeX + sizeX-1;
            xm[1] = (j + k*sizeY)*sizeX + sizeX-2;
        }
        else if (i == 1)
        {
            xm[1] = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 2)
        {
            xp[1] = (j + k*sizeY)*sizeX;
        }
        else if (i == sizeX - 1)
        {
            xp[0] = (j + k*sizeY)*sizeX;
            xp[1] = (j + k*sizeY)*sizeX + 1;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp[0] = (j+1 + k*sizeY)*sizeX + i;
            yp[1] = (j+2 + k*sizeY)*sizeX + i;
            ym[0] = (j-1 + k*sizeY)*sizeX + i;
            ym[1] = (j-2 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym[0] = (sizeY-1 + k*sizeY)*sizeX + i;
                ym[1] = (sizeY-2 + k*sizeY)*sizeX + i;
            }
            else if (j == 1)
            {
                ym[1] = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 2)
            {
                yp[1] = (k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp[0] = (k*sizeY)*sizeX + i;
                yp[1] = (1 + k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp[0] = (j + (k+1)*sizeY)*sizeX + i;
            zp[1] = (j + (k+2)*sizeY)*sizeX + i;
            zm[0] = (j + (k-1)*sizeY)*sizeX + i;
            zm[1] = (j + (k-2)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm[0] = (j + (sizeZ-1)*sizeY)*sizeX + i;
                zm[1] = (j + (sizeZ-2)*sizeY)*sizeX + i;
            }
            else if (k == 1)
            {
                zm[1] = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 2)
            {
                zp[1] = j*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp[0] = j*sizeX + i;
                zp[1] = (j + sizeY)*sizeX + i;
            }
        }

        dfdphi[1][idx] -= 2.0*kappaPhi*(phi[1][xp[0]] - 2.0*phi[1][idx] + phi[1][xm[0]])/(DELTA_X*DELTA_X);
        if (sizeY > 1)
            dfdphi[1][idx] -= 2.0*kappaPhi*(phi[1][yp[0]] - 2.0*phi[1][idx] + phi[1][ym[0]])/(DELTA_Y*DELTA_Y);
        if (sizeZ > 1)
            dfdphi[1][idx] -= 2.0*kappaPhi*(phi[1][zp[0]] - 2.0*phi[1][idx] + phi[1][zm[0]])/(DELTA_Z*DELTA_Z);

        phiNew[1][idx] = phi[1][idx] - (DELTA_t*relaxCoeff[1])*dfdphi[1][idx];

        phiNew[0][idx] = 1.0 - phiNew[1][idx];

    }
    __syncthreads();
}

void updatePhi(double **phi, double **dfdphi, double **phiNew,
               domainInfo* simDomain, controls* simControls,
               simParameters* simParams, subdomainInfo* subdomain,
               dim3 gridSize, dim3 blockSize)
{
    if (simControls->multiphase == 1 || simDomain->numPhases > 2 || simDomain->numComponents > 2)
    {
        __updatePhi__<<<gridSize, blockSize>>>(phi, dfdphi, phiNew,
                                               simParams->relax_coeff_dev,
                                               simDomain->numPhases, simDomain->numComponents,
                                               subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                               simControls->DELTA_t);
    }
    else if (simDomain->numPhases == 2 && simDomain->numComponents == 2)
    {
        __updatePhiBinary__<<<gridSize, blockSize>>>(phi, dfdphi, phiNew,
                                                     simParams->relax_coeff_dev, simParams->kappaPhi_host[0][1],
                                                     simDomain->numPhases, simDomain->numComponents,
                                                     subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                     simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                     simControls->DELTA_t);
    }
}
