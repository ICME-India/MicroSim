#include "updatePhi.cuh"

__global__
void __updatePhi__(double **phi, double **dfdphi, double **phiNew,
                   double *relaxCoeff, double *kappaPhi,
                   long NUMPHASES, long NUMCOMPONENTS,
                   long sizeX, long sizeY, long sizeZ,
                   double DELTA_X, double DELTA_Y, double DELTA_Z,
                   double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + k*sizeY)*sizeX + i;

    long xp, xm, yp, ym, zp, zm;

    double laplacian_p = 0.0, laplacian_phase = 0.0;

    long phase, p;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp = (j + k*sizeY)*sizeX + i+1;
        xm = (j + k*sizeY)*sizeX + i-1;

        if (i == 0)
        {
            xm = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 1)
        {
            xp = (j + k*sizeY)*sizeX;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp = (j+1 + k*sizeY)*sizeX + i;
            ym = (j-1 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp = (k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp = (j + (k+1)*sizeY)*sizeX + i;
            zm = (j + (k-1)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp = j*sizeX + i;
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            phiNew[phase][idx] = phi[phase][idx];

            laplacian_phase = (phi[phase][xp] - 2.0*phi[phase][idx] + phi[phase][xm])/(DELTA_X*DELTA_X);
            if (sizeY > 1)
                laplacian_phase += (phi[phase][yp] - 2.0*phi[phase][idx] + phi[phase][ym])/(DELTA_Y*DELTA_Y);
            if (sizeZ > 1)
                laplacian_phase += (phi[phase][zp] - 2.0*phi[phase][idx] + phi[phase][zm])/(DELTA_Z*DELTA_Z);

            if (fabs(laplacian_phase) > 1e-12/(DELTA_X*DELTA_X))
            {
                for (p = 0; p < NUMPHASES; p++)
                {
                    if (p == phase)
                        continue;

                    laplacian_p = (phi[p][xp] - 2.0*phi[p][idx] + phi[p][xm])/(DELTA_X*DELTA_X);
                    if (sizeY > 1)
                        laplacian_p += (phi[p][yp] - 2.0*phi[p][idx] + phi[p][ym])/(DELTA_Y*DELTA_Y);
                    if (sizeZ > 1)
                        laplacian_p += (phi[p][zp] - 2.0*phi[p][idx] + phi[p][zm])/(DELTA_Z*DELTA_Z);

                    phiNew[phase][idx] -= (DELTA_t*FunctionTau(phi, relaxCoeff, idx, NUMPHASES)/(double)NUMPHASES)*(dfdphi[phase][idx] - dfdphi[p][idx] - 2.0*kappaPhi[phase*NUMPHASES + p]*(laplacian_phase - laplacian_p));
                }
            }
        }
    }
    __syncthreads();
}

void updatePhi(double **phi, double **dfdphi, double **phiNew, double **phaseComp,
               domainInfo* simDomain, controls* simControls,
               simParameters* simParams, subdomainInfo* subdomain,
               dim3 gridSize, dim3 blockSize)
{
    __updatePhi__<<<gridSize, blockSize>>>(phi, dfdphi, phiNew,
                                           simParams->relax_coeff_dev, simParams->kappaPhi_dev,
                                           simDomain->numPhases, simDomain->numComponents,
                                           subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                           simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                           simControls->DELTA_t);
}
