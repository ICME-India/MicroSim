#include "boundary.cuh"

__global__
void applyNeumann(double **field,
                  long face, long numFields, long DIMENSION,
                  long sizeX, long sizeY, long sizeZ,
                  long yStep, long zStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + k*sizeY)*sizeX + i, idx2 = -1;

    if ((i < padding || i >= sizeX-padding || ((j < padding || j >= sizeY-padding) && DIMENSION >= 2) || ((k < padding || k >= sizeZ-padding) && DIMENSION == 3)) && i < sizeX && j < sizeY && k < sizeZ)
    {
        // X+
        if (face == 0 && i >= sizeX-padding)
        {
            idx2 = k*zStep + j*yStep + (2*(sizeX-padding) - i - 1);
        }

        // X-
        else if (face == 1 && i < padding)
        {
            idx2 = k*zStep + j*yStep + (2*padding - i - 1);
        }

        // Y+
        else if (face == 2 && j >= sizeY-padding && DIMENSION >= 2)
        {
            idx2 = k*zStep + (2*(sizeY-padding) - j - 1)*yStep + i;
        }

        // Y-
        else if (face == 3 && j < padding && DIMENSION >= 2)
        {
            idx2 = k*zStep + (2*padding - j - 1)*yStep + i;
        }

        // Z+
        else if (face == 4 && k >= sizeZ-padding && DIMENSION == 3)
        {
            idx2 = ((2*sizeZ - padding) - k - 1)*zStep + j*yStep + i;
        }

        // Z-
        else if (face == 5 && k < padding && DIMENSION == 3)
        {
            idx2 = (2*padding - k - 1)*zStep + j*yStep + i;
        }

        if (idx2 > -1)
        {
            for (long iter1 = 0; iter1 < numFields; iter1++)
            {
                field[iter1][idx] = field[iter1][idx2];
            }
        }
    }

    __syncthreads();
}

__global__
void applyPeriodic(double **field,
                   long face, long numFields, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ,
                   long yStep, long zStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = k*zStep + j*yStep + i, idx2 = -1;

    if ((i == 0 || i == sizeX-1 || ((j == 0 || j == sizeY-1) && DIMENSION >= 2) || ((k == 0 || k == sizeZ-1) && DIMENSION == 3)) && i < sizeX && j < sizeY && k < sizeZ)
    {
        // X+
        if (face == 0 && i == sizeX-1)
        {
            idx2 = k*zStep + j*yStep + 1;
        }

        // X-
        else if (face == 1 && i == 0)
        {
            idx2 = k*zStep + j*yStep + sizeX-2;
        }

        // Y+
        else if (face == 2 && j == sizeY-1 && DIMENSION >= 2)
        {
            idx2 = k*zStep + yStep + i;
        }

        // Y-
        else if (face == 3 && j == 0 && DIMENSION >= 2)
        {
            idx2 = k*zStep + (sizeY-2)*yStep + i;
        }

        // Z+
        else if (face == 4 && k == sizeZ-1 && DIMENSION == 3)
        {
            idx2 = zStep + j*yStep + i;
        }

        // Z-
        else if (face == 5 && k == 0 && DIMENSION == 3)
        {
            idx2 = (sizeZ-2)*zStep + j*yStep + i;
        }

        if (idx2 != -1)
        {
            for (long iter1 = 0; iter1 < numFields; iter1++)
            {
                field[iter1][idx] = field[iter1][idx2];
            }
        }
    }

    __syncthreads();
}

void applyBoundaryCondition(double **field, long fieldCode, long numFields,
                            domainInfo* simDomain, controls* simControls,
                            simParameters* simParams, subdomainInfo* subdomain,
                            dim3 gridSize, dim3 blockSize)
{
    if (simDomain->DIMENSION == 2)
    {
        // Y- at domain edge
        if (subdomain->rank == 0)
        {
            if (simControls->boundary[3][fieldCode].type == 1)
            {
                applyNeumann<<<gridSize, blockSize>>>(field, 3, numFields, simDomain->DIMENSION,
                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                      subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
            else if (simControls->boundary[3][fieldCode].type == 3)
            {
                applyPeriodic<<<gridSize, blockSize>>>(field, 3, numFields, simDomain->DIMENSION,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
        }

        // Y+ at domain edge
        else if (subdomain->rank == simControls->numworkers-1)
        {
            if (simControls->boundary[2][fieldCode].type == 1)
            {
                applyNeumann<<<gridSize, blockSize>>>(field, 2, numFields, simDomain->DIMENSION,
                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                      subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
            else if (simControls->boundary[2][fieldCode].type == 3)
            {
                applyPeriodic<<<gridSize, blockSize>>>(field, 2, numFields, simDomain->DIMENSION,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
        }
    }
    else if (simDomain->DIMENSION == 3)
    {
        // Z- at domain edge
        if (subdomain->rank == 0)
        {
            if (simControls->boundary[5][fieldCode].type == 1)
            {
                applyNeumann<<<gridSize, blockSize>>>(field, 5, numFields, simDomain->DIMENSION,
                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                      subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
            else if (simControls->boundary[5][fieldCode].type == 3)
            {
                applyPeriodic<<<gridSize, blockSize>>>(field, 5, numFields, simDomain->DIMENSION,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
        }

        // Z+ at domain edge
        else if (subdomain->rank == simControls->numworkers-1)
        {
            if (simControls->boundary[4][fieldCode].type == 1)
            {
                applyNeumann<<<gridSize, blockSize>>>(field, 4, numFields, simDomain->DIMENSION,
                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                      subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
            else if (simControls->boundary[4][fieldCode].type == 3)
            {
                applyPeriodic<<<gridSize, blockSize>>>(field, 4, numFields, simDomain->DIMENSION,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       subdomain->yStep, subdomain->zStep, subdomain->padding);
            }
        }


        // Y-
        if (simControls->boundary[3][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 3, numFields, simDomain->DIMENSION,
                                                  subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                  subdomain->yStep, subdomain->zStep, subdomain->padding);
        }
        else if (simControls->boundary[3][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 3, numFields, simDomain->DIMENSION,
                                                   subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                   subdomain->yStep, subdomain->zStep, subdomain->padding);
        }

        // Y+
        if (simControls->boundary[2][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 2, numFields, simDomain->DIMENSION,
                                                  subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                  subdomain->yStep, subdomain->zStep, subdomain->padding);
        }
        else if (simControls->boundary[2][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 2, numFields, simDomain->DIMENSION,
                                                   subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                   subdomain->yStep, subdomain->zStep, subdomain->padding);
        }
    }

    // X-
    if (simControls->boundary[1][fieldCode].type == 1)
    {
        applyNeumann<<<gridSize, blockSize>>>(field, 1, numFields, simDomain->DIMENSION,
                                              subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                              subdomain->yStep, subdomain->zStep, subdomain->padding);
    }
    else if (simControls->boundary[1][fieldCode].type == 3)
    {
        applyPeriodic<<<gridSize, blockSize>>>(field, 1, numFields, simDomain->DIMENSION,
                                               subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                               subdomain->yStep, subdomain->zStep, subdomain->padding);
    }

    // X+
    if (simControls->boundary[0][fieldCode].type == 1)
    {
        applyNeumann<<<gridSize, blockSize>>>(field, 0, numFields, simDomain->DIMENSION,
                                              subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                              subdomain->yStep, subdomain->zStep, subdomain->padding);
    }
    if (simControls->boundary[0][fieldCode].type == 3)
    {
        applyPeriodic<<<gridSize, blockSize>>>(field, 0, numFields, simDomain->DIMENSION,
                                               subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                               subdomain->yStep, subdomain->zStep, subdomain->padding);
    }
}
