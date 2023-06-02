#include "boundary.cuh"

__global__
void applyNeumann(double **field,
                  long face, long numFields, long DIMENSION,
                  long sizeX, long sizeY, long sizeZ,
                  long xStep, long yStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k, idx2 = -1;

    if ((i < padding || i >= sizeX-padding || ((j < padding || j >= sizeY-padding) && DIMENSION >= 2) || ((k < padding || k >= sizeZ-padding) && DIMENSION == 3)) && i < sizeX && j < sizeY && k < sizeZ)
    {
        // X+
        if (face == 0 && i >= sizeX-padding)
        {
            idx2 = (2*(sizeX-padding) - i - 1)*xStep + j*yStep + k;
        }

        // X-
        else if (face == 1 && i < padding)
        {
            idx2 = (2*padding - i - 1)*xStep + j*yStep + k;
        }

        // Y+
        else if (face == 2 && j >= sizeY-padding && DIMENSION >= 2)
        {
            idx2 = i*xStep + (2*(sizeY-padding) - j - 1)*yStep + k;
        }

        // Y-
        else if (face == 3 && j < padding && DIMENSION >= 2)
        {
            idx2 = i*xStep + (2*padding - j - 1)*yStep + k;
        }

        // Z+
        else if (face == 4 && k >= sizeZ-padding && DIMENSION == 3)
        {
            idx2 = i*xStep + j*yStep + (2*(sizeZ - padding) - k - 1);
        }

        // Z-
        else if (face == 5 && k < padding && DIMENSION == 3)
        {
            idx2 = i*xStep + j*yStep + (2*padding - k - 1);
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
                   long xStep, long yStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k, idx2 = -1;

    if ((i == 0 || i == sizeX-1 || ((j == 0 || j == sizeY-1) && DIMENSION >= 2) || ((k == 0 || k == sizeZ-1) && DIMENSION == 3)) && i < sizeX && j < sizeY && k < sizeZ)
    {
        // X+
        if (face == 0 && i == sizeX-1)
        {
            idx2 = xStep + j*yStep + k;
        }

        // X-
        else if (face == 1 && i == 0)
        {
            idx2 = (sizeX-2)*xStep + j*yStep + k;
        }

        // Y+
        else if (face == 2 && j == sizeY-1 && DIMENSION >= 2)
        {
            idx2 = i*xStep + yStep + k;
        }

        // Y-
        else if (face == 3 && j == 0 && DIMENSION >= 2)
        {
            idx2 = i*xStep + (sizeY-2)*yStep + k;
        }

        // Z+
        else if (face == 4 && k == sizeZ-1 && DIMENSION == 3)
        {
            idx2 = i*xStep + j*yStep + 1;
        }

        // Z-
        else if (face == 5 && k == 0 && DIMENSION == 3)
        {
            idx2 = i*xStep + j*yStep + (sizeZ-2);
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
                            domainInfo simDomain, controls simControls,
                            simParameters simParams, subdomainInfo subdomain,
                            dim3 gridSize, dim3 blockSize)
{
    if (subdomain.rank == 0)
    {
        if (simControls.boundary[1][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 1, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[1][fieldCode].type == 3 && subdomain.size == 1)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 1, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
    }

    if (subdomain.rank == subdomain.size-1)
    {
        if (simControls.boundary[0][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 0, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[0][fieldCode].type == 3  && subdomain.size == 1)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 0, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
    }

    if (simDomain.DIMENSION >= 2)
    {
        // Y-
        if (simControls.boundary[3][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 3, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[3][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 3, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }

        // Y+
        if (simControls.boundary[2][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 2, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[2][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 2, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
    }

    if (simDomain.DIMENSION == 3)
    {
        // Z-
        if (simControls.boundary[5][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 5, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[5][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 5, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }

        // Z+
        if (simControls.boundary[4][fieldCode].type == 1)
        {
            applyNeumann<<<gridSize, blockSize>>>(field, 4, numFields, simDomain.DIMENSION,
                                                  subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                  subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
        else if (simControls.boundary[4][fieldCode].type == 3)
        {
            applyPeriodic<<<gridSize, blockSize>>>(field, 4, numFields, simDomain.DIMENSION,
                                                   subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                   subdomain.xStep, subdomain.yStep, subdomain.padding);
        }
    }
}
