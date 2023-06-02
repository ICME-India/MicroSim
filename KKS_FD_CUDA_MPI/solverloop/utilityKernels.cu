#include "utilityKernels.cuh"
#include <cub/device/device_reduce.cuh>

/*
 *  This kernel will compute the change in every cell
 *  as the absolute value of [B - A]
 *
 *  The differences are stored in A.
 *
 *  In practice, A is the old field, B is the new field.
 */
__global__
void __computeChange__(double *A, double *B, long DIMENSION,
                       long sizeX, long sizeY, long sizeZ)
{
    long i = threadIdx.x + blockIdx.x*blockDim.x;
    long j = threadIdx.y + blockIdx.y*blockDim.y;
    long k = threadIdx.z + blockIdx.z*blockDim.z;

    long idx = (j + i*sizeY)*sizeZ + k;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        A[idx] = fabs(B[idx] - A[idx]);
    }

    __syncthreads();
}

__global__
void __resetArray__(double **arr, long numArr,
                    long DIMENSION,
                    long sizeX, long sizeY, long sizeZ)
{
    long i = threadIdx.x + blockIdx.x*blockDim.x;
    long j = threadIdx.y + blockIdx.y*blockDim.y;
    long k = threadIdx.z + blockIdx.z*blockDim.z;

    long idx = (j + i*sizeY)*sizeZ + k;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        for (long iter = 0; iter < numArr; iter++)
            arr[iter][idx] = 0.0;
    }
}

/*
 *  Aggregate max, min, relative change for each subdomain
 *
 *  This function can be modified to include mu if required.
 *  The size of maxerr, maxVal, minVal will have to be expanded to provide space for these.
 *
 */
void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize)
{
    long i, j = 0;

    void    *t_storage = NULL;
    size_t  t_storage_bytes = 0;

    // Temp is reused wherever cub returns a scalar value
    // A memcpy is called to move the returned value to the appropriate variable
    // Max, min, error are passed by reference to this kernel
    // Note that each of the above variables are arrays of size (NUMCOMPONENTS - 1 + NUMPHASES)
    // A counter var j is used to iterate through these arrays

    double *temp;
    checkCudaErrors(cudaMalloc((void**)&temp, sizeof(double)));

    // Initializing cub to some arbitrary value
    checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, phiNew[1], temp, subdomain.numCompCells));
    checkCudaErrors(cub::DeviceReduce::Sum(t_storage, t_storage_bytes, phiNew[1], temp, subdomain.numCompCells));
    checkCudaErrors(cudaMalloc((void**)&t_storage, t_storage_bytes));

    // Get all stats for compositions
    for (i = 0; i < simDomain.numComponents-1; i++)
    {
        checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, compNew[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(maxVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));
        checkCudaErrors(cub::DeviceReduce::Min(t_storage, t_storage_bytes, compNew[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(minVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));

        __computeChange__<<<gridSize, blockSize>>>(comp[i], compNew[i], simDomain.DIMENSION, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, comp[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(maxerr+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));

        j++;
    }

    // Get all stats for phi
    for (i = 0; i < simDomain.numPhases; i++)
    {
        checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, phiNew[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(maxVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));
        checkCudaErrors(cub::DeviceReduce::Min(t_storage, t_storage_bytes, phiNew[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(minVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));

        __computeChange__<<<gridSize, blockSize>>>(phi[i], phiNew[i], simDomain.DIMENSION, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, phi[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(maxerr+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));

        j++;
    }

    checkCudaErrors(cudaFree(temp));
    checkCudaErrors(cudaFree(t_storage));
}

void resetArray(double **arr, long numArr,
                long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                dim3 gridSize, dim3 blockSize)
{
    __resetArray__<<<gridSize, blockSize>>>(arr, numArr, DIMENSION, sizeX, sizeY, sizeZ);
}
