#include "utilityKernels.cuh"
#include <cub.cuh>

__device__ __host__
double evalFunc(void f(double, double*, double*), double x, double temperature)
{
    double c[2];
    c[0] = x;
    c[1] = 1.0 - c[0];

    double ans;

    f(temperature, c, &ans);

    // Non-dimensionalise
    ans /= (1.602*1e8);

    return ans;
}

__device__ __host__
double spline_eval(double x, double *controlPoints,
                   double *a, double *b, double *c, double *d,
                   int numControlPoints)
{
    double ans = 0.0;

    for (int i = 0; i < numControlPoints-1; i++)
    {
        if (x >= controlPoints[i] && x <= controlPoints[i+1])
        {
            ans = x - controlPoints[i];
            ans = d[i] + ans*(c[i] + ans*(b[i] + ans*a[i]));
        }
    }

    return ans;
}

__global__
void computeChange(double *A, double *B,
                   int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        A[idx] = fabs(B[idx] - A[idx]);
    }
    __syncthreads();
}

void printStats(double **phi, double **comp,
                double **phiNew, double **compNew,
                double *maxerr, double *maxVal, double *minVal,
                domainInfo simDomain, subdomainInfo subdomain,
                dim3 gridSize, dim3 blockSize)
{
    int i, j = 0;

    void    *t_storage = NULL;
    size_t  t_storage_bytes = 0;

    double *temp;
    cudaMalloc((void**)&temp, sizeof(double));

    cub::DeviceReduce::Max(t_storage, t_storage_bytes, phiNew[1], temp, subdomain.numCells);
    cub::DeviceReduce::Sum(t_storage, t_storage_bytes, phiNew[1], temp, subdomain.numCells);
    cudaMalloc((void**)&t_storage, t_storage_bytes);

    for (i = 0; i < simDomain.numComponents-1; i++)
    {
        cub::DeviceReduce::Max(t_storage, t_storage_bytes, compNew[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(maxVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);
        cub::DeviceReduce::Min(t_storage, t_storage_bytes, compNew[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(minVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);

        computeChange<<<gridSize, blockSize>>>(comp[i], compNew[i], subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        cub::DeviceReduce::Max(t_storage, t_storage_bytes, comp[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(maxerr+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);

        j++;
    }

    for (i = 0; i < simDomain.numPhases; i++)
    {
        cub::DeviceReduce::Max(t_storage, t_storage_bytes, phiNew[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(maxVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);
        cub::DeviceReduce::Min(t_storage, t_storage_bytes, phiNew[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(minVal+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);

        computeChange<<<gridSize, blockSize>>>(phi[i], phiNew[i], subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        cub::DeviceReduce::Max(t_storage, t_storage_bytes, phi[i] + subdomain.shiftPointer, temp, subdomain.numCells);
        cudaMemcpy(maxerr+j, temp, sizeof(double), cudaMemcpyDeviceToDevice);

        j++;
    }

    cudaFree(temp);
    cudaFree(t_storage);
}
