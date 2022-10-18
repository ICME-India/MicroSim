#include "utilityKernels.cuh"
#include <cub.cuh>

/*
 * __device__ double calcPhaseEnergy
 *
 * Calculate f_{p}(c^{p}) = \sum_{i = 1}^{K-1}\sum_{j = 1}^{K-1} A_{ij}^{p}c_{i}^{p}c_{j}^{p}
 *                          + \sum_{i = 1}^{K-1} c_{i}^{p}
 *                          + C^{p}
 *
 * Arguments:
 *              1. double **phaseComp -> all the phase compositions, ordered phase-by-phase
 *              2. long phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. long idx -> position of cell in 1D
 *              7. long NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
extern __device__
double calcPhaseEnergy(double **phaseComp, long phase,
                       double *F0_A, double *F0_B, double *F0_C,
                       long idx,
                       long NUMPHASES, long NUMCOMPONENTS)
{

    // Constant C_{phase}
    double ans = F0_C[phase];

    long index1, index2;

    for (long component1 = 0; component1 < NUMCOMPONENTS-1; component1++)
    {
        index1 = component1*NUMPHASES + phase;

        // Linear terms of the parabolic free energy
        ans += F0_B[component1 + phase*(NUMCOMPONENTS-1)]*phaseComp[index1][idx];

        for (long component2 = 0; component2 <= component1; component2++)
        {
            index2 = component2*NUMPHASES + phase;

            // Quadratic terms
            ans += F0_A[(component1 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component2]*phaseComp[index1][idx]*phaseComp[index2][idx];
        }
    }

    return ans;
}

/*
 * __device__ double calcDiffusionPotential
 *
 * Calculate \frac{df_{p}}{dc^{p}_{c}}
 *
 * Arguments:
 *              1. double **phaseComp -> all the phase compositions, ordered phase-by-phase
 *              2. long phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. long idx -> position of cell in 1D
 *              7. long NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
extern __device__
double calcDiffusionPotential(double **phaseComp,
                              long phase, long component,
                              double *F0_A, double *F0_B,
                              long idx,
                              long NUMPHASES, long NUMCOMPONENTS)
{

    // Derivative of A_{component, component}^{phase}*c_{component}^{phase}*c_{component}^{phase}
    // Rest of the quadratic terms are mixed
    double ans = 2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*phaseComp[(phase + component*NUMPHASES)][idx];

    for (long i = 0; i < NUMCOMPONENTS-1; i++)
    {

        // Derivative of A_{component, i}^{phase}*c_{component}^{phase}*c_{i}^{phase}
        if (i != component)
        {
            ans += F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)+i]*phaseComp[(phase + i*NUMPHASES)][idx];
        }
    }

    // Derivative of linear terms
    ans += F0_B[component + phase*(NUMCOMPONENTS-1)];

    return ans;
}

/*
 *  This kernel will compute the change in every cell
 *  as the absolute value of [B - A]
 *
 *  The differences are stored in A.
 *
 *  In practice, A is the old field, B is the new field.
 */
__global__
void computeChange(double *A, double *B, long DIMENSION,
                   long sizeX, long sizeY, long sizeZ)
{
    long i = threadIdx.x + blockIdx.x*blockDim.x;
    long j = threadIdx.y + blockIdx.y*blockDim.y;
    long k = threadIdx.z + blockIdx.z*blockDim.z;

    long idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        A[idx] = fabs(B[idx] - A[idx]);
    }

    __syncthreads();
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

        computeChange<<<gridSize, blockSize>>>(comp[i], compNew[i], simDomain.DIMENSION, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

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

        computeChange<<<gridSize, blockSize>>>(phi[i], phiNew[i], simDomain.DIMENSION, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        checkCudaErrors(cub::DeviceReduce::Max(t_storage, t_storage_bytes, phi[i], temp, subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(maxerr+j, temp, sizeof(double), cudaMemcpyDeviceToDevice));

        j++;
    }

    checkCudaErrors(cudaFree(temp));
    checkCudaErrors(cudaFree(t_storage));
}
