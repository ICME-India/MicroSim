#include "utilityKernels.cuh"
#include <cub.cuh>

extern __device__
int LUPDecomposeC1(double A[][MAX_NUM_COMP], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i][j];
                    A[i][j] = A[imax][j];
                    A[imax][j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
int LUPDecomposeC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k*N+i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i*N+j];
                    A[i*N+j] = A[imax*N+j];
                    A[imax*N+j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j*N+i] /= A[i*N+i];

                for (k = i + 1; k < N; k++)
                    A[j*N+k] -= A[j*N+i] * A[i*N+k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
void LUPSolveC1(double A[][MAX_NUM_COMP], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

extern __device__
void LUPSolveC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i*N+k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i*N+k] * x[k];

        x[i] /= A[i*N+i];
    }
}

extern __device__
void LUPInvertC1(double A[][MAX_NUM_COMP], int *P, long N, double IA[][MAX_NUM_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

extern __device__
void LUPInvertC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, long N, double IA[][MAX_NUM_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];

            IA[i][j] /= A[i*N+i];
        }
    }
}

extern __device__
int LUPDecomposePC1(double A[][MAX_NUM_PHASE_COMP], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i][j];
                    A[i][j] = A[imax][j];
                    A[imax][j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
int LUPDecomposePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k*N+i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i*N+j];
                    A[i*N+j] = A[imax*N+j];
                    A[imax*N+j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j*N+i] /= A[i*N+i];

                for (k = i + 1; k < N; k++)
                    A[j*N+k] -= A[j*N+i] * A[i*N+k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
void LUPSolvePC1(double A[][MAX_NUM_PHASE_COMP], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

extern __device__
void LUPSolvePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i*N+k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i*N+k] * x[k];

        x[i] /= A[i*N+i];
    }
}

extern __device__
void LUPInvertPC1(double A[][MAX_NUM_PHASE_COMP], int *P, long N, double IA[][MAX_NUM_PHASE_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

extern __device__
void LUPInvertPC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, long N, double IA[][MAX_NUM_PHASE_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];

            IA[i][j] /= A[i*N+i];
        }
    }
}

/*
 * __device__ double calcPhaseEnergy
 *
 * Calculate f_{p}(c^{p})
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
    double ans = F0_C[phase];

    long index1, index2;

    for (long component1 = 0; component1 < NUMCOMPONENTS-1; component1++)
    {
        index1 = component1*NUMPHASES + phase;

        ans += F0_B[component1 + phase*(NUMCOMPONENTS-1)]*phaseComp[index1][idx];

        for (long component2 = 0; component2 <= component1; component2++)
        {
            index2 = component2*NUMPHASES + phase;
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
    double ans = F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*phaseComp[(phase + component*NUMPHASES)][idx];

    for (long i = 0; i < NUMCOMPONENTS-1; i++)
    {
        ans += F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)+i]*phaseComp[(phase + i*NUMPHASES)][idx];
    }

    ans += F0_B[component + phase*(NUMCOMPONENTS-1)];

    return ans;
}

extern __device__
double FunctionTau(double **phi, double *relaxCoeff, long idx, long NUMPHASES)
{
  double sum=0.0, sum1=0.0;
  long a, b;
  for (a=0; a<NUMPHASES; a++) {
    for (b=0; b<NUMPHASES; b++) {
      if (a<b) {
       sum  += relaxCoeff[a*NUMPHASES + b]*phi[a][idx]*phi[b][idx];
       sum1 += phi[a][idx]*phi[b][idx];
      }
    }
  }
  if (sum1) {
    return sum/sum1;
  } else {
    return relaxCoeff[1];
  }
}

__global__
void computeChange(double *A, double *B,
                   long sizeX, long sizeY, long sizeZ)
{
    long i = threadIdx.x + blockIdx.x*blockDim.x;
    long j = threadIdx.y + blockIdx.y*blockDim.y;
    long k = threadIdx.z + blockIdx.z*blockDim.z;

    long idx = (j + k*sizeY)*sizeX + i;

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
    long i, j = 0;

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
