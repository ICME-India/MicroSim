#if ENABLE_CUFFTMP == 1

#include "box_iterator.hpp"
#include "computeElastic.cuh"

__global__
void scale_output(BoxIterator<cufftDoubleComplex> begin, BoxIterator<cufftDoubleComplex> end,
                  double scalingFactor)
{
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    begin += tid;

    if (begin < end)
    {
        *begin = {begin->x/scalingFactor, begin->y/scalingFactor};
    }
}

__global__
void setFieldZero(BoxIterator<cufftDoubleComplex> begin, BoxIterator<cufftDoubleComplex> end)
{
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    begin += tid;

    if (begin < end)
    {
        *begin = {0.0, 0.0};
    }
}

__global__
void sumBoxIterator(BoxIterator<cufftDoubleComplex> begin, BoxIterator<cufftDoubleComplex> end, double *g_odata)
{

    __shared__ double sdata[256];

    // each thread loads one element from global to shared mem
    // note use of 1D thread indices (only) in this kernel
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    begin += i;

    if (begin < end)
        sdata[threadIdx.x] = begin->x;
    else
        sdata[threadIdx.x] = 0.0;

    __syncthreads();
    // do reduction in shared mem
    for (int s = 1; s < blockDim.x; s *= 2)
    {
        int index = 2*s*threadIdx.x;;

        if (index < blockDim.x)
        {
            sdata[index] += sdata[index + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (threadIdx.x == 0)
    {
        atomicAdd(g_odata, sdata[0]);
    }
}

__global__
void compute_hphi(double **phi, BoxIterator<cufftDoubleComplex> phiXtBegin, BoxIterator<cufftDoubleComplex> phiXtEnd,
                  long phase, long NUMPHASES, long xStep, long yStep, long sizeX, long sizeY, long sizeZ)
{
    const long x = threadIdx.x + blockIdx.x*blockDim.x;
    const long y = threadIdx.y + blockIdx.y*blockDim.y;
    const long z = threadIdx.z + blockIdx.z*blockDim.z;

    const long idx = x*xStep + y*yStep + z;

    if (x != 0 && x != sizeX-1 && y != 0 && y != sizeY-1)
    {
        double phiValue = calcInterp5th(phi, phase, idx, NUMPHASES);

        if (sizeZ == 1)
            phiXtBegin += (x-1)*(sizeY-2) + (y-1);
        else
            phiXtBegin += (x-1)*(sizeY-2)*(sizeZ-2) + (y-1)*(sizeZ-2) + (z-1);

        if (phiXtBegin < phiXtEnd)
            *phiXtBegin = {phiValue, 0.0};
        else
            *phiXtBegin = {0.0, 0.0};
    }
}

__global__
void compute_hphiAvg(BoxIterator<cufftDoubleComplex> begin, BoxIterator<cufftDoubleComplex> end, double avg)
{
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;
    begin += tid;

    if (begin < end)
    {
        *begin = {begin->x - avg, 0.0};
    }
}

__global__
void copydfdphi(double **dfdphi, BoxIterator<cufftDoubleComplex> dfdphiXtBegin, BoxIterator<cufftDoubleComplex> dfdphiXtEnd,
                long phase, long xStep, long yStep, long sizeX, long sizeY, long sizeZ)
{
    const long x = threadIdx.x + blockIdx.x*blockDim.x;
    const long y = threadIdx.y + blockIdx.y*blockDim.y;
    const long z = threadIdx.z + blockIdx.z*blockDim.z;

    const long idx = x*xStep + y*yStep + z;

    if (x != 0 && x != sizeX-1 && y != 0 && y != sizeY-1)
    {
        if (sizeZ == 1)
            dfdphiXtBegin += (x-1)*(sizeY-2) + (y-1);
        else
            dfdphiXtBegin += (x-1)*(sizeY-2)*(sizeZ-2) + (y-1)*(sizeZ-2) + (z-1);

        if (dfdphiXtBegin < dfdphiXtEnd)
            dfdphi[phase][idx] = dfdphiXtBegin->x;
    }
}

// __global__
// void calcdfEldphi(BoxIterator<cufftDoubleComplex> phiBegin, BoxIterator<cufftDoubleComplex> phiEnd,
//                   BoxIterator<cufftDoubleComplex> dfEldphiBegin, BoxIterator<cufftDoubleComplex> dfEldphiEnd,
//                   double epsxx_p, double epsyy_p, double epszz_p, double epsxy_p, double epsxz_p, double epsyz_p,
//                   double epsxx_q, double epsyy_q, double epszz_q, double epsxy_q, double epsxz_q, double epsyz_q,
//                   double C11_p, double C12_p, double C44_p,
//                   double C11_q, double C12_q, double C44_q,
//                   long MESH_X, long MESH_Y, long MESH_Z,
//                   double DELTA_X, double DELTA_Y, double DELTA_Z,
//                   long DIMENSION)
// {
//     const int tid = threadIdx.x + blockIdx.x*blockDim.x;
//
//     phiBegin += tid;
//     dfEldphiBegin += tid;
//
//     if (phiBegin < phiEnd)
//     {
//         double Bpq = 0.0;
//         double C[6][6], lambda[3][3][3][3];
//         double eigstr1[3][3], eigstr2[3][3];
//         double stress1[3][3], stress2[3][3];
//
//         int i, j, k, l, m, n1;
//
//         C[0][0] = 0.5*(C11_p + C11_q);
//         C[1][1] = 0.5*(C11_p + C11_q);
//         C[2][2] = 0.5*(C11_p + C11_q);
//
//         C[0][1] = 0.5*(C12_p + C12_q);
//         C[1][0] = 0.5*(C12_p + C12_q);
//         C[0][2] = 0.5*(C12_p + C12_q);
//         C[2][0] = 0.5*(C12_p + C12_q);
//         C[1][2] = 0.5*(C12_p + C12_q);
//         C[2][1] = 0.5*(C12_p + C12_q);
//
//         C[3][3] = 0.5*(C44_p + C44_q);
//         C[4][4] = 0.5*(C44_p + C44_q);
//         C[5][5] = 0.5*(C44_p + C44_q);
//
//         for (i = 0; i < 3; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 for (k = 0; k < 3; ++k) {
//                     for (l = 0; l < 3; ++l) {
//                         lambda[i][j][k][l] = 0.0;
//                     }
//                 }
//             }
//         }
//
//         m = 0;
//         n1 = 0;
//
//         for (i = 0; i < 3; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 for (k = 0; k < 3; ++k) {
//                     for (l = 0; l < 3; ++l) {
//                         if (i == j)
//                             m = i;
//                         if (k == l)
//                             n1 = k;
//                         if ((i == 1 && j == 2) || (i == 2 && j == 1))
//                             m = 3;
//                         if ((k == 1 && l == 2) || (k == 2 && l == 1))
//                             n1 = 3;
//                         if ((i == 0 && j == 2) || (i == 2 && j == 0))
//                             m = 4;
//                         if ((k == 0 && l == 2) || (k == 2 && l == 0))
//                             n1 = 4;
//                         if ((i == 0 && j == 1) || (i == 1 && j == 0))
//                             m = 5;
//                         if ((k == 0 && l == 1) || (k == 1 && l == 0))
//                             n1= 5;
//                         lambda[i][j][k][l] = C[m][n1];
//                     }
//                 }
//             }
//         }
//
//         eigstr1[0][0] = epsxx_p;
//         eigstr1[1][1] = epsyy_p;
//         eigstr1[2][2] = epszz_p;
//
//         eigstr1[0][1] = eigstr1[1][0] = epsxy_p;
//         eigstr1[0][2] = eigstr1[2][0] = epsxz_p;
//         eigstr1[1][2] = eigstr1[2][1] = epsyz_p;
//
//         eigstr2[0][0] = epsxx_q;
//         eigstr2[1][1] = epsyy_q;
//         eigstr2[2][2] = epszz_q;
//
//         eigstr2[0][1] = eigstr2[1][0] = epsxy_q;
//         eigstr2[0][2] = eigstr2[2][0] = epsxz_q;
//         eigstr2[1][2] = eigstr2[2][1] = epsyz_q;
//
//         for (i = 0; i < 3; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 stress1[i][j] = 0.0;
//                 stress2[i][j] = 0.0;
//             }
//         }
//
//         for (i = 0; i < 3; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 for (k = 0; k < 3; ++k) {
//                     for (l = 0; l < 3; ++l) {
//                         stress1[i][j] += lambda[i][j][k][l]*eigstr1[k][l];
//                         stress2[i][j] += lambda[i][j][k][l]*eigstr2[k][l];
//                     }
//                 }
//             }
//         }
//
//         double inv_omega[3][3], n[3];
//
//         if (phiBegin.x() <= MESH_X/2)
//             n[0] = (double)phiBegin.x() * (2.0*M_PI/((double)MESH_X*DELTA_X));
//         else
//             n[0] =  (2.0*M_PI/((double)MESH_X*DELTA_X)) * ((double)phiBegin.x() - MESH_X);
//
//         if (DIMENSION >= 2)
//         {
//             if (phiBegin.y() <= MESH_Y/2)
//                 n[1] = (double)phiBegin.y() * (2.0*M_PI/((double)MESH_Y*DELTA_Y));
//             else
//                 n[1] = (2.0*M_PI/((double)MESH_Y*DELTA_Y)) * ((double)phiBegin.y() - MESH_Y);
//         }
//         else
//         {
//             n[1] = 0.0;
//         }
//
//         if (DIMENSION == 3)
//         {
//             if (phiBegin.z() <= MESH_Z/2)
//                 n[2] = (double)phiBegin.z() * (2.0*M_PI/((double)MESH_Z*DELTA_Z));
//             else
//                 n[2] = (2.0*M_PI/((double)MESH_Z*DELTA_Z)) * ((double)phiBegin.z() - MESH_Z);
//         }
//         else
//         {
//             n[2] = 0.0;
//         }
//
//         for (i = 0; i < 3; ++i)
//             for (j = 0; j < 3; ++j)
//                 inv_omega[i][j] = 0.0;
//
//         for (i = 0; i < 3; ++i)
//             for (j = 0; j < 3; ++j)
//                 for (k = 0; k < 3; ++k)
//                     for (l = 0; l < 3; ++l)
//                         inv_omega[i][j] += lambda[i][k][l][j]*n[k]*n[l];
//
//         double omega[3][3];
//
//         for (i = 0; i < 3; ++i)
//             for( j = 0; j < 3; ++j)
//                 omega[i][j] = 0.0;
//
//         double det = inv_omega[0][0]*(inv_omega[1][1]*inv_omega[2][2] - inv_omega[1][2]*inv_omega[2][1]) -
//         inv_omega[0][1]*(inv_omega[1][0]*inv_omega[2][2] - inv_omega[1][2]*inv_omega[2][0]) +
//         inv_omega[0][2]*(inv_omega[1][0]*inv_omega[2][1] - inv_omega[1][1]*inv_omega[2][0]);
//
//         if (fabs(det) > 1.0e-30)
//         {
//             omega[0][0] = (inv_omega[1][1]*inv_omega[2][2] - inv_omega[1][2]*inv_omega[2][1])/det;
//             omega[0][1] = (inv_omega[2][1]*inv_omega[0][2] - inv_omega[0][1]*inv_omega[2][2])/det;
//             omega[0][2] = (inv_omega[0][1]*inv_omega[1][2] - inv_omega[0][2]*inv_omega[1][1])/det;
//             omega[1][0] = (inv_omega[2][0]*inv_omega[1][2] - inv_omega[1][0]*inv_omega[2][2])/det;
//             omega[1][1] = (inv_omega[0][0]*inv_omega[2][2] - inv_omega[2][0]*inv_omega[0][2])/det;
//             omega[1][2] = (inv_omega[1][0]*inv_omega[0][2] - inv_omega[0][0]*inv_omega[1][2])/det;
//             omega[2][0] = (inv_omega[1][0]*inv_omega[2][1] - inv_omega[2][0]*inv_omega[1][1])/det;
//             omega[2][1] = (inv_omega[2][0]*inv_omega[0][1] - inv_omega[0][0]*inv_omega[2][1])/det;
//             omega[2][2] = (inv_omega[0][0]*inv_omega[1][1] - inv_omega[1][0]*inv_omega[0][1])/det;
//         }
//
//         for (i = 0; i < 3; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 for (k = 0; k < 3; ++k) {
//                     for (l = 0; l < 3; ++l) {
//                         Bpq += lambda[i][j][k][l]*eigstr1[i][j]*eigstr2[k][l]
//                         - n[i]*stress1[i][j]*omega[j][k]*stress2[k][l]*n[l];
//                     }
//                 }
//             }
//         }
//
//         *dfEldphiBegin = {dfEldphiBegin->x + Bpq*phiBegin->x, dfEldphiBegin->y + Bpq*phiBegin->y};
//     }
// }

__global__
void calcdfEldphi(BoxIterator<cufftDoubleComplex> phiBegin, BoxIterator<cufftDoubleComplex> phiEnd,
                  BoxIterator<cufftDoubleComplex> dfEldphiBegin, BoxIterator<cufftDoubleComplex> dfEldphiEnd, double *Bpq)
{
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;

    phiBegin += tid;
    dfEldphiBegin += tid;

    if (phiBegin < phiEnd)
    {
        *dfEldphiBegin = {dfEldphiBegin->x + Bpq[tid]*phiBegin->x, dfEldphiBegin->y + Bpq[tid]*phiBegin->y};
    }
}

__global__
void __calc_k__(double *kx, double *ky, double *kz, BoxIterator<cufftDoubleComplex> phiBegin, BoxIterator<cufftDoubleComplex> phiEnd,
                double dx, double dy, double dz,
                long MESH_X, long MESH_Y, long MESH_Z, long numCells)
{
    const int tid = threadIdx.x + blockIdx.x*blockDim.x;

    phiBegin += tid;

    if (phiBegin < phiEnd && tid < numCells)
    {
        long x = phiBegin.x();
        long y = phiBegin.y();
        long z = phiBegin.z();

        if (x <= MESH_X/2)
            kx[tid] = 2.0*M_PI*x/((double)MESH_X*dx);
        else
            kx[tid] = 2.0*M_PI*(x - MESH_X)/((double)MESH_X*dx);

        if (y <= MESH_Y/2)
            ky[tid] = 2.0*M_PI*y/((double)MESH_Y*dy);
        else
            ky[tid] = 2.0*M_PI*(y - MESH_Y)/((double)MESH_Y*dy);

        if (z <= MESH_Z/2)
            kz[tid] = 2.0*M_PI*z/((double)MESH_Z*dz);
        else
            kz[tid] = 2.0*M_PI*(z - MESH_Z)/((double)MESH_Z*dz);

    }
}

void calc_k(double *kx, double *ky, double *kz, cudaLibXtDesc *phiXt,
            long my_nx, long my_ny, long my_nz,
            domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
            cudaStream_t stream)
{
    double *kx_d, *ky_d, *kz_d;
    cudaMalloc((void**)&kx_d, sizeof(double)*my_nx*my_ny*my_nz);
    cudaMalloc((void**)&ky_d, sizeof(double)*my_nx*my_ny*my_nz);
    cudaMalloc((void**)&kz_d, sizeof(double)*my_nx*my_ny*my_nz);

    auto[phiBegin, phiEnd] = BoxIterators(CUFFT_XT_FORMAT_INPLACE_SHUFFLED, CUFFT_Z2Z,
                                          subdomain.rank, subdomain.size, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, (cufftDoubleComplex*)phiXt->descriptor->data[0]);

    static const size_t num_elements = std::distance(phiBegin, phiEnd);
    static const size_t num_threads  = 256;
    static const size_t num_blocks   = (num_elements + num_threads - 1) / num_threads;

    __calc_k__<<<num_blocks, num_threads, 0, stream>>>(kx_d, ky_d, kz_d, phiBegin, phiEnd,
                                                       simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z,
                                                       simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, my_nx*my_ny*my_nz);

    cudaMemcpy(kx, kx_d, sizeof(double)*my_nx*my_ny*my_nz, cudaMemcpyDeviceToHost);
    cudaMemcpy(ky, ky_d, sizeof(double)*my_nx*my_ny*my_nz, cudaMemcpyDeviceToHost);
    cudaMemcpy(kz, kz_d, sizeof(double)*my_nx*my_ny*my_nz, cudaMemcpyDeviceToHost);

    cudaFree(kx_d);
    cudaFree(ky_d);
    cudaFree(kz_d);
}

void moveTocudaLibXtDesc(double **phi, cudaLibXtDesc *phiXt[MAX_NUM_PHASES],
                         domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                         cudaStream_t stream, dim3 gridSize, dim3 blockSize, MPI_Comm comm)
{
    double *hphiSum, hphiLocal, hphiAvg;
    cudaMalloc((void**)&hphiSum, sizeof(double));
    cudaMemset(hphiSum, 0, sizeof(double));

    for (long p = 0; p < simDomain.numPhases; p++)
    {
        if (simControls.eigenSwitch[p] == 0)
            continue;

        auto[phiBegin, phiEnd] = BoxIterators(CUFFT_XT_FORMAT_INPLACE, CUFFT_Z2Z,
                                              subdomain.rank, subdomain.size, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, (cufftDoubleComplex*)phiXt[p]->descriptor->data[0]);

        compute_hphi<<<gridSize, blockSize, 0, stream>>>(phi, phiBegin, phiEnd,
                                                         p, simDomain.numPhases, subdomain.xStep, subdomain.yStep, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);

        static const size_t num_elements = std::distance(phiBegin, phiEnd);
        static const size_t num_threads  = 256;
        static const size_t num_blocks   = (num_elements + num_threads - 1) / num_threads;

        sumBoxIterator<<<num_blocks, num_threads, 0, stream>>>(phiBegin, phiEnd, hphiSum);
        cudaMemcpy(&hphiLocal, hphiSum, sizeof(double), cudaMemcpyDeviceToHost);
        MPI_Allreduce(&hphiLocal, &hphiAvg, 1, MPI_DOUBLE, MPI_SUM, comm);

        hphiAvg /= (double)simDomain.MESH_X*simDomain.MESH_Y*simDomain.MESH_Z;

        compute_hphiAvg<<<num_blocks, num_threads, 0, stream>>>(phiBegin, phiEnd, hphiAvg);
    }

    cudaFree(hphiSum);
}


void moveFromcudaLibXtDesc(double **dfdphi, cudaLibXtDesc *dfdphiXt[MAX_NUM_PHASES],
                           domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                           cudaStream_t stream, dim3 gridSize, dim3 blockSize)
{
    for (long p = 0; p < simDomain.numPhases; p++)
    {
        if (simControls.eigenSwitch[p] == 0)
            continue;

        auto[dfdphiBegin, dfdphiEnd] = BoxIterators(CUFFT_XT_FORMAT_INPLACE, CUFFT_Z2Z,
                                                    subdomain.rank, subdomain.size, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, (cufftDoubleComplex*)dfdphiXt[p]->descriptor->data[0]);

        static const size_t num_elements = std::distance(dfdphiBegin, dfdphiEnd);
        static const size_t num_threads  = 256;
        static const size_t num_blocks   = (num_elements + num_threads - 1) / num_threads;

        scale_output<<<num_blocks, num_threads, 0, stream>>>(dfdphiBegin, dfdphiEnd, simDomain.MESH_X*simDomain.MESH_Y*simDomain.MESH_Z);
        copydfdphi<<<gridSize, blockSize, 0, stream>>>(dfdphi, dfdphiBegin, dfdphiEnd,
                                                       p, subdomain.xStep, subdomain.yStep, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ);
    }
}

void computeDrivingForce_Elastic(cudaLibXtDesc *phi[MAX_NUM_PHASES], cudaLibXtDesc *dfeldphi[MAX_NUM_PHASES], double **Bpq,
                                 domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                                 cudaStream_t stream, MPI_Comm comm)
{
    for (long p = 0; p < simDomain.numPhases; p++)
    {
        auto[dfEldphiBegin, dfEldphiEnd] = BoxIterators(CUFFT_XT_FORMAT_INPLACE_SHUFFLED, CUFFT_Z2Z,
                                                        subdomain.rank, subdomain.size, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, (cufftDoubleComplex*)dfeldphi[p]->descriptor->data[0]);

        const size_t num_elements = std::distance(dfEldphiBegin, dfEldphiEnd);
        const size_t num_threads  = 256;
        const size_t num_blocks   = (num_elements + num_threads - 1) / num_threads;

        setFieldZero<<<num_blocks, num_threads, 0, stream>>>(dfEldphiBegin, dfEldphiEnd);

        if (simControls.eigenSwitch[p] == 0)
            continue;

        for (long q = 0; q < simDomain.numPhases; q++)
        {
            if (simControls.eigenSwitch[q] == 0)
                continue;

            auto[phiBegin, phiEnd] = BoxIterators(CUFFT_XT_FORMAT_INPLACE_SHUFFLED, CUFFT_Z2Z,
                                                  subdomain.rank, subdomain.size, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, (cufftDoubleComplex*)phi[q]->descriptor->data[0]);
/*
            calcdfEldphi<<<num_blocks, num_threads, 0, stream>>>(phiBegin, phiEnd,
                                                                 dfEldphiBegin, dfEldphiEnd,
                                                                 simParams.eigen_strain[p].xx, simParams.eigen_strain[p].yy, simParams.eigen_strain[p].zz,
                                                                 simParams.eigen_strain[p].xy, simParams.eigen_strain[p].xz, simParams.eigen_strain[p].yz,
                                                                 simParams.eigen_strain[q].xx, simParams.eigen_strain[q].yy, simParams.eigen_strain[q].zz,
                                                                 simParams.eigen_strain[q].xy, simParams.eigen_strain[q].xz, simParams.eigen_strain[q].yz,
                                                                 simParams.Stiffness_c[p].C11, simParams.Stiffness_c[p].C12, simParams.Stiffness_c[p].C44,
                                                                 simParams.Stiffness_c[q].C11, simParams.Stiffness_c[q].C12, simParams.Stiffness_c[q].C44,
                                                                 simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z,
                                                                 simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z,
                                                                 simDomain.DIMENSION);*/

           calcdfEldphi<<<num_blocks, num_threads, 0, stream>>>(phiBegin, phiEnd, dfEldphiBegin, dfEldphiEnd, Bpq[p*simDomain.numPhases + q]);
        }
    }
}

#endif
