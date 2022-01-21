#ifndef KERNELS_H_
#define KERNELS_H_

__global__ void Update_comp_phi(cufftDoubleComplex *comp, cufftDoubleComplex *phi,
                                cufftDoubleComplex *dfdc, cufftDoubleComplex *dfdphi,
                                double *B, int elast_int,
                                double dt, double diffuse,
                                double kappa_phi, double relax_coeff,
                                double* kx, double* ky, double *kz,
                                long MESH_X, long MESH_Y, long MESH_Z)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    double kpow2, lhs, lhse, n[3];
    cufftDoubleComplex  rhs, rhse;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        n[0] = kx[i];
        n[1] = ky[j];
        n[2] = kz[k];

        kpow2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];

        lhs = 1.0 + diffuse * kpow2 * dt;

        rhs.x = comp[idx].x + dt * dfdc[idx].x;
        rhs.y = comp[idx].y + dt * dfdc[idx].y;

        comp[idx].x = rhs.x/lhs;
        comp[idx].y = rhs.y/lhs;

        if (elast_int)
            lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt + relax_coeff*dt*B[idx];
        else
            lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt;

        rhse.x  = phi[idx].x - relax_coeff * dt * dfdphi[idx].x;
        rhse.y  = phi[idx].y - relax_coeff * dt * dfdphi[idx].y;

        phi[idx].x = rhse.x/lhse;
        phi[idx].y = rhse.y/lhse;

        dfdphi[idx].x = phi[idx].x;
        dfdphi[idx].y = phi[idx].y;
    }

    __syncthreads();
}

__global__ void Normalize(cufftDoubleComplex *x, double sizescale, long MESH_X, long MESH_Y, long MESH_Z)
{

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        x[idx].x = x[idx].x * sizescale;
        x[idx].y = x[idx].y * sizescale;
    }

    __syncthreads();
}

__global__ void SaveReal(double *temp, cufftDoubleComplex *x, long MESH_X, long MESH_Y, long MESH_Z)
{

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
        temp[idx] = x[idx].x;

    __syncthreads();
}

__global__ void Find_err_matrix(double *temp, cufftDoubleComplex *comp, long MESH_X, long MESH_Y, long MESH_Z)
{

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
        temp[idx] = fabs(comp[idx].x - temp[idx]);

    __syncthreads();
}

#endif
