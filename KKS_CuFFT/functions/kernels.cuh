#ifndef KERNELS_H_
#define KERNELS_H_

__global__ void  ComputeGradphi_3D(double* kx, double* ky, double *kz, int MESH_X, int MESH_Y, int MESH_Z, cufftDoubleComplex *phi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    double  n[3];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        n[0] = kx[i];
        n[1] = ky[j];
        n[2] = kz[k];

        gradphix[idx].x = -1.0*n[0]*phi[idx].y;
        gradphix[idx].y = n[0]*phi[idx].x;
        gradphiy[idx].x = -1.0*n[1]*phi[idx].y;
        gradphiy[idx].y = n[1]*phi[idx].x;
        gradphiz[idx].x = -1.0*n[2]*phi[idx].y;
        gradphiz[idx].y = n[2]*phi[idx].x;
    }

    __syncthreads();

}

__global__ void ComputeSqGradphi_3D(int MESH_X, int MESH_Y, int MESH_Z, cufftDoubleComplex *phi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz, double kappa_phi) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        phi[idx].x = 2.0 * (kappa_phi) * (gradphix[idx].x * gradphix[idx].x +
        gradphiy[idx].x * gradphiy[idx].x + gradphiz[idx].x * gradphiz[idx].x);
        phi[idx].y = 2.0 * (kappa_phi) * (gradphix[idx].y * gradphix[idx].y +
        gradphiy[idx].y * gradphiy[idx].y + gradphiz[idx].x * gradphiz[idx].x);
    }
}

__global__ void ComputeDrivForce_3D(cufftDoubleComplex *comp, cufftDoubleComplex *dfdphi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz, double f0AVminv, double f0BVminv, double c_beta_eq, double c_alpha_eq, double diff, double w, int MESH_X, int MESH_Y, int MESH_Z) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    double  interp_phi, interp_prime, g_prime;
    double  ctemp, etemp;
    double  f_alpha, f_beta, mubar;
    double  A_by_B, B_by_A;
    double  calpha, cbeta;

    A_by_B = f0AVminv/f0BVminv;
    B_by_A = f0BVminv/f0AVminv;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        ctemp  = comp[idx].x;
        etemp  = dfdphi[idx].x;

        interp_phi   = etemp * etemp * etemp *
        (6.0 * etemp * etemp - 15.0 * etemp + 10.0);

        interp_prime = 30.0 * etemp * etemp * pow((1.0 - etemp), 2.0);

        g_prime      = 2.0 * etemp * (1.0 - etemp) * (1.0 - 2.0 * etemp);

        calpha       = (ctemp - interp_phi * (c_beta_eq - c_alpha_eq * A_by_B))/
        (interp_phi * A_by_B + (1.0 - interp_phi));

        cbeta        = (ctemp + (1.0 - interp_phi) * (B_by_A * c_beta_eq - c_alpha_eq))/
        (interp_phi + B_by_A*(1.0 - interp_phi));

        comp[idx].x  = calpha * (1.0 - interp_phi) + cbeta * interp_phi;
        comp[idx].y  = 0.0;

        f_alpha      = f0AVminv * (calpha - c_alpha_eq) * (calpha - c_alpha_eq);
        f_beta       = f0BVminv * (cbeta - c_beta_eq) * (cbeta - c_beta_eq);

        gradphix[idx].x = diff * interp_prime * (calpha - cbeta) * gradphix[idx].x;
        gradphix[idx].y = diff * interp_prime * (calpha - cbeta) * gradphix[idx].y;

        gradphiy[idx].x = diff * interp_prime * (calpha - cbeta) * gradphiy[idx].x;
        gradphiy[idx].y = diff * interp_prime * (calpha - cbeta) * gradphiy[idx].y;

        gradphiz[idx].x = diff * interp_prime * (calpha - cbeta) * gradphiz[idx].x;
        gradphiz[idx].y = diff * interp_prime * (calpha - cbeta) * gradphiz[idx].y;

        mubar = 2.0 * f0BVminv * (calpha - c_alpha_eq);

        dfdphi[idx].x = interp_prime * (f_beta - f_alpha + (calpha - cbeta) * mubar) + w * g_prime;
        dfdphi[idx].y = 0.0;
    }
    __syncthreads();

}

__global__ void  ComputeDfdc_3D(cufftDoubleComplex *dfdc, cufftDoubleComplex *varmobx, cufftDoubleComplex *varmoby, cufftDoubleComplex *varmobz, int MESH_X, int MESH_Y, int MESH_Z, double* kx, double* ky, double *kz) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    double  n[3];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        n[0] = kx[i];
        n[1] = ky[j];
        n[2] = kz[k];

        dfdc[idx].x = -1.0 * (n[0] * varmobx[idx].y + n[1] * varmoby[idx].y + n[2] * varmobz[idx].y);

        dfdc[idx].y = (n[0] * varmobx[idx].x + n[1] * varmoby[idx].x + n[2] * varmobz[idx].x);
    }
    __syncthreads();
}

__global__ void Update_comp_phi_3D(cufftDoubleComplex *comp, cufftDoubleComplex *dfdc, cufftDoubleComplex *phi, cufftDoubleComplex *dfdphi, double* kx, double* ky, double *kz, double dt, double diff, double kappa_phi, double relax_coeff, int MESH_X, int MESH_Y, int MESH_Z) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    double           kpow2, lhs, lhse, n[3];
    cufftDoubleComplex  rhs, rhse;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        n[0] = kx[i];
        n[1] = ky[j];
        n[2] = kz[k];

        kpow2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];

        lhs = 1.0 + diff * kpow2 * dt;

        rhs.x = comp[idx].x + dt * dfdc[idx].x;
        rhs.y = comp[idx].y + dt * dfdc[idx].y;

        comp[idx].x = rhs.x/lhs;
        comp[idx].y = rhs.y/lhs;

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

__global__ void Normalize_3D(cufftDoubleComplex *x, double sizescale, int MESH_X, int MESH_Y, int MESH_Z)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        x[idx].x = x[idx].x * sizescale;
        x[idx].y = x[idx].y * sizescale;
    }

    __syncthreads();
}

__global__ void SaveReal_3D(double *temp, cufftDoubleComplex *x, int MESH_X, int MESH_Y, int MESH_Z)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
        temp[idx] = x[idx].x;
}

__global__ void Find_err_matrix_3D(double *temp, cufftDoubleComplex *comp_d, int MESH_X, int MESH_Y, int MESH_Z)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
        temp[idx] = fabs(comp_d[idx].x - temp[idx]);

}

__global__ void  ComputeGradphi_2D(double* kx, double* ky, int MESH_X, int MESH_Y, cufftDoubleComplex *phi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    int idx = (j + i * MESH_Y);

    double  n[2];

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
        {
        n[0] = kx[i];
        n[1] = ky[j];

        gradphix[idx].x = -1.0*n[0]*phi[idx].y;
        gradphix[idx].y = n[0]*phi[idx].x;
        gradphiy[idx].x = -1.0*n[1]*phi[idx].y;
        gradphiy[idx].y = n[1]*phi[idx].x;
    }

    __syncthreads();

}

__global__ void ComputeSqGradphi_2D(int MESH_X, int MESH_Y, cufftDoubleComplex *phi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, double kappa_phi)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    int idx = (j + i * MESH_Y);

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
    {
        phi[idx].x = 2.0 * (kappa_phi) * (gradphix[idx].x * gradphix[idx].x +
        gradphiy[idx].x * gradphiy[idx].x);
        phi[idx].y = 2.0 * (kappa_phi) * (gradphix[idx].y * gradphix[idx].y +
        gradphiy[idx].y * gradphiy[idx].y);
    }
}

__global__ void ComputeDrivForce_2D(cufftDoubleComplex *comp, cufftDoubleComplex *dfdphi, cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, double f0AVminv, double f0BVminv, double c_beta_eq, double c_alpha_eq, double diff, double w, int MESH_X, int MESH_Y)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int idx = (j + i * MESH_Y);

    double  interp_phi, interp_prime, g_prime;
    double  ctemp, etemp;
    double  f_alpha, f_beta, mubar;
    double  A_by_B, B_by_A;
    double  calpha, cbeta;

    A_by_B = f0AVminv/f0BVminv;
    B_by_A = f0BVminv/f0AVminv;

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
    {
        ctemp  = comp[idx].x;
        etemp  = dfdphi[idx].x;

        interp_phi   = etemp * etemp * etemp *
        (6.0 * etemp * etemp - 15.0 * etemp + 10.0);

        interp_prime = 30.0 * etemp * etemp * pow((1.0 - etemp), 2.0);

        g_prime      = 2.0 * etemp * (1.0 - etemp) * (1.0 - 2.0 * etemp);

        calpha       = (ctemp - interp_phi * (c_beta_eq - c_alpha_eq * A_by_B))/
        (interp_phi * A_by_B + (1.0 - interp_phi));

        cbeta        = (ctemp + (1.0 - interp_phi) * (B_by_A * c_beta_eq - c_alpha_eq))/
        (interp_phi + B_by_A*(1.0 - interp_phi));

        comp[idx].x  = calpha * (1.0 - interp_phi) + cbeta * interp_phi;
        comp[idx].y  = 0.0;

        f_alpha      = f0AVminv * (calpha - c_alpha_eq) * (calpha - c_alpha_eq);
        f_beta       = f0BVminv * (cbeta - c_beta_eq) * (cbeta - c_beta_eq);

        gradphix[idx].x = diff * interp_prime * (calpha - cbeta) * gradphix[idx].x;
        gradphix[idx].y = diff * interp_prime * (calpha - cbeta) * gradphix[idx].y;

        gradphiy[idx].x = diff * interp_prime * (calpha - cbeta) * gradphiy[idx].x;
        gradphiy[idx].y = diff * interp_prime * (calpha - cbeta) * gradphiy[idx].y;

        mubar = 2.0 * f0BVminv * (calpha - c_alpha_eq);

        dfdphi[idx].x = interp_prime * (f_beta - f_alpha + (calpha - cbeta) * mubar) + w * g_prime;
        dfdphi[idx].y = 0.0;
    }
    __syncthreads();

}

__global__ void  ComputeDfdc_2D(cufftDoubleComplex *dfdc, cufftDoubleComplex *varmobx, cufftDoubleComplex *varmoby, int MESH_X, int MESH_Y, double* kx, double* ky)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int idx = (j + i * MESH_Y);

    double  n[2];

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
    {
        n[0] = kx[i];
        n[1] = ky[j];

        dfdc[idx].x = -1.0 * (n[0] * varmobx[idx].y + n[1] * varmoby[idx].y);

        dfdc[idx].y = (n[0] * varmobx[idx].x + n[1] * varmoby[idx].x);
    }
    __syncthreads();
}

__global__ void Update_comp_phi_2D(cufftDoubleComplex *comp, cufftDoubleComplex *dfdc, cufftDoubleComplex *phi, cufftDoubleComplex *dfdphi, double* kx, double* ky, double dt, double diff, double kappa_phi, double relax_coeff, int MESH_X, int MESH_Y, int elast_int, double* B)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int idx = (j + i * MESH_Y);

    double kpow2, lhs, lhse, n[2];
    cufftDoubleComplex  rhs, rhse;

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
    {
        n[0] = kx[i];
        n[1] = ky[j];

        kpow2 = n[0] * n[0] + n[1] * n[1] ;

        lhs = 1.0 + diff * kpow2 * dt;

        rhs.x = comp[idx].x + dt * dfdc[idx].x;
        rhs.y = comp[idx].y + dt * dfdc[idx].y;

        comp[idx].x = rhs.x/lhs;
        comp[idx].y = rhs.y/lhs;

        if (elast_int == 1)
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

__global__ void Normalize_2D(cufftDoubleComplex *x, double sizescale, int MESH_X, int MESH_Y)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int idx = (j + i * MESH_Y);

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
    {
        x[idx].x = x[idx].x * sizescale;
        x[idx].y = x[idx].y * sizescale;
    }
}

__global__ void SaveReal_2D(double *temp, cufftDoubleComplex *x, int MESH_X, int MESH_Y)
{

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;

    int idx = (j + i * MESH_Y);

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
        temp[idx] = x[idx].x;
}

__global__ void Find_err_matrix_2D(double *temp,cufftDoubleComplex *comp_d, int MESH_X, int MESH_Y)
{

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int idx = (j + i * MESH_Y);

    if (i >= 0 && i < MESH_X && j >= 0 && j < MESH_Y)
        temp[idx] = fabs(comp_d[idx].x - temp[idx]);

}

__global__ void checkOverlap_cuda(int MESH_X, int MESH_Y, int MESH_Z, double centx, double centy, double centz,
                             double rad, int *occupancy, int *overlap, int shield_dist)
{
    int i = threadIdx.x + blockDim.x*blockIdx.x;
    int j = threadIdx.y + blockDim.y*blockIdx.y;
    int k = threadIdx.z + blockDim.z*blockIdx.z;

    if ((((double)i - centx)*((double)i - centx) + ((double)j - centy)*((double)j - centy) + (double)(k - centz)*(double)(k - centz))
        <= ((shield_dist*rad)*(shield_dist*rad)))
        if (occupancy[k+j*MESH_Z+i*MESH_Y*MESH_Z] == 1)
            *overlap = 1;
}

#endif
