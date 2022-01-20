#ifndef KERNELS_FFT_CUH_
#define KERNELS_FFT_CUH_

__global__ void  ComputeGradphi_FFT(cufftDoubleComplex *phi,
                                    cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz,
                                    double* kx, double* ky, double *kz,
                                    long MESH_X, long MESH_Y, long MESH_Z)
{

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

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

__global__ void ComputeDrivForce_FFT(cufftDoubleComplex *comp, cufftDoubleComplex *phi,
                                     cufftDoubleComplex *dfdphi,
                                     cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz,
                                     double f0AVminv, double f0BVminv,
                                     double cAlphaEq, double cBetaEq,
                                     double diffuse, double w,
                                     long MESH_X, long MESH_Y, long MESH_Z)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    double  interp_phi, interp_prime, g_prime;
    double  ctemp, etemp;
    double  f_alpha, f_beta, mubar;
    double  A_by_B, B_by_A;
    double  calpha, cbeta;

    A_by_B = f0AVminv/f0BVminv;
    B_by_A = f0BVminv/f0AVminv;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        ctemp  = comp[idx].x;
        etemp  = phi[idx].x;

        interp_phi   = etemp * etemp * etemp *
        (6.0 * etemp * etemp - 15.0 * etemp + 10.0);

        interp_prime = 30.0 * etemp * etemp * pow((1.0 - etemp), 2.0);

        g_prime      = 2.0 * etemp * (1.0 - etemp) * (1.0 - 2.0 * etemp);

        calpha       = (ctemp - interp_phi * (cBetaEq - cAlphaEq * A_by_B))/
        (interp_phi * A_by_B + (1.0 - interp_phi));

        cbeta        = (ctemp + (1.0 - interp_phi) * (B_by_A * cBetaEq - cAlphaEq))/
        (interp_phi + B_by_A*(1.0 - interp_phi));

        comp[idx].x  = calpha * (1.0 - interp_phi) + cbeta * interp_phi;
        comp[idx].y  = 0.0;

        f_alpha      = f0AVminv * (calpha - cAlphaEq) * (calpha - cAlphaEq);
        f_beta       = f0BVminv * (cbeta - cBetaEq) * (cbeta - cBetaEq);

        gradphix[idx].x = diffuse * interp_prime * (calpha - cbeta) * gradphix[idx].x;
        gradphix[idx].y = diffuse * interp_prime * (calpha - cbeta) * gradphix[idx].y;

        gradphiy[idx].x = diffuse * interp_prime * (calpha - cbeta) * gradphiy[idx].x;
        gradphiy[idx].y = diffuse * interp_prime * (calpha - cbeta) * gradphiy[idx].y;

        gradphiz[idx].x = diffuse * interp_prime * (calpha - cbeta) * gradphiz[idx].x;
        gradphiz[idx].y = diffuse * interp_prime * (calpha - cbeta) * gradphiz[idx].y;

        mubar = 2.0 * f0AVminv * (calpha - cAlphaEq);

        dfdphi[idx].x = interp_prime * (f_beta - f_alpha + (calpha - cbeta) * mubar) + w * g_prime;
        dfdphi[idx].y = 0.0;
    }

    __syncthreads();
}

// __global__ void ComputeDrivForce_FFT_TW(cufftDoubleComplex *comp, cufftDoubleComplex *phi, cufftDoubleComplex *dfdphi,
//                                        cufftDoubleComplex *gradphix, cufftDoubleComplex *gradphiy, cufftDoubleComplex *gradphiz,
//                                        double cAlphaEq, double cBetaEq,
//                                        double diffuse, double w, double temperature,
//                                        long MESH_X, long MESH_Y, long MESH_Z)
// {
//     long i = threadIdx.x + blockIdx.x * blockDim.x;
//     long j = threadIdx.y + blockIdx.y * blockDim.y;
//     long k = threadIdx.z + blockIdx.z * blockDim.z;
//
//     long idx = (j + i * MESH_Y)*MESH_Z + k;
//
//     double  interp_phi, interp_prime, g_prime;
//     double  ctemp, etemp;
//     double  f_alpha, f_beta;
//     double  tol_alpha = 1.0, tol_beta = 1.0;
//     double  c_alpha[2], c_beta[2];
//     double  c_alpha_new, c_beta_new;
//
//     double  dgA_dcA, dgA2_d2cA;
//     double  dgB_dcB, dgB2_d2cB;
//
//     double func[2], Jac[4], detJac;
//
//     if (i < MESH_X && j < MESH_Y && k < MESH_Z)
//     {
//         ctemp  = comp[idx].x;
//         etemp  = phi[idx].x;
//
//         interp_phi   = etemp * etemp * etemp *
//         (6.0 * etemp * etemp - 15.0 * etemp + 10.0);
//
//         interp_prime = 30.0 * etemp * etemp * pow((1.0 - etemp), 2.0);
//
//         g_prime      = 2.0 * etemp * (1.0 - etemp) * (1.0 - 2.0 * etemp);
//
//         c_alpha[0] = cAlphaEq;
//         c_alpha[1] = 1.0 - c_alpha[0];
//         c_beta[0] = cBetaEq;
//         c_beta[1] = 1.0 - c_beta[0];
//
//         c_alpha_new = 1;
//         c_beta_new = 1;
//
//         int iter = 0;
//
//         if (etemp > 0.001 &&  etemp < 0.999)
//         {
//             do
//             {
//                 func[0] = c_beta[0]*interp_phi + c_alpha[0]*(1.0 - interp_phi) - ctemp;
//                 dGES(temperature, c_alpha, &dgA_dcA);
//                 dGEL(temperature, c_beta, &dgB_dcB);
//                 func[1] = (dgA_dcA - dgB_dcB);
//
//                 Jac[0] = 1.0 - interp_phi;
//                 Jac[1] = interp_phi;
//                 ddGES(temperature, c_alpha, &dgA2_d2cA);
//                 ddGEL(temperature, c_beta, &dgB2_d2cB);
//                 Jac[2] = dgA2_d2cA;
//                 Jac[3] = -dgB2_d2cB;
//
//                 detJac = (Jac[0]*Jac[3] - Jac[1]*Jac[2]);
//
//                 c_alpha_new = c_alpha[0] - (func[0]*Jac[3] - func[1]*Jac[1])/detJac;
//                 c_beta_new  = c_beta[0]  + (func[0]*Jac[2] - func[1]*Jac[0])/detJac;
//
//                 tol_alpha = c_alpha_new - c_alpha[0];
//                 if (tol_alpha < 0.0)
//                     tol_alpha *= -1.0;
//
//                 tol_beta = c_beta_new - c_beta[0];
//                 if (tol_beta < 0.0)
//                     tol_beta *= -1.0;
//
//                 c_alpha[0] = c_alpha_new;
//                 c_alpha[1] = 1.0 - c_alpha[0];
//                 c_beta[0]  = c_beta_new;
//                 c_beta[1] = 1.0 - c_beta[0];
//
//                 iter++;
//
//                 if (iter > 25)
//                     break;
//
//             } while (tol_alpha > 1e-3 || tol_beta > 1e-3);
//         }
//
//         if (c_alpha[0] < cAlphaEq)
//         {
//             c_alpha[0] = cAlphaEq;
//             c_alpha[1] = 1.0 - cAlphaEq;
//         }
//         if (c_beta[0] > cBetaEq)
//         {
//             c_beta[0] = cBetaEq;
//             c_beta[1] = 1.0 - cBetaEq;
//         }
//
//         GES(temperature, c_alpha, &f_alpha);
//         GEL(temperature, c_beta, &f_beta);
//
//         dGES(temperature, c_alpha, &dgA_dcA);;
//
//         dfdphi[idx].x = (interp_prime * (f_beta - f_alpha + (c_alpha[0] - c_beta[0]) * dgA_dcA)/(8.314*temperature) + w * g_prime);
//         if (i == MESH_X/2 && j == MESH_Y/2 && k == MESH_Z/2)
//             printf("interp_prime = %lf\ng_prime = %lf\n", interp_prime, g_prime);
//         dfdphi[idx].y = 0.0;
//
//         gradphix[idx].x = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphix[idx].x;
//         gradphix[idx].y = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphix[idx].y;
//
//         gradphiy[idx].x = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphiy[idx].x;
//         gradphiy[idx].y = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphiy[idx].y;
//
//         gradphiz[idx].x = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphiz[idx].x;
//         gradphiz[idx].y = diffuse * interp_prime * (c_alpha[0]- c_beta[0]) * gradphiz[idx].y;
//
//     }
//     __syncthreads();
// }

__global__ void  ComputeDfdc_FFT(cufftDoubleComplex *dfdc,
                                 cufftDoubleComplex *gradPhi_x, cufftDoubleComplex *gradPhi_y, cufftDoubleComplex *gradPhi_z,
                                 double* kx, double* ky, double *kz,
                                 long MESH_X, long MESH_Y, long MESH_Z)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    double  n[3];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z) {
        n[0] = kx[i];
        n[1] = ky[j];
        n[2] = kz[k];

        dfdc[idx].x = -1.0 * (n[0] * gradPhi_x[idx].y + n[1] * gradPhi_y[idx].y + n[2] * gradPhi_z[idx].y);

        dfdc[idx].y = (n[0] * gradPhi_x[idx].x + n[1] * gradPhi_y[idx].x + n[2] * gradPhi_z[idx].x);
    }

    __syncthreads();
}

#endif
