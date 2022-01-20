#ifndef KERNELS_FD_CUH_
#define KERNELS_FD_CUH_

__global__ void ComputeDrivForce_FD(cufftDoubleComplex *comp, cufftDoubleComplex *phi,
                                    cufftDoubleComplex *dfdc, cufftDoubleComplex *dfdphi,
                                    double f0AVminv, double f0BVminv, double c_alpha_eq, double c_beta_eq,
                                    double diffuse, double w, long MESH_X, long MESH_Y, long MESH_Z)
{

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    double  interp_phi, interp_prime, g_prime;
    double  ctemp, etemp;
    double  f_alpha, f_beta;
    double  A_by_B, B_by_A;
    double  calpha, cbeta;

    A_by_B = f0AVminv/f0BVminv;
    B_by_A = f0BVminv/f0AVminv;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        ctemp  = comp[idx].x;
        etemp  = phi[idx].x;

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

        dfdc[idx].x  = 0.0;
        dfdc[idx].y  = diffuse * interp_prime * (calpha - cbeta);

        dfdphi[idx].x = interp_prime * (f_beta - f_alpha + (calpha - cbeta) * 2.0 * f0AVminv * (calpha - c_alpha_eq)) + w * g_prime;
        dfdphi[idx].y = 0.0;
    }
    __syncthreads();
}

// __global__ void ComputeDrivForce_FD_TW(cufftDoubleComplex *comp, cufftDoubleComplex *phi,
//                                        cufftDoubleComplex *dfdc, cufftDoubleComplex *dfdphi,
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
//         dfdc[idx].x  = 0.0;
//         dfdc[idx].y  = diffuse * interp_prime * (c_alpha[0] - c_beta[0]);
//
//         GES(temperature, c_alpha, &f_alpha);
//         GEL(temperature, c_beta, &f_beta);
//
//         dGES(temperature, c_alpha, &dgA_dcA);;
//
//         dfdphi[idx].x = (interp_prime * (f_beta - f_alpha + (c_alpha[0] - c_beta[0]) * dgA_dcA)/(8.314*temperature) + w * g_prime);
//         if (i == MESH_X/2 && j == MESH_Y/2 && k == MESH_Z/2)
//             printf("dfdphi = %lf\tdfdc = %lf\n", dfdphi[idx].x, dfdc[idx].y);
//         dfdphi[idx].y = 0.0;
//     }
//     __syncthreads();
// }

// __global__ void ComputeGradphi_FD(double *tempPhi,
//                                   double *gradPhi_x, double *gradPhi_y, double *gradPhi_z,
//                                   long MESH_X, long MESH_Y, long MESH_Z,
//                                   double DELTA_X, double DELTA_Y, double DELTA_Z, int DIMENSION)
// {
//     long idx;
//
//     long i = threadIdx.x + blockIdx.x * blockDim.x;
//     long j = threadIdx.y + blockIdx.y * blockDim.y;
//     long k = threadIdx.z + blockIdx.z * blockDim.z;
//
//     long xp[2], xm[2], yp[2], ym[2], zp[2], zm[2];
//
//     if (i < MESH_X && j < MESH_Y && k < MESH_Z)
//     {
//         idx = (j + i*MESH_Y)*MESH_Z + k;
//
//         if (i == 0)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (MESH_X-2)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == 1)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = j*MESH_Z + k;
//             xm[1] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == MESH_X - 2)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = j*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == MESH_X - 1)
//         {
//             xp[0] = j*MESH_Z + k;
//             xp[1] = (j + MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//         else
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//
//         if (j == 0)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (MESH_Y-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == 1)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = i*MESH_Y*MESH_Z + k;
//             ym[1] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == MESH_Y - 2)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = i*MESH_Y*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == MESH_Y - 1)
//         {
//             yp[0] = i*MESH_Y*MESH_Z + k;
//             yp[1] = (1 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//
//         if (k == 0)
//         {
//             zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//             zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//             zm[0] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
//             zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-2;
//         }
//         else if (k == 1)
//         {
//             zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//             zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//             zm[0] = (j + i*MESH_Y)*MESH_Z;
//             zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
//         }
//         else if (k == MESH_Z - 2)
//         {
//             zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//             zp[1] = (j + i*MESH_Y)*MESH_Z;
//             zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//             zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//         }
//         else if (k == MESH_Z - 1)
//         {
//             zp[0] = (j + i*MESH_Y)*MESH_Z;
//             zp[1] = (j + i*MESH_Y)*MESH_Z + 1;
//             zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//             zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//         }
//         else
//         {
//             zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//             zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//             zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//             zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//         }
//
//         gradPhi_x[idx] = (-1*tempPhi[xp[1]] + 8*tempPhi[xp[0]] - 8*tempPhi[xm[0]] + tempPhi[xm[1]])/(12*DELTA_X);
//         gradPhi_y[idx] = (-1*tempPhi[yp[1]] + 8*tempPhi[yp[0]] - 8*tempPhi[ym[0]] + tempPhi[ym[1]])/(12*DELTA_Y);
//         if (DIMENSION == 3)
//             gradPhi_z[idx] = (-1*tempPhi[zp[1]] + 8*tempPhi[zp[0]] - 8*tempPhi[zm[0]] + tempPhi[zm[1]])/(12*DELTA_Z);
//         else
//             gradPhi_z[idx] = 0.0;
//     }
//     __syncthreads();
// }

__global__ void ComputeGradphi_x(double *tempPhi,
                                 double *gradPhi_x,
                                 long MESH_X, long MESH_Y, long MESH_Z,
                                 double DELTA_X)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i*MESH_Y)*MESH_Z + k;

//     __shared__ double sPhi[36][32][1];
//
//     sPhi[threadIdx.x+2][threadIdx.y][threadIdx.z] = tempPhi[idx];
//
//     __syncthreads();
//
//     if (threadIdx.x == 0)
//     {
//         sPhi[threadIdx.x][threadIdx.y][threadIdx.z] = tempPhi[(j + (MESH_X-2)*MESH_Y)*MESH_Z + k];
//         sPhi[threadIdx.x+1][threadIdx.y][threadIdx.z] = tempPhi[(j + (MESH_X-1)*MESH_Y)*MESH_Z + k];
//     }
//     else if (threadIdx.x == 31)
//     {
//         sPhi[threadIdx.x+3][threadIdx.y][threadIdx.z] = tempPhi[j*MESH_Z + k];
//         sPhi[threadIdx.x+4][threadIdx.y][threadIdx.z] = tempPhi[(j + MESH_Y)*MESH_Z + k];
//     }
//
//     __syncthreads();

    long xp[2], xm[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (i == 0)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (MESH_X-2)*MESH_Y)*MESH_Z + k;
        }
        else if (i == 1)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = j*MESH_Z + k;
            xm[1] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
        }
        else if (i == MESH_X - 2)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = j*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }
        else if (i == MESH_X - 1)
        {
            xp[0] = j*MESH_Z + k;
            xp[1] = (j + MESH_Y)*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }
        else
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }

        gradPhi_x[idx] = (-1*tempPhi[xp[1]] + 8*tempPhi[xp[0]] - 8*tempPhi[xm[0]] + tempPhi[xm[1]])/(12*DELTA_X);
    }
    __syncthreads();
}

__global__ void ComputeGradphi_y(double *tempPhi,
                                  double *gradPhi_y,
                                  long MESH_X, long MESH_Y, long MESH_Z,
                                  double DELTA_Y)
{
    long idx;

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long yp[2], ym[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        idx = (j + i*MESH_Y)*MESH_Z + k;

        if (j == 0)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (MESH_Y-2 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == 1)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = i*MESH_Y*MESH_Z + k;
            ym[1] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == MESH_Y - 2)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = i*MESH_Y*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == MESH_Y - 1)
        {
            yp[0] = i*MESH_Y*MESH_Z + k;
            yp[1] = (1 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }
        else
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }

        gradPhi_y[idx] = (-1*tempPhi[yp[1]] + 8*tempPhi[yp[0]] - 8*tempPhi[ym[0]] + tempPhi[ym[1]])/(12*DELTA_Y);
    }
    __syncthreads();
}

__global__ void ComputeGradphi_z(double *tempPhi,
                                  double *gradPhi_z,
                                  long MESH_X, long MESH_Y, long MESH_Z,
                                  double DELTA_Z)
{
    long idx;

    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long zp[2], zm[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        idx = (j + i*MESH_Y)*MESH_Z + k;
        if (k == 0)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-2;
        }
        else if (k == 1)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z;
            zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
        }
        else if (k == MESH_Z - 2)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }
        else if (k == MESH_Z - 1)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z;
            zp[1] = (j + i*MESH_Y)*MESH_Z + 1;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }
        else
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }

        gradPhi_z[idx] = (-1*tempPhi[zp[1]] + 8*tempPhi[zp[0]] - 8*tempPhi[zm[0]] + tempPhi[zm[1]])/(12*DELTA_Z);
    }
    __syncthreads();
}

// __global__ void  ComputeDfdc_FD(cufftDoubleComplex *dfdc, double *gradPhi_x, double *gradPhi_y, double *gradPhi_z,
//                                 long MESH_X, long MESH_Y, long MESH_Z,
//                                 double DELTA_X, double DELTA_Y, double DELTA_Z, int DIMENSION)
// {
//     long i = threadIdx.x + blockIdx.x * blockDim.x;
//     long j = threadIdx.y + blockIdx.y * blockDim.y;
//     long k = threadIdx.z + blockIdx.z * blockDim.z;
//
//     long idx = (j + i * MESH_Y)*MESH_Z + k;
//
//     long xp[2] = {0}, xm[2] = {0}, yp[2] = {0}, ym[2] = {0}, zp[2] = {0}, zm[2] = {0};
//
//     if (i < MESH_X && j < MESH_Y && k < MESH_Z)
//     {
//         if (i == 0)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (MESH_X-2)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == 1)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = j*MESH_Z + k;
//             xm[1] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == MESH_X - 2)
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = j*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//         else if (i == MESH_X - 1)
//         {
//             xp[0] = j*MESH_Z + k;
//             xp[1] = (j + MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//         else
//         {
//             xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
//             xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
//             xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
//             xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
//         }
//
//         if (j == 0)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (MESH_Y-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == 1)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = i*MESH_Y*MESH_Z + k;
//             ym[1] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == MESH_Y - 2)
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = i*MESH_Y*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else if (j == MESH_Y - 1)
//         {
//             yp[0] = i*MESH_Y*MESH_Z + k;
//             yp[1] = (1 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//         else
//         {
//             yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
//             yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
//             ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
//             ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
//         }
//
//         if (DIMENSION == 3)
//             if (k == 0)
//             {
//                 zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//                 zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//                 zm[0] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
//                 zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-2;
//             }
//             else if (k == 1)
//             {
//                 zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//                 zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//                 zm[0] = (j + i*MESH_Y)*MESH_Z;
//                 zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
//             }
//             else if (k == MESH_Z - 2)
//             {
//                 zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//                 zp[1] = (j + i*MESH_Y)*MESH_Z;
//                 zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//                 zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//             }
//             else if (k == MESH_Z - 1)
//             {
//                 zp[0] = (j + i*MESH_Y)*MESH_Z;
//                 zp[1] = (j + i*MESH_Y)*MESH_Z + 1;
//                 zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//                 zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//             }
//             else
//             {
//                 zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
//                 zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
//                 zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
//                 zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
//             }
//     }
//
//     dfdc[idx].x = ((-1*gradPhi_x[xp[1]] + 8*gradPhi_x[xp[0]] - 8*gradPhi_x[xm[0]] + gradPhi_x[xm[1]])/DELTA_X
//     + (-1*gradPhi_y[yp[1]] + 8*gradPhi_y[yp[0]] - 8*gradPhi_y[ym[0]] + gradPhi_y[ym[1]])/DELTA_Y
//     + (-1*gradPhi_z[zp[1]] + 8*gradPhi_z[zp[0]] - 8*gradPhi_z[zm[0]] + gradPhi_z[zm[1]])/DELTA_Z)/12.0;
//     dfdc[idx].y = 0.0;
//     __syncthreads();
// }

__global__ void  ComputeDfdc_x(cufftDoubleComplex *dfdc, double *gradPhi_x, long MESH_X, long MESH_Y, long MESH_Z, double DELTA_X)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    long xp[2], xm[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (i == 0)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (MESH_X-2)*MESH_Y)*MESH_Z + k;
        }
        else if (i == 1)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = j*MESH_Z + k;
            xm[1] = (j + (MESH_X-1)*MESH_Y)*MESH_Z + k;
        }
        else if (i == MESH_X - 2)
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = j*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }
        else if (i == MESH_X - 1)
        {
            xp[0] = j*MESH_Z + k;
            xp[1] = (j + MESH_Y)*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }
        else
        {
            xp[0] = (j + (i+1)*MESH_Y)*MESH_Z + k;
            xp[1] = (j + (i+2)*MESH_Y)*MESH_Z + k;
            xm[0] = (j + (i-1)*MESH_Y)*MESH_Z + k;
            xm[1] = (j + (i-2)*MESH_Y)*MESH_Z + k;
        }

        dfdc[idx].x += (-1*dfdc[xp[1]].y*gradPhi_x[xp[1]] + 8*dfdc[xp[0]].y*gradPhi_x[xp[0]]
                    - 8*dfdc[xm[0]].y*gradPhi_x[xm[0]] + dfdc[xm[1]].y*gradPhi_x[xm[1]])/(12.0*DELTA_X);
    }
    __syncthreads();
}

__global__ void  ComputeDfdc_y(cufftDoubleComplex *dfdc, double *gradPhi_y, long MESH_X, long MESH_Y, long MESH_Z, double DELTA_Y)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    long yp[2], ym[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (j == 0)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (MESH_Y-2 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == 1)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = i*MESH_Y*MESH_Z + k;
            ym[1] = (MESH_Y-1 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == MESH_Y - 2)
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = i*MESH_Y*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }
        else if (j == MESH_Y - 1)
        {
            yp[0] = i*MESH_Y*MESH_Z + k;
            yp[1] = (1 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }
        else
        {
            yp[0] = (j+1 + i*MESH_Y)*MESH_Z + k;
            yp[1] = (j+2 + i*MESH_Y)*MESH_Z + k;
            ym[0] = (j-1 + i*MESH_Y)*MESH_Z + k;
            ym[1] = (j-2 + i*MESH_Y)*MESH_Z + k;
        }

        dfdc[idx].x += (-1*dfdc[yp[1]].y*gradPhi_y[yp[1]] + 8*dfdc[yp[0]].y*gradPhi_y[yp[0]]
                    - 8*dfdc[ym[0]].y*gradPhi_y[ym[0]] + dfdc[ym[1]].y*gradPhi_y[ym[1]])/(12.0*DELTA_Y);
    }
    __syncthreads();
}

__global__ void  ComputeDfdc_z(cufftDoubleComplex *dfdc, double *gradPhi_z, long MESH_X, long MESH_Y, long MESH_Z, double DELTA_Z)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    long zp[2], zm[2];

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
    {
        if (k == 0)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-2;
        }
        else if (k == 1)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z;
            zm[1] = (j + i*MESH_Y)*MESH_Z + MESH_Z-1;
        }
        else if (k == MESH_Z - 2)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }
        else if (k == MESH_Z - 1)
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z;
            zp[1] = (j + i*MESH_Y)*MESH_Z + 1;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }
        else
        {
            zp[0] = (j + i*MESH_Y)*MESH_Z + k+1;
            zp[1] = (j + i*MESH_Y)*MESH_Z + k+2;
            zm[0] = (j + i*MESH_Y)*MESH_Z + k-1;
            zm[1] = (j + i*MESH_Y)*MESH_Z + k-2;
        }

        dfdc[idx].x += (-1*dfdc[zp[1]].y*gradPhi_z[zp[1]] + 8*dfdc[zp[0]].y*gradPhi_z[zp[0]]
        - 8*dfdc[zm[0]].y*gradPhi_z[zm[0]] + dfdc[zm[1]].y*gradPhi_z[zm[1]])/(12.0*DELTA_Z);
    }
    __syncthreads();
}

__global__ void reset_dfdc(cufftDoubleComplex *dfdc, long MESH_X, long MESH_Y, long MESH_Z)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = (j + i * MESH_Y)*MESH_Z + k;

    if (i < MESH_X && j < MESH_Y && k < MESH_Z)
        dfdc[idx].y = 0.0;
}
#endif
