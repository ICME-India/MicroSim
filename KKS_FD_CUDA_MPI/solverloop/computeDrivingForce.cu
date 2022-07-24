#include "computeDrivingForce.cuh"

__device__
double calcInterp5thDiff(double **phi, int a, int b, int idx, int NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double ans = 0.0;
    double temp = 0.0;

    double phiValue = phi[a][idx];
    double const1 = 7.5*((double)NUMPHASES-2.0)/((double)NUMPHASES-1.0);

    if (a == b)
    {
        ans  = 5.0*pow(phiValue, 4)*(6.0 - const1);
        ans += 4.0*pow(phiValue, 3)*(-15.0 + 3.0*const1);

        for (int i = 1; i < NUMPHASES; i++)
        {
            if (i != a)
            {
                for (int j = 0; j < i; j++)
                {
                    if (j != a)
                    {
                        temp += (phi[j][idx] - phi[i][idx])*(phi[j][idx] - phi[i][idx]);
                    }
                }
            }
        }

        temp *= 7.5*(3.0*phiValue - 2.0)/((double)NUMPHASES-1.0);
        temp += 3.0*phiValue*(10.0 - 3.0*const1) + 2.0*const1;
        ans  += phiValue*temp;

        temp  = 0.0;
        if (NUMPHASES > 3)
        {
            for (int i = 2; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (int j = 1; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (int k = 0; k < j; k++)
                            {
                                if (k != a)
                                {
                                    temp += phi[i][idx]*phi[j][idx]*phi[k][idx];
                                }
                            }
                        }
                    }
                }
            }
            ans += 30.0*phiValue*temp;
            temp = 0.0;
        }

        if (NUMPHASES > 4)
        {
            for (int i = 3; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (int j = 2; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (int k = 1; k < j; k++)
                            {
                                if (k != a)
                                {
                                    for (int l = 0; l < k; l++)
                                    {
                                        if (l != a)
                                        {
                                            temp += phi[i][idx]*phi[j][idx]*phi[k][idx]*phi[l][idx];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ans += 24.0*temp;
        }
    }
    else
    {
        if (NUMPHASES == 2)
        {
            ans = -30.0*phiValue*phiValue*(1.0-phiValue)*(1.0-phiValue);
        }

        if (NUMPHASES > 2)
        {
            for (int i = 1; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (int j = 0; j < i; j++)
                    {
                        if (j != a)
                        {
                            if (i == b)
                                temp += -2.0*(phi[j][idx] - phi[i][idx]);
                            else if (j == b)
                                temp += 2.0*(phi[j][idx] - phi[i][idx]);
                        }
                    }
                }
            }
            ans = (phiValue-1.0)*phiValue*phiValue*7.5*temp/((double)NUMPHASES-1.0);
            temp = 0.0;
        }
        if (NUMPHASES > 3)
        {
            for (int i = 2; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (int j = 1; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (int k = 0; k < j; k++)
                            {
                                if (k != a)
                                {
                                    if (i == b)
                                        temp += phi[j][idx]*phi[k][idx];
                                    else if (j == b)
                                        temp += phi[i][idx]*phi[k][idx];
                                    else if (k == b)
                                        temp += phi[i][idx]*phi[j][idx];
                                }
                            }
                        }
                    }
                }
            }
            ans += 15*phiValue*phiValue*temp;
            temp = 0.0;
        }

        if (NUMPHASES > 4)
        {
            for (int i = 3; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (int j = 2; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (int k = 1; k < j; k++)
                            {
                                if (k != a)
                                {
                                    for (int l = 0; l < k; l++)
                                    {
                                        if (l != a)
                                        {
                                            if (i == b)
                                                temp += phi[j][idx]*phi[k][idx]*phi[l][idx];
                                            else if (j == b)
                                                temp += phi[i][idx]*phi[k][idx]*phi[l][idx];
                                            else if (k == b)
                                                temp += phi[i][idx]*phi[j][idx]*phi[l][idx];
                                            else if (l == b)
                                                temp += phi[i][idx]*phi[j][idx]*phi[k][idx];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ans += 24.0*phiValue*temp;
        }
    }

    return ans;
}

/*
 * __device__ double calcDoubleWellDerivative
 *
 * Calculate g'(\phi)
 *
 * Arguments:
 *              1. cufftDoubleComplex **phi -> all the phase volume-fraction values
 *              2. int phase -> differentiate wrt to this phase
 *              3. double *theta_i   -> coefficients for theta_i, one per phase
 *              4. double *theta_ij  -> coefficients for theta_ij, one per pair of phases
 *              5. double *theta_ijk -> coefficients for theta_ijk, one per triplet of phases
 *              6. int idx -> position of cell in 1D
 *              7. int NUMPHASES -> number of phases
 * Return:
 *              numerical evaluation of derivative of interpolation polynomial, as a double datatype
 */
__device__
double calcDoubleWellDerivative(double **phi, int phase,
                                double *theta_i, double *theta_ij, double *theta_ijk,
                                int idx, int NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double ans = 0.0, temp1 = 0.0, temp2 = 0.0;

    for (int i = 0; i < NUMPHASES; i++)
    {
        temp1 = phi[i][idx]*phi[i][idx];

        // Derivative of \phi^{2}(1 - \phi)^{2}
        if (i == phase)
            ans += theta_i[i]*2.0*phi[i][idx]*(1.0 - phi[i][idx])*(1.0 - 2.0*phi[i][idx]);

        for (int j = 0; j < NUMPHASES; j++)
        {
            if (j == i)
                continue;

            temp2 = phi[j][idx]*phi[j][idx];

            // Derivative of \sum_{i=1}^{N} \sum_{j=1, j!= i}^{N} \phi^{2}_{i}\phi^{2}_{j}
            if (i == phase)
                ans += 2.0*theta_ij[j + i*NUMPHASES]*phi[i][idx]*temp2;
            else if (j == phase)
                ans += 2.0*theta_ij[j + i*NUMPHASES]*temp1*phi[j][idx];

            for (int k = 0; k < NUMPHASES; k++)
            {
                if (k == i || k == j)
                    continue;
                // Derivate of (ijk)^2
                if (i == phase)
                    ans += 2.0*theta_ijk[(j + i*NUMPHASES)*NUMPHASES + k]*phi[i][idx]*temp2*phi[k][idx]*phi[k][idx];
                else if (j == phase)
                    ans += 2.0*theta_ijk[(j + i*NUMPHASES)*NUMPHASES + k]*temp1*phi[j][idx]*phi[k][idx]*phi[k][idx];
                else if (k == phase)
                    ans += 2.0*theta_ijk[(j + i*NUMPHASES)*NUMPHASES + k]*temp1*temp2*phi[k][idx];
            }
        }
    }

    return ans;
}

/*
 * __device__ double calcPhaseEnergy
 *
 * Calculate f_{p}(c^{p})
 *
 * Arguments:
 *              1. double **phaseComp -> all the phase compositions, ordered phase-by-phase
 *              2. int phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. int idx -> position of cell in 1D
 *              7. int NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
__device__
double calcPhaseEnergy(double **phaseComp, int phase,
                       double *F0_A, double *F0_B, double *F0_C,
                       int idx,
                       int NUMPHASES, int NUMCOMPONENTS)
{
    double ans = F0_C[phase];

    int index1, index2;

    for (int component1 = 0; component1 < NUMCOMPONENTS-1; component1++)
    {
        index1 = component1*NUMPHASES + phase;

        ans += F0_B[component1 + phase*(NUMCOMPONENTS-1)]*phaseComp[index1][idx];

        for (int component2 = 0; component2 <= component1; component2++)
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
 *              2. int phase -> phase for which energy is being calculated
 *              3. double *F0_A -> coefficients for F0_A (quadratic)
 *              4. double *F0_B -> coefficients for F0_B (linear)
 *              5. double *F0_C -> coefficients for F0_C (constant in a phase)
 *              6. int idx -> position of cell in 1D
 *              7. int NUMCOMPONENTS -> number of components
 * Return:
 *              numerical evaluation of bulk energy of the phase, as a double datatype
 */
__device__
double calcDiffusionPotential(double **phaseComp,
                              int phase, int component,
                              double *F0_A, double *F0_B,
                              int idx,
                              int NUMPHASES, int NUMCOMPONENTS)
{
    double ans = F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*phaseComp[(phase + component*NUMPHASES)][idx];

    for (int i = 0; i < NUMCOMPONENTS-1; i++)
    {
        ans += F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)+i]*phaseComp[(phase + i*NUMPHASES)][idx];
    }

    ans += F0_B[component + phase*(NUMCOMPONENTS-1)];

    return ans;
}

__global__ void __computeDrivingForce_01_03__(double **phi, double **comp,
                                              double **dfdphi,
                                              double **phaseComp,
                                              double *F0_A, double *F0_B, double *F0_C,
                                              double *diffusivity, double *kappaPhi,
                                              double *theta_i, double *theta_ij, double *theta_ijk,
                                              int NUMPHASES, int NUMCOMPONENTS,
                                              int sizeX, int sizeY, int sizeZ,
                                              double DELTA_X, double DELTA_Y, double DELTA_Z)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp[2], xm[2], yp[2], ym[2], zp[2], zm[2];

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp[0] = (j + k*sizeY)*sizeX + i+1;
        xp[1] = (j + k*sizeY)*sizeX + i+2;
        xm[0] = (j + k*sizeY)*sizeX + i-1;
        xm[1] = (j + k*sizeY)*sizeX + i-2;

        if (i == 0)
        {
            xm[0] = (j + k*sizeY)*sizeX + sizeX-1;
            xm[1] = (j + k*sizeY)*sizeX + sizeX-2;
        }
        else if (i == 1)
        {
            xm[1] = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 2)
        {
            xp[1] = (j + k*sizeY)*sizeX;
        }
        else if (i == sizeX - 1)
        {
            xp[0] = (j + k*sizeY)*sizeX;
            xp[1] = (j + k*sizeY)*sizeX + 1;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp[0] = (j+1 + k*sizeY)*sizeX + i;
            yp[1] = (j+2 + k*sizeY)*sizeX + i;
            ym[0] = (j-1 + k*sizeY)*sizeX + i;
            ym[1] = (j-2 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym[0] = (sizeY-1 + k*sizeY)*sizeX + i;
                ym[1] = (sizeY-2 + k*sizeY)*sizeX + i;
            }
            else if (j == 1)
            {
                ym[1] = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 2)
            {
                yp[1] = (k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp[0] = (k*sizeY)*sizeX + i;
                yp[1] = (1 + k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp[0] = (j + (k+1)*sizeY)*sizeX + i;
            zp[1] = (j + (k+2)*sizeY)*sizeX + i;
            zm[0] = (j + (k-1)*sizeY)*sizeX + i;
            zm[1] = (j + (k-2)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm[0] = (j + (sizeZ-1)*sizeY)*sizeX + i;
                zm[1] = (j + (sizeZ-2)*sizeY)*sizeX + i;
            }
            else if (k == 1)
            {
                zm[1] = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 2)
            {
                zp[1] = j*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp[0] = j*sizeX + i;
                zp[1] = (j + sizeY)*sizeX + i;
            }
        }

        for (int phase = 0; phase < NUMPHASES; phase++)
        {
            // f'_{bulk}
            dfdphi[phase][idx] = 0.0;
            for (int p = 0; p < NUMPHASES; p++)
            {
                dfdphi[phase][idx] += calcInterp5thDiff(phi, p, phase, idx, NUMPHASES)*calcPhaseEnergy(phaseComp, p, F0_A, F0_B, F0_C, idx, NUMPHASES, NUMCOMPONENTS);

                for (int component = 0; component < NUMCOMPONENTS-1; component++)
                    dfdphi[phase][idx] -= calcInterp5thDiff(phi, p, phase, idx, NUMPHASES)
                    *calcDiffusionPotential(phaseComp, phase, component, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS)
                    *phaseComp[(component*NUMPHASES + p)][idx];

                if (phase == p)
                    continue;

                dfdphi[phase][idx] -= 2.0*kappaPhi[phase*NUMPHASES + p]*(-phi[phase][xp[1]] + 16.0*phi[phase][xp[0]] - 30.0*phi[phase][idx] + 16.0*phi[phase][xm[0]] - phi[phase][xm[1]])/(12.0*DELTA_X*DELTA_X);
                if (sizeY > 1)
                    dfdphi[phase][idx] -= 2.0*kappaPhi[phase*NUMPHASES + p]*(-phi[phase][yp[1]] + 16.0*phi[phase][yp[0]] - 30.0*phi[phase][idx] + 16.0*phi[phase][ym[0]] - phi[phase][ym[1]])/(12.0*DELTA_Y*DELTA_Y);
                if (sizeZ > 1)
                    dfdphi[phase][idx] -= 2.0*kappaPhi[phase*NUMPHASES + p]*(-phi[phase][zp[1]] + 16.0*phi[phase][zp[0]] - 30.0*phi[phase][idx] + 16.0*phi[phase][zm[0]] - phi[phase][zm[1]])/(12.0*DELTA_Z*DELTA_Z);
            }

            // Potential function
            dfdphi[phase][idx] += calcDoubleWellDerivative(phi, phase,
                                                           theta_i, theta_ij, theta_ijk,
                                                           idx, NUMPHASES);


        }
    }
    __syncthreads();
}

__global__
void __computeDrivingForceBinary_02__(double **phi, double **comp,
                                      double **dfdphi,
                                      double **phaseComp,
                                      double *F0_A, double *F0_B, double *F0_C,
                                      double *diffusivity,
                                      double *theta_i, double *theta_ij, double *theta_ijk,
                                      double temperature, int *thermo_phase,
                                      int NUMPHASES, int NUMCOMPONENTS,
                                      int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    double mu;

    double interp_phi, interp_prime, g_prime;
    double etemp;
    double calpha, cbeta, f_alpha, f_beta;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        etemp  = phi[1][idx];

        calpha = phaseComp[0][idx];
        cbeta  = phaseComp[1][idx];

        interp_phi   = etemp*etemp*etemp*(6.0*etemp*etemp - 15.0*etemp + 10.0);
        interp_prime = 30.0*etemp*etemp*(1.0 - etemp)*(1.0 - etemp);
        g_prime      = 2.0*etemp*(1.0 - etemp)*(1.0 - 2.0*etemp);

        comp[0][idx]  = calpha*(1.0 - interp_phi) + cbeta*interp_phi;

        f_alpha = evalFunc(free_energy_tdb_dev[thermo_phase[0]], phaseComp[0][idx], temperature);
        f_beta = evalFunc(free_energy_tdb_dev[thermo_phase[1]], phaseComp[1][idx], temperature);

        mu = evalFunc(Mu_tdb_dev[thermo_phase[0]], phaseComp[0][idx], temperature);

        dfdphi[1][idx] = interp_prime * (f_beta - f_alpha - (cbeta - calpha)*mu) + theta_i[1]*g_prime;
    }
    __syncthreads();
}

__global__ void __computeDrivingForceBinary_01_03__(double **phi, double **comp,
                                                    double **dfdphi,
                                                    double **phaseComp,
                                                    double *F0_A, double *F0_B, double *F0_C,
                                                    double *diffusivity,
                                                    double *theta_i, double *theta_ij, double *theta_ijk,
                                                    int NUMPHASES, int NUMCOMPONENTS,
                                                    int sizeX, int sizeY, int sizeZ)
{

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;


    double interp_phi, interp_prime, g_prime;
    double etemp;
    double calpha, cbeta, f_alpha, f_beta;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        etemp  = phi[1][idx];

        calpha = phaseComp[0][idx];
        cbeta  = phaseComp[1][idx];

        interp_phi   = etemp*etemp*etemp*(6.0*etemp*etemp - 15.0*etemp + 10.0);
        interp_prime = 30.0*etemp*etemp*(1.0 - etemp)*(1.0 - etemp);
        g_prime      = 2.0*etemp*(1.0 - etemp)*(1.0 - 2.0*etemp);

        comp[0][idx]  = calpha*(1.0 - interp_phi) + cbeta*interp_phi;

        f_alpha      = F0_A[0]*calpha*calpha + F0_B[0]*calpha + F0_C[0];
        f_beta       = F0_A[1]*cbeta*cbeta + F0_B[1]*cbeta + F0_C[1];

        dfdphi[1][idx] = interp_prime * (f_beta - f_alpha - (cbeta - calpha)*(2.0*F0_A[1]*cbeta + F0_B[1])) + theta_ij[1]*g_prime;
    }
    __syncthreads();
}

void computeDrivingForce(double **phi, double **comp,
                         double **dfdphi,
                         double **phaseComp,
                         domainInfo* simDomain, controls* simControls,
                         simParameters* simParams, subdomainInfo* subdomain,
                         dim3 gridSize, dim3 blockSize)
{
    if (simControls->multiphase == 1 || simDomain->numPhases > 2 || simDomain->numComponents > 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3  || simControls->FUNCTION_F == 4)
        {
            __computeDrivingForce_01_03__<<<gridSize, blockSize>>>(phi, comp,
                                                                   dfdphi,
                                                                   phaseComp,
                                                                   simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                                   simParams->diffusivity_dev, simParams->kappaPhi_dev,
                                                                   simParams->theta_i_dev, simParams->theta_ij_dev, simParams->theta_ijk_dev,
                                                                   simDomain->numPhases, simDomain->numComponents,
                                                                   subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                                   simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z);
        }
        else if (simControls->FUNCTION_F == 2)
            printf("MPMC Exact not implemented\n");

    }
    else if (simDomain->numPhases == 2 && simDomain->numComponents == 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3  || simControls->FUNCTION_F == 4)
        {
            __computeDrivingForceBinary_01_03__<<<gridSize, blockSize>>>(phi, comp,
                                                                         dfdphi,
                                                                         phaseComp,
                                                                         simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                                         simParams->diffusivity_dev,
                                                                         simParams->theta_i_dev, simParams->theta_ij_dev, simParams->theta_ijk_dev,
                                                                         simDomain->numPhases, simDomain->numComponents,
                                                                         subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);
        }
        else if (simControls->FUNCTION_F == 2)
        {
            __computeDrivingForceBinary_02__<<<gridSize, blockSize>>>(phi, comp,
                                                                      dfdphi,
                                                                      phaseComp,
                                                                      simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                                      simParams->diffusivity_dev,
                                                                      simParams->theta_i_dev, simParams->theta_ij_dev, simParams->theta_ijk_dev,
                                                                      simParams->T, simDomain->thermo_phase_dev,
                                                                      simDomain->numPhases, simDomain->numComponents,
                                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);

        }
    }

    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }
}
