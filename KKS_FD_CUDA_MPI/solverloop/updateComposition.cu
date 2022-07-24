#include "updateComposition.cuh"

__device__
double calcInterp5th2(double **phi, int a, int idx, int NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0;

    double ans = 0.0, temp = 0.0;
    double phiValue = phi[a][idx];
    double const1 = 7.5*(NUMPHASES-2.0)/(NUMPHASES-1.0);

    ans  = pow(phiValue, 5)*(6.0 - const1);
    ans += pow(phiValue, 4)*(-15.0 + 3.0*const1);

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
    temp *= 7.5*(phiValue - 1.0)/((double)NUMPHASES-1.0);
    temp += phiValue*(10.0 - 3.0*const1) + const1;
    ans  += pow(phiValue, 2)*temp;

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
        ans += 15.0*phiValue*phiValue*temp;
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
        ans += 24.0*phiValue*temp;
    }

    return ans;
}

__device__
double calcDiffusionPotential2(double **phaseComp,
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

__global__
void __updateComposition__(double **phi,
                           double **comp, double **compNew,
                           double **phaseComp,
                           double *F0_A, double *F0_B,
                           double *diffusivity,
                           int NUMPHASES, int NUMCOMPONENTS,
                           int sizeX, int sizeY, int sizeZ,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp[2], xm[2], yp[2], ym[2], zp[2], zm[2];

    double mu[13];
    double mobility[7];
    double RHS = 0.0;

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

        for (int component = 0; component < NUMCOMPONENTS-1; component++)
        {
            mu[0] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS);

            mu[1] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, xp[0], NUMPHASES, NUMCOMPONENTS);
            mu[2] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, xp[1], NUMPHASES, NUMCOMPONENTS);
            mu[3] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, xm[0], NUMPHASES, NUMCOMPONENTS);
            mu[4] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, xm[1], NUMPHASES, NUMCOMPONENTS);

            mobility[0] = 0.0;

            mobility[1] = 0.0;
            mobility[2] = 0.0;

            for (int phase = 0; phase < NUMPHASES; phase++)
            {
                mobility[0] += calcInterp5th2(phi, phase, idx, NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);

                mobility[1] += calcInterp5th2(phi, phase, xp[0], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
                mobility[2] += calcInterp5th2(phi, phase, xp[1], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
            }

            if (sizeY > 1)
            {
                mobility[3] = 0.0;
                mobility[4] = 0.0;

                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    mobility[3] += calcInterp5th2(phi, phase, yp[0], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                    /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
                    mobility[4] += calcInterp5th2(phi, phase, yp[1], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                    /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
                }

                mu[5] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, yp[0], NUMPHASES, NUMCOMPONENTS);
                mu[6] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, yp[1], NUMPHASES, NUMCOMPONENTS);
                mu[7] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, ym[0], NUMPHASES, NUMCOMPONENTS);
                mu[8] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, ym[1], NUMPHASES, NUMCOMPONENTS);
            }

            if (sizeZ > 1)
            {
                mobility[5] = 0.0;
                mobility[6] = 0.0;

                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    mobility[5] += calcInterp5th2(phi, phase, zp[0], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                    /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
                    mobility[6] += calcInterp5th2(phi, phase, zp[1], NUMPHASES)*diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]
                    /(2.0*F0_A[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]);
                }

                mu[9] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, zp[0], NUMPHASES, NUMCOMPONENTS);
                mu[10] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, zp[1], NUMPHASES, NUMCOMPONENTS);
                mu[11] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, zm[0], NUMPHASES, NUMCOMPONENTS);
                mu[12] = calcDiffusionPotential2(phaseComp, 0, component, F0_A, F0_B, zm[1], NUMPHASES, NUMCOMPONENTS);
            }

            RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[3] + 1.0*mu[4])/(4.0*DELTA_X*DELTA_X);
            RHS +=  4.0*mobility[1]*(3.0*mu[1] - 4.0*mu[0] + 1.0*mu[3])/(4.0*DELTA_X*DELTA_X);
            RHS += -1.0*mobility[2]*(3.0*mu[2] - 4.0*mu[1] + 1.0*mu[0])/(4.0*DELTA_X*DELTA_X);

            if (sizeY > 1)
            {
                RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[7] + 1.0*mu[8])/(4.0*DELTA_Y*DELTA_Y);
                RHS +=  4.0*mobility[3]*(3.0*mu[5] - 4.0*mu[0] + 1.0*mu[7])/(4.0*DELTA_Y*DELTA_Y);
                RHS += -1.0*mobility[4]*(3.0*mu[6] - 4.0*mu[5] + 1.0*mu[0])/(4.0*DELTA_Y*DELTA_Y);
            }

            if (sizeZ > 1)
            {
                RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[11] + 1.0*mu[12])/(4.0*DELTA_Z*DELTA_Z);
                RHS +=  4.0*mobility[5]*(3.0*mu[9] - 4.0*mu[0] + 1.0*mu[11])/(4.0*DELTA_Z*DELTA_Z);
                RHS += -1.0*mobility[6]*(3.0*mu[10] - 4.0*mu[9] + 1.0*mu[0])/(4.0*DELTA_Z*DELTA_Z);
            }

            compNew[component][idx] = comp[component][idx] + RHS*DELTA_t;
        }
    }
    __syncthreads();
}

__global__
void __updateCompositionBinary__(double **phi,
                                 double **comp, double **compNew,
                                 double **phaseComp,
                                 double *F0_A, double *F0_B,
                                 double *diffusivity,
                                 int NUMPHASES, int NUMCOMPONENTS,
                                 int sizeX, int sizeY, int sizeZ,
                                 double DELTA_X, double DELTA_Y, double DELTA_Z,
                                 double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp[2], xm[2], yp[2], ym[2], zp[2], zm[2];

    double mu[13];
    double mobility[7];
    double RHS = 0.0;

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

        mu[0] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS);

        mu[1] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, xp[0], NUMPHASES, NUMCOMPONENTS);
        mu[2] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, xp[1], NUMPHASES, NUMCOMPONENTS);
        mu[3] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, xm[0], NUMPHASES, NUMCOMPONENTS);
        mu[4] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, xm[1], NUMPHASES, NUMCOMPONENTS);

        mobility[0] = 0.0;

        mobility[1] = 0.0;
        mobility[2] = 0.0;

        for (int phase = 0; phase < NUMPHASES; phase++)
        {
            mobility[0] += calcInterp5th2(phi, phase, idx, NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
            /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);

            mobility[1] += calcInterp5th2(phi, phase, xp[0], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
            /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
            mobility[2] += calcInterp5th2(phi, phase, xp[1], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
            /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
        }

        if (sizeY > 1)
        {
            mobility[3] = 0.0;
            mobility[4] = 0.0;

            for (int phase = 0; phase < NUMPHASES; phase++)
            {
                mobility[3] += calcInterp5th2(phi, phase, yp[0], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
                /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
                mobility[4] += calcInterp5th2(phi, phase, yp[1], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
                /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
            }

            mu[5] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, yp[0], NUMPHASES, NUMCOMPONENTS);
            mu[6] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, yp[1], NUMPHASES, NUMCOMPONENTS);
            mu[7] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, ym[0], NUMPHASES, NUMCOMPONENTS);
            mu[8] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, ym[1], NUMPHASES, NUMCOMPONENTS);
        }

        if (sizeZ > 1)
        {
            mobility[5] = 0.0;
            mobility[6] = 0.0;

            for (int phase = 0; phase < NUMPHASES; phase++)
            {
                mobility[5] += calcInterp5th2(phi, phase, zp[0], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
                /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
                mobility[6] += calcInterp5th2(phi, phase, zp[1], NUMPHASES)*diffusivity[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]
                /(2.0*F0_A[(phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1)]);
            }

            mu[9] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, zp[0], NUMPHASES, NUMCOMPONENTS);
            mu[10] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, zp[1], NUMPHASES, NUMCOMPONENTS);
            mu[11] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, zm[0], NUMPHASES, NUMCOMPONENTS);
            mu[12] = calcDiffusionPotential2(phaseComp, 0, 0, F0_A, F0_B, zm[1], NUMPHASES, NUMCOMPONENTS);
        }

        RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[3] + 1.0*mu[4])/(4.0*DELTA_X*DELTA_X);
        RHS +=  4.0*mobility[1]*(3.0*mu[1] - 4.0*mu[0] + 1.0*mu[3])/(4.0*DELTA_X*DELTA_X);
        RHS += -1.0*mobility[2]*(3.0*mu[2] - 4.0*mu[1] + 1.0*mu[0])/(4.0*DELTA_X*DELTA_X);

        if (sizeY > 1)
        {
            RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[7] + 1.0*mu[8])/(4.0*DELTA_Y*DELTA_Y);
            RHS +=  4.0*mobility[3]*(3.0*mu[5] - 4.0*mu[0] + 1.0*mu[7])/(4.0*DELTA_Y*DELTA_Y);
            RHS += -1.0*mobility[4]*(3.0*mu[6] - 4.0*mu[5] + 1.0*mu[0])/(4.0*DELTA_Y*DELTA_Y);
        }

        if (sizeZ > 1)
        {
            RHS += -3.0*mobility[0]*(3.0*mu[0] - 4.0*mu[11] + 1.0*mu[12])/(4.0*DELTA_Z*DELTA_Z);
            RHS +=  4.0*mobility[5]*(3.0*mu[9] - 4.0*mu[0] + 1.0*mu[11])/(4.0*DELTA_Z*DELTA_Z);
            RHS += -1.0*mobility[6]*(3.0*mu[10] - 4.0*mu[9] + 1.0*mu[0])/(4.0*DELTA_Z*DELTA_Z);
        }

        compNew[0][idx] = comp[0][idx] + RHS*DELTA_t;
    }
    __syncthreads();
}

void updateComposition(double **phi, double **comp, double **compNew,
                       double **phaseComp,
                       domainInfo* simDomain, controls* simControls,
                       simParameters* simParams, subdomainInfo* subdomain,
                       dim3 gridSize, dim3 blockSize)
{
    if (simControls->multiphase == 1 || simDomain->numPhases > 2 || simDomain->numComponents > 2)
    {
        __updateComposition__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                       phaseComp,
                                                       simParams->F0_A_dev, simParams->F0_B_dev,
                                                       simParams->diffusivity_dev,
                                                       simDomain->numPhases, simDomain->numComponents,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                       simControls->DELTA_t);
    }
    else if (simDomain->numPhases == 2 && simDomain->numComponents == 2)
    {
        __updateCompositionBinary__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                             phaseComp,
                                                             simParams->F0_A_dev, simParams->F0_B_dev,
                                                             simParams->diffusivity_dev,
                                                             simDomain->numPhases, simDomain->numComponents,
                                                             subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                             simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                             simControls->DELTA_t);
    }
}
