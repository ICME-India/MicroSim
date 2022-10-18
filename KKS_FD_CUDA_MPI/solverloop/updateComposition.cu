#include "updateComposition.cuh"

__global__
void __updateComposition__(double **phi,
                           double **comp, double **compNew,
                           double **phaseComp,
                           double *F0_A, double *F0_B,
                           double *mobility,
                           long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                           long sizeX, long sizeY, long sizeZ,
                           long yStep, long zStep, long padding,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx[7] = {-1};

    idx[0] = k*zStep + j*yStep + i;

    double mu[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        // x-direction
        idx[1] = k*zStep + j*yStep + i+1;
        idx[2] = k*zStep + j*yStep + i-1;

        // y-direction
        if (DIMENSION >= 2)
        {
            idx[3] = k*zStep + (j+1)*yStep + i;
            idx[4] = k*zStep + (j-1)*yStep + i;
        }

        // z-direction
        if (DIMENSION == 3)
        {
            idx[5] = (k+1)*zStep + j*yStep + i;
            idx[6] = (k-1)*zStep + j*yStep + i;
        }

        for (long component = 0; component < NUMCOMPONENTS-1; component++)
        {
            J_xp = 0.0;
            J_xm = 0.0;
            J_yp = 0.0;
            J_ym = 0.0;
            J_zp = 0.0;
            J_zm = 0.0;

            for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
            {

                for (long iter = 0; iter < 7; iter++)
                    effMobility[iter] = 0.0;

                for (long phase = 0; phase < NUMPHASES; phase++)
                {
                    effMobility[0] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[0], NUMPHASES);

                    effMobility[1] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[1], NUMPHASES);
                    effMobility[2] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[2], NUMPHASES);

                    if (DIMENSION >= 2)
                    {
                        effMobility[3] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[3], NUMPHASES);
                        effMobility[4] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[4], NUMPHASES);
                    }

                    if (DIMENSION == 3)
                    {
                        effMobility[5] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[5], NUMPHASES);
                        effMobility[6] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th(phi, phase, idx[6], NUMPHASES);
                    }
                }


                mu[0] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[0], NUMPHASES, NUMCOMPONENTS);

                mu[1] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[1], NUMPHASES, NUMCOMPONENTS);
                mu[2] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[2], NUMPHASES, NUMCOMPONENTS);

                if (DIMENSION >= 2)
                {
                    mu[3] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[3], NUMPHASES, NUMCOMPONENTS);
                    mu[4] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[4], NUMPHASES, NUMCOMPONENTS);
                }

                if (DIMENSION == 3)
                {
                    mu[5] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[5], NUMPHASES, NUMCOMPONENTS);
                    mu[6] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, idx[6], NUMPHASES, NUMCOMPONENTS);
                }

                J_xp += ((effMobility[1] + effMobility[0])/2.0)*(mu[1] - mu[0])/DELTA_X;
                J_xm += ((effMobility[0] + effMobility[2])/2.0)*(mu[0] - mu[2])/DELTA_X;

                if (DIMENSION >= 2)
                {
                    J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[3] - mu[0])/DELTA_Y;
                    J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0] - mu[4])/DELTA_Y;
                }

                if (DIMENSION == 3)
                {
                    J_zp += ((effMobility[5] + effMobility[0])/2.0)*(mu[5] - mu[0])/DELTA_Z;
                    J_zm += ((effMobility[0] + effMobility[6])/2.0)*(mu[0] - mu[6])/DELTA_Z;
                }
            }

            compNew[component][idx[0]] = comp[component][idx[0]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);
        }
    }
}

__global__
void __updateComposition_02__(double **phi,
                              double **comp, double **compNew, double **mu,
                              double **phaseComp, long *thermo_phase,
                              double *diffusivity, double temperature, double molarVolume,
                              long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                              long sizeX, long sizeY, long sizeZ,
                              long yStep, long zStep, long padding,
                              double DELTA_X, double DELTA_Y, double DELTA_Z,
                              double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx[7] = {-1};

    idx[0] = k*zStep + j*yStep + i;

    double muLocal[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;
    double tol = 1e-6;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        // x-direction
        idx[1] = k*zStep + j*yStep + i+1;
        idx[2] = k*zStep + j*yStep + i-1;

        // y-direction
        if (DIMENSION >= 2)
        {
            idx[3] = k*zStep + (j+1)*yStep + i;
            idx[4] = k*zStep + (j-1)*yStep + i;
        }

        // z-direction
        if (DIMENSION == 3)
        {
            idx[5] = (k+1)*zStep + j*yStep + i;
            idx[6] = (k-1)*zStep + j*yStep + i;
        }

        double dmudc[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
        double y[MAX_NUM_COMP];
        double dmudcInv[MAX_NUM_COMP][MAX_NUM_COMP];
        int P[MAX_NUM_COMP];
        double mobility[MAX_NUM_COMP][MAX_NUM_COMP];

        long maxPos;

        long interface = 1, bulkphase = 0;

        for (long is = 0; is < NUMPHASES; is++)
        {
            if (phi[is][idx[0]] > 0.99999)
            {
                bulkphase = is;
                interface = 0;
                break;
            }
        }

        if (interface)
        {
            for (int component = 0; component < NUMCOMPONENTS-1; component++)
            {
                // Fluxes
                J_xp = 0.0;
                J_xm = 0.0;
                J_yp = 0.0;
                J_ym = 0.0;
                J_zp = 0.0;
                J_zm = 0.0;

                // Computing the inner derivative and mobilities to get the fluxes
                for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    for (long iter = 0; iter < 7; iter++)
                        effMobility[iter] = 0.0;

                    if (DIMENSION == 3)
                        maxPos = 7;
                    else if (DIMENSION == 2)
                        maxPos = 5;
                    else
                        maxPos = 3;

                    for (long pos = 0; pos < maxPos; pos++)
                    {
                        // M_{ij} = \sum_{\phi} M(\phi) = \sum_{\phi} D*dcdmu
                        for (long phase = 0; phase < NUMPHASES; phase++)
                        {
                            double tmp0 = 0.0;

                            for (long is = 0; is < NUMCOMPONENTS-1; is++)
                            {
                                y[is] = phaseComp[is*NUMPHASES + phase][idx[pos]];
                                tmp0  += y[is];
                            }

                            y[NUMCOMPONENTS-1] = 1.0 - tmp0;


                            // Get dmudc for the current phase
                            (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, dmudc);

                            // Invert dmudc to get dcdmu for the current phase
                            LUPDecomposeC2(dmudc, NUMCOMPONENTS-1, tol, P);
                            LUPInvertC2(dmudc, P, NUMCOMPONENTS-1, dmudcInv);

                            // multiply diffusivity with dcdmu
                            for (long iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
                            {
                                for (long iter2 = 0; iter2 < NUMCOMPONENTS-1; iter2++)
                                {
                                    mobility[iter1][iter2] = 0.0;

                                    for (long iter3 = 0; iter3 < NUMCOMPONENTS-1; iter3++)
                                    {
                                        mobility[iter1][iter2] += diffusivity[(iter1 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter3]*dmudcInv[iter3][iter2];
                                    }
                                }
                            }

                            // Summing over all phases, weighting with the interpolation fn.
                            effMobility[pos] += mobility[component][component2]*calcInterp5th(phi, phase, idx[pos], NUMPHASES);
                        }
                    }


                    muLocal[0] = mu[component2][idx[0]];

                    muLocal[1] = mu[component2][idx[1]];
                    muLocal[2] = mu[component2][idx[2]];

                    if (DIMENSION >= 2)
                    {
                        muLocal[3] = mu[component2][idx[3]];
                        muLocal[4] = mu[component2][idx[4]];
                    }
                    if (DIMENSION == 3)
                    {
                        muLocal[5] = mu[component2][idx[5]];
                        muLocal[6] = mu[component2][idx[6]];
                    }

                    J_xp += ((effMobility[1] + effMobility[0])/2.0)*(muLocal[1] - muLocal[0])/DELTA_X;
                    J_xm += ((effMobility[0] + effMobility[2])/2.0)*(muLocal[0] - muLocal[2])/DELTA_X;

                    if (DIMENSION >= 2)
                    {
                        J_yp += ((effMobility[3] + effMobility[0])/2.0)*(muLocal[3] - muLocal[0])/DELTA_Y;
                        J_ym += ((effMobility[0] + effMobility[4])/2.0)*(muLocal[0] - muLocal[4])/DELTA_Y;
                    }

                    if (DIMENSION == 3)
                    {
                        J_zp += ((effMobility[5] + effMobility[0])/2.0)*(muLocal[5] - muLocal[0])/DELTA_Z;
                        J_zm += ((effMobility[0] + effMobility[6])/2.0)*(muLocal[0] - muLocal[6])/DELTA_Z;
                    }
                }

                compNew[component][idx[0]] = comp[component][idx[0]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);
            }
        }
        else
        {
            for (int component = 0; component < NUMCOMPONENTS-1; component++)
            {
                // Fluxes
                J_xp = 0.0;
                J_xm = 0.0;
                J_yp = 0.0;
                J_ym = 0.0;
                J_zp = 0.0;
                J_zm = 0.0;

                // Computing the inner derivative and mobilities to get the fluxes
                for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    for (long iter = 0; iter < 7; iter++)
                        effMobility[iter] = 0.0;

                    if (DIMENSION == 3)
                        maxPos = 7;
                    else if (DIMENSION == 2)
                        maxPos = 5;
                    else
                        maxPos = 3;

                    for (long pos = 0; pos < maxPos; pos++)
                    {
                        // M_{ij} = \sum_{\phi} M(\phi) = \sum_{\phi} D*dcdmu
                        double tmp0 = 0.0;

                        for (long is = 0; is < NUMCOMPONENTS-1; is++)
                        {
                            y[is] = phaseComp[is*NUMPHASES + bulkphase][idx[pos]];
                            tmp0  += y[is];
                        }

                        y[NUMCOMPONENTS-1] = 1.0 - tmp0;


                        // Get dmudc for the current phase
                        (*dmudc_tdb_dev[thermo_phase[bulkphase]])(temperature, y, dmudc);

                        // Invert dmudc to get dcdmu for the current phase
                        LUPDecomposeC2(dmudc, NUMCOMPONENTS-1, tol, P);
                        LUPInvertC2(dmudc, P, NUMCOMPONENTS-1, dmudcInv);

                        // multiply diffusivity with dcdmu
                        for (long iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
                        {
                            for (long iter2 = 0; iter2 < NUMCOMPONENTS-1; iter2++)
                            {
                                mobility[iter1][iter2] = 0.0;

                                for (long iter3 = 0; iter3 < NUMCOMPONENTS-1; iter3++)
                                {
                                    mobility[iter1][iter2] += diffusivity[(iter1 + bulkphase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter3]*dmudcInv[iter3][iter2];
                                }
                            }
                        }

                        // Summing over all phases, weighting with the interpolation fn.
                        effMobility[pos] += mobility[component][component2]*calcInterp5th(phi, bulkphase, idx[pos], NUMPHASES);
                    }


                    muLocal[0] = mu[component2][idx[0]];

                    muLocal[1] = mu[component2][idx[1]];
                    muLocal[2] = mu[component2][idx[2]];

                    if (DIMENSION >= 2)
                    {
                        muLocal[3] = mu[component2][idx[3]];
                        muLocal[4] = mu[component2][idx[4]];
                    }
                    if (DIMENSION == 3)
                    {
                        muLocal[5] = mu[component2][idx[5]];
                        muLocal[6] = mu[component2][idx[6]];
                    }

                    J_xp += ((effMobility[1] + effMobility[0])/2.0)*(muLocal[1] - muLocal[0])/DELTA_X;
                    J_xm += ((effMobility[0] + effMobility[2])/2.0)*(muLocal[0] - muLocal[2])/DELTA_X;

                    if (DIMENSION >= 2)
                    {
                        J_yp += ((effMobility[3] + effMobility[0])/2.0)*(muLocal[3] - muLocal[0])/DELTA_Y;
                        J_ym += ((effMobility[0] + effMobility[4])/2.0)*(muLocal[0] - muLocal[4])/DELTA_Y;
                    }

                    if (DIMENSION == 3)
                    {
                        J_zp += ((effMobility[5] + effMobility[0])/2.0)*(muLocal[5] - muLocal[0])/DELTA_Z;
                        J_zm += ((effMobility[0] + effMobility[6])/2.0)*(muLocal[0] - muLocal[6])/DELTA_Z;
                    }
                }

                compNew[component][idx[0]] = comp[component][idx[0]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);
            }
        }
    }
}

__global__
void __updateMu_02__(double **phi, double **comp,
                     double **phiNew, double **compNew,
                     double **phaseComp, double **mu,
                     long *thermo_phase, double temperature, double molarVolume,
                     long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                     long sizeX, long sizeY, long sizeZ,
                     long yStep, long zStep, long padding,
                     double DELTA_X, double DELTA_Y, double DELTA_Z,
                     double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = k*zStep + j*yStep + i;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        double RHS[MAX_NUM_COMP] = {0.0}, sum = 0.0;
        double tol = 1e-6;

        long bulkphase = 0, interface = 1;

        for (long is = 0; is < NUMPHASES; is++)
        {
            if (phi[is][idx] > 0.99999)
            {
                bulkphase = is;
                interface = 0;
                break;
            }
        }

        if (interface)
        {
            double dmudc[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
            double dcdmu[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
            double y[MAX_NUM_COMP];
            double Inv[MAX_NUM_COMP][MAX_NUM_COMP];
            int P[MAX_NUM_COMP];

            for (long component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    dcdmu[component*(NUMCOMPONENTS-1) + component2] = 0.0;
                }
            }

            for (long phase = 0; phase < NUMPHASES; phase++)
            {
                sum = 0.0;

                for (long component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    y[component] = phaseComp[component*NUMPHASES + phase][idx];
                    sum += y[component];
                }

                y[NUMCOMPONENTS-1] = 1.0 - sum;

                (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, dmudc);

                LUPDecomposeC2(dmudc, NUMCOMPONENTS-1, tol, P);
                LUPInvertC2(dmudc, P, NUMCOMPONENTS-1, Inv);

                for (long component = 0; component < NUMCOMPONENTS-1; component++)
                    for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                        dcdmu[component*(NUMCOMPONENTS-1) + component2] += calcInterp5th(phi, phase, idx, NUMPHASES)*Inv[component][component2];
            }

            LUPDecomposeC2(dcdmu, NUMCOMPONENTS-1, tol, P);
            LUPInvertC2(dcdmu, P, NUMCOMPONENTS-1, Inv);

            for (long component = 0; component < NUMCOMPONENTS-1; component++)
            {
                RHS[component] = (compNew[component][idx] - comp[component][idx]);

                for (long phase = 0; phase < NUMPHASES; phase++)
                {
                    sum = 0.0;

                    for (long phase2 = 0; phase2 < NUMPHASES; phase2++)
                    {
                        sum += calcInterp5thDiff(phi, phase, phase2, idx, NUMPHASES)*(phiNew[phase2][idx] - phi[phase2][idx]);
                    }

                    RHS[component] -= phaseComp[phase + component*NUMPHASES][idx]*sum;
                }
            }

            for (long component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (long component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    mu[component][idx] += Inv[component][component2]*RHS[component2];
                }
            }
        }
        else
        {
            double y[MAX_NUM_COMP];
            double mu1[MAX_NUM_COMP];

            sum = 0.0;

            for (long is = 0; is < NUMCOMPONENTS-1; is++)
            {
                y[is] = compNew[is][idx];
                sum += y[is];
            }

            y[NUMCOMPONENTS-1] = 1.0 - sum;

            (*Mu_tdb_dev[thermo_phase[bulkphase]])(temperature, y, mu1);

            for (long is = 0; is < NUMCOMPONENTS-1; is++)
                mu[is][idx] = mu1[is];
        }
    }
}

void updateComposition(double **phi, double **comp, double **phiNew, double **compNew,
                       double **phaseComp, double **mu,
                       domainInfo* simDomain, controls* simControls,
                       simParameters* simParams, subdomainInfo* subdomain,
                       dim3 gridSize, dim3 blockSize)
{
    if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3 || simControls->FUNCTION_F == 4)
    {
        __updateComposition__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                       phaseComp,
                                                       simParams->F0_A_dev, simParams->F0_B_dev,
                                                       simParams->mobility_dev,
                                                       simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                       subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                       subdomain->yStep, subdomain->zStep, subdomain->padding,
                                                       simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                       simControls->DELTA_t);

        applyBoundaryCondition(compNew, 2, simDomain->numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);
    }
    else if (simControls->FUNCTION_F == 2)
    {
        __updateComposition_02__<<<gridSize, blockSize>>>(phi, comp,
                                                          compNew, mu,
                                                          phaseComp, simDomain->thermo_phase_dev,
                                                          simParams->diffusivity_dev, simParams->T, simParams->molarVolume,
                                                          simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                          subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                          subdomain->yStep, subdomain->zStep, subdomain->padding,
                                                          simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                          simControls->DELTA_t);

        applyBoundaryCondition(compNew, 2, simDomain->numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        __updateMu_02__<<<gridSize, blockSize>>>(phi, comp,
                                                 phiNew, compNew,
                                                 phaseComp, mu,
                                                 simDomain->thermo_phase_dev, simParams->T, simParams->molarVolume,
                                                 simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                 subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                 subdomain->yStep, subdomain->zStep, subdomain->padding,
                                                 simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                 simControls->DELTA_t);

        applyBoundaryCondition(mu, 1, simDomain->numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

    }
}
