#include "calcPhaseComp.cuh"

__global__
void __initMu__(double **phi, double **comp, double **phaseComp, double **mu,
                long *thermo_phase, double temperature,
                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                long sizeX, long sizeY, long sizeZ,
                long xStep, long yStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        double y[MAX_NUM_COMP], mu0[MAX_NUM_COMP];
        double sum = 0.0;

        for (long phase = 0; phase < NUMPHASES; phase++)
        {
            // Bulk
            if (phi[phase][idx] == 1.0)
            {
                sum = 0.0;

                for (long i = 0; i < NUMCOMPONENTS-1; i++)
                {
                    phaseComp[phase + NUMPHASES*i][idx] = comp[i][idx];
                    y[i] = comp[i][idx];
                    sum += y[i];
                }

                y[NUMCOMPONENTS-1] = 1.0 - sum;

                (*Mu_tdb_dev[thermo_phase[phase]])(temperature, y, mu0);

                for (long i = 0; i < NUMCOMPONENTS-1; i++)
                    mu[i][idx] = mu0[i];
            }
            else
            {
                for (long i = 0; i < NUMCOMPONENTS-1; i++)
                    phaseComp[phase + NUMPHASES*i][idx] = 0.0;
            }
        }
    }
}

__global__
void __calcPhaseComp__(double **phi, double **comp,
                       double **phaseComp,
                       double *F0_A, double *F0_B, double *F0_C,
                       long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                       long sizeX, long sizeY, long sizeZ,
                       long xStep, long yStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k;

    // Number of phase compositions
    long N = NUMPHASES*(NUMCOMPONENTS-1);

    // Iterate to fill jacobian and function vector
    long iter1, iter2, iter3, iter4;
    long index1, index2;

    // Delta vector, function vector, jacobian matrix
    double sol[MAX_NUM_PHASE_COMP], B[MAX_NUM_PHASE_COMP];
    double A[MAX_NUM_PHASE_COMP][MAX_NUM_PHASE_COMP];

    // Permutation matrix required for LUP linear system solver
    int P[MAX_NUM_PHASE_COMP+1];

    // Tolerance for LU solver
    double tol = 1e-10;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        // Calculate function vector using x^n
        for (iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
        {
            for (iter2 = 0; iter2 < NUMPHASES; iter2++)
            {
                index1 = iter2 + iter1*NUMPHASES;
                if (iter2)
                {
                    // Constant part (linear part of free-energy) taken to the RHS of AX = B
                    // B^{K}_{P+1} - B^{K}_{P}
                    B[index1] = F0_B[iter1 + iter2*(NUMCOMPONENTS-1)] - F0_B[iter1 + (iter2-1)*(NUMCOMPONENTS-1)];
                }
                else
                {
                    // Composition conservation equation
                    B[index1] = comp[iter1][idx];
                }
            }
        }

        // Calculate jacobian using x^n
        for (iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
        {
            for (iter2 = 0; iter2 < NUMPHASES; iter2++)
            {
                index1 = iter2 + iter1*NUMPHASES;
                if (iter2)
                {
                    for (iter3 = 0; iter3 < NUMCOMPONENTS-1; iter3++)
                    {
                        for (iter4 = 0; iter4 < NUMPHASES; iter4++)
                        {
                            index2 = iter4 + iter3*NUMPHASES;

                            if (iter4 == iter2-1)
                            {
                                if (iter3 == iter1)
                                    A[index1][index2] = 2.0*F0_A[(iter3 + iter4*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter1];
                                else
                                    A[index1][index2] = F0_A[(iter3 + iter4*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter1];
                            }
                            else if (iter4 == iter2)
                            {
                                if (iter3 == iter1)
                                    A[index1][index2] = -2.0*F0_A[(iter3 + iter4*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter1];
                                else
                                    A[index1][index2] = -1.0*F0_A[(iter3 + iter4*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter1];
                            }
                            else
                                A[index1][index2] = 0.0;
                        }
                    }
                } // if (iter2)
                else
                {
                    for (iter3 = 0; iter3 < NUMCOMPONENTS-1; iter3++)
                    {
                        for (iter4 = 0; iter4 < NUMPHASES; iter4++)
                        {
                            index2 = iter4 + iter3*NUMPHASES;

                            if (iter3 == iter1)
                                A[index1][index2] = calcInterp5th(phi, iter4, idx, NUMPHASES);
                            else
                                A[index1][index2] = 0.0;
                        }
                    }
                }
            } // for (iter2 = 0 ... )
        } // for (iter1 = 0 ... )

        // Get x^(n+1) - x^(n)
        LUPDecomposePC1(A, N, tol, P);
        LUPSolvePC1(A, P, B, N, sol);

        for (iter1 = 0; iter1 < N; iter1++)
        {
            // Update phase composition
            phaseComp[iter1][idx] = sol[iter1];
        }
    }
}

__global__
void __calcPhaseComp_02__(double **phi, double **comp,
                          double **phaseComp, double **mu, double *cguess,
                          double temperature, long *thermo_phase,
                          long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                          long sizeX, long sizeY, long sizeZ,
                          long xStep, long yStep, long padding)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        double fun[MAX_NUM_COMP], jacInv[MAX_NUM_COMP][MAX_NUM_COMP], cn[MAX_NUM_COMP], co[MAX_NUM_COMP];
        double tmp0, norm;
        double retdmuphase[MAX_NUM_COMP*MAX_NUM_COMP], retdmuphase2[MAX_NUM_COMP][MAX_NUM_COMP], y[MAX_NUM_COMP], mu0[MAX_NUM_COMP];

        double tol = 1e-6;

        long interface = 1;
        long bulkphase;

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
            // Number of iterations for Newton-Raphson
            long count = 0;
            // Number of iterations for diffusion-potential correction
            long count2 = 0;

            long is, is1, is2;
            long maxCount = 10000;

            // Permutation matrix required by LU decomposition routine
            int P[MAX_NUM_COMP];

            double dmudc[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
            double dcdmu[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
            double Inv[MAX_NUM_COMP][MAX_NUM_COMP];

            double deltac[MAX_NUM_COMP] = {0.0};
            long deltac_flag = 0;

            do
            {
                count2++;

                for (long phase = 0; phase < NUMPHASES; phase++)
                {
                    for (is = 0; is < NUMCOMPONENTS-1; is++)
                    {
                        cn[is] = cguess[(phase + phase*(NUMPHASES))*(NUMCOMPONENTS-1) + is];
                    }

                    do
                    {
                        count++;

                        tmp0 = 0.0;

                        // Getting phase-compositions at the node
                        for (is = 0; is < NUMCOMPONENTS-1; is++)
                        {
                            co[is] = cn[is];
                            y[is]  = co[is];
                            tmp0  += co[is];
                        }
                        y[NUMCOMPONENTS-1] = 1.0 - tmp0;

                        // Getting local diffusion potential from tdb function
                        (*Mu_tdb_dev[thermo_phase[phase]])(temperature, y, mu0);

                        // Deviation of mu obtained from evolution from mu obtained from tdb
                        for (is = 0; is < NUMCOMPONENTS-1; is++)
                            fun[is] = (mu0[is] - mu[is][idx]);

                        // Second derivative of free-energy
                        (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, retdmuphase);

                        // Translating 2D array to 1D
                        for (is1 = 0; is1 < NUMCOMPONENTS-1; is1++)
                        {
                            for (is2 = 0; is2 < NUMCOMPONENTS-1; is2++)
                            {
                                retdmuphase2[is1][is2] = retdmuphase[is1*(NUMCOMPONENTS-1) + is2];
                            }
                        }

                        // Inverting dmudc to get dcdmu
                        LUPDecomposeC1(retdmuphase2, NUMCOMPONENTS-1, tol, P);
                        LUPInvertC1(retdmuphase2, P, NUMCOMPONENTS-1, jacInv);

                        // Newton-Raphson (-J^{-1}F)
                        for (is1 = 0; is1 < NUMCOMPONENTS-1; is1++)
                        {
                            tmp0 = 0.0;
                            for (is2 = 0; is2 < NUMCOMPONENTS-1; is2++)
                            {
                                tmp0 += jacInv[is1][is2] * fun[is2];
                            }

                            cn[is1] = co[is1] - tmp0;
                        }

                        // L-inf norm
                        norm = 0.0;
                        for (is = 0; is < NUMCOMPONENTS-1; is++)
                            if (fabs(cn[is] - co[is]) > 1e-6)
                                norm = 1.0;
                    } while (count < maxCount && norm > 0.0);

                    for (is = 0; is < NUMCOMPONENTS-1; is++)
                        phaseComp[is*NUMPHASES + phase][idx] = cn[is];
                }

                // Check conservation of comp
                deltac_flag = 0;
                for (is = 0; is < NUMCOMPONENTS-1; is++)
                {
                    deltac[is] = 0.0;

                    for (int phase = 0; phase < NUMPHASES; phase++)
                    {
                        deltac[is] += phaseComp[is*NUMPHASES + phase][idx]*calcInterp5th(phi, phase, idx, NUMPHASES);
                    }

                    deltac[is] = comp[is][idx] - deltac[is];

                    if (fabs(deltac[is]) > 1e-6)
                        deltac_flag = 1;
                }

                // deltac_flag will be 1 if not conserved
                // mu-correction will be carried out consequently, and the Newton-Raphson routine will be repeated

                if (deltac_flag)
                {
                    for (int component = 0; component < NUMCOMPONENTS-1; component++)
                    {
                        for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                        {
                            dcdmu[component*(NUMCOMPONENTS-1) + component2] = 0.0;
                        }
                    }

                    for (long phase = 0; phase < NUMPHASES; phase++)
                    {
                        double sum = 0.0;

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

                    for (int component = 0; component < NUMCOMPONENTS-1; component++)
                    {
                        for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                        {
                            mu[component][idx] += Inv[component][component2]*deltac[component2];
                        }
                    }
                }
            } while (count2 < 1000 && deltac_flag);
        }
        else
        {
            for (long component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (long phase = 0; phase < NUMPHASES; phase++)
                {
                    if (phase == bulkphase)
                        phaseComp[bulkphase + NUMPHASES*component][idx] = comp[component][idx];
                    else
                        phaseComp[phase + NUMPHASES*component][idx] = 0.0;
                }
            }
        }
    }
}

void calcPhaseComp(double **phi, double **comp,
                   double **phaseComp, double **mu,
                   domainInfo* simDomain, controls* simControls,
                   simParameters* simParams, subdomainInfo* subdomain,
                   dim3 gridSize, dim3 blockSize)
{
    if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3 || simControls->FUNCTION_F == 4)
    {
        __calcPhaseComp__<<<gridSize, blockSize>>>(phi, comp,
                                                   phaseComp,
                                                   simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                   simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                   subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                   subdomain->xStep, subdomain->yStep, subdomain->padding);
    }
    else if (simControls->FUNCTION_F == 2)
    {
        if (simControls->startTime == simControls->count && (simControls->restart == 0 && simControls->startTime == 0))
        {
            __initMu__<<<gridSize, blockSize>>>(phi, comp, phaseComp, mu,
                                                simDomain->thermo_phase_dev, simParams->Teq,
                                                simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                subdomain->xStep, subdomain->yStep, subdomain->padding);
        }
        __calcPhaseComp_02__<<<gridSize, blockSize>>>(phi, comp,
                                                      phaseComp, mu, simParams->cguess_dev,
                                                      simParams->T, simDomain->thermo_phase_dev,
                                                      simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                      subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                      subdomain->xStep, subdomain->yStep, subdomain->padding);
    }
}
