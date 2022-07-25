#include "calcPhaseComp.cuh"

__device__
int LUPDecompose(double A[][NUM_PHASE_COMP], int N, double Tol, int *P)
{
    int i, j, k, imax;
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

__device__
void LUPSolve(double A[][NUM_PHASE_COMP], int *P, double *b, int N, double *x)
{
    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

__device__
double calcInterp5th(double **phi, int a, int idx, int NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double ans = 0.0, temp = 0.0;
    double phiValue = phi[a][idx];
    double const1 = 7.5*((double)NUMPHASES-2.0)/((double)NUMPHASES-1.0);

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

__global__
void __calcPhaseComp_01_03__(double **phi, double **comp,
                             double **phaseComp,
                             double *F0_A, double *F0_B, double *F0_C,
                             int NUMPHASES, int NUMCOMPONENTS,
                             int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    // Number of phase compositions
    int N = NUMPHASES*(NUMCOMPONENTS-1);

    // Iterate to fill jacobian and function vector
    int iter1, iter2, iter3, iter4;
    int index1, index2;

    // Delta vector, function vector, jacobian matrix
    double sol[NUM_PHASE_COMP], B[NUM_PHASE_COMP];
    double A[NUM_PHASE_COMP][NUM_PHASE_COMP];

    // Permutation matrix required for LUP linear system solver
    int P[NUM_PHASE_COMP+1];

    // Tolerance for LU solver
    double tol = 1e-10;

    if (i < sizeX && j < sizeY && k < sizeZ)
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
        LUPDecompose(A, N, tol, P);
        LUPSolve(A, P, B, N, sol);

        for (iter1 = 0; iter1 < N; iter1++)
        {
            // Update phase composition
            phaseComp[iter1][idx] = sol[iter1];
        }

        for (int component = 0; component < NUMCOMPONENTS-1; component++)
        {
            comp[component][idx] = 0.0;

            for (int phase = 0; phase < NUMPHASES; phase++)
                comp[component][idx] += phaseComp[(phase + component*NUMPHASES)][idx]*calcInterp5th(phi, phase, idx, NUMPHASES);
        }
    }

    __syncthreads();
}

__device__
double newtonRaphson(void f(double, double*, double*), double f_const,
                     void df(double , double *, double *),
                     double temperature, double molarVolume,
                     double initialGuess)
{
    return 0;
    int iter = 0, maxIter = 1000;
    double tol = 1.0e-6;

    double ans = initialGuess;
    double delta = 1.0;

    double URF = 1.0e-0;

    double fValue, dfValue;

    do
    {
        iter++;

        if (ans < 1.0e-3)
            fValue = evalFunc(df, 1.0e-3, temperature)*ans/molarVolume + evalFunc(f, 1.0e-3, temperature)/molarVolume - f_const;
        else if (ans > 1.0 - 1.0e-3)
            fValue = evalFunc(df, 1.0 - 1.0e-3, temperature)*ans/molarVolume + evalFunc(f, 1.0 - 1.0e-3, temperature)/molarVolume - f_const;
        else
            fValue = evalFunc(f, ans, temperature)/molarVolume - f_const;

        if (ans < 1.0e-3)
            dfValue = evalFunc(df, 1.0e-3, temperature)/molarVolume;
        else if (ans > 1.0 - 1.0e-3)
            dfValue = evalFunc(df, 1.0 - 1.0e-3, temperature)/molarVolume;
        else
            dfValue = evalFunc(df, ans, temperature)/molarVolume;

        delta = -fValue/dfValue;

        ans += URF*delta;

    } while (fabs(fValue) > tol && iter < maxIter);

    return ans;
}

__global__
void initPhaseComp_02(double **phi, double **comp,
                      double **phaseComp, double *cguess,
                      double temperature, double molarVolume, int *thermo_phase,
                      int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        for (int phase = 0; phase < 2; phase++)
        {
                phaseComp[phase][idx] = comp[0][idx];

            for (long phase2 = 0; phase2 < 2; phase2++)
            {
                if (phase2 == phase)
                    continue;

                phaseComp[phase2][idx] = newtonRaphson(Mu_tdb_dev[thermo_phase[phase]],
                                                       evalFunc(Mu_tdb_dev[thermo_phase[phase]], phaseComp[phase][idx], temperature),
                                                       dmudc_tdb_dev[thermo_phase[phase]],
                                                       temperature, molarVolume,
                                                       cguess[(phase + phase2*2)]);
            }
        }
    }
}

__global__
void __calcPhaseCompBinary_02__(double **phi, double **comp,
                                double **phaseComp,
                                double temperature, double molarVolume, int *thermo_phase,
                                int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        double  interp_phi;
        double  ctemp, etemp;

        int iter    = 0;
        int maxIter = 1000;
        double tol  = 1e-8;

        double f[2], g[4];
        double z[2], z0;
        double h[3];
        double alpha[4];
        double ans[2];

        ans[0] = phaseComp[0][idx];
        ans[1] = phaseComp[1][idx];

        ctemp = comp[0][idx];
        etemp = phi[1][idx];

        interp_phi = etemp*etemp*etemp*(6.0*etemp*etemp - 15.0*etemp + 10.0);

        if (etemp > 0.01 && etemp < 0.99)
        {
            while (iter++ < maxIter)
            {
                /*** Step 3 **/
                /*** f[0] ***/
                f[0] = (1.0 - interp_phi)*ans[0] + interp_phi*ans[1] - ctemp;

                /*** f[1] ***/
                f[1] = (evalFunc(Mu_tdb_dev[thermo_phase[0]], ans[0], temperature) - evalFunc(Mu_tdb_dev[thermo_phase[1]], ans[1], temperature))/molarVolume;

                g[1] = f[0]*f[0] + f[1]*f[1];

                z[0] = 2.0*((1.0 - interp_phi)*f[0] + evalFunc(dmudc_tdb_dev[thermo_phase[0]], ans[0], temperature)/molarVolume*f[1]);
                z[1] = 2.0*(interp_phi*f[0] - evalFunc(dmudc_tdb_dev[thermo_phase[1]], ans[1], temperature)/molarVolume*f[1]);

                z0 = sqrt(z[0]*z[0] + z[1]*z[1]);

                /*** Step 4 ***/
                if (z0 == 0.0)
                {
                    printf("z0 = 0\n");
                    break;
                }

                /*** Step 5 ***/
                z[0] /= z0;
                z[1] /= z0;
                alpha[1] = 0.0;
                alpha[3] = 1.0;

                /*** f[0] ***/
                f[0] = (1.0 - interp_phi)*(ans[0] - alpha[3]*z[0]) + interp_phi*(ans[1] - alpha[3]*z[1]) - ctemp;

                /*** f[1] ***/
                f[1] = (evalFunc(Mu_tdb_dev[thermo_phase[0]], ans[0] - alpha[3]*z[0], temperature) - evalFunc(Mu_tdb_dev[thermo_phase[1]], ans[1] - alpha[3]*z[1], temperature))/molarVolume;

                /*** g(3) ***/
                g[3] = f[0]*f[0] + f[1]*f[1];

                /** Step 6, 7, 8 ***/
                while (g[3] >= g[1])
                {
                    /*** Step 7 ***/
                    alpha[3] /= 2.0;

                    /*** f[0] ***/
                    f[0] = (1.0 - interp_phi)*(ans[0] - alpha[3]*z[0]) + interp_phi*(ans[1] - alpha[3]*z[1]) - ctemp;

                    /*** f[1] ***/
                    f[1] = (evalFunc(Mu_tdb_dev[thermo_phase[0]], ans[0] - alpha[3]*z[0], temperature) - evalFunc(Mu_tdb_dev[thermo_phase[1]], ans[1] - alpha[3]*z[1], temperature))/molarVolume;

                    /*** g(3) ***/
                    g[3] = f[0]*f[0] + f[1]*f[1];

                    /*** Step 8 ***/
                    if (alpha[3] < tol/2.0)
                    {
                        break;
                    }
                }

                if (alpha[3] < tol/2.0)
                {
                    break;
                }

                /*** Step 9 ***/
                alpha[2] = alpha[3]/2.0;

                /*** f[0] ***/
                f[0] = (1.0 - interp_phi)*(ans[0] - alpha[2]*z[0]) + interp_phi*(ans[1] - alpha[2]*z[1]) - ctemp;

                /*** f[1] ***/
                f[1] = (evalFunc(Mu_tdb_dev[thermo_phase[0]], ans[0] - alpha[2]*z[0], temperature) - evalFunc(Mu_tdb_dev[thermo_phase[1]], ans[1] - alpha[2]*z[1], temperature))/molarVolume;

                /*** g(2) ***/
                g[2] = f[0]*f[0] + f[1]*f[1];

                /*** Step 10 ***/
                h[0] = (g[2] - g[1])/alpha[2];
                h[1] = (g[3] - g[2])/(alpha[3] - alpha[2]);
                h[2] = (h[1] - h[0])/alpha[3];

                /*** Step 11 ***/
                alpha[0] = 0.5*(alpha[2] - h[0]/h[2]);

                /*** f[0] ***/
                f[0] = (1.0 - interp_phi)*(ans[0] - alpha[0]*z[0]) + interp_phi*(ans[1] - alpha[0]*z[1]) - ctemp;

                /*** f[1] ***/
                f[1] = (evalFunc(Mu_tdb_dev[thermo_phase[0]], ans[0] - alpha[0]*z[0], temperature) - evalFunc(Mu_tdb_dev[thermo_phase[1]], ans[1] - alpha[0]*z[1], temperature))/molarVolume;

                /*** g(2) ***/
                g[0] = f[0]*f[0] + f[1]*f[1];

                /*** Step 12, 13 ***/
                if (g[3] < g[0])
                {
                    g[0] = g[3];
                    ans[0] -= alpha[3]*z[0];
                    ans[1] -= alpha[3]*z[1];
                }
                else
                {
                    ans[0] -= alpha[0]*z[0];
                    ans[1] -= alpha[0]*z[1];
                }

                /*** Step 14 ***/
                if (fabs(g[0] - g[1]) < tol || (fabs(f[0]) < 1e-5 && fabs(f[1]) < 1e-5))
                {
                    break;
                }
            } // while (iter++ < maxIter)

            phaseComp[0][idx] = ans[0];
            phaseComp[1][idx] = ans[1];

            comp[0][idx] = ans[0]*(1.0 - interp_phi) + ans[1]*interp_phi;
        }
        else
        {
            if (etemp <= 0.01)
            {
                phaseComp[0][idx] = ctemp;
                phaseComp[1][idx] = newtonRaphson(Mu_tdb_dev[thermo_phase[1]],
                                                  evalFunc(Mu_tdb_dev[thermo_phase[0]], phaseComp[0][idx], temperature)/molarVolume,
                                                  dmudc_tdb_dev[thermo_phase[1]],
                                                  temperature, molarVolume,
                                                  phaseComp[1][idx]);
            }
            else
            {
                phaseComp[1][idx] = ctemp;
                phaseComp[0][idx] = newtonRaphson(Mu_tdb_dev[thermo_phase[0]],
                                                  evalFunc(Mu_tdb_dev[thermo_phase[1]], phaseComp[1][idx], temperature)/molarVolume,
                                                  dmudc_tdb_dev[thermo_phase[0]],
                                                  temperature, molarVolume,
                                                  phaseComp[0][idx]);
            }
        }
    } // if (i < MESH_X ...)
    __syncthreads();
}

__global__ void __calcPhaseCompBinary_01_03__(double **phi, double **comp,
                                              double **phaseComp,
                                              double *F0_A, double *F0_B, double *F0_C,
                                              int NUMPHASES, int NUMCOMPONENTS,
                                              int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        double ctemp  = comp[0][idx];
        double etemp  = phi[1][idx];

        double interp_phi   = etemp * etemp * etemp * (6.0 * etemp * etemp - 15.0 * etemp + 10.0);

        phaseComp[0][idx] = (2.0*F0_A[1]*ctemp + (F0_B[1] - F0_B[0])*interp_phi)/(2.0*(F0_A[1] + interp_phi*(F0_A[0] - F0_A[1])));
        phaseComp[1][idx] = (2.0*F0_A[0]*phaseComp[0][idx] + F0_B[0] - F0_B[1])/(2.0*F0_A[1]);
    }
}

void calcPhaseComp(double **phi, double **comp,
                   double **phaseComp,
                   domainInfo* simDomain, controls* simControls,
                   simParameters* simParams, subdomainInfo* subdomain,
                   dim3 gridSize, dim3 blockSize)
{
    if (simControls->multiphase == 1 || simDomain->numPhases > 2 || simDomain->numComponents > 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3 || simControls->FUNCTION_F == 4)
        {
            __calcPhaseComp_01_03__<<<gridSize, blockSize>>>(phi, comp,
                                                             phaseComp,
                                                             simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                             simDomain->numPhases, simDomain->numComponents,
                                                             subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);
        }
        else if (simControls->FUNCTION_F == 2)
            printf("MPMC Exact not implemented\n");
    }
    else if (simDomain->numPhases == 2 && simDomain->numComponents == 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3  || simControls->FUNCTION_F == 4)
        {
            __calcPhaseCompBinary_01_03__<<<gridSize, blockSize>>>(phi, comp,
                                                                   phaseComp,
                                                                   simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                                   simDomain->numPhases, simDomain->numComponents,
                                                                   subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);
        }
        else if (simControls->FUNCTION_F == 2)
        {
            if (simControls->count == 0)
                initPhaseComp_02<<<gridSize, blockSize>>>(phi, comp,
                                                          phaseComp, simParams->cguess_dev,
                                                          simParams->T, simParams->molarVolume, simDomain->thermo_phase_dev,
                                                          subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);

            __calcPhaseCompBinary_02__<<<gridSize, blockSize>>>(phi, comp,
                                                                phaseComp,
                                                                simParams->T, simParams->molarVolume, simDomain->thermo_phase_dev,
                                                                subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);
        }
    }
    cudaError_t error = cudaGetLastError();;

    if(error != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }
}
