#include "updateComposition.cuh"

__device__
int LUPDecompose3(double A[(MAX_NUM_COMP-1)*(MAX_NUM_COMP-1)], int N, double Tol, int *P)
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
            if ((absA = fabs(A[k*N+i])) > maxA)
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
                    ptr = A[i*N+j];
                    A[i*N+j] = A[imax*N+j];
                    A[imax*N+j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j*N+i] /= A[i*N+i];

                for (k = i + 1; k < N; k++)
                    A[j*N+k] -= A[j*N+i] * A[i*N+k];
            }
    }

    return 1;  //decomposition done
}

__device__
void LUPInvert3(double A[(MAX_NUM_COMP-1)*(MAX_NUM_COMP-1)], int *P, int N, double IA[][MAX_NUM_COMP-1])
{
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];

            IA[i][j] /= A[i*N+i];
        }
    }
}

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
double calcInterp5thDiff2(double **phi, int a, int b, int idx, int NUMPHASES)
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
                           double *mobility,
                           int NUMPHASES, int NUMCOMPONENTS,
                           int sizeX, int sizeY, int sizeZ,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp, xm, yp, ym, zp, zm;

    double mu[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp = (j + k*sizeY)*sizeX + i+1;
        xm = (j + k*sizeY)*sizeX + i-1;

        if (i == 0)
        {
            xm = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 1)
        {
            xp = (j + k*sizeY)*sizeX;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp = (j+1 + k*sizeY)*sizeX + i;
            ym = (j-1 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp = (k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp = (j + (k+1)*sizeY)*sizeX + i;
            zm = (j + (k-1)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp = j*sizeX + i;
            }
        }

        for (int component = 0; component < NUMCOMPONENTS-1; component++)
        {
            compNew[component][idx] = comp[component][idx];

            J_xp = 0.0;
            J_xm = 0.0;
            J_yp = 0.0;
            J_ym = 0.0;
            J_zp = 0.0;
            J_zm = 0.0;

            for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
            {

                for (int iter = 0; iter < 7; iter++)
                    effMobility[iter] = 0.0;

                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    effMobility[0] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, idx, NUMPHASES);

                    effMobility[1] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xp, NUMPHASES);
                    effMobility[2] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xm, NUMPHASES);

                    if (sizeY > 1)
                    {
                        effMobility[3] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, yp, NUMPHASES);
                        effMobility[4] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, ym, NUMPHASES);
                    }

                    if (sizeZ > 1)
                    {
                        effMobility[5] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zp, NUMPHASES);
                        effMobility[6] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zm, NUMPHASES);
                    }
                }

                mu[0] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS);

                mu[1] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, xp, NUMPHASES, NUMCOMPONENTS);
                mu[2] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, xm, NUMPHASES, NUMCOMPONENTS);

                if (sizeY > 1)
                {
                    mu[3] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, yp, NUMPHASES, NUMCOMPONENTS);
                    mu[4] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, ym, NUMPHASES, NUMCOMPONENTS);
                }

                if (sizeZ > 1)
                {
                    mu[5]  = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, zp, NUMPHASES, NUMCOMPONENTS);
                    mu[6] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, zm, NUMPHASES, NUMCOMPONENTS);
                }

                J_xp += ((effMobility[1] + effMobility[0])/2.0)*(mu[1] - mu[0])/DELTA_X;
                J_xm += ((effMobility[0] + effMobility[2])/2.0)*(mu[0] - mu[2])/DELTA_X;

                if (sizeY > 1)
                {
                    J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[3] - mu[0])/DELTA_Y;
                    J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0] - mu[4])/DELTA_Y;
                }

                if (sizeZ > 1)
                {
                    J_zp += ((effMobility[5] + effMobility[0])/2.0)*(mu[5] - mu[0])/DELTA_Z;
                    J_zm += ((effMobility[0] + effMobility[6])/2.0)*(mu[0] - mu[6])/DELTA_Z;
                }
            }

            compNew[component][idx] +=  DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);
        }
    }
    __syncthreads();
}

__global__
void __updateComposition_02__(double **phi,
                              double **comp, double **compNew, double **mu,
                              double **phaseComp, int *thermo_phase,
                              double *diffusivity, double temperature, double molarVolume,
                              int NUMPHASES, int NUMCOMPONENTS,
                              int sizeX, int sizeY, int sizeZ,
                              double DELTA_X, double DELTA_Y, double DELTA_Z,
                              double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp, xm, yp, ym, zp, zm;

    double muLocal[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;
    double tol = 1e-6;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp = (j + k*sizeY)*sizeX + i+1;
        xm = (j + k*sizeY)*sizeX + i-1;

        if (i == 0)
        {
            xm = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 1)
        {
            xp = (j + k*sizeY)*sizeX;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp = (j+1 + k*sizeY)*sizeX + i;
            ym = (j-1 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp = (k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp = (j + (k+1)*sizeY)*sizeX + i;
            zm = (j + (k-1)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp = j*sizeX + i;
            }
        }

        double dmudc[(MAX_NUM_COMP-1)*(MAX_NUM_COMP-1)];
        double y[MAX_NUM_COMP];
        double dmudcInv[MAX_NUM_COMP-1][MAX_NUM_COMP-1];
        int P[MAX_NUM_COMP];
        double mobility[MAX_NUM_COMP-1][MAX_NUM_COMP-1];

        for (int component = 0; component < NUMCOMPONENTS-1; component++)
        {
            compNew[component][idx] = comp[component][idx];

            J_xp = 0.0;
            J_xm = 0.0;
            J_yp = 0.0;
            J_ym = 0.0;
            J_zp = 0.0;
            J_zm = 0.0;

            for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
            {

                for (int iter = 0; iter < 7; iter++)
                    effMobility[iter] = 0.0;

                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    double tmp0 = 0.0;

                    for (int is = 0; is < NUMCOMPONENTS-1; is++)
                    {
                        y[is] = phaseComp[is*NUMPHASES + phase][idx];
                        tmp0  += y[is];
                    }

                    y[NUMCOMPONENTS-1] = 1.0 - tmp0;

                    (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, dmudc);


                    LUPDecompose3(dmudc, NUMCOMPONENTS-1, tol, P);
                    LUPInvert3(dmudc, P, NUMCOMPONENTS-1, dmudcInv);

                    for (int iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
                    {
                        for (int iter2 = 0; iter2 < NUMCOMPONENTS-1; iter2++)
                        {
                            mobility[iter1][iter2] = 0.0;

                            for (int iter3 = 0; iter3 < NUMCOMPONENTS-1; iter3++)
                            {
                                mobility[iter1][iter2] += diffusivity[(iter3 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + iter1]*dmudcInv[iter3][iter2];
                            }
                        }
                    }

                    effMobility[0] += mobility[component][component2]*calcInterp5th2(phi, phase, idx, NUMPHASES);

                    effMobility[1] += mobility[component][component2]*calcInterp5th2(phi, phase, xp, NUMPHASES);
                    effMobility[2] += mobility[component][component2]*calcInterp5th2(phi, phase, xm, NUMPHASES);

                    if (sizeY > 1)
                    {
                        effMobility[3] += mobility[component][component2]*calcInterp5th2(phi, phase, yp, NUMPHASES);
                        effMobility[4] += mobility[component][component2]*calcInterp5th2(phi, phase, ym, NUMPHASES);
                    }

                    if (sizeZ > 1)
                    {
                        effMobility[5] += mobility[component][component2]*calcInterp5th2(phi, phase, zp, NUMPHASES);
                        effMobility[6] += mobility[component][component2]*calcInterp5th2(phi, phase, zm, NUMPHASES);
                    }
                }

                muLocal[0] = mu[component2][idx];

                muLocal[1] = mu[component2][xp];
                muLocal[2] = mu[component2][xm];

                if (sizeY > 1)
                {
                    muLocal[3] = mu[component2][yp];
                    muLocal[4] = mu[component2][ym];
                }
                if (sizeZ > 1)
                {
                    muLocal[5] = mu[component2][zp];
                    muLocal[6] = mu[component2][zm];
                }

                J_xp += ((effMobility[1] + effMobility[0])/2.0)*(muLocal[1] - muLocal[0])/DELTA_X;
                J_xm += ((effMobility[0] + effMobility[2])/2.0)*(muLocal[0] - muLocal[2])/DELTA_X;

                if (sizeY > 1)
                {
                    J_yp += ((effMobility[3] + effMobility[0])/2.0)*(muLocal[3] - muLocal[0])/DELTA_Y;
                    J_ym += ((effMobility[0] + effMobility[4])/2.0)*(muLocal[0] - muLocal[4])/DELTA_Y;
                }

                if (sizeZ > 1)
                {
                    J_zp += ((effMobility[5] + effMobility[0])/2.0)*(muLocal[5] - muLocal[0])/DELTA_Z;
                    J_zm += ((effMobility[0] + effMobility[6])/2.0)*(muLocal[0] - muLocal[6])/DELTA_Z;
                }
            }

            compNew[component][idx] +=  DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);
        }
    }
    __syncthreads();
}

__global__
void __updateMu_02__(double **phi, double **comp,
                     double **phiNew, double **compNew,
                     double **phaseComp, double **mu,
                     int *thermo_phase, double temperature, double molarVolume,
                     int NUMPHASES, int NUMCOMPONENTS,
                     int sizeX, int sizeY, int sizeZ,
                     double DELTA_X, double DELTA_Y, double DELTA_Z,
                     double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        double RHS = 0.0, sum = 0.0;
        double tol = 1e-6;

        int bulkphase = 0, interface = 1;

        for (int is = 0; is < NUMPHASES; is++)
        {
            if (phi[is][idx] > 0.999)
            {
                bulkphase = is;
                interface = 0;
                break;
            }
        }

        if (interface)
        {
            double dmudc[(MAX_NUM_COMP-1)*(MAX_NUM_COMP-1)];
            double dcdmu[(MAX_NUM_COMP-1)*(MAX_NUM_COMP-1)];
            double y[MAX_NUM_COMP];
            double Inv[MAX_NUM_COMP-1][MAX_NUM_COMP-1];
            int P[MAX_NUM_COMP];

            for (int component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    dcdmu[component*(NUMCOMPONENTS-1) + component2] = 0.0;
                }
            }

            for (int component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    sum = 0.0;

                    for (int is = 0; is < NUMCOMPONENTS-1; is++)
                    {
                        y[is] = phaseComp[is*NUMPHASES + phase][idx];
                        sum += y[is];
                    }

                    y[NUMCOMPONENTS-1] = 1.0 - sum;

                    sum = 0.0;

                    (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, dmudc);

                    //                 for (int iter1 = 0; iter1 < NUMCOMPONENTS-1; iter1++)
                    //                     dmudc[iter1] /= (1.602*1e8*molarVolume);

                    LUPDecompose3(dmudc, NUMCOMPONENTS-1, tol, P);
                    LUPInvert3(dmudc, P, NUMCOMPONENTS-1, Inv);

                    for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                        dcdmu[component*(NUMCOMPONENTS-1) + component2] += calcInterp5th2(phi, phase, idx, NUMPHASES)*(Inv[component][component2]);
                }

                LUPDecompose3(dcdmu, NUMCOMPONENTS-1, tol, P);
                LUPInvert3(dcdmu, P, NUMCOMPONENTS-1, Inv);
            }

            for (int component = 0; component < NUMCOMPONENTS-1; component++)
            {
                RHS = (compNew[component][idx] - comp[component][idx]);

                for (int phase = 0; phase < NUMPHASES; phase++)
                {
                    for (int phase2 = 0; phase2 < NUMPHASES; phase2++)
                    {
                        RHS += phaseComp[phase + component*NUMPHASES][idx]*calcInterp5thDiff2(phi, phase, phase2, idx, NUMPHASES)*(phiNew[phase2][idx] - phi[phase2][idx]);
                    }
                }

                for (int component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                    mu[component][idx] += RHS*dcdmu[component*(NUMCOMPONENTS-1) + component2];
            }
        }
        else
        {
            double y[MAX_NUM_COMP];
            double mu1[MAX_NUM_COMP];

            sum = 0.0;

            for (int is = 0; is < NUMCOMPONENTS-1; is++)
            {
                y[is] = comp[is][idx];
                sum += y[is];
            }

            y[NUMCOMPONENTS-1] = 1.0 - sum;

            (*Mu_tdb_dev[thermo_phase[bulkphase]])(temperature, y, mu1);

            for (int is = 0; is < NUMCOMPONENTS-1; is++)
                mu[is][idx] = mu1[is];
        }
    }
}

__global__
void __updateCompositionBinary__(double **phi,
                                 double **comp, double **compNew,
                                 double **phaseComp,
                                 double *F0_A, double *F0_B,
                                 double *mobility,
                                 int NUMPHASES, int NUMCOMPONENTS,
                                 int sizeX, int sizeY, int sizeZ,
                                 double DELTA_X, double DELTA_Y, double DELTA_Z,
                                 double DELTA_t)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int idx = (j + k*sizeY)*sizeX + i;

    int xp, xm, yp, ym, zp, zm;

    double mu[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp = (j + k*sizeY)*sizeX + i+1;
        xm = (j + k*sizeY)*sizeX + i-1;

        if (i == 0)
        {
            xm = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 1)
        {
            xp = (j + k*sizeY)*sizeX;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp = (j+1 + k*sizeY)*sizeX + i;
            ym = (j-1 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp = (k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp = (j + (k+1)*sizeY)*sizeX + i;
            zm = (j + (k-1)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp = j*sizeX + i;
            }
        }

        int component = 0;
        int component2 = 0;

        J_xp = 0.0;
        J_xm = 0.0;
        J_yp = 0.0;
        J_ym = 0.0;
        J_zp = 0.0;
        J_zm = 0.0;

        for (int iter = 0; iter < 7; iter++)
            effMobility[iter] = 0.0;

        for (int phase = 0; phase < NUMPHASES; phase++)
        {
            effMobility[0] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, idx, NUMPHASES);

            effMobility[1] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xp, NUMPHASES);
            effMobility[2] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xm, NUMPHASES);

            if (sizeY > 1)
            {
                effMobility[3] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, yp, NUMPHASES);
                effMobility[4] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, ym, NUMPHASES);
            }

            if (sizeZ > 1)
            {
                effMobility[5] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zp, NUMPHASES);
                effMobility[6] += mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zm, NUMPHASES);
            }
        }

        mu[0] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS);

        mu[1] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, xp, NUMPHASES, NUMCOMPONENTS);
        mu[2] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, xm, NUMPHASES, NUMCOMPONENTS);

        if (sizeY > 1)
        {
            mu[3] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, yp, NUMPHASES, NUMCOMPONENTS);
            mu[4] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, ym, NUMPHASES, NUMCOMPONENTS);
        }

        if (sizeZ > 1)
        {
            mu[5]  = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, zp, NUMPHASES, NUMCOMPONENTS);
            mu[6] = calcDiffusionPotential2(phaseComp, 0, component2, F0_A, F0_B, zm, NUMPHASES, NUMCOMPONENTS);
        }

        J_xp += ((effMobility[1] + effMobility[0])/2.0)*(mu[1] - mu[0])/DELTA_X;
        J_xm += ((effMobility[0] + effMobility[2])/2.0)*(mu[0] - mu[2])/DELTA_X;

        if (sizeY > 1)
        {
            J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[3] - mu[0])/DELTA_Y;
            J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0] - mu[4])/DELTA_Y;
        }

        if (sizeZ > 1)
        {
            J_zp += ((effMobility[5] + effMobility[0])/2.0)*(mu[5] - mu[0])/DELTA_Z;
            J_zm += ((effMobility[0] + effMobility[6])/2.0)*(mu[0] - mu[6])/DELTA_Z;
        }
    }

    compNew[0][idx] = comp[0][idx] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z);

    __syncthreads();
}

__global__
void __updateCompositionBinary_02__(double **phi,
                                    double **comp, double **compNew,
                                    double **phaseComp, double **mu,
                                    int *thermo_phase, double molarVolume, double temperature,
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

    int xp, xm, yp, ym, zp, zm;

    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    if (i < sizeX && j < sizeY && k < sizeZ)
    {
        // x-direction
        xp = (j + k*sizeY)*sizeX + i+1;
        xm = (j + k*sizeY)*sizeX + i-1;

        if (i == 0)
        {
            xm = (j + k*sizeY)*sizeX + sizeX-1;
        }
        else if (i == sizeX - 1)
        {
            xp = (j + k*sizeY)*sizeX;
        }

        // y-direction
        if (sizeY > 1)
        {
            yp = (j+1 + k*sizeY)*sizeX + i;
            ym = (j-1 + k*sizeY)*sizeX + i;

            if (j == 0)
            {
                ym = (sizeY-1 + k*sizeY)*sizeX + i;
            }
            else if (j == sizeY - 1)
            {
                yp = (k*sizeY)*sizeX + i;
            }
        }

        // z-direction
        if (sizeZ > 1)
        {
            zp = (j + (k+1)*sizeY)*sizeX + i;
            zm = (j + (k-1)*sizeY)*sizeX + i;

            if (k == 0)
            {
                zm = (j + (sizeZ-1)*sizeY)*sizeX + i;
            }
            else if (k == sizeZ - 1)
            {
                zp = j*sizeX + i;
            }
        }

        int component = 0;
        int component2 = 0;

        J_xp = 0.0;
        J_xm = 0.0;
        J_yp = 0.0;
        J_ym = 0.0;
        J_zp = 0.0;
        J_zm = 0.0;

        for (int iter = 0; iter < 7; iter++)
            effMobility[iter] = 0.0;

        for (int phase = 1; phase < NUMPHASES; phase++)
        {
/*
            (*dmudc_tdb_dev[thermo_phase[phase]])(temperature, y, dmudc);

            effMobility[0] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, idx, NUMPHASES)/dmudc[phase];

            effMobility[1] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xp, NUMPHASES)/dmudc[phase];
            effMobility[2] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, xm, NUMPHASES)/dmudc[phase];*/

            effMobility[0] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, idx, NUMPHASES)*(phaseComp[0][idx] - phaseComp[1][idx]);

            effMobility[1] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, xp, NUMPHASES)*(phaseComp[0][xp] - phaseComp[1][xp]);
            effMobility[2] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, xm, NUMPHASES)*(phaseComp[0][xm] - phaseComp[1][xm]);

            if (sizeY > 1)
            {
//                 effMobility[3] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, yp, NUMPHASES)/dmudc[phase];
//                 effMobility[4] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, ym, NUMPHASES)/dmudc[phase];

                effMobility[3] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, yp, NUMPHASES)*(phaseComp[0][yp] - phaseComp[1][yp]);
                effMobility[4] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, ym, NUMPHASES)*(phaseComp[0][ym] - phaseComp[1][ym]);
            }

            if (sizeZ > 1)
            {
//                 effMobility[5] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zp, NUMPHASES)/dmudc[phase];
//                 effMobility[6] += diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5th2(phi, phase, zm, NUMPHASES)/dmudc[phase];

                effMobility[5] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, zp, NUMPHASES)*(phaseComp[0][zp] - phaseComp[1][zp]);
                effMobility[6] = diffusivity[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]*calcInterp5thDiff2(phi, phase, phase, zm, NUMPHASES)*(phaseComp[0][zm] - phaseComp[1][zm]);
            }
        }

//         (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//         mu[0] = dmudc[0];
//
//         y[0] = phaseComp[0][xp];
//         y[1] = phaseComp[1][xp];
//         (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//         mu[1] = dmudc[0];
//
//         y[0] = phaseComp[0][xm];
//         y[1] = phaseComp[1][xm];
//         (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//         mu[2] = dmudc[0];
//
//         if (sizeY > 1)
//         {
//             y[0] = phaseComp[0][yp];
//             y[1] = phaseComp[1][ym];
//             (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//             mu[3] = dmudc[0];
//
//             y[0] = phaseComp[0][ym];
//             y[1] = phaseComp[1][ym];
//             (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//             mu[4] = dmudc[0];
//
//         }
//
//         if (sizeZ > 1)
//         {
//             y[0] = phaseComp[0][zp];
//             y[1] = phaseComp[1][zm];
//             (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//             mu[5] = dmudc[0];
//
//             y[0] = phaseComp[0][zm];
//             y[1] = phaseComp[1][zm];
//             (*Mu_tdb_dev[thermo_phase[0]])(temperature, y, dmudc);
//             mu[6] = dmudc[0];
//         }

        J_xp += ((effMobility[1] + effMobility[0])/2.0)*(mu[0][xp] - mu[0][idx])/DELTA_X;
        J_xm += ((effMobility[0] + effMobility[2])/2.0)*(mu[0][idx] - mu[0][xm])/DELTA_X;

        if (sizeY > 1)
        {
            J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[0][yp] - mu[0][idx])/DELTA_Y;
            J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0][idx] - mu[0][ym])/DELTA_Y;
        }

        if (sizeZ > 1)
        {
            J_zp += ((effMobility[5] + effMobility[0])/2.0)*(mu[0][zp] - mu[0][idx])/DELTA_Z;
            J_zm += ((effMobility[0] + effMobility[6])/2.0)*(mu[0][idx] - mu[0][zm])/DELTA_Z;
        }
    }

    compNew[0][idx] = comp[0][idx] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z)/molarVolume;

    __syncthreads();
}

void updateComposition(double **phi, double **comp, double **phiNew, double **compNew,
                       double **phaseComp, double **mu,
                       domainInfo* simDomain, controls* simControls,
                       simParameters* simParams, subdomainInfo* subdomain,
                       dim3 gridSize, dim3 blockSize)
{
    if (simControls->multiphase == 1 || simDomain->numPhases > 2 || simDomain->numComponents > 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3 || simControls->FUNCTION_F == 4)
        {
            __updateComposition__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                           phaseComp,
                                                           simParams->F0_A_dev, simParams->F0_B_dev,
                                                           simParams->mobility_dev,
                                                           simDomain->numPhases, simDomain->numComponents,
                                                           subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                           simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                           simControls->DELTA_t);
        }
        else if (simControls->FUNCTION_F == 2)
        {
            __updateComposition_02__<<<gridSize, blockSize>>>(phi, comp,
                                                              compNew, mu,
                                                              phaseComp, simDomain->thermo_phase_dev,
                                                              simParams->diffusivity_dev, simParams->T, simParams->molarVolume,
                                                              simDomain->numPhases, simDomain->numComponents,
                                                              subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                              simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                              simControls->DELTA_t);


            __updateMu_02__<<<gridSize, blockSize>>>(phi, comp,
                                                     phiNew, compNew,
                                                     phaseComp, mu,
                                                     simDomain->thermo_phase_dev, simParams->T, simParams->molarVolume,
                                                     simDomain->numPhases, simDomain->numComponents,
                                                     subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                     simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                     simControls->DELTA_t);
        }
    }
    else if (simDomain->numPhases == 2 && simDomain->numComponents == 2)
    {
        if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3 || simControls->FUNCTION_F == 4)
        {
            __updateCompositionBinary__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                                 phaseComp,
                                                                 simParams->F0_A_dev, simParams->F0_B_dev,
                                                                 simParams->mobility_dev,
                                                                 simDomain->numPhases, simDomain->numComponents,
                                                                 subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                                 simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                                 simControls->DELTA_t);
        }
        else if (simControls->FUNCTION_F == 2)
        {
            __updateCompositionBinary_02__<<<gridSize, blockSize>>>(phi, comp, compNew,
                                                                    phaseComp, mu,
                                                                    simDomain->thermo_phase_dev, simParams->molarVolume, simParams->T,
                                                                    simParams->diffusivity_dev,
                                                                    simDomain->numPhases, simDomain->numComponents,
                                                                    subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                                    simDomain->DELTA_X, simDomain->DELTA_Y, simDomain->DELTA_Z,
                                                                    simControls->DELTA_t);
        }

    }
}
