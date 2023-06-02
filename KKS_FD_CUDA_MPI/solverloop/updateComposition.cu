//#include <cub/cub.cuh>
#include "updateComposition.cuh"

__global__
void __updateComposition__(double **phi, double **phiNew,
                           double **comp, double **compNew,
                           double **phaseComp,
                           double *F0_A, double *F0_B,
                           double *mobility, double *diffusivity, double *kappaPhi, double *theta_ij,
                           long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                           long sizeX, long sizeY, long sizeZ,
                           long xStep, long yStep, long padding, int antiTrapping,
                           double DELTA_X, double DELTA_Y, double DELTA_Z,
                           double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long index[3][3][3] = {-1}, maxPos = 5;

    double mu[7];
    double effMobility[7];
    double mobilityLocal;
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    // Antitrapping variables
    double alpha[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};
    double modgradphi[MAX_NUM_PHASES][7] = {0.0};
    double dphidt[MAX_NUM_PHASES][7] = {0.0};

    double gradx_phi[MAX_NUM_PHASES][5] = {0.0};
    double grady_phi[MAX_NUM_PHASES][5] = {0.0};
    double gradz_phi[MAX_NUM_PHASES][5] = {0.0};

    double phix[MAX_NUM_PHASES][7] = {0.0}, phiy[MAX_NUM_PHASES][7] = {0.0}, phiz[MAX_NUM_PHASES][7] = {0.0};

    double alphidot[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};
    double scalprodct[MAX_NUM_PHASES-1][7] = {0.0};

    double jatr[MAX_NUM_COMP-1] = {0.0};
    double jatc[MAX_NUM_PHASES-1][MAX_NUM_COMP-1] = {0.0};
    double jat[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};

    long component, component2, phase;

    int interface = 1;
    double divphi1 = 0.0, divphi2 = 0.0;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        /*
         * index is a 3D matrix that stores location indices for the ...
         * 26 gridpoints around the gridpoint at which computation is being carried out.
         *
         * Eg:
         * index[0][1][2] is (i-1, j, k+1)
         *
         * for 2D simulations, the third index will be set to 1 (index[x][y][1])
         *
         */
        for (long x = 0; x < 3; x++)
        {
            for (long y = 0; y < 3; y++)
            {
                for (long z = 0; z < 3; z++)
                {
                    index[x][y][z] = (k+z-1) + (j+y-1)*yStep + (i+x-1)*xStep;
                }
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            // Determine if interfacial point

            divphi1 = (phi[phase][index[2][1][1]] - 2.0*phi[phase][index[1][1][1]] + phi[phase][index[0][1][1]])/(DELTA_X*DELTA_X);
            if (DIMENSION >= 2)
                divphi1 += (phi[phase][index[1][2][1]] - 2.0*phi[phase][index[1][1][1]] + phi[phase][index[1][0][1]])/(DELTA_Y*DELTA_Y);
            if (DIMENSION == 3)
                divphi1 += (phi[phase][index[1][1][2]] - 2.0*phi[phase][index[1][1][1]] + phi[phase][index[1][1][0]])/(DELTA_Z*DELTA_Z);

            divphi2 = (phiNew[phase][index[2][1][1]] - 2.0*phiNew[phase][index[1][1][1]] + phiNew[phase][index[0][1][1]])/(DELTA_X*DELTA_X);
            if (DIMENSION >= 2)
                divphi2 += (phiNew[phase][index[1][2][1]] - 2.0*phiNew[phase][index[1][1][1]] + phiNew[phase][index[1][0][1]])/(DELTA_Y*DELTA_Y);
            if (DIMENSION == 3)
                divphi2 += (phiNew[phase][index[1][1][2]] - 2.0*phiNew[phase][index[1][1][1]] + phiNew[phase][index[1][1][0]])/(DELTA_Z*DELTA_Z);

            if (fabs(divphi1) < 1e-3/DELTA_X && fabs(divphi2) < 1e-3/DELTA_X)
                interface = 0;
        }

        if (DIMENSION == 3)
            maxPos = 7;
        else if (DIMENSION == 2)
            maxPos = 5;
        else
            maxPos = 3;

        if (interface && antiTrapping)
        {
            //Antitrapping Calculations
            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                double A1 = sqrt(2.0*kappaPhi[phase*NUMPHASES + NUMPHASES-1]/theta_ij[phase*NUMPHASES + NUMPHASES-1]);

                for (component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    alpha[phase][component][0] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][1]] - phaseComp[component*NUMPHASES + phase][index[1][1][1]]);

                    alpha[phase][component][1] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[2][1][1]] - phaseComp[component*NUMPHASES + phase][index[1][1][1]]);
                    alpha[phase][component][2] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[0][1][1]] - phaseComp[component*NUMPHASES + phase][index[0][1][1]]);

                    if (DIMENSION == 2)
                    {
                        alpha[phase][component][3] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][2][1]] - phaseComp[component*NUMPHASES + phase][index[1][2][1]]);
                        alpha[phase][component][4] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][0][1]] - phaseComp[component*NUMPHASES + phase][index[1][0][1]]);
                    }
                    else if (DIMENSION == 3)
                    {
                        alpha[phase][component][3] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][2][1]] - phaseComp[component*NUMPHASES + phase][index[1][2][1]]);
                        alpha[phase][component][4] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][0][1]] - phaseComp[component*NUMPHASES + phase][index[1][0][1]]);

                        alpha[phase][component][5] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][2]] - phaseComp[component*NUMPHASES + phase][index[1][1][2]]);
                        alpha[phase][component][6] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][0]] - phaseComp[component*NUMPHASES + phase][index[1][1][0]]);
                    }
                }

                dphidt[phase][0] = (phiNew[phase][index[1][1][1]] - phi[phase][index[1][1][1]])/DELTA_t;

                dphidt[phase][1] = (phiNew[phase][index[2][1][1]] - phi[phase][index[2][1][1]])/DELTA_t;
                dphidt[phase][2] = (phiNew[phase][index[0][1][1]] - phi[phase][index[0][1][1]])/DELTA_t;

                if (DIMENSION == 2)
                {
                    dphidt[phase][3] = (phiNew[phase][index[1][2][1]] - phi[phase][index[1][2][1]])/DELTA_t;
                    dphidt[phase][4] = (phiNew[phase][index[1][0][1]] - phi[phase][index[1][0][1]])/DELTA_t;
                }
                else if (DIMENSION == 3)
                {
                    dphidt[phase][3] = (phiNew[phase][index[1][2][1]] - phi[phase][index[1][2][1]])/DELTA_t;
                    dphidt[phase][4] = (phiNew[phase][index[1][0][1]] - phi[phase][index[1][0][1]])/DELTA_t;

                    dphidt[phase][5] = (phiNew[phase][index[1][1][2]] - phi[phase][index[1][1][2]])/DELTA_t;
                    dphidt[phase][6] = (phiNew[phase][index[1][1][0]] - phi[phase][index[1][1][0]])/DELTA_t;
                }
            }

            for (phase = 0; phase < NUMPHASES; phase++)
            {
                gradx_phi[phase][0] = (phi[phase][index[2][1][1]] - phi[phase][index[0][1][1]])/(2.0*DELTA_X);

                if (DIMENSION == 2)
                {
                    gradx_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[0][2][1]])/(2.0*DELTA_X);
                    gradx_phi[phase][2] = (phi[phase][index[2][0][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_X);

                    grady_phi[phase][0] = (phi[phase][index[1][2][1]] - phi[phase][index[1][0][1]])/(2.0*DELTA_Y);
                    grady_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[2][0][1]])/(2.0*DELTA_Y);
                    grady_phi[phase][2] = (phi[phase][index[0][2][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_Y);
                }
                else if (DIMENSION == 3)
                {
                    gradx_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[0][2][1]])/(2.0*DELTA_X);
                    gradx_phi[phase][2] = (phi[phase][index[2][0][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_X);
                    gradx_phi[phase][3] = (phi[phase][index[2][1][2]] - phi[phase][index[0][1][2]])/(2.0*DELTA_X);
                    gradx_phi[phase][4] = (phi[phase][index[2][1][0]] - phi[phase][index[0][1][0]])/(2.0*DELTA_X);

                    grady_phi[phase][0] = (phi[phase][index[1][2][1]] - phi[phase][index[1][0][1]])/(2.0*DELTA_Y);
                    grady_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[2][0][1]])/(2.0*DELTA_Y);
                    grady_phi[phase][2] = (phi[phase][index[0][2][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_Y);
                    grady_phi[phase][3] = (phi[phase][index[1][2][2]] - phi[phase][index[1][0][2]])/(2.0*DELTA_Y);
                    grady_phi[phase][4] = (phi[phase][index[1][2][0]] - phi[phase][index[1][0][0]])/(2.0*DELTA_Y);

                    gradz_phi[phase][0] = (phi[phase][index[1][1][2]] - phi[phase][index[1][1][0]])/(2.0*DELTA_Z);
                    gradz_phi[phase][1] = (phi[phase][index[2][1][2]] - phi[phase][index[2][1][0]])/(2.0*DELTA_Z);
                    gradz_phi[phase][2] = (phi[phase][index[0][1][2]] - phi[phase][index[0][1][0]])/(2.0*DELTA_Z);
                    gradz_phi[phase][3] = (phi[phase][index[1][2][2]] - phi[phase][index[1][2][0]])/(2.0*DELTA_Z);
                    gradz_phi[phase][4] = (phi[phase][index[1][0][2]] - phi[phase][index[1][0][0]])/(2.0*DELTA_Z);
                }
            }

            /*
             * 0 -> central node
             * 1 -> (i+0.5, j, k)
             * 2 -> (i-0.5, j, k)
             * 3 -> (i, j+0.5, k)
             * 4 -> (i, j-0.5, k)
             * 5 -> (i, j, k+0.5)
             * 6 -> (i, j, k-0.5)
             */

            for (phase = 0; phase < NUMPHASES; phase++)
            {
                phix[phase][0] = gradx_phi[phase][0];
                phix[phase][1] = (phi[phase][index[2][1][1]] - phi[phase][index[1][1][1]])/(DELTA_X);
                phix[phase][2] = (phi[phase][index[1][1][1]] - phi[phase][index[0][1][1]])/(DELTA_X);

                if (DIMENSION == 2)
                {
                    phix[phase][3] = (gradx_phi[phase][0] + gradx_phi[phase][1])/2.0;
                    phix[phase][4] = (gradx_phi[phase][0] + gradx_phi[phase][2])/2.0;
                }
                else if (DIMENSION == 3)
                {
                    phix[phase][3] = (gradx_phi[phase][0] + gradx_phi[phase][1])/2.0;
                    phix[phase][4] = (gradx_phi[phase][0] + gradx_phi[phase][2])/2.0;
                    phix[phase][5] = (gradx_phi[phase][0] + gradx_phi[phase][3])/2.0;
                    phix[phase][6] = (gradx_phi[phase][0] + gradx_phi[phase][4])/2.0;
                }

                if (DIMENSION >= 2)
                {
                    phiy[phase][0] = grady_phi[phase][0];
                    phiy[phase][1] = (grady_phi[phase][0] + grady_phi[phase][1])/2.0;
                    phiy[phase][2] = (grady_phi[phase][0] + grady_phi[phase][2])/2.0;
                    phiy[phase][3] = (phi[phase][index[1][2][1]] - phi[phase][index[1][1][1]])/(DELTA_Y);
                    phiy[phase][4] = (phi[phase][index[1][1][1]] - phi[phase][index[1][0][1]])/(DELTA_Y);

                    if (DIMENSION ==  3)
                    {
                        phiy[phase][5] = (grady_phi[phase][0] + grady_phi[phase][3])/2.0;
                        phiy[phase][6] = (grady_phi[phase][0] + grady_phi[phase][4])/2.0;

                        phiz[phase][0] = gradz_phi[phase][0];
                        phiz[phase][1] = (gradz_phi[phase][0] + gradz_phi[phase][1])/2.0;
                        phiz[phase][2] = (gradz_phi[phase][0] + gradz_phi[phase][2])/2.0;
                        phiz[phase][3] = (gradz_phi[phase][0] + gradz_phi[phase][3])/2.0;
                        phiz[phase][4] = (gradz_phi[phase][0] + gradz_phi[phase][4])/2.0;
                        phiz[phase][5] = (phi[phase][index[1][1][2]] - phi[phase][index[1][1][1]])/(DELTA_Z);
                        phiz[phase][6] = (phi[phase][index[1][1][1]] - phi[phase][index[1][1][0]])/(DELTA_Z);
                    }
                }
            }

            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                for (component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    alphidot[phase][component][1] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][1]*dphidt[phase][1]))/2.0;
                    alphidot[phase][component][2] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][2]*dphidt[phase][2]))/2.0;

                    if (DIMENSION == 2)
                    {
                        alphidot[phase][component][3] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][3]*dphidt[phase][3]))/2.0;
                        alphidot[phase][component][4] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][4]*dphidt[phase][4]))/2.0;
                    }
                    else if (DIMENSION == 3)
                    {
                        alphidot[phase][component][3] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][3]*dphidt[phase][3]))/2.0;
                        alphidot[phase][component][4] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][4]*dphidt[phase][4]))/2.0;

                        alphidot[phase][component][5] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][5]*dphidt[phase][5]))/2.0;
                        alphidot[phase][component][6] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][6]*dphidt[phase][6]))/2.0;
                    }
                }
            }

            for (phase = 0; phase < NUMPHASES; phase++)
            {
                for (long iter = 0; iter < maxPos; iter++)
                {
                    modgradphi[phase][iter] = phix[phase][iter]*phix[phase][iter];
                    if (DIMENSION == 2)
                    {
                        modgradphi[phase][iter] += phiy[phase][iter]*phiy[phase][iter];
                    }
                    else if (DIMENSION == 3)
                    {
                        modgradphi[phase][iter] += phiy[phase][iter]*phiy[phase][iter];
                        modgradphi[phase][iter] += phiz[phase][iter]*phiz[phase][iter];
                    }

                    modgradphi[phase][iter] = sqrt(modgradphi[phase][iter]);
                }
            }

            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                for (long iter = 0; iter < maxPos; iter++)
                {
                    scalprodct[phase][iter] = -1.0*(phix[phase][iter]*phix[NUMPHASES-1][iter] + phiy[phase][iter]*phiy[NUMPHASES-1][iter] + phiz[phase][iter]*phiz[NUMPHASES-1][iter]);

                    if (modgradphi[NUMPHASES-1][iter] > 0.0)
                    {
                        scalprodct[phase][iter] /= (modgradphi[phase][iter]*modgradphi[NUMPHASES-1][iter]);
                    }
                }
            }

            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                for (component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    double diffLocal = 1.0 - diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]/diffusivity[(component + (NUMPHASES-1)*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component];

                    jat[phase][component][1] = (alphidot[phase][component][1]*phix[phase][1]/modgradphi[phase][1])*diffLocal*fabs(scalprodct[phase][1]);
                    jat[phase][component][2] = (alphidot[phase][component][2]*phix[phase][2]/modgradphi[phase][2])*diffLocal*fabs(scalprodct[phase][2]);

                    if (DIMENSION == 2)
                    {
                        jat[phase][component][3] = (alphidot[phase][component][3]*phiy[phase][3]/modgradphi[phase][3])*diffLocal*fabs(scalprodct[phase][3]);
                        jat[phase][component][4] = (alphidot[phase][component][4]*phiy[phase][4]/modgradphi[phase][4])*diffLocal*fabs(scalprodct[phase][4]);
                    }
                    else if (DIMENSION == 3)
                    {
                        jat[phase][component][3] = ((alphidot[phase][component][3]*phiy[phase][3])/(modgradphi[phase][3]))*diffLocal*fabs(scalprodct[phase][3]);
                        jat[phase][component][4] = ((alphidot[phase][component][4]*phiy[phase][4])/(modgradphi[phase][4]))*diffLocal*fabs(scalprodct[phase][4]);

                        jat[phase][component][5] = ((alphidot[phase][component][5]*phiz[phase][5])/(modgradphi[phase][5]))*diffLocal*fabs(scalprodct[phase][5]);
                        jat[phase][component][6] = ((alphidot[phase][component][6]*phiz[phase][6])/(modgradphi[phase][6]))*diffLocal*fabs(scalprodct[phase][6]);
                    }
                }
            }

            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                for (component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    for (long iter = 0; iter < maxPos; iter++)
                    {
                        if (modgradphi[phase][iter] == 0.0)
                        {
                            jat[phase][component][iter] = 0.0;
                        }
                    }
                }
            }

            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                for (component = 0; component < NUMCOMPONENTS-1; component++)
                {
                    jatc[phase][component] = (jat[phase][component][1] - jat[phase][component][2])/DELTA_X;
                    if (DIMENSION >= 2)
                        jatc[phase][component] += (jat[phase][component][3] - jat[phase][component][4])/DELTA_Y;
                    if (DIMENSION == 3)
                        jatc[phase][component] += (jat[phase][component][5] - jat[phase][component][6])/DELTA_Z;
                }
            }

            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (phase = 0; phase < NUMPHASES-1; phase++)
                {
                    jatr[component] += jatc[phase][component];
                }
            }
        }
        else
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                jatr[component] = 0.0;
            }
        }

        for (component = 0; component < NUMCOMPONENTS-1; component++)
        {
            J_xp = 0.0;
            J_xm = 0.0;
            J_yp = 0.0;
            J_ym = 0.0;
            J_zp = 0.0;
            J_zm = 0.0;

            for (component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
            {
                effMobility[0] = 0.0;
                effMobility[1] = 0.0;
                effMobility[2] = 0.0;
                effMobility[3] = 0.0;
                effMobility[4] = 0.0;
                effMobility[5] = 0.0;
                effMobility[6] = 0.0;

                for (phase = 0; phase < NUMPHASES; phase++)
                {
                    mobilityLocal = mobility[(component2 + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component];

                    if (phi[phase][index[1][1][1]] > 0.999)
                        effMobility[0] = mobilityLocal;
                    else if (phi[phase][index[1][1][1]] > 0.001)
                        effMobility[0] += mobilityLocal*calcInterp5th(phi, phase, index[1][1][1], NUMPHASES);

                    if (phi[phase][index[2][1][1]] > 0.999)
                        effMobility[1] = mobilityLocal;
                    else if (phi[phase][index[2][1][1]] > 0.001)
                        effMobility[1] += mobilityLocal*calcInterp5th(phi, phase, index[2][1][1], NUMPHASES);

                    if (phi[phase][index[0][1][1]] > 0.999)
                        effMobility[2] = mobilityLocal;
                    else if (phi[phase][index[0][1][1]] > 0.001)
                        effMobility[2] += mobilityLocal*calcInterp5th(phi, phase, index[0][1][1], NUMPHASES);

                    if (DIMENSION == 2)
                    {
                        if (phi[phase][index[1][2][1]] > 0.999)
                            effMobility[3] = mobilityLocal;
                        else if (phi[phase][index[1][2][1]] > 0.001)
                            effMobility[3] += mobilityLocal*calcInterp5th(phi, phase, index[1][2][1], NUMPHASES);

                        if (phi[phase][index[1][0][1]] > 0.999)
                            effMobility[4] = mobilityLocal;
                        else if (phi[phase][index[1][0][1]] > 0.001)
                            effMobility[4] += mobilityLocal*calcInterp5th(phi, phase, index[1][0][1], NUMPHASES);
                    }
                    else if (DIMENSION == 3)
                    {
                        if (phi[phase][index[1][2][1]] > 0.999)
                            effMobility[3] = mobilityLocal;
                        else if (phi[phase][index[1][2][1]] > 0.001)
                            effMobility[3] += mobilityLocal*calcInterp5th(phi, phase, index[1][2][1], NUMPHASES);

                        if (phi[phase][index[1][0][1]] > 0.999)
                            effMobility[4] = mobilityLocal;
                        else if (phi[phase][index[1][0][1]] > 0.001)
                            effMobility[4] += mobilityLocal*calcInterp5th(phi, phase, index[1][0][1], NUMPHASES);

                        if (phi[phase][index[1][1][2]] > 0.999)
                            effMobility[5] = mobilityLocal;
                        else if (phi[phase][index[1][1][2]] > 0.001)
                            effMobility[5] += mobilityLocal*calcInterp5th(phi, phase, index[1][1][2], NUMPHASES);

                        if (phi[phase][index[1][1][0]] > 0.999)
                            effMobility[6] = mobilityLocal;
                        else if (phi[phase][index[1][1][0]] > 0.001)
                            effMobility[6] += mobilityLocal*calcInterp5th(phi, phase, index[1][1][0], NUMPHASES);
                    }
                }

                mu[0] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][1][1], NUMPHASES, NUMCOMPONENTS);

                mu[1] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[2][1][1], NUMPHASES, NUMCOMPONENTS);
                mu[2] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[0][1][1], NUMPHASES, NUMCOMPONENTS);

                if (DIMENSION == 2)
                {
                    mu[3] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][2][1], NUMPHASES, NUMCOMPONENTS);
                    mu[4] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][0][1], NUMPHASES, NUMCOMPONENTS);
                }
                else if (DIMENSION == 3)
                {
                    mu[3] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][2][1], NUMPHASES, NUMCOMPONENTS);
                    mu[4] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][0][1], NUMPHASES, NUMCOMPONENTS);

                    mu[5] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][1][2], NUMPHASES, NUMCOMPONENTS);
                    mu[6] = calcDiffusionPotential(phaseComp, NUMPHASES-1, component2, F0_A, F0_B, index[1][1][0], NUMPHASES, NUMCOMPONENTS);
                }

                /*
                 * Calculate flux by averaging effective mobility at the face centres, ...
                 * and using the difference in diffusion potential between nearest point to the face centre and cell centre
                 */

                J_xp += ((effMobility[1] + effMobility[0])/2.0)*(mu[1] - mu[0])/DELTA_X;
                J_xm += ((effMobility[0] + effMobility[2])/2.0)*(mu[0] - mu[2])/DELTA_X;

                if (DIMENSION == 2)
                {
                    J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[3] - mu[0])/DELTA_Y;
                    J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0] - mu[4])/DELTA_Y;

                    J_zp = J_zm = 0.0;
                }
                else if (DIMENSION == 3)
                {
                    J_yp += ((effMobility[3] + effMobility[0])/2.0)*(mu[3] - mu[0])/DELTA_Y;
                    J_ym += ((effMobility[0] + effMobility[4])/2.0)*(mu[0] - mu[4])/DELTA_Y;

                    J_zp += ((effMobility[5] + effMobility[0])/2.0)*(mu[5] - mu[0])/DELTA_Z;
                    J_zm += ((effMobility[0] + effMobility[6])/2.0)*(mu[0] - mu[6])/DELTA_Z;
                }
            }

            compNew[component][index[1][1][1]] = comp[component][index[1][1][1]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z + jatr[component]);
        }

    }
}

__global__
void __updateComposition_02__(double **phi, double **phiNew,
                              double **comp, double **compNew, double **mu,
                              double **phaseComp, long *thermo_phase,
                              double *diffusivity, double *kappaPhi, double *theta_ij,
                              double temperature, double molarVolume,
                              long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                              long sizeX, long sizeY, long sizeZ,
                              long xStep, long yStep, long padding,
                              double DELTA_X, double DELTA_Y, double DELTA_Z,
                              double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long index[3][3][3] = {-1};

    long idx[7] = {-1}, maxPos = 5;

    idx[0] = i*xStep + j*yStep + k;

    double muLocal[7];
    double effMobility[7];
    double J_xp = 0.0, J_xm = 0.0, J_yp = 0.0, J_ym = 0.0, J_zp = 0.0, J_zm = 0.0;

    // Antitrapping variables
    double alpha[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};
    double modgradphi[MAX_NUM_PHASES][7] = {0.0};
    double dphidt[MAX_NUM_PHASES][7] = {0.0};

    double gradx_phi[MAX_NUM_PHASES][5] = {0.0};
    double grady_phi[MAX_NUM_PHASES][5] = {0.0};
    double gradz_phi[MAX_NUM_PHASES][5] = {0.0};
    double phix[MAX_NUM_PHASES][7] = {0.0}, phiy[MAX_NUM_PHASES][7] = {0.0}, phiz[MAX_NUM_PHASES][7] = {0.0};

    double alphidot[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};
    double scalprodct[MAX_NUM_PHASES-1][7] = {0.0};

    double jatr[MAX_NUM_COMP-1] = {0.0};
    double jatc[MAX_NUM_PHASES-1][MAX_NUM_COMP-1] = {0.0};
    double jat[MAX_NUM_PHASES-1][MAX_NUM_COMP-1][7] = {0.0};

    long component, component2, phase;

    double tol = 1e-6;

    if (i >= padding && i < sizeX-padding && ((j >= padding && j < sizeY-padding && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k >= padding && k < sizeZ-padding && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        for (long x = 0; x < 3; x++)
        {
            for (long y = 0; y < 3; y++)
            {
                for (long z = 0; z < 3; z++)
                {
                    index[x][y][z] = (k+z-1) + (j+y-1)*yStep + (i+x-1)*xStep;
                }
            }
        }

        if (DIMENSION == 3)
            maxPos = 7;
        else if (DIMENSION == 2)
            maxPos = 5;
        else
            maxPos = 3;

        // Antitrapping calculations
        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            double A1 = sqrt(2.0*kappaPhi[phase*NUMPHASES + NUMPHASES-1]/theta_ij[phase*NUMPHASES + NUMPHASES-1]);

            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                alpha[phase][component][0] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][1]] - phaseComp[component*NUMPHASES + phase][index[1][1][1]]);

                alpha[phase][component][1] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[2][1][1]] - phaseComp[component*NUMPHASES + phase][index[1][1][1]]);
                alpha[phase][component][2] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[0][1][1]] - phaseComp[component*NUMPHASES + phase][index[0][1][1]]);

                if (DIMENSION == 2)
                {
                    alpha[phase][component][3] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][2][1]] - phaseComp[component*NUMPHASES + phase][index[1][2][1]]);
                    alpha[phase][component][4] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][0][1]] - phaseComp[component*NUMPHASES + phase][index[1][0][1]]);
                }
                else if (DIMENSION == 3)
                {
                    alpha[phase][component][3] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][2][1]] - phaseComp[component*NUMPHASES + phase][index[1][2][1]]);
                    alpha[phase][component][4] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][0][1]] - phaseComp[component*NUMPHASES + phase][index[1][0][1]]);

                    alpha[phase][component][5] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][2]] - phaseComp[component*NUMPHASES + phase][index[1][1][2]]);
                    alpha[phase][component][6] = A1*(phaseComp[component*NUMPHASES + NUMPHASES-1][index[1][1][0]] - phaseComp[component*NUMPHASES + phase][index[1][1][0]]);
                }
            }

            dphidt[phase][0] = (phiNew[phase][index[1][1][1]] - phi[phase][index[1][1][1]])/DELTA_t;

            dphidt[phase][1] = (phiNew[phase][index[2][1][1]] - phi[phase][index[2][1][1]])/DELTA_t;
            dphidt[phase][2] = (phiNew[phase][index[0][1][1]] - phi[phase][index[0][1][1]])/DELTA_t;

            if (DIMENSION == 2)
            {
                dphidt[phase][3] = (phiNew[phase][index[1][2][1]] - phi[phase][index[1][2][1]])/DELTA_t;
                dphidt[phase][4] = (phiNew[phase][index[1][0][1]] - phi[phase][index[1][0][1]])/DELTA_t;
            }
            else if (DIMENSION == 3)
            {
                dphidt[phase][3] = (phiNew[phase][index[1][2][1]] - phi[phase][index[1][2][1]])/DELTA_t;
                dphidt[phase][4] = (phiNew[phase][index[1][0][1]] - phi[phase][index[1][0][1]])/DELTA_t;

                dphidt[phase][5] = (phiNew[phase][index[1][1][2]] - phi[phase][index[1][1][2]])/DELTA_t;
                dphidt[phase][6] = (phiNew[phase][index[1][1][0]] - phi[phase][index[1][1][0]])/DELTA_t;
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            gradx_phi[phase][0] = (phi[phase][index[2][1][1]] - phi[phase][index[0][1][1]])/(2.0*DELTA_X);

            if (DIMENSION == 2)
            {
                gradx_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[0][2][1]])/(2.0*DELTA_X);
                gradx_phi[phase][2] = (phi[phase][index[2][0][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_X);

                grady_phi[phase][0] = (phi[phase][index[1][2][1]] - phi[phase][index[1][0][1]])/(2.0*DELTA_Y);
                grady_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[2][0][1]])/(2.0*DELTA_Y);
                grady_phi[phase][2] = (phi[phase][index[0][2][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_Y);
            }
            else if (DIMENSION == 3)
            {
                gradx_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[0][2][1]])/(2.0*DELTA_X);
                gradx_phi[phase][2] = (phi[phase][index[2][0][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_X);
                gradx_phi[phase][3] = (phi[phase][index[2][1][2]] - phi[phase][index[0][1][2]])/(2.0*DELTA_X);
                gradx_phi[phase][4] = (phi[phase][index[2][1][0]] - phi[phase][index[0][1][0]])/(2.0*DELTA_X);

                grady_phi[phase][0] = (phi[phase][index[1][2][1]] - phi[phase][index[1][0][1]])/(2.0*DELTA_Y);
                grady_phi[phase][1] = (phi[phase][index[2][2][1]] - phi[phase][index[2][0][1]])/(2.0*DELTA_Y);
                grady_phi[phase][2] = (phi[phase][index[0][2][1]] - phi[phase][index[0][0][1]])/(2.0*DELTA_Y);
                grady_phi[phase][3] = (phi[phase][index[1][2][2]] - phi[phase][index[1][0][2]])/(2.0*DELTA_Y);
                grady_phi[phase][4] = (phi[phase][index[1][2][0]] - phi[phase][index[1][0][0]])/(2.0*DELTA_Y);

                gradz_phi[phase][0] = (phi[phase][index[1][1][2]] - phi[phase][index[1][1][0]])/(2.0*DELTA_Z);
                gradz_phi[phase][1] = (phi[phase][index[2][1][2]] - phi[phase][index[2][1][0]])/(2.0*DELTA_Z);
                gradz_phi[phase][2] = (phi[phase][index[0][1][2]] - phi[phase][index[0][1][0]])/(2.0*DELTA_Z);
                gradz_phi[phase][3] = (phi[phase][index[1][2][2]] - phi[phase][index[1][2][0]])/(2.0*DELTA_Z);
                gradz_phi[phase][4] = (phi[phase][index[1][0][2]] - phi[phase][index[1][0][0]])/(2.0*DELTA_Z);
            }
        }

        /*
         * 0 -> central node
         * 1 -> (i+0.5, j, k)
         * 2 -> (i-0.5, j, k)
         * 3 -> (i, j+0.5, k)
         * 4 -> (i, j-0.5, k)
         * 5 -> (i, j, k+0.5)
         * 6 -> (i, j, k-0.5)
         */

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            phix[phase][0] = gradx_phi[phase][0];
            phix[phase][1] = (phi[phase][index[2][1][1]] - phi[phase][index[1][1][1]])/(DELTA_X);
            phix[phase][2] = (phi[phase][index[1][1][1]] - phi[phase][index[0][1][1]])/(DELTA_X);

            if (DIMENSION == 2)
            {
                phix[phase][3] = (gradx_phi[phase][0] + gradx_phi[phase][1])/2.0;
                phix[phase][4] = (gradx_phi[phase][0] + gradx_phi[phase][2])/2.0;
            }
            else if (DIMENSION == 3)
            {
                phix[phase][3] = (gradx_phi[phase][0] + gradx_phi[phase][1])/2.0;
                phix[phase][4] = (gradx_phi[phase][0] + gradx_phi[phase][2])/2.0;
                phix[phase][5] = (gradx_phi[phase][0] + gradx_phi[phase][3])/2.0;
                phix[phase][6] = (gradx_phi[phase][0] + gradx_phi[phase][4])/2.0;
            }

            if (DIMENSION >= 2)
            {
                phiy[phase][0] = grady_phi[phase][0];
                phiy[phase][1] = (grady_phi[phase][0] + grady_phi[phase][1])/2.0;
                phiy[phase][2] = (grady_phi[phase][0] + grady_phi[phase][2])/2.0;
                phiy[phase][3] = (phi[phase][index[1][2][1]] - phi[phase][index[1][1][1]])/(DELTA_Y);
                phiy[phase][4] = (phi[phase][index[1][1][1]] - phi[phase][index[1][0][1]])/(DELTA_Y);

                if (DIMENSION ==  3)
                {
                    phiy[phase][5] = (grady_phi[phase][0] + grady_phi[phase][3])/2.0;
                    phiy[phase][6] = (grady_phi[phase][0] + grady_phi[phase][4])/2.0;

                    phiz[phase][0] = gradz_phi[phase][0];
                    phiz[phase][1] = (gradz_phi[phase][0] + gradz_phi[phase][1])/2.0;
                    phiz[phase][2] = (gradz_phi[phase][0] + gradz_phi[phase][2])/2.0;
                    phiz[phase][3] = (gradz_phi[phase][0] + gradz_phi[phase][3])/2.0;
                    phiz[phase][4] = (gradz_phi[phase][0] + gradz_phi[phase][4])/2.0;
                    phiz[phase][5] = (phi[phase][index[1][1][2]] - phi[phase][index[1][1][1]])/(DELTA_Z);
                    phiz[phase][6] = (phi[phase][index[1][1][1]] - phi[phase][index[1][1][0]])/(DELTA_Z);
                }
            }
        }

        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                alphidot[phase][component][1] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][1]*dphidt[phase][1]))/2.0;
                alphidot[phase][component][2] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][2]*dphidt[phase][2]))/2.0;

                if (DIMENSION == 2)
                {
                    alphidot[phase][component][3] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][3]*dphidt[phase][3]))/2.0;
                    alphidot[phase][component][4] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][4]*dphidt[phase][4]))/2.0;
                }
                else if (DIMENSION == 3)
                {
                    alphidot[phase][component][3] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][3]*dphidt[phase][3]))/2.0;
                    alphidot[phase][component][4] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][4]*dphidt[phase][4]))/2.0;

                    alphidot[phase][component][5] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][5]*dphidt[phase][5]))/2.0;
                    alphidot[phase][component][6] = ((alpha[phase][component][0]*dphidt[phase][0]) + (alpha[phase][component][6]*dphidt[phase][6]))/2.0;
                }
            }
        }

        for (phase = 0; phase < NUMPHASES; phase++)
        {
            for (long iter = 0; iter < maxPos; iter++)
            {
                modgradphi[phase][iter] = phix[phase][iter]*phix[phase][iter];
                if (DIMENSION == 2)
                {
                    modgradphi[phase][iter] += phiy[phase][iter]*phiy[phase][iter];
                }
                else if (DIMENSION == 3)
                {
                    modgradphi[phase][iter] += phiy[phase][iter]*phiy[phase][iter];
                    modgradphi[phase][iter] += phiz[phase][iter]*phiz[phase][iter];
                }

                modgradphi[phase][iter] = sqrt(modgradphi[phase][iter]);
            }
        }

        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            for (long iter = 0; iter < maxPos; iter++)
            {
                scalprodct[phase][iter] = -1.0*(phix[phase][iter]*phix[NUMPHASES-1][iter] + phiy[phase][iter]*phiy[NUMPHASES-1][iter] + phiz[phase][iter]*phiz[NUMPHASES-1][iter]);

                if (modgradphi[NUMPHASES-1][iter] > 0.0)
                {
                    scalprodct[phase][iter] /= (modgradphi[phase][iter]*modgradphi[NUMPHASES-1][iter]);
                }
            }
        }

        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                double diffLocal = 1.0 - diffusivity[(component + phase*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component]/diffusivity[(component + (NUMPHASES-1)*(NUMCOMPONENTS-1))*(NUMCOMPONENTS-1) + component];

                jat[phase][component][1] = (alphidot[phase][component][1]*phix[phase][1]/modgradphi[phase][1])*diffLocal*fabs(scalprodct[phase][1]);
                jat[phase][component][2] = (alphidot[phase][component][2]*phix[phase][2]/modgradphi[phase][2])*diffLocal*fabs(scalprodct[phase][2]);

                if (DIMENSION == 2)
                {
                    jat[phase][component][3] = (alphidot[phase][component][3]*phiy[phase][3]/modgradphi[phase][3])*diffLocal*fabs(scalprodct[phase][3]);
                    jat[phase][component][4] = (alphidot[phase][component][4]*phiy[phase][4]/modgradphi[phase][4])*diffLocal*fabs(scalprodct[phase][4]);
                }
                else if (DIMENSION == 3)
                {
                    jat[phase][component][3] = ((alphidot[phase][component][3]*phiy[phase][3])/(modgradphi[phase][3]))*diffLocal*fabs(scalprodct[phase][3]);
                    jat[phase][component][4] = ((alphidot[phase][component][4]*phiy[phase][4])/(modgradphi[phase][4]))*diffLocal*fabs(scalprodct[phase][4]);

                    jat[phase][component][5] = ((alphidot[phase][component][5]*phiz[phase][5])/(modgradphi[phase][5]))*diffLocal*fabs(scalprodct[phase][5]);
                    jat[phase][component][6] = ((alphidot[phase][component][6]*phiz[phase][6])/(modgradphi[phase][6]))*diffLocal*fabs(scalprodct[phase][6]);
                }
            }
        }

        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                for (long iter = 0; iter < maxPos; iter++)
                {
                    if (modgradphi[phase][iter] == 0.0)
                    {
                        jat[phase][component][iter] = 0.0;
                    }
                }
            }
        }

        for (phase = 0; phase < NUMPHASES-1; phase++)
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                jatc[phase][component] = (jat[phase][component][1] - jat[phase][component][2])/DELTA_X;
                if (DIMENSION >= 2)
                    jatc[phase][component] += (jat[phase][component][3] - jat[phase][component][4])/DELTA_Y;
                if (DIMENSION == 3)
                    jatc[phase][component] += (jat[phase][component][5] - jat[phase][component][6])/DELTA_Z;
            }
        }

        for (component = 0; component < NUMCOMPONENTS-1; component++)
        {
            for (phase = 0; phase < NUMPHASES-1; phase++)
            {
                jatr[component] += jatc[phase][component];
            }
        }

        /******End of antitrapping calcs*********/
        /*
         *
         *
         */
        /******Start Function_F02 calcs**********/

        // x-direction
        idx[1] = (i+1)*xStep + j*yStep + k;
        idx[2] = (i-1)*xStep + j*yStep + k;

        // y-direction
        if (DIMENSION >= 2)
        {
            idx[3] = i*xStep + (j+1)*yStep + k;
            idx[4] = i*xStep + (j-1)*yStep + k;
        }

        // z-direction
        if (DIMENSION == 3)
        {
            idx[5] = i*xStep + j*yStep + k+1;
            idx[6] = i*xStep + j*yStep + k-1;
        }

        double dmudc[(MAX_NUM_COMP)*(MAX_NUM_COMP)];
        double y[MAX_NUM_COMP];
        double dmudcInv[MAX_NUM_COMP][MAX_NUM_COMP];
        int P[MAX_NUM_COMP];
        double mobility[MAX_NUM_COMP][MAX_NUM_COMP];

        long interface = 1, bulkphase = 0;

        for (long is = 0; is < NUMPHASES; is++)
        {
            if (phi[is][index[1][1][1]] > 0.99999)
            {
                bulkphase = is;
                interface = 0;
                break;
            }
        }

        if (interface)
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                // Fluxes
                J_xp = 0.0;
                J_xm = 0.0;
                J_yp = 0.0;
                J_ym = 0.0;
                J_zp = 0.0;
                J_zm = 0.0;

                // Computing the inner derivative and mobilities to get the fluxes
                for (component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
                {
                    for (long iter = 0; iter < 7; iter++)
                        effMobility[iter] = 0.0;

                    for (long pos = 0; pos < maxPos; pos++)
                    {
                        // M_{ij} = \sum_{\phi} M(\phi) = \sum_{\phi} D*dcdmu
                        for (phase = 0; phase < NUMPHASES; phase++)
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

                compNew[component][idx[0]] = comp[component][idx[0]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z + jatr[component]);
            }
        }
        else
        {
            for (component = 0; component < NUMCOMPONENTS-1; component++)
            {
                // Fluxes
                J_xp = 0.0;
                J_xm = 0.0;
                J_yp = 0.0;
                J_ym = 0.0;
                J_zp = 0.0;
                J_zm = 0.0;

                // Computing the inner derivative and mobilities to get the fluxes
                for (component2 = 0; component2 < NUMCOMPONENTS-1; component2++)
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

                compNew[component][idx[0]] = comp[component][idx[0]] + DELTA_t*((J_xp - J_xm)/DELTA_X + (J_yp - J_ym)/DELTA_Y + (J_zp - J_zm)/DELTA_Z + jatr[component]);
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
                     long xStep, long yStep, long padding,
                     double DELTA_X, double DELTA_Y, double DELTA_Z,
                     double DELTA_t)
{
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = i*xStep + j*yStep + k;

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
                       domainInfo simDomain, controls simControls,
                       simParameters simParams, subdomainInfo subdomain,
                       dim3 gridSize, dim3 blockSize)
{
    if (simControls.FUNCTION_F == 1 || simControls.FUNCTION_F == 3 || simControls.FUNCTION_F == 4)
    {
        __updateComposition__<<<gridSize, blockSize>>>(phi, phiNew,
                                                       comp, compNew,
                                                       phaseComp,
                                                       simParams.F0_A_dev, simParams.F0_B_dev,
                                                       simParams.mobility_dev, simParams.diffusivity_dev, simParams.kappaPhi_dev, simParams.theta_ij_dev,
                                                       simDomain.numPhases, simDomain.numComponents, simDomain.DIMENSION,
                                                       subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                       subdomain.xStep, subdomain.yStep, subdomain.padding, simControls.antiTrapping,
                                                       simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z,
                                                       simControls.DELTA_t);

        // applyBoundaryCondition(compNew, 2, simDomain.numComponents-1,
        //                        simDomain, simControls,
        //                        simParams, subdomain,
        //                        gridSize, blockSize);
    }
    else if (simControls.FUNCTION_F == 2)
    {
        __updateComposition_02__<<<gridSize, blockSize>>>(phi, phiNew,
                                                          comp, compNew,
                                                          mu, phaseComp, simDomain.thermo_phase_dev,
                                                          simParams.diffusivity_dev, simParams.kappaPhi_dev, simParams.theta_ij_dev,
                                                          simParams.T, simParams.molarVolume,
                                                          simDomain.numPhases, simDomain.numComponents, simDomain.DIMENSION,
                                                          subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                          subdomain.xStep, subdomain.yStep, subdomain.padding,
                                                          simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z,
                                                          simControls.DELTA_t);

        applyBoundaryCondition(compNew, 2, simDomain.numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        __updateMu_02__<<<gridSize, blockSize>>>(phi, comp,
                                                 phiNew, compNew,
                                                 phaseComp, mu,
                                                 simDomain.thermo_phase_dev, simParams.T, simParams.molarVolume,
                                                 simDomain.numPhases, simDomain.numComponents, simDomain.DIMENSION,
                                                 subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ,
                                                 subdomain.xStep, subdomain.yStep, subdomain.padding,
                                                 simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z,
                                                 simControls.DELTA_t);

        applyBoundaryCondition(mu, 1, simDomain.numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

    }
}
