#include "computeDrivingForce.cuh"

__global__
void __computeDrivingForce__(double **phi, double **comp,
                             double **dfdphi,
                             double **phaseComp,
                             double *F0_A, double *F0_B, double *F0_C,
                             double molarVolume,
                             double *theta_i, double *theta_ij, double *theta_ijk,
                             long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                             long sizeX, long sizeY, long sizeZ,
                             long yStep, long zStep, long padding)
{
    /*
     * Get thread coordinates
     */
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = k*zStep + j*yStep + i;

    /*
     * Calculate grand potential density for every phase
     */
    double psi = 0.0;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        for (long phase = 0; phase < NUMPHASES; phase++)
        {
            // f'_{bulk}
            dfdphi[phase][idx] = 0.0;

            for (long p = 0; p < NUMPHASES; p++)
            {
                /*
                 * \psi_{p} = f^{p} - \sum_{i=1}^{K-1} (c^{p}_{i}\mu_{i})
                 */
                psi = 0.0;

                psi += calcPhaseEnergy(phaseComp, p, F0_A, F0_B, F0_C, idx, NUMPHASES, NUMCOMPONENTS);
                for (long component = 0; component < NUMCOMPONENTS-1; component++)
                    psi -= calcDiffusionPotential(phaseComp, p, component, F0_A, F0_B, idx, NUMPHASES, NUMCOMPONENTS)*phaseComp[(component*NUMPHASES) + p][idx];

                /*
                 * \frac{\delta F}{\delta\phi_{phase}} += \sum_{p=1}^{N} (\frac{\partial h(\phi_{p})}{\partial \phi_{phase}} \frac{\psi_{p}}{V_{m}})
                 */
                psi *= calcInterp5thDiff(phi, p, phase, idx, NUMPHASES);

                dfdphi[phase][idx] += psi/molarVolume;
            }

            /*
             * Potential function
             * \frac{\delta F}{\delta\phi_{phase}} += \frac{\partial g(\phi_{phase})}{\partial\phi_{phase}}
             */
            dfdphi[phase][idx] += calcDoubleWellDerivative(phi, phase,
                                                           theta_i, theta_ij, theta_ijk,
                                                           idx, NUMPHASES);
        }
    }
}

__global__
void __computeDrivingForce_02__(double **phi, double **comp,
                                double **dfdphi, double **phaseComp,
                                double **mu,
                                double molarVolume,
                                double *theta_i, double *theta_ij, double *theta_ijk,
                                double temperature, long *thermo_phase,
                                long NUMPHASES, long NUMCOMPONENTS, long DIMENSION,
                                long sizeX, long sizeY, long sizeZ,
                                long yStep, long zStep, long padding)
{
    /*
     * Get thread coordinates
     */
    long i = threadIdx.x + blockIdx.x * blockDim.x;
    long j = threadIdx.y + blockIdx.y * blockDim.y;
    long k = threadIdx.z + blockIdx.z * blockDim.z;

    long idx = k*zStep + j*yStep + i;

    if (i < sizeX && ((j < sizeY && DIMENSION >= 2) || (DIMENSION == 1 && j == 0)) && ((k < sizeZ && DIMENSION == 3) || (DIMENSION < 3 && k == 0)))
    {
        /*
         * Calculate grand potential density for every phase
         */
        double y[MAX_NUM_COMP];
        double phaseEnergy = 0.0;
        double sum = 0.0;
        double psi = 0.0;

        for (long phase = 0; phase < NUMPHASES; phase++)
        {
            // f'_{bulk}
            dfdphi[phase][idx] = 0.0;

            for (long p = 0; p < NUMPHASES; p++)
            {
                sum = 0.0;

                for (long is = 0; is < NUMCOMPONENTS-1; is++)
                {
                    y[is] = phaseComp[is*NUMPHASES + p][idx];
                    sum += y[is];
                }

                y[NUMCOMPONENTS-1] = 1.0 - sum;

                (*free_energy_tdb_dev[thermo_phase[p]])(temperature, y, &phaseEnergy);

                /*
                 * \psi_{p} = f^{p} - \sum_{i=1}^{K-1} (c^{p}_{i}\mu_{i})
                 */
                psi = 0.0;

                psi += phaseEnergy;
                for (long component = 0; component < NUMCOMPONENTS-1; component++)
                    psi -= mu[component][idx]*phaseComp[(component*NUMPHASES + p)][idx];


                /*
                 * \frac{\delta F}{\delta\phi_{phase}} += \sum_{p=1}^{N} (\frac{\partial h(\phi_{p})}{\partial \phi_{phase}} \frac{\psi_{p}}{V_{m}})
                 */
                psi *= calcInterp5thDiff(phi, p, phase, idx, NUMPHASES);
                dfdphi[phase][idx] += psi/molarVolume;
            }

            /*
             * Potential function
             * \frac{\delta F}{\delta\phi_{phase}} += \frac{\partial g(\phi_{phase})}{\partial\phi_{phase}}
             */
            dfdphi[phase][idx] += calcDoubleWellDerivative(phi, phase,
                                                           theta_i, theta_ij, theta_ijk,
                                                           idx, NUMPHASES);
        }
    }
}

void computeDrivingForce(double **phi, double **comp,
                         double **dfdphi,
                         double **phaseComp, double **mu,
                         domainInfo* simDomain, controls* simControls,
                         simParameters* simParams, subdomainInfo* subdomain,
                         dim3 gridSize, dim3 blockSize)
{

    if (simControls->FUNCTION_F == 1 || simControls->FUNCTION_F == 3  || simControls->FUNCTION_F == 4)
    {
        __computeDrivingForce__<<<gridSize, blockSize>>>(phi, comp,
                                                         dfdphi,
                                                         phaseComp,
                                                         simParams->F0_A_dev, simParams->F0_B_dev, simParams->F0_C_dev,
                                                         simParams->molarVolume,
                                                         simParams->theta_i_dev, simParams->theta_ij_dev, simParams->theta_ijk_dev,
                                                         simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                         subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                         subdomain->yStep, subdomain->zStep, subdomain->padding);
    }
    else if (simControls->FUNCTION_F == 2)
    {
        __computeDrivingForce_02__<<<gridSize, blockSize>>>(phi, comp,
                                                            dfdphi, phaseComp,
                                                            mu,
                                                            simParams->molarVolume,
                                                            simParams->theta_i_dev, simParams->theta_ij_dev, simParams->theta_ijk_dev,
                                                            simParams->T, simDomain->thermo_phase_dev,
                                                            simDomain->numPhases, simDomain->numComponents, simDomain->DIMENSION,
                                                            subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ,
                                                            subdomain->yStep, subdomain->zStep, subdomain->padding);
    }
}
