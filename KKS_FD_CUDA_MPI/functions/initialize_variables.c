#include "initialize_variables.h"

/*
 *  Distributing the domain to all the MPI processes
 *  The subdomain local and global coordinates will be initialized here
 *  Each worker will be assigned its respective neighbours in this function
 *  Decomposition is presently only done in 1 dimension (slab)
 */
void decomposeDomain(domainInfo simDomain, controls *simControls, subdomainInfo *subdomain,
                     int rank, int size)
{
    subdomain->rank = rank;

    // No decomposition along the x-axis regardless of domain dimension
    subdomain->xS = 0;
    subdomain->xE = simDomain.MESH_X - 1;
    subdomain->numCells = (subdomain->xE - subdomain->xS + 1);

    if (simDomain.DIMENSION == 3 && simDomain.MESH_Z >= size)
    {
        // Decomposing along the y-axis
        subdomain->yS = 0;
        subdomain->yE = simDomain.MESH_Y - 1;
        subdomain->numCells *= (subdomain->yE - subdomain->yS + 1);

        if (simDomain.MESH_Z % size == 0)
        {
            // Decomposing along the z-axis
            subdomain->zS = simDomain.MESH_Z*rank/size;
            subdomain->zE = simDomain.MESH_Z*(rank + 1)/size - 1;
        }
        else
        {
            if (rank < simDomain.MESH_Z % size)
            {
                subdomain->zS = simDomain.MESH_Z*rank/size + rank;
                subdomain->zE = simDomain.MESH_Z*(rank + 1)/size + rank;
            }
            else
            {
                subdomain->zS = simDomain.MESH_Z*rank/size + (simDomain.MESH_Z % size);
                subdomain->zE = simDomain.MESH_Z*(rank + 1)/size + (simDomain.MESH_Z % size) - 1;
            }
        }

        subdomain->numCells *= (subdomain->zE - subdomain->zS) + 1;

        // Setting neighbours for every block
        if (rank == 0)
        {
            subdomain->nbBack = size - 1;
            subdomain->nbFront = rank + 1;
        }
        else if (rank == size - 1)
        {
            subdomain->nbBack = rank - 1;
            subdomain->nbFront = 0;
        }
        else
        {
            subdomain->nbBack = rank - 1;
            subdomain->nbFront = rank + 1;
        }

    }

    else if (simDomain.DIMENSION == 2)
    {
        if (simDomain.MESH_Y % size == 0)
        {
            // Decomposing along the y-axis
            subdomain->yS = simDomain.MESH_Y*rank/size;
            subdomain->yE = simDomain.MESH_Y*(rank + 1)/size - 1;
        }
        else
        {
            if (rank < simDomain.MESH_Y % size)
            {
                subdomain->yS = simDomain.MESH_Y*rank/size + rank;
                subdomain->yE = simDomain.MESH_Y*(rank + 1)/size + rank;
            }
            else
            {
                subdomain->yS = simDomain.MESH_Y*rank/size + (simDomain.MESH_Y % size);
                subdomain->yE = simDomain.MESH_Y*(rank + 1)/size + (simDomain.MESH_Y % size) - 1;
            }
        }
        subdomain->zS = 0;
        subdomain->zE = 0;

        subdomain->numCells *= (subdomain->yE - subdomain->yS) + 1;

        // Setting neighbours for every block
        if (rank == 0)
        {
            subdomain->nbBack = size - 1;
            subdomain->nbFront = rank + 1;
        }
        else if (rank == size - 1)
        {
            subdomain->nbBack = rank - 1;
            subdomain->nbFront = 0;
        }
        else
        {
            subdomain->nbBack = rank - 1;
            subdomain->nbFront = rank + 1;
        }
    }

    subdomain->sizeX = subdomain->xE - subdomain->xS + 1;
    subdomain->sizeY = subdomain->yE - subdomain->yS + 1;
    subdomain->sizeZ = subdomain->zE - subdomain->zS + 1;

    if (simDomain.DIMENSION == 2)
    {
        subdomain->shiftPointer = simControls->padding*subdomain->sizeX;
        subdomain->sizeY += 2*simControls->padding;
        subdomain->numCompCells = subdomain->numCells + 2*subdomain->shiftPointer;
        subdomain->padding = simControls->padding;

        subdomain->sizeZ = 1;
    }
    else if (simDomain.DIMENSION == 3)
    {
        subdomain->shiftPointer = simControls->padding*subdomain->sizeX*subdomain->sizeY;
        subdomain->sizeZ += 2*simControls->padding;
        subdomain->numCompCells = subdomain->numCells + 2*subdomain->shiftPointer;
        subdomain->padding = simControls->padding;
    }

}

/*
 *  All simulation constants are calculated here using parameters read from the input file
 *
 *  Since device-side data is best stored in 1-D contiguous arrays, the data transfer is also done here
 *  in an element-by-element manner.
 */
void moveParamsToGPU(domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    // Calculating kappa, theta
    for (int i = 0; i < simDomain->numPhases; i++)
    {
        simParams->theta_i_host[i] = 0.0;

        for (int j = 0; j < simDomain->numPhases; j++)
        {
            if (i != j)
            {
                simParams->kappaPhi_host[i][j] =  (3.0*simParams->gamma_host[i][j]*simParams->epsilon)/(2.0*simParams->alpha);
                simParams->relax_coeff_host[i][j] = 1.0/(simParams->Tau_host[i][j]*simParams->epsilon);
                simParams->theta_ij_host[i][j] = 6.0 * simParams->alpha * simParams->gamma_host[i][j] / simParams->epsilon;
            }
            else
            {
                simParams->kappaPhi_host[i][j] = 0.0;
                simParams->theta_ij_host[i][j] = 0.0;
                simParams->relax_coeff_host[i][j] = 0.0;
            }
        }
    }

    double Tol = 1e-6;
    int P[simDomain->numComponents];

    if (simControls->FUNCTION_F != 2)
    {
        double ***dmudc = malloc3M(simDomain->numPhases, simDomain->numComponents-1, simDomain->numComponents-1);
        double **inverted = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);

        for (long i = 0; i < simDomain->numPhases; i++)
        {
            for (long j = 0; j < simDomain->numComponents-1; j++)
            {
                for (long k = 0; k < simDomain->numComponents-1; k++)
                {
                    if (j == k)
                        dmudc[i][j][k] = 2.0*simParams->F0_A_host[i][j][k];
                    else
                        dmudc[i][j][k] = simParams->F0_A_host[i][j][k];
                }
            }
        }

        // Mobility matrix
        for (long i = 0; i < simDomain->numPhases; i++)
        {
            LUPDecompose(dmudc[i], simDomain->numComponents-1, Tol, P);
            LUPInvert(dmudc[i], P, simDomain->numComponents-1, inverted);

            matrixMultiply(simParams->diffusivity_host[i], inverted, simParams->mobility_host[i], simDomain->numComponents-1);
        }

        // Relaxation Coefficients
        double **mobilityInv = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);
        double deltac[simDomain->numComponents-1], deltamu[simDomain->numComponents-1];
        double sum = 0.0;
        double minTau = 1e12;

        // Reusing inverted to hold mobility of the reference phase
        for (long i = 0; i < simDomain->numComponents-1; i++)
        {
            for (long j = 0; j < simDomain->numComponents-1; j++)
            {
                inverted[i][j] = simParams->mobility_host[simDomain->numPhases-1][i][j];
            }
        }

        LUPDecompose(inverted, simDomain->numComponents-1, Tol, P);
        LUPInvert(inverted, P, simDomain->numComponents-1, mobilityInv);

        for (long a = 0; a < simDomain->numPhases-1; a++)
        {
            for (long k = 0; k < simDomain->numComponents-1; k++)
                deltac[k] = simParams->ceq_host[a][simDomain->numPhases-1][k] - simParams->ceq_host[a][a][k];

            for (long i = 0; i < simDomain->numComponents-1; i++)
            {
                deltamu[i] = 0.0;

                for (long j = 0; j < simDomain->numComponents-1; j++)
                {
                    deltamu[i] += deltac[j]*mobilityInv[i][j];
                }
            }

            sum = 0.0;
            for (long i = 0; i < simDomain->numComponents-1; i++)
            {
                sum += deltac[i]*deltamu[i];
            }

            simParams->Tau_host[a][simDomain->numPhases-1] = (6.0*0.783333*simParams->kappaPhi_host[a][simDomain->numPhases-1]*sum)/(simParams->theta_ij_host[a][simDomain->numPhases-1]*simParams->molarVolume);
            simParams->Tau_host[simDomain->numPhases-1][a] = simParams->Tau_host[a][simDomain->numPhases-1];

            if (a == 0)
                minTau = simParams->Tau_host[a][simDomain->numPhases-1];

            if (simParams->Tau_host[a][simDomain->numPhases-1] < minTau)
                minTau = simParams->Tau_host[a][simDomain->numPhases-1];
        }

        for (long a = 0; a < simDomain->numPhases; a++)
        {
            for (long b = 0; b < simDomain->numPhases; b++)
            {
                simParams->Tau_host[a][b] = minTau;
                simParams->relax_coeff_host[a][b] = 1.0/minTau;
            }
        }

        free2M(mobilityInv, simDomain->numComponents-1);
        free2M(inverted, simDomain->numComponents-1);
        free3M(dmudc, simDomain->numPhases, simDomain->numComponents-1);
    }
    else if (simControls->FUNCTION_F == 2)
    {
        // Relaxation Coefficients
        double **diffInv = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);
        double **diff = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);
        double dmudc[(simDomain->numComponents-1)*(simDomain->numComponents-1)];
        double **mobilityInv = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);
        double deltac[simDomain->numComponents-1], deltamu[simDomain->numComponents-1];
        double y[simDomain->numComponents-1];
        double sum = 0.0;
        double minTau = 1e12;

        // Get diffusivity of matrix phase
        for (long i = 0; i < simDomain->numComponents-1; i++)
        {
            for (long j = 0; j < simDomain->numComponents-1; j++)
            {
                diff[i][j] = simParams->diffusivity_host[simDomain->numPhases-1][i][j];
            }
        }

        // Get D^{-1}
        LUPDecompose(diff, simDomain->numComponents-1, Tol, P);
        LUPInvert(diff, P, simDomain->numComponents-1, diffInv);

        // Get equilibrium compositions of matrix phase
        for (long i = 0; i < simDomain->numComponents-1; i++)
        {
            y[i] = simParams->ceq_host[simDomain->numPhases-1][simDomain->numPhases-1][i];
            sum += y[i];
        }
        y[simDomain->numComponents-1] = 1.0 - sum;

        // Get dmudc in matrix
        (*dmudc_tdb[simDomain->thermo_phase_host[simDomain->numPhases-1]])(simParams->T, y, dmudc);

        // Get mobility inverse as dmudc*D^{-1}
        for (long i = 0; i < simDomain->numComponents-1; i++)
        {
            for (long j = 0; j < simDomain->numComponents-1; j++)
            {
                mobilityInv[i][j] = 0.0;

                for (long k = 0; k < simDomain->numComponents-1; k++)
                {
                    mobilityInv[i][j] += dmudc[i*(simDomain->numComponents-1) + k]*diffInv[k][j];
                }
            }
        }

        for (long a = 0; a < simDomain->numPhases-1; a++)
        {
            for (long k = 0; k < simDomain->numComponents-1; k++)
                deltac[k] = simParams->ceq_host[a][simDomain->numPhases-1][k] - simParams->ceq_host[a][a][k];

            for (long i = 0; i < simDomain->numComponents-1; i++)
            {
                deltamu[i] = 0.0;

                for (long j = 0; j < simDomain->numComponents-1; j++)
                {
                    deltamu[i] += deltac[j]*mobilityInv[i][j];
                }
            }

            sum = 0.0;
            for (long i = 0; i < simDomain->numComponents-1; i++)
            {
                sum += deltac[i]*deltamu[i];
            }

            simParams->Tau_host[a][simDomain->numPhases-1] = (6.0*0.783333*simParams->kappaPhi_host[a][simDomain->numPhases-1]*sum)/(simParams->theta_ij_host[a][simDomain->numPhases-1]*simParams->molarVolume);
            simParams->Tau_host[simDomain->numPhases-1][a] = simParams->Tau_host[a][simDomain->numPhases-1];

            if (a == 0)
                minTau = simParams->Tau_host[a][simDomain->numPhases-1];

            if (simParams->Tau_host[a][simDomain->numPhases-1] < minTau)
                minTau = simParams->Tau_host[a][simDomain->numPhases-1];
        }

        for (long a = 0; a < simDomain->numPhases; a++)
        {
            for (long b = 0; b < simDomain->numPhases; b++)
            {
                simParams->Tau_host[a][b] = minTau;
                simParams->relax_coeff_host[a][b] = 1.0/minTau;
            }
        }

        free2M(diff, simDomain->numComponents-1);
        free2M(mobilityInv, simDomain->numComponents-1);
        free2M(diffInv, simDomain->numComponents-1);
    }

    //printf("Relax coeff = %le\n", simParams->relax_coeff_host[0][1]);


    // Move to GPU
    cudaMemcpy(simDomain->thermo_phase_dev, simDomain->thermo_phase_host, sizeof(long)*simDomain->numPhases, cudaMemcpyHostToDevice);

    for (int i = 0; i < simDomain->numPhases; i++)
    {
        cudaMemcpy(&simParams->F0_C_dev[i], &simParams->F0_C_host[i], sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&simParams->theta_i_dev[i], &simParams->theta_i_host[i], sizeof(double), cudaMemcpyHostToDevice);

        for (int j = 0; j < simDomain->numPhases; j++)
        {
            cudaMemcpy(&simParams->gamma_dev[i*simDomain->numPhases + j], &simParams->gamma_host[i][j], sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(&simParams->relax_coeff_dev[i*simDomain->numPhases + j], &simParams->relax_coeff_host[i][j], sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(&simParams->theta_ij_dev[i*simDomain->numPhases + j], &simParams->theta_ij_host[i][j], sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(&simParams->kappaPhi_dev[i*simDomain->numPhases + j], &simParams->kappaPhi_host[i][j], sizeof(double),cudaMemcpyHostToDevice);

            for (int k = 0; k < simDomain->numPhases; k++)
            {
                cudaMemcpy(&simParams->theta_ijk_dev[(i*simDomain->numPhases + j)*simDomain->numPhases + k], &simParams->theta_ijk_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
            }

            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                cudaMemcpy(&simParams->ceq_dev[(i*simDomain->numPhases + j)*(simDomain->numComponents-1) + k], &simParams->ceq_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(&simParams->cfill_dev[(i*simDomain->numPhases + j)*(simDomain->numComponents-1) + k], &simParams->cfill_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(&simParams->cguess_dev[(i*simDomain->numPhases + j)*(simDomain->numComponents-1) + k], &simParams->cguess_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
            }
        }

        for (int j = 0; j < simDomain->numComponents-1; j++)
        {
            cudaMemcpy(&simParams->F0_B_dev[i*(simDomain->numComponents-1) + j], &simParams->F0_B_host[i][j], sizeof(double), cudaMemcpyHostToDevice);

            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                cudaMemcpy(&simParams->mobility_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->mobility_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(&simParams->diffusivity_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->diffusivity_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(&simParams->F0_A_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->F0_A_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
            }
        }
    }
}

/*
 *  Kernel launch parameters
 *  Number of blocks, and size of each block is calculated here
 */
void calcKernelParams(dim3 *gridSize, dim3 *blockSize, domainInfo simDomain, controls simControls, subdomainInfo *subdomain)
{
    if (simDomain.DIMENSION == 2)
    {
        blockSize->x = 16;
        blockSize->y = 16;
        blockSize->z = 1;
    }
    else if (simDomain.DIMENSION == 3)
    {
        blockSize->x = 8;
        blockSize->y = 8;
        blockSize->z = 4;
    }

    gridSize->x = ceil(subdomain->sizeX/blockSize->x);
    if (subdomain->sizeX % blockSize->x)
        gridSize->x += 1;

    gridSize->y = ceil(subdomain->sizeY/blockSize->y);
    if (subdomain->sizeY % blockSize->y)
        gridSize->y += 1;

    gridSize->z = ceil(subdomain->sizeZ/blockSize->z);
    if (subdomain->sizeZ % blockSize->z)
        gridSize->z += 1;
}
