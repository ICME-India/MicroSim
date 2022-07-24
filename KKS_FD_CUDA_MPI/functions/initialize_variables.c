#include "initialize_variables.h"

void decomposeDomain(domainInfo simDomain, subdomainInfo *subdomain,
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
}

void moveParamsToGPU(domainInfo *simDomain, simParameters *simParams)
{
    cudaMemcpy(simDomain->thermo_phase_dev, simDomain->thermo_phase_host, sizeof(int)*simDomain->numPhases, cudaMemcpyHostToDevice);

    for (int i = 0; i < simDomain->numPhases; i++)
    {
        simParams->theta_i_host[i] = 0.0;

        for (int j = 0; j < simDomain->numPhases; j++)
        {
            if (i != j)
            {
                simParams->kappaPhi_host[i][j] = 3.0/(2.0*simParams->alpha) * (simParams->gamma_host[i][j]*simParams->epsilon);
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
                cudaMemcpy(&simParams->diffusivity_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->diffusivity_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
                cudaMemcpy(&simParams->F0_A_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->F0_A_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
            }
        }
    }
}

void calcKernelParams(dim3 *gridSize, dim3 *blockSize, domainInfo simDomain, controls simControls, subdomainInfo *subdomain)
{
    if (simDomain.DIMENSION == 2)
    {
        blockSize->x = 16;
        blockSize->y = 16;
        blockSize->z = 1;

        subdomain->shiftPointer = simControls.padding*subdomain->sizeX;
        subdomain->sizeY += 2*simControls.padding;
        subdomain->numCompCells = subdomain->numCells + 2*subdomain->shiftPointer;
        subdomain->padding = simControls.padding;

        subdomain->sizeZ = 1;
    }
    else if (simDomain.DIMENSION == 3)
    {
        blockSize->x = 8;
        blockSize->y = 8;
        blockSize->z = 4;

        subdomain->shiftPointer = simControls.padding*subdomain->sizeX*subdomain->sizeY;
        subdomain->sizeZ += 2*simControls.padding;
        subdomain->numCompCells = subdomain->numCells + 2*subdomain->shiftPointer;
        subdomain->padding = simControls.padding;
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
