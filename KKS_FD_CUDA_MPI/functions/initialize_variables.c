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
    // Assign subdomain its rank
    subdomain->rank = rank;
    subdomain->padding = simControls->padding;

    /*
     * In 3-D, decomposition is done along the z-axis
     */
    if (simDomain.DIMENSION == 3 && simDomain.MESH_Z >= size)
    {
        // x-axis
        subdomain->xS = 0;
        subdomain->xE = simDomain.MESH_X-1;
        subdomain->sizeX = subdomain->xE - subdomain->xS + 1 + 2*simControls->padding;
        subdomain->numCells = subdomain->sizeX - 2*simControls->padding;
        subdomain->numCompCells = subdomain->sizeX;

        // Decomposing along the y-axis
        subdomain->yS = 0;
        subdomain->yE = simDomain.MESH_Y-1;
        subdomain->sizeY = subdomain->yE - subdomain->yS + 1 + 2*simControls->padding;
        subdomain->numCells *= subdomain->sizeY - 2*simControls->padding;
        subdomain->numCompCells *= subdomain->sizeY;

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

        subdomain->sizeZ = subdomain->zE - subdomain->zS + 1 + 2*simControls->padding;
        subdomain->numCells *= subdomain->sizeZ - 2*simControls->padding;
        subdomain->numCompCells *= subdomain->sizeZ;

        subdomain->shiftPointer = simControls->padding*subdomain->sizeX*subdomain->sizeY;
        subdomain->yStep = subdomain->sizeX;
        subdomain->zStep = subdomain->sizeX*subdomain->sizeY;

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

    else if (simDomain.DIMENSION == 2 && simDomain.MESH_Y >= size)
    {
        // x-axis
        subdomain->xS = 0;
        subdomain->xE = simDomain.MESH_X-1;
        subdomain->sizeX = subdomain->xE - subdomain->xS + 1 + 2*simControls->padding;
        subdomain->numCells = subdomain->sizeX - 2*simControls->padding;
        subdomain->numCompCells = subdomain->sizeX;

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

        subdomain->sizeY = subdomain->yE - subdomain->yS + 1 + 2*simControls->padding;
        subdomain->numCells *= subdomain->sizeY - 2*simControls->padding;
        subdomain->numCompCells *= subdomain->sizeY;

        simDomain.MESH_Z = 1;
        subdomain->zS = 0;
        subdomain->zE = 0;
        subdomain->sizeZ = 3;

        subdomain->shiftPointer = simControls->padding*subdomain->sizeX;
        subdomain->yStep = subdomain->sizeX;
        subdomain->zStep = 0;

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
    else if (simDomain.DIMENSION == 1)
    {
        if (simDomain.MESH_X % size == 0)
        {
            // Decomposing along the y-axis
            subdomain->xS = simDomain.MESH_X*rank/size;
            subdomain->xE = simDomain.MESH_X*(rank + 1)/size - 1;
        }
        else
        {
            if (rank < simDomain.MESH_X % size)
            {
                subdomain->xS = simDomain.MESH_X*rank/size + rank;
                subdomain->xE = simDomain.MESH_X*(rank + 1)/size + rank;
            }
            else
            {
                subdomain->xS = simDomain.MESH_X*rank/size + (simDomain.MESH_X % size);
                subdomain->xE = simDomain.MESH_X*(rank + 1)/size + (simDomain.MESH_X % size) - 1;
            }
        }

        subdomain->sizeX = subdomain->xE - subdomain->xS + 1 + 2*simControls->padding;
        subdomain->numCells = subdomain->sizeX - 2*simControls->padding;
        subdomain->numCompCells = subdomain->sizeX;

        simDomain.MESH_Y = 1;
        subdomain->yS = 0;
        subdomain->yE = 0;
        subdomain->sizeY = 3;

        simDomain.MESH_Z = 1;
        subdomain->zS = 0;
        subdomain->zE = 0;
        subdomain->sizeZ = 3;

        subdomain->shiftPointer = simControls->padding;
        subdomain->yStep = 0;
        subdomain->zStep = 0;

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
}

/*
 *  All simulation constants are calculated here using parameters read from the input file
 *
 *  Since device-side data is best stored in 1-D contiguous arrays, the data transfer is also done here
 *  in an element-by-element manner.
 */
void moveParamsToGPU(domainInfo *simDomain, controls *simControls, simParameters *simParams)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;

    MPI_Comm_rank(comm, &rank);

    // Calculate kappa, theta
    for (int i = 0; i < simDomain->numPhases; i++)
    {
        simParams->theta_i_host[i] = 0.0;

        for (int j = 0; j < simDomain->numPhases; j++)
        {
            if (i != j)
            {
                simParams->kappaPhi_host[i][j] =  (3.0*simParams->gamma_host[i][j]*simParams->epsilon)/(2.0*simParams->alpha);
                simParams->theta_ij_host[i][j] = 6.0 * simParams->alpha * simParams->gamma_host[i][j] / simParams->epsilon;
            }
            else
            {
                simParams->kappaPhi_host[i][j] = 0.0;
                simParams->theta_ij_host[i][j] = 0.0;
            }
        }
    }

    // Calculate phase-field mobility
    calculateTau(simDomain, simControls, simParams);

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
            cudaMemcpy(&simParams->dab_dev[i*simDomain->numPhases + j], &simParams->dab_host[i][j], sizeof(double),cudaMemcpyHostToDevice);

            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    cudaMemcpy(&simParams->Rotation_matrix_dev[((i*simDomain->numPhases + j)*simDomain->numPhases + k)*3 + l], &simParams->Rotation_matrix_host[i][j][k][l], sizeof(double),cudaMemcpyHostToDevice);
                    cudaMemcpy(&simParams->Inv_Rotation_matrix_dev[((i*simDomain->numPhases + j)*simDomain->numPhases + k)*3 + l], &simParams->Inv_Rotation_matrix_host[i][j][k][l], sizeof(double),cudaMemcpyHostToDevice);
                }
            }

            for (int k = 0; k < simDomain->numPhases; k++)
            {
                simParams->theta_ijk_host[i][j][k] *= (6.0*simParams->alpha)/simParams->epsilon;
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
    if (simDomain.DIMENSION == 1)
    {
        blockSize->x = 256;
        blockSize->y = 1;
        blockSize->z = 1;

        gridSize->x = ceil(subdomain->sizeX/blockSize->x);
        if (subdomain->sizeX % blockSize->x)
            gridSize->x += 1;

        gridSize->y = 1;
        gridSize->z = 1;

    }
    else if (simDomain.DIMENSION == 2)
    {
        blockSize->x = 16;
        blockSize->y = 16;
        blockSize->z = 1;

        gridSize->x = ceil(subdomain->sizeX/blockSize->x);
        if (subdomain->sizeX % blockSize->x)
            gridSize->x += 1;

        gridSize->y = ceil(subdomain->sizeY/blockSize->y);
        if (subdomain->sizeY % blockSize->y)
            gridSize->y += 1;

        gridSize->z = 1;
    }
    else if (simDomain.DIMENSION == 3)
    {
        blockSize->x = 8;
        blockSize->y = 8;
        blockSize->z = 4;

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

    if (subdomain->rank == 0)
        printf("\nGrid size: (%ld, %ld, %ld)\nBlock size: (%ld, %ld, %ld)\n\n", gridSize->x, gridSize->y, gridSize->z, blockSize->x, blockSize->y, blockSize->z);
}
