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
    subdomain->size = size;
    subdomain->padding = simControls->padding;

    /*
     * In 3-D, decomposition is done along the z-axis
     */
    if (simDomain.DIMENSION == 3)
    {
        // z-axis
        subdomain->zS = 0;
        subdomain->zE = simDomain.MESH_Z;

        subdomain->zS_c = 0;
        subdomain->zE_c = simDomain.MESH_Z + 2*subdomain->padding;

        subdomain->zS_r = subdomain->padding;
        subdomain->zE_r = simDomain.MESH_Z + subdomain->padding;

        subdomain->numZ = subdomain->zE - subdomain->zS ;
        subdomain->sizeZ = subdomain->zE_c - subdomain->zS_c;
        subdomain->numCells = subdomain->numZ;
        subdomain->numCompCells = subdomain->sizeZ;

        // y-axis
        subdomain->yS = 0;
        subdomain->yE = simDomain.MESH_Y;

        subdomain->yS_c = 0;
        subdomain->yE_c = simDomain.MESH_Y + 2*subdomain->padding;

        subdomain->yS_r = subdomain->padding;
        subdomain->yE_r = simDomain.MESH_Y + subdomain->padding;

        subdomain->numY = subdomain->yE - subdomain->yS;
        subdomain->sizeY = subdomain->yE_c - subdomain->yS_c;
        subdomain->numCells *= subdomain->numY;
        subdomain->numCompCells *= subdomain->sizeY;
    }
    else if (simDomain.DIMENSION == 2)
    {
        // z-axis
        subdomain->zS = 0;
        subdomain->zE = 1;

        subdomain->zS_c = 0;
        subdomain->zE_c = 1;

        subdomain->zS_r = 0;
        subdomain->zE_r = 1;

        subdomain->numZ = subdomain->zE - subdomain->zS;
        subdomain->sizeZ = subdomain->zE_c - subdomain->zS_c;
        subdomain->numCells = subdomain->numZ;
        subdomain->numCompCells = subdomain->sizeZ;

        // y-axis
        subdomain->yS = 0;
        subdomain->yE = simDomain.MESH_Y;

        subdomain->yS_c = 0;
        subdomain->yE_c = simDomain.MESH_Y + 2*subdomain->padding;

        subdomain->yS_r = subdomain->padding;
        subdomain->yE_r = simDomain.MESH_Y + subdomain->padding;

        subdomain->numY = subdomain->yE - subdomain->yS;
        subdomain->sizeY = subdomain->yE_c - subdomain->yS_c;
        subdomain->numCells *= subdomain->numY;
        subdomain->numCompCells *= subdomain->sizeY;

    }
    else if (simDomain.DIMENSION == 1)
    {
        // z-axis
        subdomain->zS = 0;
        subdomain->zE = 1;

        subdomain->zS_c = 0;
        subdomain->zE_c = 1;

        subdomain->zS_r = 0;
        subdomain->zE_r = 1;

        subdomain->numZ = subdomain->zE - subdomain->zS;
        subdomain->sizeZ = subdomain->zE_c - subdomain->zS_c;
        subdomain->numCells = subdomain->numZ;
        subdomain->numCompCells = subdomain->sizeZ;

        // y-axis
        subdomain->yS = 0;
        subdomain->yE = 1;

        subdomain->yS_c = 0;
        subdomain->yE_c = 1;

        subdomain->yS_r = 0;
        subdomain->yE_r = 1;

        subdomain->numY = subdomain->yE - subdomain->yS;
        subdomain->sizeY = subdomain->yE_c - subdomain->yS_c;
        subdomain->numCells *= subdomain->numY;
        subdomain->numCompCells *= subdomain->sizeY;
    }

    // Decomposing along the x-axis
    int ranks_cutoff = simDomain.MESH_X % size;
    size_t my_nx = (simDomain.MESH_X / size) + (rank < ranks_cutoff ? 1 : 0);

    if (size > 1)
        subdomain->xS = (rank < ranks_cutoff ? rank*my_nx : ranks_cutoff*(my_nx+1) + (rank - ranks_cutoff)*my_nx);
    else
        subdomain->xS = 0;

    subdomain->xE = subdomain->xS + my_nx;

    subdomain->xS_c = 0;
    subdomain->xE_c = my_nx + 2*subdomain->padding;

    subdomain->xS_r = subdomain->padding;
    subdomain->xE_r = my_nx + subdomain->padding;

    subdomain->numX = subdomain->xE - subdomain->xS;
    subdomain->sizeX = subdomain->xE_c - subdomain->xS_c;
    subdomain->numCells *= subdomain->numX;
    subdomain->numCompCells *= subdomain->sizeX;

    subdomain->shiftPointer = simControls->padding*subdomain->sizeY*subdomain->sizeZ;

    subdomain->xStep = subdomain->sizeY*subdomain->sizeZ;
    subdomain->yStep = subdomain->sizeZ;
    subdomain->zStep = 1;

    // Setting neighbours for every block
    if (rank == 0 && size)
    {
        subdomain->nbBack = size - 1;
        subdomain->nbFront = rank + 1;
    }
    else if (rank == size - 1 && size)
    {
        subdomain->nbBack = rank - 1;
        subdomain->nbFront = 0;
    }
    else if (size)
    {
        subdomain->nbBack = rank - 1;
        subdomain->nbFront = rank + 1;
    }
    else
    {
        subdomain->nbBack = 0;
        subdomain->nbFront = 0;
    }

    if (rank == 0)
    {
        printf("\nCompleted domain decomposition\n");

        printf("Subdomain information for rank 0:\n");
        printf("Size: (%ld, %ld, %ld)\n", subdomain->sizeX, subdomain->sizeY, subdomain->sizeZ);
        printf("Boundary: (x: (%ld, %ld), y:(%ld, %ld), z:(%ld, %ld))\n", subdomain->xS, subdomain->xE, subdomain->yS, subdomain->yE, subdomain->zS, subdomain->zE);
        printf("Computational Boundary: (x: (%ld, %ld), y:(%ld, %ld), z:(%ld, %ld))\n", subdomain->xS_c, subdomain->xE_c, subdomain->yS_c, subdomain->yE_c, subdomain->zS_c, subdomain->zE_c);
        printf("Steps: (%ld, %ld, %ld)\n", subdomain->xStep, subdomain->yStep, subdomain->zStep);
        printf("Num. Cells: %ld, Num. Comp. Cells: %ld, Shift Pointer: %ld\n", subdomain->numCells, subdomain->numCompCells, subdomain->shiftPointer);
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

    simControls->antiTrapping = 1;

    // for (int i = 0; i < simDomain->numPhases; i++)
    // {
    //     for (int j = 0; j < simDomain->numComponents-1; j++)
    //     {
    //         for (int k = 0; k < simDomain->numComponents-1; k++)
    //         {
    //             if (simParams->diffusivity_host[i][j][k] == 0.0)
    //             {
    //                 simControls->antiTrapping = 1;
    //             }
    //         }
    //     }
    // }

    // Eigen-strain
    for (int i = 0; i < simDomain->numPhases; i++)
    {
        if (simParams->eigen_strain[i].xx == 0.0 && simParams->eigen_strain[i].yy == 0.0 && simParams->eigen_strain[i].zz == 0.0 && simParams->eigen_strain[i].xy == 0.0 && simParams->eigen_strain[i].xz == 0.0 && simParams->eigen_strain[i].yz == 0.0)
        {
            simControls->eigenSwitch[i] = 0;

        }
        else
        {
            simControls->eigenSwitch[i] = 1;
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
                    //printf("%d\t%d\t%d\t%d\t%le\n", i, j, k, l, simParams->Rotation_matrix_host[i][j][k][l]);
                    cudaMemcpy(&simParams->Rotation_matrix_dev[((i*simDomain->numPhases + j)*3 + k)*3 + l], &simParams->Rotation_matrix_host[i][j][k][l], sizeof(double),cudaMemcpyHostToDevice);
                    cudaMemcpy(&simParams->Inv_Rotation_matrix_dev[((i*simDomain->numPhases + j)*3 + k)*3 + l], &simParams->Inv_Rotation_matrix_host[i][j][k][l], sizeof(double),cudaMemcpyHostToDevice);
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
        blockSize->x = 2;
        blockSize->y = 128;
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
        blockSize->x = 4;
        blockSize->y = 8;
        blockSize->z = 8;

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
        printf("\nGrid size: (%ld, %ld, %ld)\nBlock size: (%ld, %ld, %ld)\n\n", (long)gridSize->x, (long)gridSize->y, (long)gridSize->z, (long)blockSize->x, (long)blockSize->y, (long)blockSize->z);
}
