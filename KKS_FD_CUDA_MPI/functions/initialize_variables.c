#include "initialize_variables.h"

int LUPDecompose(double **A, int N, double Tol, int *P) {

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done
}

void LUPInvert(double **A, int *P, int N, double **IA) {

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

void matrixMultiply(double **A, double **B, double **C, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = 0.0;

            for (int k = 0; k < N; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

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
    double Tol = 1e-6;
    int P[simDomain->numComponents];

    double ***dmudc = malloc3M(simDomain->numPhases, simDomain->numComponents-1, simDomain->numComponents-1);
    double **inverted = malloc2M(simDomain->numComponents-1, simDomain->numComponents-1);

    for (int i = 0; i < simDomain->numPhases; i++)
    {
        for (int j = 0; j < simDomain->numComponents-1; j++)
        {
            for (int k = 0; k < simDomain->numComponents-1; k++)
            {
                if (j == k)
                    dmudc[i][j][k] = 2.0*simParams->F0_A_host[i][j][k];
                else
                    dmudc[i][j][k] = simParams->F0_A_host[i][j][k];
            }
        }
    }

    for (int i = 0; i < simDomain->numPhases; i++)
    {
        LUPDecompose(dmudc[i], simDomain->numComponents-1, Tol, P);
        LUPInvert(dmudc[i], P, simDomain->numComponents-1, inverted);

        matrixMultiply(simParams->diffusivity_host[i], inverted, simParams->mobility_host[i], simDomain->numComponents-1);
    }

    free2M(inverted, simDomain->numComponents-1);
    free3M(dmudc, simDomain->numPhases, simDomain->numComponents-1);

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

                //printf("%le\t%d\t%d\n", simParams->relax_coeff_host[i][j], i, j);
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
                cudaMemcpy(&simParams->mobility_dev[(i*(simDomain->numComponents-1) + j)*(simDomain->numComponents-1) + k], &simParams->mobility_host[i][j][k], sizeof(double), cudaMemcpyHostToDevice);
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
