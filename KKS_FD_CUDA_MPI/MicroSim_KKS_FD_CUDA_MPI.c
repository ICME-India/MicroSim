// Standard libraries
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Constant definitions
#define MASTER 0

#ifndef MAX_NUM_PHASES
#define MAX_NUM_PHASES 5
#endif

#ifndef MAX_NUM_COMP
#define MAX_NUM_COMP 5
#endif

#ifndef MAX_NUM_PHASE_COMP
#define MAX_NUM_PHASE_COMP 16
#endif

// C header files
// From ./functions/
#include "structures.h"
#include "inputReader.h"
#include "initialize_variables.h"
#include "filling.h"
#include "utilityFunctions.h"
// From ./solverloop/
#include "fileWriter.h"

// CUDA header files
// From ./functions
#include "functionF.cuh"
// From ./solverloop/
#include "calcPhaseComp.cuh"
#include "computeDrivingForce.cuh"
#include "updateComposition.cuh"
#include "updatePhi.cuh"
#include "utilityKernels.cuh"


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);                 // Initialise MPI
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;                         // MPI current proc. and num. of procs.
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    domainInfo      simDomain;              // Global mesh size and cell size
    controls        simControls;            // Timestep, num. of iterations, etc.
    simParameters   simParams;              // All the simulation parameters and details
    fillParameters  simFill;                // Parameters for filling the domain
    subdomainInfo   subdomain;              // MPI subdomain info

    int numDevices;
    cudaGetDeviceCount(&numDevices);

    cudaError_t error;

    // Phase-field variables on host-side
    double *phi;
    double *comp;

    // Phase-field variables on the GPUs
    // The blocks hold all the phases' and components' data contiguously,
    // phase-by-phase or component-by-component
    double **phiHost, **phiDev;
    double **compHost, **compDev;
    double **dfdphiHost, **dfdphiDev;
    double **dfdcHost, **dfdcDev;
    double **phiNewHost, **phiNewDev;
    double **compNewHost, **compNewDev;
    double **phaseCompHost, **phaseCompDev;
    double **muHost, **muDev;

    // for-loop iterators
    int i, j;

    // Kernel launch parameters
    dim3 gridSize, blockSize;

    double *maxerr;
    double *maxerr_h;
    double *maxVal, *minVal;
    double *maxVal_h, *minVal_h;

    // Run with at least two processes
    // Even if there is only 1 GPU, it will run.
    if (size < 2)
    {
        printf("\n\nNot enough processes. Try again\n\n");
        exit(0);
    }

    // Divide GPUs amongst the MPI processes
    cudaSetDevice(rank % numDevices);

    // Set up simulation using only one process
    // Using proc 0 here
    if (rank == MASTER)
    {
        // Create directory for output storage
        mkdir("DATA", 0777);
    }

    MPI_Barrier(comm);
    cudaDeviceSynchronize();

    // Read input from specified input file
    readInput_MPI(&simDomain, &simControls, &simParams, size, argv);
    MPI_Barrier(comm);

    calcFreeEnergyCoeffs(&simDomain, &simControls, &simParams);
    moveParamsToGPU(&simDomain, &simParams);

    // Fixed padding size for now
    simControls.padding = 2;

    // Create subdirectories for every processor's share of the output
    char directory[1000];
    sprintf(directory, "DATA/Processor_%d", rank);
    mkdir(directory, 0777);
    MPI_Barrier(comm);

    // Distribute data to all the processors
    // Only 1-D slab decomposition is available currently
    decomposeDomain(simDomain, &subdomain, rank, size);
    MPI_Barrier(comm);

    // Allocate memory on CPUs
    // Equivalent to gridinfo in the other MicroSim modules
    phi = (double*)malloc(sizeof(double) * simDomain.numPhases * subdomain.numCells);
    comp = (double*)malloc(sizeof(double) * (simDomain.numComponents-1) * subdomain.numCells);

    // Decide appropriate kernel grid size and block size
    // Subdomain padding is done here
    calcKernelParams(&gridSize, &blockSize, simDomain, simControls, &subdomain);
    MPI_Barrier(comm);
    cudaDeviceSynchronize();

    if (simControls.restart == 0 && simControls.startTime == 0)
    {
        // Read geometry from specified filling file
        readFill(&simFill, argv, rank);

        // Initialise the domain on the CPU using the read geometry
        fillDomain(simDomain, subdomain, simParams, phi, comp, &simFill);
    }

    if (simControls.restart != 0 || simControls.startTime != 0)
    {
        read_domain(phi, comp, simDomain, subdomain, simControls.startTime, rank, comm, argv);
    }

    MPI_Barrier(comm);


    // Allocate memory for values required to be calculated during solution,
    // like residuals and explicit calculations of derivatives

    // These pointers can only be dereferenced by device-side code
    // Allows dereferencing in kernels using [phase/component][location]
    cudaMalloc((void**)&compDev, sizeof(double*)*(simDomain.numComponents-1));
    cudaMalloc((void**)&dfdcDev, sizeof(double*)*(simDomain.numComponents-1));
    cudaMalloc((void**)&compNewDev, sizeof(double*)*(simDomain.numComponents-1));
    cudaMalloc((void**)&phiDev, sizeof(double*)*(simDomain.numPhases));
    cudaMalloc((void**)&dfdphiDev, sizeof(double*)*(simDomain.numPhases));
    cudaMalloc((void**)&phiNewDev, sizeof(double*)*(simDomain.numPhases));
    cudaMalloc((void**)&phaseCompDev, sizeof(double*)*simDomain.numPhases*(simDomain.numComponents-1));

    if (simControls.FUNCTION_F == 2)
        cudaMalloc((void**)&muDev, sizeof(double*)*(simDomain.numComponents-1));

    // These pointers can be dereferenced at the phase or component level on the host-side
    // Useful for data transfer and CUB compatibility
    compHost    = (double**)malloc((simDomain.numComponents-1)*sizeof(double*));
    dfdcHost    = (double**)malloc((simDomain.numComponents-1)*sizeof(double*));
    compNewHost = (double**)malloc((simDomain.numComponents-1)*sizeof(double*));
    phiHost     = (double**)malloc(simDomain.numPhases*sizeof(double*));
    dfdphiHost  = (double**)malloc(simDomain.numPhases*sizeof(double*));
    phiNewHost  = (double**)malloc(simDomain.numPhases*sizeof(double*));
    phaseCompHost = (double**)malloc(sizeof(double*)*simDomain.numPhases*(simDomain.numComponents-1));

    if (simControls.FUNCTION_F == 2)
        muHost = (double**)malloc((simDomain.numComponents-1)*sizeof(double*));

    // Memory on the device is allocated using the host-side pointers
    // The pointers are then copied to the device-side pointers so that they point to the same data
    for (j = 0; j < simDomain.numComponents-1; j++)
    {
        cudaMalloc((void**)&compHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&compDev[j], &compHost[j], sizeof(double*), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&dfdcHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&dfdcDev[j], &dfdcHost[j], sizeof(double*), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&compNewHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&compNewDev[j], &compNewHost[j], sizeof(double*), cudaMemcpyHostToDevice);

        if (simControls.FUNCTION_F == 2)
        {
            cudaMalloc((void**)&muHost[j], sizeof(double)*subdomain.numCompCells);
            cudaMemcpy(&muDev[j], &muHost[j], sizeof(double*), cudaMemcpyHostToDevice);
        }
    }

    for (j = 0; j < simDomain.numPhases; j++)
    {
        cudaMalloc((void**)&phiHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&phiDev[j], &phiHost[j], sizeof(double*), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&dfdphiHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&dfdphiDev[j], &dfdphiHost[j], sizeof(double*), cudaMemcpyHostToDevice);

        cudaMalloc((void**)&phiNewHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&phiNewDev[j], &phiNewHost[j], sizeof(double*), cudaMemcpyHostToDevice);
    }

    for (j = 0; j < simDomain.numPhases*(simDomain.numComponents-1); j++)
    {
        cudaMalloc((void**)&phaseCompHost[j], sizeof(double)*subdomain.numCompCells);
        cudaMemcpy(&phaseCompDev[j], &phaseCompHost[j], sizeof(double*), cudaMemcpyHostToDevice);
    }

    // Required for computing and printing statistics
    cudaMalloc((void**)&maxerr, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    cudaMalloc((void**)&maxVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    cudaMalloc((void**)&minVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));

    maxerr_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    maxVal_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    minVal_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));

    if (rank == MASTER)
        printf("\nAllocated memory on GPUs\n\n");

    MPI_Barrier(comm);
    cudaDeviceSynchronize();

    // Move fields from host to device
    for (i = 0; i < simDomain.numPhases; i++)
        cudaMemcpy(phiHost[i] + subdomain.shiftPointer, phi + i*subdomain.numCells, sizeof(double)*subdomain.numCells, cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    for (i = 0; i < simDomain.numComponents-1; i++)
        cudaMemcpy(compHost[i] + subdomain.shiftPointer, comp + i*subdomain.numCells, sizeof(double)*subdomain.numCells, cudaMemcpyHostToDevice);

    MPI_Barrier(comm);
    cudaDeviceSynchronize();

    // Copy old field to new field
    for (i = 0; i < simDomain.numPhases; i++)
        cudaMemcpy(phiNewHost[i], phiHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);

    for (i = 0; i < simDomain.numComponents-1; i++)
        cudaMemcpy(compNewHost[i], compHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);

    // Solution timeloop
    for (simControls.count = simControls.startTime; simControls.count <= simControls.startTime + simControls.numSteps; simControls.count++)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();

        /* Print statistics */
        if (simControls.count % simControls.trackProgress == 0)
        {
            if (rank == MASTER)
                printf("\nTime: %le\n", (double)(simControls.count)*simControls.DELTA_t);

            printStats(phiHost, compHost,
                       phiNewHost, compNewHost,
                       maxerr, maxVal, minVal,
                       simDomain, subdomain,
                       gridSize, blockSize);

            cudaMemcpy(maxerr_h, maxerr, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost);
            cudaMemcpy(maxVal_h, maxVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost);
            cudaMemcpy(minVal_h, minVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost);

            double ans1, ans2, ans3;

            for (i = 0; i < simDomain.numPhases; i++)
            {
                if (simControls.multiphase != 1 && i == 0 && simDomain.numPhases == 2 && simDomain.numComponents == 2)
                    continue;

                MPI_Reduce(maxerr_h+i+simDomain.numComponents-1, &ans1, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(maxVal_h+i+simDomain.numComponents-1, &ans2, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(minVal_h+i+simDomain.numComponents-1, &ans3, 1, MPI_DOUBLE, MPI_MIN, MASTER, comm);

                MPI_Barrier(comm);

                if (rank == MASTER)
                    printf("%*s, Max = %lf, Min = %lf, Relative_Change = %lf\n", 5, simDomain.phaseNames[i], ans2, ans3, ans1);
            }

            for (i = 0; i < simDomain.numComponents-1; i++)
            {
                MPI_Reduce(maxerr_h+i, &ans1, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(maxVal_h+i, &ans2, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(minVal_h+i, &ans3, 1, MPI_DOUBLE, MPI_MIN, MASTER, comm);

                MPI_Barrier(comm);

                if (rank == MASTER)
                    printf("%*s, Max = %lf, Min = %lf, Relative_Change = %lf\n", 5, simDomain.componentNames[i], ans2, ans3, ans1);
            }
        }

        // Copy new field to old field
        for (i = 0; i < simDomain.numPhases; i++)
            cudaMemcpy(phiHost[i], phiNewHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);

        for (i = 0; i < simDomain.numComponents-1; i++)
            cudaMemcpy(compHost[i], compNewHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);


        //Moving buffers to neighbours
        /*
         *  If CUDA-aware MPI is not found, the data will be staged through the host
         *  This is, of course, highly inefficient. Installing CUDA-aware MPI with UCX and GDRcopy is highly recommended.
         *  Instructions to get it running can be found in this module's manual, and on the OpenMPI website.
         */
        for (i = 0; i < simDomain.numPhases; i++)
        {
            MPI_Sendrecv(phiHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 0,
                         phiHost[i]+subdomain.shiftPointer+subdomain.numCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 0,
                         comm, MPI_STATUS_IGNORE);

            MPI_Sendrecv(phiHost[i]+subdomain.numCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 1,
                         phiHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 1,
                         comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(comm);
        cudaDeviceSynchronize();

        for (i = 0; i < simDomain.numComponents-1; i++)
        {
            MPI_Sendrecv(compHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                         compHost[i]+subdomain.shiftPointer+subdomain.numCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                         comm, MPI_STATUS_IGNORE);

            MPI_Sendrecv(compHost[i]+subdomain.numCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                         compHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                         comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(comm);
        cudaDeviceSynchronize();

        // Writing to file
        if (simControls.count % simControls.saveInterval == 0 || simControls.count == simControls.startTime + simControls.numSteps)
        {
            for (i = 0; i < simDomain.numPhases; i++)
                cudaMemcpy(phi + i*subdomain.numCells, phiHost[i] + subdomain.shiftPointer, sizeof(double)*subdomain.numCells, cudaMemcpyDeviceToHost);
            for (i = 0; i < simDomain.numComponents-1; i++)
                cudaMemcpy(comp + i*subdomain.numCells, compHost[i] + subdomain.shiftPointer, sizeof(double)*subdomain.numCells, cudaMemcpyDeviceToHost);

            if (simControls.writeFormat == 0)
                writeVTK_BINARY(phi, comp, simDomain, subdomain, simControls.count, rank, comm, argv);
            else
                writeVTK_ASCII(phi, comp, simDomain, subdomain, simControls.count, rank, comm, argv);


            if (rank == MASTER)
                printf("Wrote to file\n");

            if (simControls.count == simControls.startTime + simControls.numSteps)
                break;
        }

        calcPhaseComp(phiDev, compDev,
                      phaseCompDev, muDev,
                      &simDomain, &simControls,
                      &simParams, &subdomain,
                      gridSize, blockSize);

        computeDrivingForce(phiDev, compDev,
                            dfdphiDev, muDev, phaseCompDev,
                            &simDomain, &simControls,
                            &simParams, &subdomain,
                            gridSize, blockSize);

        updatePhi(phiDev, dfdphiDev, phiNewDev,
                  &simDomain, &simControls,
                  &simParams, &subdomain,
                  gridSize, blockSize);

        updateComposition(phiDev, compDev, phiNewDev, compNewDev,
                          phaseCompDev, muDev,
                          &simDomain, &simControls,
                          &simParams, &subdomain,
                          gridSize, blockSize);

        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    } //End solution timeloop

    free(minVal_h);
    free(maxVal_h);
    free(maxerr_h);

    cudaFree(minVal);
    cudaFree(maxVal);
    cudaFree(maxerr);

    for (j = 0; j < simDomain.numComponents-1; j++)
    {
        cudaFree(compHost[j]);
        cudaFree(dfdcHost[j]);
        cudaFree(compNewHost[j]);

        if (simControls.FUNCTION_F == 2)
            cudaFree(muHost[j]);
    }

    cudaFree(compDev);
    cudaFree(dfdcDev);
    cudaFree(compNewDev);

    if (simControls.FUNCTION_F == 2)
        cudaFree(muDev);

    cudaFree(phiDev);
    cudaFree(dfdphiDev);
    cudaFree(phiNewDev);

    cudaFree(phaseCompDev);

    free(compHost);
    free(dfdcHost);
    free(compNewHost);

    if (simControls.FUNCTION_F == 2)
        free(muHost);

    free(phiHost);
    free(dfdphiHost);
    free(phiNewHost);

    free(phaseCompHost);

    freeVars(&simDomain, &simParams);

    free(phi);
    free(comp);

    MPI_Finalize();

    return 0;
}
