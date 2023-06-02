// Standard libraries
#include <numeric>
#include <vector>
#include <complex>
#include <random>
#include <cstdlib>
#include <cstdio>

#include <mpi.h>
#include <mpi-ext.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cuda.h>
#include <cuda_runtime.h>

// Flags and constants
#define MASTER 0

#ifndef ENABLE_CUFFTMP
#define ENABLE_CUFFTMP 0
#endif

#ifndef ENABLE_HDF5
#define ENABLE_HDF5 0
#endif

// cuFFTMp
#if ENABLE_CUFFTMP == 1
//#include "box_iterator.hpp"
#include "error_checks.hpp"
#include <cufft.h>
#include <cufftMp.h>
#include "calc_bn.h"
#endif


// C header files
// From ./functions/
#include "structures.h"
#include "inputReader.h"
#include "initialize_variables.h"
#include "filling.h"
#include "utilityFunctions.h"
#include "helper_string.h"
#include "helper_cuda.h"
// From ./solverloop/
#include "fileWriter.h"

// CUDA header files
// From ./functions
#include "functionF.cuh"
// From ./solverloop/
#include "boundary.cuh"
#include "smooth.cuh"
#include "calcPhaseComp.cuh"
#include "computeDrivingForce.cuh"
#include "updateComposition.cuh"
#include "updatePhi.cuh"
#include "computeElastic.cuh"
#include "utilityKernels.cuh"


int main(int argc, char *argv[])
{
    cudaDeviceReset();

    MPI_Init(&argc, &argv);                 // Initialise MPI
    MPI_Comm comm = MPI_COMM_WORLD;

    int MPI_Enabled = 1;

    int rank, size;                         // MPI current proc. and num. of procs.
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (size == 1)
        MPI_Enabled = 0;

    if (MPI_Enabled && rank == MASTER)
        printf("MPI enabled\n");

    domainInfo      simDomain;              // Global mesh size and cell size
    controls        simControls;            // Timestep, num. of iterations, etc.
    simParameters   simParams;              // All the simulation parameters and details
    fillParameters  simFill;                // Parameters for filling the domain
    subdomainInfo   subdomain;              // MPI subdomain info

    int numDevices;                         // Number of visible GPUs in the node
    cudaGetDeviceCount(&numDevices);

    cudaError_t error;

    // Phase-field variables on host-side
    double *phi;
    double *comp;
    double *mu;

    // Phase-field variables on the GPUs
    // The blocks hold all the phases' and components' data contiguously,
    // phase-by-phase or component-by-component
    double **phiHost, **phiDev;
    double **compHost, **compDev;
    double **dfdphiHost, **dfdphiDev;
    double **phiNewHost, **phiNewDev;
    double **compNewHost, **compNewDev;
    double **phaseCompHost, **phaseCompDev;
    double **muHost, **muDev;

    // for-loop iterators
    long i, j;

    // Kernel launch parameters
    dim3 gridSize, blockSize;

    // Variables for finding and storing statistics
    double *maxerr;
    double *maxerr_h;
    double *maxVal, *minVal;
    double *maxVal_h, *minVal_h;

    // Divide GPUs amongst the MPI processes
    cudaSetDevice(rank % numDevices);

    // Set up simulation using only one process
    // Using proc 0 here
    if (rank == MASTER)
    {
        // Create directory for output storage
        mkdir("DATA", 0777);
    }

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    // Read input from specified input file
    if (readInput_MPI(&simDomain, &simControls, &simParams, rank, argv))
    {
        if (rank == MASTER)
        {
            printf("\n----------------------------------------------------------------------------ERROR---------------------------------------------------------------------------\n");
            printf("\nSolver not compiled to handle required number of phases and/or components.\n");
            printf("Please select the appropriate number for both, using the MicroSim GUI or by running GEdata_writer.py with the desired input file through the command line\n");
            if (MAX_NUM_PHASES < simDomain.numPhases)
                printf("Currently, the maximum number of phases supported is %ld, but attempted to run with %d phases\n", MAX_NUM_PHASES, simDomain.numPhases);
            if (MAX_NUM_COMP < simDomain.numComponents)
                printf("Currently, the maximum number of components supported is %ld, but you have attempted to run with %d components\n", MAX_NUM_COMP, simDomain.numComponents);
            printf("\n------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
        }

        cudaDeviceReset();
        MPI_Abort(comm, 1);
    }

    read_boundary_conditions(&simDomain, &simControls, &simParams, rank, argv);

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
    }

    // Fix padding size
    simControls.padding = 1;
    simControls.numworkers = size;

    // Initialize free energy, simulation parameters
    calcFreeEnergyCoeffs(&simDomain, &simControls, &simParams);
    moveParamsToGPU(&simDomain, &simControls, &simParams);

    // Create subdirectories for every processor's share of the output
    char directory[1000];
    sprintf(directory, "DATA/Processor_%d", rank);
    mkdir(directory, 0777);

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
    }

    // Distribute data to all the processors
    // Only 1-D slab decomposition is available currently
    decomposeDomain(simDomain, &simControls, &subdomain, rank, size);

    // Decide appropriate kernel grid size and block size
    calcKernelParams(&gridSize, &blockSize, simDomain, simControls, &subdomain);

    // Allocate memory on CPUs
    // Equivalent to gridinfo in the other MicroSim modules
    phi = (double*)malloc(sizeof(double) * simDomain.numPhases * subdomain.numCompCells);
    comp = (double*)malloc(sizeof(double) * (simDomain.numComponents-1) * subdomain.numCompCells);
    if (simControls.FUNCTION_F == 2)
        mu = (double*)malloc(sizeof(double) * (simDomain.numComponents-1) * subdomain.numCompCells);

    for (i = 0; i < simDomain.numPhases*subdomain.numCompCells; i++)
        phi[i] = 0.0;

    for (i = 0; i < (simDomain.numComponents-1)*subdomain.numCompCells; i++)
        comp[i] = 0.0;

    if (simControls.FUNCTION_F == 2)
        for (i = 0; i < (simDomain.numComponents-1)*subdomain.numCompCells; i++)
            mu[i] = 0.0;

    #if ENABLE_CUFFTMP == 0
    simControls.ELASTICITY = 0;
    #elif ENABLE_CUFFTMP == 1

    int FFT_Alternate= 0;

    double **BpqDev, **BpqHost, *B_calc, *kx, *ky, *kz;

    const size_t my_nx = subdomain.numX;
    const size_t my_ny = subdomain.numY;
    const size_t my_nz = subdomain.numZ;

    printf("cuFFTMp is enabled\n");
    printf("For processor %d: %d, %d, %d\n", rank, (int)my_nx, (int)my_ny, (int)my_nz);

    // Create a plan
    cufftHandle plan = 0;

    // Attach to a stream
    cudaStream_t stream = nullptr;
    CUDA_CHECK(cudaStreamCreate(&stream));

    if (simControls.ELASTICITY)
    {
        CUFFT_CHECK(cufftCreate(&plan));
        CUFFT_CHECK(cufftMpAttachComm(plan, CUFFT_COMM_MPI, &comm));
        CUFFT_CHECK(cufftSetStream(plan, stream));
    }

    // Make the plan
    size_t workspace;
    if (simControls.ELASTICITY)
        CUFFT_CHECK(cufftMakePlan3d(plan, simDomain.MESH_X, simDomain.MESH_Y, simDomain.MESH_Z, CUFFT_Z2Z, &workspace));

    // cuFFTMp variables
    cudaLibXtDesc *phiElDev[simDomain.numPhases], *dfEldphi[simDomain.numPhases], *tempFFT_Shuffled;

    if (simControls.ELASTICITY)
    {
        for (i = 0; i < simDomain.numPhases; i++)
        {
            CUFFT_CHECK(cufftXtMalloc(plan, &phiElDev[i], CUFFT_XT_FORMAT_INPLACE));
            CUFFT_CHECK(cufftXtMalloc(plan, &dfEldphi[i], CUFFT_XT_FORMAT_INPLACE));
        }

        CUFFT_CHECK(cufftXtMalloc(plan, &tempFFT_Shuffled, CUFFT_XT_FORMAT_INPLACE_SHUFFLED));
    }

    if (simControls.ELASTICITY)
    {
        cudaMalloc((void**)&BpqDev, sizeof(double*)*simDomain.numPhases*simDomain.numPhases);
        BpqHost = (double**)malloc(sizeof(double*)*simDomain.numPhases*simDomain.numPhases);

        for (i = 0; i < simDomain.numPhases*simDomain.numPhases; i++)
        {
            checkCudaErrors(cudaMalloc((void**)&BpqHost[i], sizeof(double)*my_nx*my_ny*my_nz));
            checkCudaErrors(cudaMemcpy(&BpqDev[i], &BpqHost[i], sizeof(double*), cudaMemcpyHostToDevice));
        }

        B_calc = (double*)malloc(sizeof(double)*my_nx*my_ny*my_nz);
        kx     = (double*)malloc(sizeof(double)*my_nx*my_ny*my_nz);
        ky     = (double*)malloc(sizeof(double)*my_nx*my_ny*my_nz);
        kz     = (double*)malloc(sizeof(double)*my_nx*my_ny*my_nz);

        calc_k(kx, ky, kz, tempFFT_Shuffled,
               my_nx, my_ny, my_nz,
               simDomain, simControls, simParams, subdomain,
               stream);

        for (i = 0; i < simDomain.numPhases; i++)
        {
            for (j = 0; j < simDomain.numPhases; j++)
            {
                calculate_Bn(B_calc, kx, ky, kz, simDomain, simParams, subdomain, i, j);
                cudaMemcpy(BpqHost[i*simDomain.numPhases + j], B_calc, sizeof(double)*my_nx*my_ny*my_nz, cudaMemcpyHostToDevice);
            }
        }

        CUFFT_CHECK(cufftXtFree(tempFFT_Shuffled));
        free(B_calc);
        free(kx);
        free(ky);
        free(kz);
    }
    #endif

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    // Not restarting
    if (simControls.restart == 0 && simControls.startTime == 0)
    {
        // Read geometry from specified filling file
        readFill(&simFill, argv, rank);

        // Initialise the domain on the CPU using the read geometry
        fillDomain(simDomain, subdomain, simParams, phi, comp, &simFill);

        if (rank == MASTER)
            printf("Finished filling\n");
    }

    // If restarting
    if (!(simControls.restart == 0) || !(simControls.startTime == 0))
    {
        if (rank == MASTER)
            printf("Reading from disk\n");

        #if ENABLE_HDF5 == 1
        if (simControls.writeHDF5)
        {
            readHDF5(phi, comp, mu,
                     simDomain, subdomain,
                     simControls, rank, comm, argv);
        }
        else
            #endif
        {
            read_domain(phi, comp, mu, simDomain, subdomain, simControls, rank, comm, argv);
        }
    }

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    if (rank == MASTER)
        printf("\nAllocating memory on GPUs\n");

    // Allocate memory for values required to be calculated during solution,
    // like residuals and explicit calculations of derivatives

    // These pointers can only be dereferenced by device-side code
    // Allows dereferencing in kernels using [phase/component][location]
    checkCudaErrors(cudaMalloc((void**)&compDev, sizeof(double*)*(simDomain.numComponents-1)));
    checkCudaErrors(cudaMalloc((void**)&compNewDev, sizeof(double*)*(simDomain.numComponents-1)));
    checkCudaErrors(cudaMalloc((void**)&phiDev, sizeof(double*)*(simDomain.numPhases)));
    checkCudaErrors(cudaMalloc((void**)&dfdphiDev, sizeof(double*)*(simDomain.numPhases)));
    checkCudaErrors(cudaMalloc((void**)&phiNewDev, sizeof(double*)*(simDomain.numPhases)));
    checkCudaErrors(cudaMalloc((void**)&phaseCompDev, sizeof(double*)*simDomain.numPhases*(simDomain.numComponents-1)));

    if (simControls.FUNCTION_F == 2)
        cudaMalloc((void**)&muDev, sizeof(double*)*(simDomain.numComponents-1));

    // These pointers can be dereferenced at the phase or component level on the host-side
    // Useful for data transfer and CUB compatibility
    compHost    = (double**)malloc((simDomain.numComponents-1)*sizeof(double*));
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
        checkCudaErrors(cudaMalloc((void**)&compHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&compDev[j], &compHost[j], sizeof(double*), cudaMemcpyHostToDevice));

        checkCudaErrors(cudaMalloc((void**)&compNewHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&compNewDev[j], &compNewHost[j], sizeof(double*), cudaMemcpyHostToDevice));

        if (simControls.FUNCTION_F == 2)
        {
            checkCudaErrors(cudaMalloc((void**)&muHost[j], sizeof(double)*subdomain.numCompCells));
            checkCudaErrors(cudaMemcpy(&muDev[j], &muHost[j], sizeof(double*), cudaMemcpyHostToDevice));
        }
    }

    for (j = 0; j < simDomain.numPhases; j++)
    {
        checkCudaErrors(cudaMalloc((void**)&phiHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&phiDev[j], &phiHost[j], sizeof(double*), cudaMemcpyHostToDevice));

        checkCudaErrors(cudaMalloc((void**)&dfdphiHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&dfdphiDev[j], &dfdphiHost[j], sizeof(double*), cudaMemcpyHostToDevice));

        checkCudaErrors(cudaMalloc((void**)&phiNewHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&phiNewDev[j], &phiNewHost[j], sizeof(double*), cudaMemcpyHostToDevice));
    }

    for (j = 0; j < simDomain.numPhases*(simDomain.numComponents-1); j++)
    {
        checkCudaErrors(cudaMalloc((void**)&phaseCompHost[j], sizeof(double)*subdomain.numCompCells));
        checkCudaErrors(cudaMemcpy(&phaseCompDev[j], &phaseCompHost[j], sizeof(double*), cudaMemcpyHostToDevice));
    }

    // Required for computing and printing statistics
    checkCudaErrors(cudaMalloc((void**)&maxerr, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1)));
    checkCudaErrors(cudaMalloc((void**)&maxVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1)));
    checkCudaErrors(cudaMalloc((void**)&minVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1)));

    maxerr_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    maxVal_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));
    minVal_h = (double*)malloc(sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1));

    if (rank == MASTER)
        printf("Allocated memory on GPUs\n");

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    // Move fields from host to device
    for (i = 0; i < simDomain.numPhases; i++)
        checkCudaErrors(cudaMemcpy(phiHost[i], phi + i*subdomain.numCompCells, sizeof(double)*subdomain.numCompCells, cudaMemcpyHostToDevice));

    for (i = 0; i < simDomain.numComponents-1; i++)
        checkCudaErrors(cudaMemcpy(compHost[i], comp + i*subdomain.numCompCells, sizeof(double)*subdomain.numCompCells, cudaMemcpyHostToDevice));

    if (simControls.FUNCTION_F == 2)
        for (i = 0; i < simDomain.numComponents-1; i++)
            checkCudaErrors(cudaMemcpy(muHost[i], mu + i*subdomain.numCompCells, sizeof(double)*subdomain.numCompCells, cudaMemcpyHostToDevice));

    //Moving buffers to neighbours
    /*
     *  If CUDA-aware MPI is not found, the data will be staged through the host
     *  This is, of course, highly inefficient. Installing CUDA-aware MPI with UCX and GDRcopy is highly recommended.
     *  Instructions to get it running can be found in this module's manual, and on the OpenMPI website.
     *
     */

    if (MPI_Enabled)
    {
        for (i = 0; i < simDomain.numPhases; i++)
        {
            MPI_Sendrecv(phiHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 0,
                         phiHost[i]+subdomain.numCompCells-subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 0,
                         comm, MPI_STATUS_IGNORE);

            MPI_Sendrecv(phiHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 1,
                         phiHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 1,
                         comm, MPI_STATUS_IGNORE);
        }

        for (i = 0; i < simDomain.numComponents-1; i++)
        {
            MPI_Sendrecv(compHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                         compHost[i]+subdomain.numCompCells-subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                         comm, MPI_STATUS_IGNORE);

            MPI_Sendrecv(compHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                         compHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                         comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    applyBoundaryCondition(phiDev, 0, simDomain.numPhases,
                           simDomain, simControls,
                           simParams, subdomain,
                           gridSize, blockSize);

    applyBoundaryCondition(compDev, 2, simDomain.numComponents-1,
                           simDomain, simControls,
                           simParams, subdomain,
                           gridSize, blockSize);


    // Copy old field to new field
    for (i = 0; i < simDomain.numPhases; i++)
        checkCudaErrors(cudaMemcpy(phiNewHost[i], phiHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice));

    for (i = 0; i < simDomain.numComponents-1; i++)
        checkCudaErrors(cudaMemcpy(compNewHost[i], compHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice));

    if (MPI_Enabled)
    {
        MPI_Barrier(comm);
        cudaDeviceSynchronize();
    }

    // Smooth
    for (j = 1; j <= simControls.nsmooth; j++)
    {
        smooth(phiDev, phiNewDev,
               &simDomain, &simControls,
               &simParams, &subdomain,
               gridSize, blockSize);

        for (i = 0; i < simDomain.numPhases; i++)
        {
            checkCudaErrors(cudaMemcpy(phiHost[i], phiNewHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice));
        }

        if (MPI_Enabled)
        {
            for (i = 0; i < simDomain.numPhases; i++)
            {
                MPI_Sendrecv(phiHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 0,
                             phiHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 0,
                             comm, MPI_STATUS_IGNORE);

                MPI_Sendrecv(phiHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 1,
                             phiHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 1,
                             comm, MPI_STATUS_IGNORE);
            }
        }

        applyBoundaryCondition(phiDev, 0, simDomain.numPhases,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        calcPhaseComp(phiDev, compDev,
                      phaseCompDev, muDev,
                      &simDomain, &simControls,
                      &simParams, &subdomain,
                      gridSize, blockSize);

        if (j == simControls.nsmooth)
        {
            for (i = 0; i < simDomain.numPhases; i++)
            {
                checkCudaErrors(cudaMemcpy(phiNewHost[i], phiHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice));
            }

            if (rank == MASTER)
                printf("\nFinished smoothing\n");
        }
    }

    if (rank == MASTER)
        printf("\nStarting solution procedure\n");

    // start2 = MPI_Wtime();

    // Solution timeloop
    for (simControls.count = simControls.startTime; simControls.count <= simControls.startTime + simControls.numSteps; simControls.count++)
    {
        /*
         *
         * Print max, min, change
         *
         */


        if (simControls.count % simControls.trackProgress == 0)
        {
            if (rank == MASTER)
                printf("\nTime: %le\n", (double)(simControls.count)*simControls.DELTA_t);

            printStats(phiHost, compHost,
                       phiNewHost, compNewHost,
                       maxerr, maxVal, minVal,
                       simDomain, subdomain,
                       gridSize, blockSize);

            checkCudaErrors(cudaMemcpy(maxerr_h, maxerr, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(maxVal_h, maxVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(minVal_h, minVal, sizeof(double)*(simDomain.numPhases+simDomain.numComponents-1), cudaMemcpyDeviceToHost));

            double ans1, ans2, ans3;

            for (i = 0; i < simDomain.numPhases; i++)
            {
                MPI_Reduce(maxerr_h+i+simDomain.numComponents-1, &ans1, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(maxVal_h+i+simDomain.numComponents-1, &ans2, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(minVal_h+i+simDomain.numComponents-1, &ans3, 1, MPI_DOUBLE, MPI_MIN, MASTER, comm);

                if (rank == MASTER)
                    printf("%*s, Max = %le, Min = %le, Relative_Change = %le\n", 5, simDomain.phaseNames[i], ans2, ans3, ans1);

                if (fabs(ans1) > 2 && rank == MASTER)
                {
                    cudaDeviceReset();
                    MPI_Abort(comm, 1);
                }
            }

            for (i = 0; i < simDomain.numComponents-1; i++)
            {
                MPI_Reduce(maxerr_h+i, &ans1, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(maxVal_h+i, &ans2, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
                MPI_Reduce(minVal_h+i, &ans3, 1, MPI_DOUBLE, MPI_MIN, MASTER, comm);

                if (rank == MASTER)
                    printf("%*s, Max = %le, Min = %le, Relative_Change = %le\n", 5, simDomain.componentNames[i], ans2, ans3, ans1);

                if (fabs(ans1) > 2 && rank == MASTER)
                {
                    cudaDeviceReset();
                    MPI_Abort(comm, 1);
                }
            }
        }

        // Copy new field to old field
        for (i = 0; i < simDomain.numPhases; i++)
            cudaMemcpy(phiHost[i], phiNewHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);

        for (i = 0; i < simDomain.numComponents-1; i++)
            cudaMemcpy(compHost[i], compNewHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToDevice);

        /*
         *
         * Writing to file
         * Write mu only if FUNCTION_F = 2 (Exact)
         *
         */
        if (simControls.count % simControls.saveInterval == 0 || simControls.count == simControls.startTime + simControls.numSteps)
        {
            for (i = 0; i < simDomain.numPhases; i++)
                checkCudaErrors(cudaMemcpy(phi + i*subdomain.numCompCells, phiHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToHost));
            for (i = 0; i < simDomain.numComponents-1; i++)
                checkCudaErrors(cudaMemcpy(comp + i*subdomain.numCompCells, compHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToHost));
            if (simControls.FUNCTION_F == 2)
                for (i = 0; i < simDomain.numComponents-1; i++)
                    checkCudaErrors(cudaMemcpy(mu + i*subdomain.numCompCells, muHost[i], sizeof(double)*subdomain.numCompCells, cudaMemcpyDeviceToHost));

            if (simControls.writeFormat == 0)
                writeVTK_BINARY(phi, comp, mu, simDomain, subdomain, simControls, rank, comm, argv);
            else if (simControls.writeFormat == 1)
                writeVTK_ASCII(phi, comp, mu, simDomain, subdomain, simControls, rank, comm, argv);

            #if ENABLE_HDF5 == 1
            if (simControls.writeHDF5)
                writeHDF5(phi, comp, mu,
                          simDomain, subdomain,
                          simControls, rank, comm, argv);
                if (rank == MASTER && simControls.writeHDF5)
                    printf("Wrote to file\n");
            #endif

            if (rank == MASTER && simControls.writeHDF5 == 0)
                printf("Wrote to file\n");

            if (simControls.count == simControls.startTime + simControls.numSteps)
                break;
        }

        //Moving buffers to neighbours
        /*
         *  If CUDA-aware MPI is not found, the data will have to be staged through the host
         *  This is highly inefficient. Hence, installing CUDA-aware MPI with UCX and GDRcopy is preferable.
         *  Instructions to get it running can be found in this module's manual, and on the OpenMPI website.
         *
         */

        if (MPI_Enabled)
        {
            for (i = 0; i < simDomain.numPhases; i++)
            {
                MPI_Sendrecv(phiHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 0,
                             phiHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 0,
                             comm, MPI_STATUS_IGNORE);

                MPI_Sendrecv(phiHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 1,
                             phiHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 1,
                             comm, MPI_STATUS_IGNORE);
            }

            for (i = 0; i < simDomain.numComponents-1; i++)
            {

                MPI_Sendrecv(compHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                             compHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                             comm, MPI_STATUS_IGNORE);


                MPI_Sendrecv(compHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                             compHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                             comm, MPI_STATUS_IGNORE);

            }

            if (simControls.FUNCTION_F == 2)
            {
                for (i = 0; i < simDomain.numComponents-1; i++)
                {
                    MPI_Sendrecv(muHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                                 muHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                                 comm, MPI_STATUS_IGNORE);

                    MPI_Sendrecv(muHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                                 muHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                                 comm, MPI_STATUS_IGNORE);
                }
            }
        }

        applyBoundaryCondition(phiDev, 0, simDomain.numPhases,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        applyBoundaryCondition(compDev, 2, simDomain.numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        if (simControls.FUNCTION_F == 2)
        {
            applyBoundaryCondition(muDev, 1, simDomain.numComponents-1,
                                   simDomain, simControls,
                                   simParams, subdomain,
                                   gridSize, blockSize);
        }

        /*
         *
         *  Solution procedure
         *  Kernel calls are wrapped using host-side functions
         *
         */

        calcPhaseComp(phiDev, compDev,
                      phaseCompDev, muDev,
                      &simDomain, &simControls,
                      &simParams, &subdomain,
                      gridSize, blockSize);

        resetArray(dfdphiDev, simDomain.numPhases, simDomain.DIMENSION, subdomain.sizeX, subdomain.sizeY, subdomain.sizeZ, gridSize, blockSize);

        #if ENABLE_CUFFTMP == 1
        if (simControls.ELASTICITY)
        {
            if (FFT_Alternate == 1)
            {
                moveTocudaLibXtDesc(phiDev, phiElDev,
                                    simDomain, simControls, simParams, subdomain,
                                    stream, gridSize, blockSize, comm);

                for (i = 0; i < simDomain.numPhases; i++)
                {
                    if (simControls.eigenSwitch[i] == 1)
                    {
                        CUFFT_CHECK(cufftXtExecDescriptor(plan, phiElDev[i], phiElDev[i], CUFFT_FORWARD));
                        //CUFFT_CHECK(cufftXtExecDescriptor(plan, dfEldphi[i], dfEldphi[i], CUFFT_FORWARD));
                    }
                }

                computeDrivingForce_Elastic(phiElDev, dfEldphi, BpqHost,
                                            simDomain, simControls, simParams, subdomain,
                                            stream, comm);

                for (i = 0; i < simDomain.numPhases; i++)
                {
                    if (simControls.eigenSwitch[i] == 1)
                    {
                        CUFFT_CHECK(cufftXtExecDescriptor(plan, dfEldphi[i], dfEldphi[i], CUFFT_INVERSE));
                        //CUFFT_CHECK(cufftXtExecDescriptor(plan, phiElDev[i], phiElDev[i], CUFFT_INVERSE));
                    }
                }

                moveFromcudaLibXtDesc(dfdphiDev, dfEldphi,
                                      simDomain, simControls, simParams, subdomain,
                                      stream, gridSize, blockSize);
                FFT_Alternate = 0;
            }
            else if (FFT_Alternate == 0)
            {
                moveTocudaLibXtDesc(phiDev, dfEldphi,
                                    simDomain, simControls, simParams, subdomain,
                                    stream, gridSize, blockSize, comm);

                for (i = 0; i < simDomain.numPhases; i++)
                {
                    if (simControls.eigenSwitch[i] == 1)
                    {
                        //CUFFT_CHECK(cufftXtExecDescriptor(plan, phiElDev[i], phiElDev[i], CUFFT_FORWARD));
                        CUFFT_CHECK(cufftXtExecDescriptor(plan, dfEldphi[i], dfEldphi[i], CUFFT_FORWARD));
                    }
                }

                computeDrivingForce_Elastic(dfEldphi, phiElDev, BpqHost,
                                            simDomain, simControls, simParams, subdomain,
                                            stream, comm);

                for (i = 0; i < simDomain.numPhases; i++)
                {
                    if (simControls.eigenSwitch[i] == 1)
                    {
                        //CUFFT_CHECK(cufftXtExecDescriptor(plan, dfEldphi[i], dfEldphi[i], CUFFT_INVERSE));
                        CUFFT_CHECK(cufftXtExecDescriptor(plan, phiElDev[i], phiElDev[i], CUFFT_INVERSE));
                    }
                }

                moveFromcudaLibXtDesc(dfdphiDev, phiElDev,
                                      simDomain, simControls, simParams, subdomain,
                                      stream, gridSize, blockSize);

                FFT_Alternate = 1;
            }
        }
        #endif

        computeDrivingForce_Chemical(phiDev, compDev,
                                     dfdphiDev, phaseCompDev, muDev,
                                     &simDomain, &simControls,
                                     &simParams, &subdomain,
                                     gridSize, blockSize);

        updatePhi(phiDev, dfdphiDev, phiNewDev, phaseCompDev,
                  &simDomain, &simControls,
                  &simParams, &subdomain,
                  gridSize, blockSize);

        if (MPI_Enabled)
        {
            for (i = 0; i < simDomain.numPhases; i++)
            {
                MPI_Sendrecv(phiNewHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 0,
                             phiNewHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 0,
                             comm, MPI_STATUS_IGNORE);

                MPI_Sendrecv(phiNewHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 1,
                             phiNewHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 1,
                             comm, MPI_STATUS_IGNORE);
            }
        }

        applyBoundaryCondition(phiNewDev, 0, simDomain.numPhases,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        updateComposition(phiDev, compDev, phiNewDev, compNewDev,
                          phaseCompDev, muDev,
                          simDomain, simControls,
                          simParams, subdomain,
                          gridSize, blockSize);

        if (MPI_Enabled)
        {
            for (i = 0; i < simDomain.numComponents-1; i++)
            {

                MPI_Sendrecv(compNewHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                             compNewHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                             comm, MPI_STATUS_IGNORE);


                MPI_Sendrecv(compNewHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                             compNewHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                             comm, MPI_STATUS_IGNORE);

            }

            if (simControls.FUNCTION_F == 2)
            {
                for (i = 0; i < simDomain.numComponents-1; i++)
                {
                    MPI_Sendrecv(muHost[i]+subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 2,
                                 muHost[i]-subdomain.shiftPointer+subdomain.numCompCells, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 2,
                                 comm, MPI_STATUS_IGNORE);

                    MPI_Sendrecv(muHost[i]+subdomain.numCompCells-2*subdomain.shiftPointer, subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbFront, 3,
                                 muHost[i], subdomain.shiftPointer, MPI_DOUBLE, subdomain.nbBack, 3,
                                 comm, MPI_STATUS_IGNORE);
                }
            }
        }

        applyBoundaryCondition(compNewDev, 2, simDomain.numComponents-1,
                               simDomain, simControls,
                               simParams, subdomain,
                               gridSize, blockSize);

        if (simControls.FUNCTION_F == 2)
        {
            applyBoundaryCondition(muDev, 1, simDomain.numComponents-1,
                                   simDomain, simControls,
                                   simParams, subdomain,
                                   gridSize, blockSize);
        }

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
        cudaFree(compNewHost[j]);

        if (simControls.FUNCTION_F == 2)
            cudaFree(muHost[j]);
    }

    cudaFree(compDev);
    cudaFree(compNewDev);

    if (simControls.FUNCTION_F == 2)
        cudaFree(muDev);

    for (j = 0; j < simDomain.numPhases; j++)
    {
        cudaFree(phiHost[j]);
        cudaFree(phiNewHost[j]);
        cudaFree(dfdphiHost[j]);
    }

    #if ENABLE_CUFFTMP == 1
    if (simControls.ELASTICITY)
    {
        for (j = 0; j < simDomain.numPhases; j++)
        {
            CUFFT_CHECK(cufftXtFree(phiElDev[j]));
            CUFFT_CHECK(cufftXtFree(dfEldphi[j]));
        }
    }
    #endif

    cudaFree(phiDev);
    cudaFree(dfdphiDev);
    cudaFree(phiNewDev);

    for (j = 0; j < simDomain.numPhases*(simDomain.numComponents-1); j++)
    {
        cudaFree(phaseCompHost[j]);
    }

    cudaFree(phaseCompDev);

    free(compHost);
    free(compNewHost);

    if (simControls.FUNCTION_F == 2)
    {
        free(muHost);
        free(mu);
    }

    free(phiHost);
    free(dfdphiHost);
    free(phiNewHost);

    free(phaseCompHost);

    freeVars(&simDomain, &simControls, &simParams);

    free(phi);
    free(comp);

    #if ENABLE_CUFFTMP == 1
    if (simControls.ELASTICITY)
    {
        for (i = 0; i < simDomain.numPhases*simDomain.numPhases; i++)
        {
            cudaFree(BpqHost[i]);
        }

        cudaFree(BpqDev);
        free(BpqHost);

        CUFFT_CHECK(cufftDestroy(plan));
        CUDA_CHECK(cudaStreamDestroy(stream));
    }
    #endif

    MPI_Finalize();

    cudaDeviceReset();

    return 0;
}
