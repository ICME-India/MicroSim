#ifndef COMPUTEELASTIC_CUH_
#define COMPUTEELASTIC_CUH_

#ifndef ENABLE_CUFFTMP
#define ENABLE_CUFFTMP 0
#endif

#if ENABLE_CUFFTMP == 1

#include "structures.h"
#include "defines.h"
#include "functionH.cuh"

#include <numeric>
#include <vector>
#include <complex>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cufft.h>
#include <cufftMp.h>
#include <mpi.h>
#include <sys/stat.h>

/*
 * Wrapper function for a kernel that will return kx, ky, kz arrays for Bn calculations
 * This is the simplest way to ensure Bn aligns with the cuFFTXt data structure
 */
void calc_k(double *kx, double *ky, double *kz, cudaLibXtDesc *phiXt,
            long my_nx, long my_ny, long my_nz,
            domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
            cudaStream_t stream);

/*
 * Wrapper function to compute h(phi) and move data from MicroSim array to cuFFTMp data structures.
 * To be called from main function.
 */
void moveTocudaLibXtDesc(double **phi, cudaLibXtDesc *phiXt[MAX_NUM_PHASES],
                         domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                         cudaStream_t stream, dim3 gridSize, dim3 blockSize, MPI_Comm comm);

/*
 * Wrapper function to compute h(phi) and move data from MicroSim array to cuFFTMp data structures.
 * To be called from main function.
 */
void moveFromcudaLibXtDesc(double **dfdphi, cudaLibXtDesc *dfdphiXt[MAX_NUM_PHASES],
                           domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                           cudaStream_t stream, dim3 gridSize, dim3 blockSize);

/*
 * Wrapper function for all Fourier-space elastic calculations.
 * To be called from main function.
 */
void computeDrivingForce_Elastic(cudaLibXtDesc *phi[MAX_NUM_PHASES], cudaLibXtDesc *dfeldphi[MAX_NUM_PHASES], double **Bpq,
                                 domainInfo simDomain, controls simControls, simParameters simParams, subdomainInfo subdomain,
                                 cudaStream_t stream, MPI_Comm comm);

#endif // ENABLE_CUFFTMP
#endif // ifndef
