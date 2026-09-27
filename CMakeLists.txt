cmake_minimum_required(VERSION 3.18)
project(MicroSim VERSION 1.0 LANGUAGES C CXX)

# Global options
option(BUILD_GRAND_POTENTIAL_SERIAL  "Build the Grand-Potential Serial solver"    ON)
option(BUILD_GRAND_POTENTIAL_MPI     "Build the Grand-Potential MPI solver"       ON)
option(BUILD_KKS_CUFFT               "Build the KKS cuFFT (CUDA) solver"          ON)
option(BUILD_KKS_FD_CUDA_MPI         "Build the KKS FD CUDA MPI solver"           ON)
option(BUILD_KKS_OPENCL              "Build the KKS OpenCL solver"                ON)
option(BUILD_CAHN_HILLIARD_FFT_2D    "Build the Cahn-Hilliard FFT 2D solver"      ON)

# Append local cmake/ directory to the module search path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(BUILD_GRAND_POTENTIAL_SERIAL)
    add_subdirectory(Grand_Potential_Serial)
endif()

if(BUILD_GRAND_POTENTIAL_MPI)
    add_subdirectory(Grand_Potential_MPI)
endif()

if(BUILD_KKS_CUFFT)
    add_subdirectory(KKS_CuFFT)
endif()

if(BUILD_KKS_FD_CUDA_MPI)
    add_subdirectory(KKS_FD_CUDA_MPI)
endif()

if(BUILD_KKS_OPENCL)
    add_subdirectory(KKS_OpenCl)
endif()

if(BUILD_CAHN_HILLIARD_FFT_2D)
    add_subdirectory(Cahn_Hilliard_FFT_2D)
endif()
